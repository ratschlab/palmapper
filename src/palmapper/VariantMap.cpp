#include <palmapper/VariantMap.h>
#include <string>
#include <list>
#include <map>
#include <palmapper/Genome.h>
#include <palmapper/Config.h>
#include <stdlib.h> 
#include <palmapper/Util.h>
#include <pthread.h>
#include <palmapper/QPalma.h>

VariantMap::VariantMap(Genome const &genome_)
{
	genome = &genome_ ;
	unsigned int nbchr = genome->nrChromosomes();
	
	variantlist = new std::deque<Variant>[nbchr];
	next_variant_id =0;
	limit_known_variants =-1;
	
	int ret = pthread_mutex_init(&variant_mutex, NULL) ;// = PTHREAD_MUTEX_INITIALIZER ;
	assert(ret==0) ;
	
}

VariantMap::~VariantMap()
{
	delete[] variantlist;	
}


void VariantMap::report_non_variant(const Chromosome * chr, std::vector<int> & aligned_positions, std::vector<int> & exons) 
{
	if (aligned_positions.size()==0)
		return ;
	int start_pos = aligned_positions[0] ;
	int end_pos = aligned_positions[aligned_positions.size()-1] ;
	if (start_pos>end_pos)
	{
		end_pos = aligned_positions[0] ;
		start_pos = aligned_positions[aligned_positions.size()-1] ;
	}
	if (start_pos>exons[0])
		start_pos = exons[0] ;
	if (end_pos<exons[exons.size()-1])
		end_pos = exons[exons.size()-1] ;
	
	std::vector<bool> map(end_pos-start_pos+1, false) ;

	// mark aligned positions
	for (unsigned int i=0; i<aligned_positions.size(); i++)
		map[aligned_positions[i]-start_pos]=true ;
	// mark region in vicinity of splice sites as "aligned"
	/*for (unsigned int i=2; i<exons.size(); i+=2)
	{
		const int intron_region = 10 ;
		
		for (int j=exons[i-1]; j<exons[i]+1 && j<exons[i-1]+intron_region; j++)
			map[j-start_pos] = true ;
		for (int j=exons[i]; j>=exons[i-1] && j>exons[i]-intron_region; j--)
			map[j-start_pos] = true ;
			}*/

	std::deque<Variant>::iterator it = my_lower_bound(variantlist[chr->nr()].begin(), variantlist[chr->nr()].end(), start_pos) ;
	
	for (; it!=variantlist[chr->nr()].end(); it++)
	{
		if ((*it).position>=end_pos)
			break ;
		if ((*it).position<start_pos)
			continue ;
		assert((*it).position>=start_pos && (*it).position<end_pos) ;
		
		// mark SNPs contraticting the alignment
		if ((*it).type==pt_SNP && map[(*it).position-start_pos])
			(*it).non_conf_count++ ;

		// mark deletions and substitution contraticting the alignment
		// we say its contraticting, the whole deleted/substituted region was "aligned"
		if ((*it).type==pt_deletion || (*it).type==pt_substitution) 
		{
			bool all_aligned=true ;
			
			for (int i=(*it).position-start_pos; i<(*it).end_position-start_pos; i++)
			{
				if (i<0) 
					continue ;
				if (i >= end_pos-start_pos || !map[i])
				{
					all_aligned=false ;
					break ;
				}
			}
			//fprintf(stdout, "delsub %i: %i-%i (%i,%i), %i\n", (*it).type, (*it).position, (*it).end_position, (*it).ref_len, (*it).variant_len, all_aligned) ;
			
			if (all_aligned)
				(*it).non_conf_count++ ;
		}

		// mark insertions contradicting the alignment
		// we say that it is contradicting, if the adjacent positions are "aligned"
		if ((*it).type==pt_insertion) 
		{
			if (map[(*it).position-start_pos] && (*it).position-start_pos+1<end_pos && map[(*it).position-start_pos+1])
				(*it).non_conf_count++ ;
				}
	}
}

void VariantMap::insert_variant(int chr, int pos, int ref_len, int variant_len, const std::string & ref_str, const std::string & variant_str, int conf_count, int used_count, const std::string & read_id, int read_pos)
{
	enum polytype pt = pt_unknown ;
	if (ref_len==1 && variant_len==1)
		pt = pt_SNP ;
	else if (ref_len==0 && variant_len>0)
		pt = pt_insertion ;
	else if (variant_len==0 && ref_len>0)
		pt = pt_deletion ;
	else if (variant_len>0 && ref_len>0)
		pt = pt_substitution ;

	Variant j;
	int end = pos+ref_len ;

	j.type = pt ;
	j.position = pos ;
	j.end_position=end ;
	j.ref_len = ref_len ;
	j.variant_len = variant_len ;
	j.ref_str=ref_str ;
	j.variant_str = variant_str ;
	j.conf_count = conf_count ;
	j.non_conf_count = 0 ;
	j.used_count = used_count ;
	j.read_id = read_id ;
	j.read_pos= read_pos;
	
	insert_variant(j, chr) ;
}

void VariantMap::insert_variant(Variant & j, int chr)
{
	lock() ;

	if (variantlist[chr].empty())
	{
		
		j.id=next_variant_id;
		next_variant_id++;
		variantlist[chr].push_back(j);

		unlock() ;
		return;
	}

	std::deque<Variant>::iterator it = my_lower_bound(variantlist[chr].begin(), variantlist[chr].end(), j.position);
	if (it !=variantlist[chr].begin())
		it--;
	
	for (; it!=variantlist[chr].end(); it++)
	{
		if (variant_cmp(j, *it)<0)
		{
			j.id=next_variant_id;
			next_variant_id++;
			variantlist[chr].insert(it, j);
			unlock() ;
			return;
		}
		if (variant_cmp(j, *it)==0)
		{
			(*it).used_count += j.used_count ;
			(*it).conf_count += j.conf_count ;
			(*it).non_conf_count += j.non_conf_count ;
			if	((*it).read_id ==""){
				(*it).read_id = j.read_id;
				(*it).read_pos = j.read_pos;
			}
			
			
			unlock() ;
			return;
		}
		continue;
	}

	j.id=next_variant_id;
	next_variant_id++;
	variantlist[chr].push_back(j);

	unlock() ;
	return ;
}

int VariantMap::init_from_sdi(std::string &sdi_fname)
{

	fprintf(stdout, "initializing genome variant list with SDI file %s\n", sdi_fname.c_str()) ;

	FILE * fd=Util::openFile(sdi_fname.c_str(), "r") ;
	if (!fd)
		return -1 ;
	int variant_lines = 0, variant_lines_checked = 0 ;

	while (!feof(fd))
	{
		char chr_name[1000], ref_str[100000], variant_str[100000], tmp[1000], buf[250000] ;
		int position, lendiff ;
		
		Util::skip_comment_lines(fd) ;
		
		if (fgets(buf, 250000, fd)==NULL)
			break ; 

		//Scan sdi line
		int num = sscanf(buf, "%1000s\t%i\t%i\t%100000s\t%100000s\t%1000s\t%1000s\n", chr_name, &position, &lendiff, ref_str, variant_str, tmp, tmp) ;  
		if (num<5)
		{
			if (feof(fd))
				break ;
			fprintf(stdout, "sdi line only contained %i columns (5 expected), aborting (%s)\n", num, chr_name) ;
		}

		int chr_idx = genome->find_desc(chr_name) ;
		if (chr_idx==-1)
		{
			fprintf(stderr, "chromosome %s not found. known chromosome names:\n", chr_name) ;
			genome->print_desc(stderr) ;
			fclose(fd) ;
			return -1 ;
		}
		if (strcmp(ref_str, "-")==0)
			strcpy(ref_str, "") ;
		if (strcmp(variant_str, "-")==0)
			strcpy(variant_str, "") ;

		int ref_len = strlen(ref_str) ;
		int variant_len = strlen(variant_str) ;
		assert(lendiff==variant_len-ref_len) ;

		// validate variants on genome sequence
		if (ref_len>0)
		{
			//fprintf(stdout, "pos=%i\tref_len=%i\tvariant_len=%i\n", pos, ref_len, variant_len) ;
			for (int i=0; i<ref_len; i++)
			{
				if (genome->chromosome(chr_idx)[position+i-1]!=ref_str[i] && 
					ref_str[i]!='N' && ref_str[i]!='Y' && ref_str[i]!='W' && ref_str[i]!='K' && ref_str[i]!='S' && ref_str[i]!='M' && ref_str[i]!='R' && ref_str[i]!='D')
				{
					fprintf(stdout, "error: variant map disagrees with genome: %i\t%i\t%c\t%c\t%s\n", i, position+i, genome->chromosome(chr_idx)[position+i-1], ref_str[i], ref_str) ;
					exit(-1) ;
				}
			}
			variant_lines_checked++ ;
		}
			
		insert_variant(chr_idx, position-1, ref_len, variant_len, ref_str, variant_str, 0, 0,"",-1);
		variant_lines++ ;
	}		

	fclose(fd) ;

	fprintf(stdout, "read %i variant lines (checked %i)\n", variant_lines, variant_lines_checked) ;

	return 0 ;

}

int VariantMap::report_to_sdi(std::string &sdi_fname)
{
	lock() ;

	int nb_variants=0;
	
	fprintf(stdout, "report genome list in SDI file %s\n", sdi_fname.c_str()) ;
	
	FILE * fd=Util::openFile(sdi_fname.c_str(), "w") ;
	if (!fd)
		return -1 ;	

	for (unsigned int i=0; i<genome->nrChromosomes(); i++){
		
		const char * chr= genome->get_desc(i);
		std::deque<Variant>::iterator it;
		
		for (it=variantlist[i].begin(); it!=variantlist[i].end(); it++)
		{			
			
			if (false && (((*it).conf_count<2 || (double)(*it).conf_count/(double)(*it).non_conf_count<0.2) && (*it).used_count<=0))
				continue ;
			
			std::string ref_str = (*it).ref_str ;
			if (ref_str.size()==0)
				ref_str+='-' ;
			std::string variant_str = (*it).variant_str ;
			if (variant_str.size()==0)
				variant_str+='-' ;
			
			fprintf(fd,"%s\t%i\t%i\t%s\t%s\t%i\t%i\t%i\t%s\t%i\t%s\n",
					chr, (*it).position+1, (*it).variant_len-(*it).ref_len, ref_str.c_str(), variant_str.c_str(), (*it).conf_count, (*it).non_conf_count,(*it).used_count, (*it).read_id.c_str(),(*it).read_pos+1,(*it).id<=limit_known_variants?"known":"discovered");
			nb_variants++;
		}
	}
	fclose(fd) ;
	fprintf(stdout, "reported %i variants\n", nb_variants) ;	

	unlock() ;

	return 0;
	
}

int  VariantMap::init_from_mafs(std::string &maf_fname,std::string &ref_name)
{
	int previousfound=0;
	int found=maf_fname.find(",");
	std::string filename;
	
	while (found >= 0)
	{
		
		filename = maf_fname.substr(previousfound, found-previousfound);
		int ret = init_from_maf(filename,ref_name);
		if (ret!=0)
			return ret;
		check_variant_order() ;
		
		previousfound=found+1;
		found=maf_fname.find(",",found+1);
	}
	
	filename=maf_fname.substr(previousfound);
	int ret=init_from_maf(filename,ref_name);
	limit_known_variants=next_variant_id-1;
	
	if (ret!=0)
		return  ret;
	
	return ret;
}

int VariantMap::init_from_sdis(std::string &sdi_fname)
{

	int previousfound=0;
	int found=sdi_fname.find(",");
	std::string filename;
	
	while (found >= 0)
	{
		
		filename = sdi_fname.substr(previousfound, found-previousfound);
		int ret = init_from_sdi(filename);
		if (ret!=0)
			return ret;
		check_variant_order() ;
		
		previousfound=found+1;
		found=sdi_fname.find(",",found+1);
	}
	
	filename=sdi_fname.substr(previousfound);
	int ret=init_from_sdi(filename);
	limit_known_variants=next_variant_id-1;
	
	if (ret!=0)
		return  ret;
	
	return ret;
	
}


int VariantMap::insert_variants_from_multiple_alignments(std::string ref_align,int ref_len, std::vector<std::string> variant_align,int start_position,int chr_idx, char strand)
{

	int start_variant=-1;
	enum polytype type_variant = pt_unknown ;
	std::string diff_seq;
	int num_variants=0;
	
	

	//Compare reference sequence with alignments from different genomes 
	for (unsigned int i=0;i<variant_align.size();i++){

		std::string variant_seq=variant_align[i];
		if (variant_seq.length() != ref_align.length())
			continue;
		int ref_position=-1;
		fprintf(stdout, "%s\n%s\n",(char*)ref_align.c_str(), (char*)variant_seq.c_str()) ;
		for (unsigned j=0;j<ref_align.length();j++){

			//Gap on both sequences
			if (ref_align[j]=='-' && variant_seq[j]=='-')
				continue;
			
			//Match
			if (ref_align[j]==variant_seq[j]){
				ref_position++;
				//Next position after an insertion or deletion
				if(start_variant!=-1){
					if (type_variant == pt_insertion){
						if (strand =='+')
							insert_variant(chr_idx, start_position -1 + ref_position, 0, diff_seq.length(), "", diff_seq, 0, 0,"",-1);
						else
							insert_variant(chr_idx, start_position -1 + ref_len-ref_position, 0,  diff_seq.length(), "", QPalma::reverse(QPalma::complement(diff_seq)), 0, 0,"",-1);
					}
					if (type_variant == pt_deletion){
						if (strand =='+')
							insert_variant(chr_idx, start_position -1 +start_variant, diff_seq.length(), 0,  diff_seq, "", 0, 0,"",-1);
						else
							insert_variant(chr_idx, start_position -1 + ref_len-1 -start_variant, diff_seq.length(), 0, QPalma::reverse(QPalma::complement(diff_seq)), "", 0, 0,"",-1);
					}
					start_variant=-1;
					diff_seq.clear();
					num_variants++;
				}
				continue;
			}
			
			//Insertion on ref
			if (ref_align[j]=='-'){
				//This position corresponds to the end of a deletion
				if (start_variant!=-1 && type_variant == pt_deletion){
					if (strand =='+')
						insert_variant(chr_idx, start_position -1 +start_variant,  diff_seq.length(), 0,  diff_seq, "", 0, 0,"",-1);
					else
						insert_variant(chr_idx, start_position -1 + ref_len-1 -start_variant, diff_seq.length(), 0, QPalma::reverse(QPalma::complement(diff_seq)), "", 0, 0,"",-1);
					diff_seq.clear();
					start_variant=ref_position;
					type_variant=pt_insertion;
					diff_seq+=variant_seq[j];
					num_variants++;
				}
				else{
					//New insertion
					if (start_variant==-1){
						start_variant=ref_position;
						type_variant=pt_insertion;
					}
					diff_seq+=variant_seq[j];
				}
				
				continue;
			}
			
			//Deletion on ref
			if (variant_seq[j]=='-'){
				ref_position++;
				//This position corresponds to the end of an insertion
				if (start_variant!=-1 && type_variant == pt_insertion){
					if (strand =='+')
						insert_variant(chr_idx, start_position -1 + ref_position, 0, diff_seq.length(), "", diff_seq, 0, 0,"",-1);
					else
						insert_variant(chr_idx, start_position -1 + ref_len-ref_position, 0,  diff_seq.length(), "", QPalma::reverse(QPalma::complement(diff_seq)), 0, 0,"",-1);
					diff_seq.clear();
					start_variant=ref_position;
					type_variant=pt_deletion;
					diff_seq+=ref_align[j];
					num_variants++;
				}
				else{
					//New deletion
					if (start_variant==-1){
						start_variant=ref_position;
						type_variant=pt_deletion;
					}
					diff_seq+=ref_align[j];
					
				}
				
				continue;
			}


			//Mismatch => SNP
			ref_position++;
			//Next position after an insertion or deletion
			if(start_variant!=-1){
				if (type_variant == pt_insertion){
					if (strand =='+')
						insert_variant(chr_idx, start_position -1 + ref_position, 0, diff_seq.length(), "", diff_seq, 0, 0,"",-1);
					else
						insert_variant(chr_idx, start_position -1 + ref_len-ref_position, 0,  diff_seq.length(), "", QPalma::reverse(QPalma::complement(diff_seq)), 0, 0,"",-1);
				}
				if (type_variant == pt_deletion){
					if (strand =='+')
						insert_variant(chr_idx, start_position -1 +start_variant,  diff_seq.length(), 0,  diff_seq, "", 0, 0,"",-1);
					else
						insert_variant(chr_idx, start_position -1 + ref_len-1 -start_variant, diff_seq.length(), 0, QPalma::reverse(QPalma::complement(diff_seq)), "", 0, 0,"",-1);
				}
				start_variant=-1;
				diff_seq.clear();
				num_variants++;
			}
			//Insert SNP
			if (strand =='+')
				insert_variant(chr_idx, start_position -1 +ref_position,  1, 1,  &ref_align[j], &variant_seq[j], 0, 0,"",-1);
			else{
				std::string ref;
				std::string var;
				char r=QPalma::complement(ref_align[j]);
				char v=QPalma::complement(variant_seq[j]);
				ref+=r;
				var+=v;
				insert_variant(chr_idx, start_position -1 + ref_len-1 -ref_position, 1,1,&ref[0], &var[0], 0, 0,"",-1);
			}
			
			num_variants++;
		}
		//Next position after an insertion or deletion
		if(start_variant!=-1){
			if (type_variant == pt_insertion){
				if (strand =='+')
					insert_variant(chr_idx, start_position -1 + ref_position, 0, diff_seq.length(), "", diff_seq, 0, 0,"",-1);
				else
					insert_variant(chr_idx, start_position -1 + ref_len-ref_position, 0,  diff_seq.length(), "", QPalma::reverse(QPalma::complement(diff_seq)), 0, 0,"",-1);
			}
			if (type_variant == pt_deletion){
				if (strand =='+')
					insert_variant(chr_idx, start_position -1 +start_variant,  diff_seq.length(), 0,  diff_seq, "", 0, 0,"",-1);
				else
					insert_variant(chr_idx, start_position -1 + ref_len-1 -start_variant, diff_seq.length(), 0, QPalma::reverse(QPalma::complement(diff_seq)), "", 0, 0,"",-1);
			}
			start_variant=-1;
			diff_seq.clear();
			num_variants++;
		}//End comparaison reference against one other genome
		

	}//End loop over non reference genomes
	
	return num_variants;
	
}


int VariantMap::init_from_maf(std::string &maf_fname,std::string &ref_genome)
{

	fprintf(stdout, "initializing genome variant list with MAF file %s\n", maf_fname.c_str()) ;

	FILE * fd=Util::openFile(maf_fname.c_str(), "r") ;
	if (!fd)
		return -1 ;
	int variant_lines = 0;
	
	std::string ref_align="";
	std::vector<std::string> variant_align;
	char strand='+';
	int start_position=-1;
	int chr_idx=-1;
	int ref_len=-1;

	
	while (!feof(fd))
	{
		
		Util::skip_comment_lines(fd) ;
		
		/**************************************/
		/* Parse MAF line                     */
		/* Initialize multiple alignments     */
		/**************************************/

		//Scan MAF line type
		char name_maf[1000], alignment_maf[100000];
		int position_maf, len_maf ;
		char strand_maf;
		char type_maf;
		int tmp;
		int num = fscanf(fd, "%c\t%1000s\t%i\t%i\t%c\t%i\t%100000s\n", &type_maf, name_maf, &position_maf, &len_maf, &strand_maf, &tmp,alignment_maf) ;  
		if (num<1)
		{
			if (feof(fd))
				break ;
			fprintf(stdout, "maf line only contained %i columns (7 expected), aborting\n", num) ;
			continue;
			
		}

		//New alignment block: get variants from the previous block
		if (type_maf == 'a'){
			//Block with reference name
			if(!ref_align.empty())
				variant_lines+=insert_variants_from_multiple_alignments(ref_align,ref_len,variant_align,start_position,chr_idx,strand);
			ref_align.clear();
			variant_align.clear();
		}
		
			

		//Line should be a sequence line and have 7 fields
		if (type_maf != 's' || num <7)
			continue;
		

		//Get genome name and chromosome
		std::string name_tmp(name_maf);
		unsigned int position_tmp = name_tmp.find_last_of('.');
		if (position_tmp == std::string::npos)
		{
			fprintf(stderr, "Name %s does not contain genome and chromosome information\n", name_maf) ;
			fclose(fd) ;
			return -1 ;
		}

		std::string genome_name=name_tmp.substr(0,position_tmp);
		std::string chr_name=name_tmp.substr(position_tmp+1);
		name_tmp.clear();
		

		//Reference case
		if (genome_name == ref_genome){
			
			int chr_tmp = genome->find_desc((char*)chr_name.c_str()) ;
			if (chr_tmp==-1)
			{
				fprintf(stderr, "chromosome %s not found. known chromosome names:\n", (char*)chr_name.c_str()) ;
				genome->print_desc(stderr) ;
				fclose(fd) ;
				return -1 ;
			}

			start_position=position_maf;
			chr_idx=chr_tmp;
			strand=strand_maf;
			ref_len=len_maf;
			ref_align.assign(alignment_maf);
			std::transform(ref_align.begin(), ref_align.end(),ref_align.begin(), ::toupper);			

			
			

		}
		else{
			std::string var_seq(alignment_maf);
			std::transform(var_seq.begin(), var_seq.end(),var_seq.begin(), ::toupper);			
			variant_align.push_back(var_seq);
			var_seq.clear();
			
		}
		

	}		

	//Last block if exists
	if (!ref_align.empty()){
		variant_lines+=insert_variants_from_multiple_alignments(ref_align,ref_len,variant_align,start_position,chr_idx,strand);
		ref_align.clear();
		variant_align.clear();
	}
	
	fclose(fd) ;

	fprintf(stdout, "Insert %i variant lines\n", variant_lines) ;

	return 0 ;

}
