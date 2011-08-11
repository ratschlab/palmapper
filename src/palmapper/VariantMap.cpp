#include <palmapper/VariantMap.h>
#include <string>
#include <list>
#include <map>
#include <palmapper/Genome.h>
#include <palmapper/Config.h>
#include <stdlib.h> 
#include <palmapper/Util.h>
#include <pthread.h>


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

void VariantMap::insert_variant(int chr, int pos, int ref_len, int variant_len, const std::string & ref_str, const std::string & variant_str, int conf_count, const std::string & read_id, int read_pos)
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
	j.read_id = read_id ;
	j.read_pos= read_pos;
	j.used_to_map=false;
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

	std::deque<Variant>::iterator it = my_lower_bound(variantlist[chr].begin(), variantlist[chr].end(), j.position) ;
	
	for (; it!=variantlist[chr].end(); it++)
	{
		if (variant_cmp(j, *it)<0)
		{
			variantlist[chr].insert(it, j);
			unlock() ;
			return;
		}
		if (variant_cmp(j, *it)==0)
		{
			if (j.used_to_map){
				j.id=next_variant_id;
				next_variant_id++;
				variantlist[chr].insert(it, j);
			}
			else{	
				(*it).conf_count += j.conf_count ;
				(*it).non_conf_count += j.non_conf_count ;
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
			
		insert_variant(chr_idx, position-1, ref_len, variant_len, ref_str, variant_str, 0, "",-1);
		variant_lines++ ;
	}		

	fclose(fd) ;

	fprintf(stdout, "read %i variant lines (checked %i)\n", variant_lines, variant_lines_checked) ;

	return 0 ;

}

int VariantMap::report_to_sdi(std::string &sdi_fname, bool used_to_map)
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
			if (used_to_map && ! (*it).used_to_map)
				continue;
			
			if (((*it).conf_count<2 || (double)(*it).conf_count/(double)(*it).non_conf_count<0.2) && !used_to_map)
				continue ;
			
			std::string ref_str = (*it).ref_str ;
			if (ref_str.size()==0)
				ref_str+='-' ;
			std::string variant_str = (*it).variant_str ;
			if (variant_str.size()==0)
				variant_str+='-' ;
			
			fprintf(fd,"%s\t%i\t%i\t%s\t%s\t%i\t%i\t%s\t%i\n",
					chr, (*it).position+1, (*it).variant_len-(*it).ref_len, ref_str.c_str(), variant_str.c_str(), (*it).conf_count, (*it).non_conf_count, (*it).read_id.c_str(),(*it).read_pos+1);
			nb_variants++;
		}
	}
	fclose(fd) ;
	fprintf(stdout, "reported %i variants\n", nb_variants) ;	

	unlock() ;

	return 0;
	
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
