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
#include <algorithm>

VariantMap::VariantMap(Genome const &genome_)
{
	genome = &genome_ ;
	unsigned int nbchr = genome->nrChromosomes();
	
	variantlist = new std::deque<Variant>[nbchr];
	next_variant_id =0;
	known_variants_limit =-1;
	
	int ret = pthread_mutex_init(&variant_mutex, NULL) ;
	assert(ret==0) ;

	validate_variants=true ;
	exit_on_validation_error=false ;
	insert_unsorted=false ;
	max_variant_len=50 ;
	
}

VariantMap::~VariantMap()
{
	for (unsigned int i=0; i<genome->nrChromosomes(); i++)
		variantlist[i].clear() ;
	
	delete[] variantlist;	
}


void VariantMap::filter_variants(int min_conf_count, double max_nonconf_ratio, std::vector<std::string> & accept_sources, int filter_by_map, const GenomeMaps & genomemaps) 
{
	fprintf(stdout, "Filtering variants, requiring\n* %i as minimum confirmation count\n* %1.2f as the ratio of confirmed vs. non-confirmed\n* Accepting %ld specific sources (independent of conditions above)\n", 
			min_conf_count, max_nonconf_ratio, accept_sources.size()) ;
	if (filter_by_map>=0)
		fprintf(stdout, "* requiring variation next to mapped read or annotated exon with distance at most %i bp\n", filter_by_map) ;

	int nbchr = genome->nrChromosomes();
	int N=0, T=0;
	 
	for (int i=0; i<nbchr; i++) 
	{
		int n=0,t=0 ;
		
		std::deque<Variant> filtered ;
		for (unsigned int j=0; j<variantlist[i].size(); j++)
		{
			n++ ;
			bool take=false ;
			if (variantlist[i][j].conf_count>=min_conf_count)
				take = true ;

			if (((double)variantlist[i][j].non_conf_count/(double)variantlist[i][j].conf_count)>max_nonconf_ratio)
				take = false ;
			
			if (!take)
			{
				for (unsigned int k=0; k<accept_sources.size(); k++)
					if (accept_sources[k]==variantlist[i][j].read_id)
						take = true ;
			}
			if (take && filter_by_map>=0)
			{
				bool map=false ;
				for (int p=-filter_by_map; p<=filter_by_map && !map; p++)
					for (int l=variantlist[i][j].position; l<=variantlist[i][j].end_position && !map; l++)
						if (l+p>=0 && l+p<(int)genome->chromosome(i).length())
							map |= genomemaps.CHR_MAP(genome->chromosome(i), l+p) ;
				if (!map)
					take=false ;
			}
			if (take)
			{
				t++ ;
				filtered.push_back(variantlist[i][j]) ;
			}
		}
		fprintf(stdout, "%s: analyzed %i variants, accepted %i variants\n", genome->chromosome(i).desc(), n, t) ;
		variantlist[i].clear() ;
		variantlist[i]=filtered ;
		N+=n ;
		T+=t ;
	}
	fprintf(stdout, "All: analyzed %i variants, accepted %i variants\n", N, T) ;
}

void VariantMap::report_non_variant(const Chromosome * chr, std::vector<int> & aligned_positions, std::vector<int> & exons, int no_gap_end) 
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
	if (no_gap_end<0)
		no_gap_end=0 ;
	
	// mark aligned positions
	for (unsigned int i=no_gap_end; i<aligned_positions.size()-no_gap_end; i++)
		map[aligned_positions[i]-start_pos]=true ;

	// mark region in vicinity of splice sites as "aligned"
	for (unsigned int i=2; i<exons.size(); i+=2)
	{
		const int intron_region = 10 ;
		
		for (int j=exons[i-1]; j<exons[i]+1 && j<exons[i-1]+intron_region; j++)
			map[j-start_pos] = true ;
		for (int j=exons[i]; j>=exons[i-1] && j>exons[i]-intron_region; j--)
			map[j-start_pos] = true ;
	}

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

void VariantMap::insert_variant(int chr, int pos, int ref_len, int variant_len, const std::string & ref_str, const std::string & variant_str, int conf_count, int non_conf_count, int used_count, int non_used_count, const std::string & read_id, int read_pos, int read_len, const char* flank)
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
	j.non_conf_count = non_conf_count ;
	j.used_count = used_count ;
	j.non_used_count = non_used_count ;
	j.read_id = read_id ;
	j.read_pos= read_pos;
	j.read_len= read_len;
	
	insert_variant(j, chr, flank) ;
}

inline int min(int a, int b)
{
	if (a<b)
		return a ;
	return b ;
}

int VariantMap::update_variant(int index, int chr, const Variant &v,const char *flank)
{

	if (validate_variants)
		if (!validate_variant(v, chr, flank))
			return 0;
	
	if (v.variant_len>max_variant_len)
		return 0;
	
	for (unsigned int i=index; i<variantlist[chr].size(); i++)
	{
		
		if (v.id == variantlist[chr][i].id){
			variantlist[chr][i].used_count += v.used_count ;
			variantlist[chr][i].non_used_count += v.non_used_count ;
			variantlist[chr][i].conf_count += v.conf_count ;
			variantlist[chr][i].non_conf_count += v.non_conf_count ;
			int old_dist = min(variantlist[chr][i].read_pos, variantlist[chr][i].read_len-variantlist[chr][i].read_pos) ;
			int new_dist = min(v.read_pos, v.read_len-v.read_pos) ;
			if	((new_dist>old_dist && variantlist[chr][i].read_len>0) || variantlist[chr][i].read_id.size()==0 || v.read_len<=0 )
			{
				variantlist[chr][i].read_id = v.read_id;
				variantlist[chr][i].read_pos = v.read_pos;
				variantlist[chr][i].read_len = v.read_len;
			}
			return 1;
		}
		
		if (v.position > variantlist[chr][i].position)
			break;
	}
	
	return 0;
	
}

bool VariantMap::validate_variant(const Variant & j, int chr, const char *flank) const
{
	if (j.ref_len>0)
	{
		std::string genome_str = "" ;
		for (int i=-6; i<j.ref_len+6; i++)
		{
			if (i==0 || i==j.ref_len)
				genome_str+='|' ;
			genome_str+=genome->chromosome(chr)[j.position+i] ;
		}

		for (int i=0; i<j.ref_len; i++)
		{
			if (genome->chromosome(chr)[j.position+i]!=j.ref_str[i] && 
				j.ref_str[i]!='N' && j.ref_str[i]!='Y' && j.ref_str[i]!='W' && j.ref_str[i]!='K' && j.ref_str[i]!='S' && j.ref_str[i]!='M' && j.ref_str[i]!='R' && j.ref_str[i]!='D')
			{
				FILE* of = stdout ;
				if (exit_on_validation_error)
				{
					of=stderr ;
					fprintf(of, "ERROR: variant map disagrees with genome: %i\t%i\tgenome=%c\tref=%c\tvariant=%s\n", i, j.position+i, genome->chromosome(chr)[j.position+i], j.ref_str[i], j.variant_str.c_str()) ;
				}
				else
					fprintf(of, "WARNING: variant map disagrees with genome: %i\t%i\tgenome=%c\tref=%c\tvariant=%s\n", i, j.position+i, genome->chromosome(chr)[j.position+i], j.ref_str[i], j.variant_str.c_str()) ;
				if (j.type==pt_SNP)
					fprintf(of, "SNP\t%s\t%i\t%c\t%c\t%s\n", genome->chromosome(chr).desc(), j.position, j.ref_str[0], j.variant_str[0], genome_str.c_str()) ;
				if (j.type==pt_deletion)
					fprintf(of, "Deletion\t%s\t%i\t%i\t%s\t%s\n", genome->chromosome(chr).desc(), j.position, j.end_position, j.ref_str.c_str(), genome_str.c_str()) ;
			
				if (exit_on_validation_error)
				{
					fprintf(of, "writing current variants to variants.tmp\n") ;
					report_to_sdi(std::string("variants.tmp")) ;
					exit(-1) ;
				}
				else
					return false ;
			}
		}
	}
	if (flank!=NULL)
	{
		std::string genome_str = "" ;
		for (int i=-4; i<4; i++)
		{
			if (i==0)
				genome_str+='|' ;
			genome_str+=genome->chromosome(chr)[j.position+i] ;
		}
		if ((genome_str[3]!=flank[0] && flank[0]!='N') || (genome_str[5]!=flank[1] && flank[1]!='N'))
		{
			if (exit_on_validation_error)
			{
				fprintf(stderr, "ERROR: flanking region does not match: %s %s\nvariant: %s \n", flank, genome_str.c_str(), j.variant_str.c_str()) ;
				fprintf(stderr, "writing current variants to variants.tmp\n") ;
				report_to_sdi(std::string("variants.tmp")) ;
				exit(-1) ;
			} else
			{
				fprintf(stdout, "WARNING: flanking region does not match: %s %s\n", flank, genome_str.c_str()) ;
				return false ;
			}
		}
	}
	//fprintf(stdout, "variant ok\n") ;
	
	return true ;
}


void VariantMap::insert_variant(Variant & j, int chr, const char* flank)
{
	if (validate_variants)
		if (!validate_variant(j, chr, flank))
			return ;

	if (j.variant_len>max_variant_len)
		return ;
	

	if (insert_unsorted || variantlist[chr].empty())
	{
		j.id=next_variant_id;

		lock() ;
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
			lock() ;
			next_variant_id++;
			variantlist[chr].insert(it, j);
			unlock() ;
			return;
		}
		if (variant_cmp(j, *it)==0)
		{
			(*it).used_count += j.used_count ;
			(*it).non_used_count += j.non_used_count ;
			(*it).conf_count += j.conf_count ;
			(*it).non_conf_count += j.non_conf_count ;
			int old_dist = min((*it).read_pos, (*it).read_len-(*it).read_pos) ;
			int new_dist = min(j.read_pos, j.read_len-j.read_pos) ;
			if	((new_dist>old_dist && (*it).read_len>0) || (*it).read_id.size()==0 || j.read_len<=0 )
			{
				(*it).read_id = j.read_id;
				(*it).read_pos = j.read_pos;
				(*it).read_len = j.read_len;
			}
			return;
		}
		continue;
	}

	j.id=next_variant_id;

	lock() ;
	next_variant_id++;
	variantlist[chr].push_back(j);
	unlock() ;

	return ;
}

int VariantMap::init_from_sdi(const std::string &sdi_fname)
{

	fprintf(stdout, "initializing genome variant list with SDI file %s\n", sdi_fname.c_str()) ;

	FILE * fd=Util::openFile(sdi_fname.c_str(), "r") ;
	if (!fd)
		return -1 ;
	int variant_lines = 0, variant_lines_checked = 0 ;

	while (!feof(fd))
	{
		char chr_name[1001]="", ref_str[100001]="", variant_str[100001]="", buf[250000]="", prop[1001]="", source_id[1001]=""  ;
		int position, lendiff, read_pos=-1, read_len=-1, conf_count=0, non_conf_count=0, used_count=0, non_used_count=0 ;
		
		Util::skip_comment_lines(fd) ;
		
		if (fgets(buf, 250000, fd)==NULL)
			break ; 

		//Scan sdi line
		//int num = sscanf(buf, "%1000s\t%i\t%i\t%100000s\t%100000s\t%1000s\t%1000s\n", chr_name, &position, &lendiff, ref_str, variant_str, tmp, tmp) ;  

		int num1 = sscanf(buf,"%1000s\t%i\t%i\t%100000s\t%100000s\t%i\t%i\t%i\t%i\t%1000s\t%i/%i\t%1000s\n",
						 chr_name, &position, &lendiff, ref_str, variant_str, &conf_count, &non_conf_count, &used_count,&non_used_count, source_id, &read_pos, &read_len, prop);
		// compatibility with old format
		int num = sscanf(buf,"%1000s\t%i\t%i\t%100000s\t%100000s\t%i\t%i\t%i\t%1000s\t%i/%i\t%1000s\n",
						 chr_name, &position, &lendiff, ref_str, variant_str, &conf_count, &non_conf_count, &used_count, source_id, &read_pos, &read_len, prop);
		if (num1>num)
		{
			num = sscanf(buf,"%1000s\t%i\t%i\t%100000s\t%100000s\t%i\t%i\t%i\t%i\t%1000s\t%i/%i\t%1000s\n",
						 chr_name, &position, &lendiff, ref_str, variant_str, &conf_count, &non_conf_count, &used_count,&non_used_count, source_id, &read_pos, &read_len, prop);
		}

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
		if (ref_len>0 && validate_variants && false)
		{
			//fprintf(stdout, "pos=%i\tref_len=%i\tvariant_len=%i\n", pos, ref_len, variant_len) ;
			for (int i=0; i<ref_len; i++)
			{
				if (genome->chromosome(chr_idx)[position+i-1]!=ref_str[i] && 
					ref_str[i]!='N' && ref_str[i]!='Y' && ref_str[i]!='W' && ref_str[i]!='K' && ref_str[i]!='S' && ref_str[i]!='M' && ref_str[i]!='R' && ref_str[i]!='D')
				{
					if (exit_on_validation_error)
					{
						fprintf(stderr, "ERROR: variant map disagrees with genome: %i\t%i\t%c\t%c\t%s\n%s\n", i, position+i, genome->chromosome(chr_idx)[position+i-1], ref_str[i], ref_str, buf) ;
						exit(-1) ;
					}
					else
						fprintf(stdout, "WARNING: variant map disagrees with genome: %i\t%i\t%c\t%c\t%s\t%s\n", i, position+i, genome->chromosome(chr_idx)[position+i-1], ref_str[i], ref_str, buf) ;
				}
			}
			variant_lines_checked++ ;
		}
		
		//SDI file does not come from PALMapper and does not provide all fields
		if (num <12){
			insert_variant(chr_idx, position-1, ref_len, variant_len, ref_str, variant_str, 0, 0, 0,0, "", -2, -1);
		}
		//SDI file from PALMapper with counter values
		else{
			insert_variant(chr_idx, position-1, ref_len, variant_len, ref_str, variant_str, conf_count, non_conf_count, used_count,non_used_count, source_id, read_pos-1, read_len);
		}
		
		variant_lines++ ;

		if (variant_lines%10000==0)
		{
			long pos = ftell(fd) ;
			fprintf(stdout, "num_variants=%i\t\t%ld Mb read\r", variant_lines, pos/1024/1024) ; 
		}
	}		

	fclose(fd) ;

	fprintf(stdout, "read %i variants (checked %i)\n", variant_lines, variant_lines_checked) ;

	return 0 ;

}

int VariantMap::init_from_snp(const std::string &sdi_fname, const std::vector<std::string> & take_lines)
{
	std::string ref_line ="" ;
	int ref_line_idx=-1 ;
	
	std::vector<std::string> take_lines_names ;
	std::vector<int> take_lines_idx ;
	
	fprintf(stdout, "initializing genome variant list with SNP file %s\n", sdi_fname.c_str()) ;
	
	FILE * fd=Util::openFile(sdi_fname.c_str(), "r") ;
	if (!fd)
		return -1 ;
	int variant_lines = 0 ;
	
	{
		char list[250001], buf[250001] ;
		if (fgets(buf, 250000, fd)==NULL)
		{
			fprintf(stdout, "Error reading file\n") ;
			exit(-1) ;
		}
		int num = sscanf(buf,"Chromosome,Position,%100000s", list);
		if (num<1)
		{
			fprintf(stderr, "SNP header only contained %i columns (>=3 expected), aborting \n%s\n", num, buf) ;
			exit(-1) ;
		}
		int column=0 ;
		std::string line="" ;
		for (unsigned int i=0; i<strlen(list); i++) 
		{
			if (list[i]==',')
			{
				for (unsigned int j=0; j<take_lines.size(); j++)
					if (take_lines[j]==line || take_lines[j]=="*")
					{
						take_lines_idx.push_back(column) ;
						take_lines_names.push_back(line) ;
					}
				
				if (ref_line==line)
					ref_line_idx=column ;
				line="" ;
				column++ ;
			}
			else
				line+=list[i] ;
		}
	}
	
	while (!feof(fd))
	{
		char chr_name[1001]="", ref_str[100001]="", variant_str[100001]="", buf[250001] ;
		int position=0 ;
		char list[250001] ;
		
		Util::skip_comment_lines(fd) ;
		
		if (fgets(buf, 250000, fd)==NULL)
			break ;  

		int num = sscanf(buf,"%1000[^,],%i,%100000s", chr_name, &position, list);
		
		if (num<3)
		{
			if (feof(fd))
				break ;
			fprintf(stderr, "SNP line only contained %i columns (>=3 expected), aborting (%s)\n%s\n", num, chr_name, buf) ;
			exit(-1) ;
		}

		int chr_idx = genome->find_desc(chr_name) ;
		if (chr_idx==-1)
		{
			chr_idx = genome->find_desc((std::string("Chr") + std::string(chr_name)).c_str()) ;
			if (chr_idx==-1)
			{
				fprintf(stderr, "chromosome %s not found. known chromosome names:\n", chr_name) ;
				genome->print_desc(stderr) ;
				fclose(fd) ;
				return -1 ;
			}
		}

		int column=0 ;
		std::string line="" ;
		char SNP[take_lines_idx.size()] ;
		char ref='N' ;
		for (unsigned int j=0; j<take_lines_idx.size(); j++)
			SNP[j]='N' ;
		
		for (unsigned int i=0; i<strlen(list); i++) 
		{
			if (list[i]==',')
			{
				for (unsigned int j=0; j<take_lines_idx.size(); j++)
					if (take_lines_idx[j]==column)
					{
						assert(line.size()==1) ;
						SNP[j] = line[0] ;
					}
				if (column==ref_line_idx)
				{
					assert(line.size()==1) ;
					ref=line[0] ;
				}
				line="" ;
				column++ ;
			}
			else
				line+=list[i] ;
		}
				
		for (unsigned int s=0; s<take_lines_idx.size(); s++)
		{
			int spos=position ;
		
			strcpy(ref_str, "N") ;
			if (ref!='N')
				ref_str[0]=ref ;
			else
				ref_str[0]=genome->chromosome(chr_idx)[spos-1] ;
				
			strcpy(variant_str, "N") ;
			variant_str[0]=SNP[s] ;

			if (ref_str[0]!=variant_str[0])
			{
				int ref_len = strlen(ref_str) ; assert(ref_len==1) ;
				int variant_len = strlen(variant_str) ; assert(variant_len==1) ;
				
				std::string source="SNP:" ; 
				source+= take_lines_names[s] ;	
				char flank[3]="NN" ;
				insert_variant(chr_idx, spos-1, ref_len, variant_len, ref_str, variant_str, 1, 0, 0, 0, source.c_str(), -1, -1, flank) ;
			
				variant_lines++ ;
				
				if (variant_lines%10000==0)
				{
					long pos = ftell(fd) ;
					fprintf(stdout, "num_variants=%i\t\t%ld Mb read\r", variant_lines, pos/1024/1024) ; 
				}
			}
		}		
	}
	
	fclose(fd) ;

	fprintf(stdout, "read %i SNPs\n", variant_lines) ;

	return 0 ;

}

int VariantMap::init_from_samtools(const std::string &sdi_fname)
{
	fprintf(stdout, "initializing genome variant list with indel/samtools file %s\n", sdi_fname.c_str()) ;

	FILE * fd=Util::openFile(sdi_fname.c_str(), "r") ;
	if (!fd)
		return -1 ;
	int variant_lines = 0, variant_lines_checked = 0 ;

	while (!feof(fd))
	{
		char chr_name[1001]="", ref_str[100001]="", variant_str[100001]="", buf[250001] ;
		int position=0 ;
		char str[2][1000001], rest[250001] ;
		
		Util::skip_comment_lines(fd) ;
		
		if (fgets(buf, 250000, fd)==NULL)
			break ; 

		int num = sscanf(buf,"%1000s\t%i\t*\t%100000s\t%s", chr_name, &position, str[0], rest);
		
		if (num<3)
		{
			if (feof(fd))
				break ;
			fprintf(stdout, "sdi line only contained %i columns (>=4 expected), aborting (%s)\n%s\n", num, chr_name, buf) ;
		}

		int chr_idx = genome->find_desc(chr_name) ;
		if (chr_idx==-1)
		{
			fprintf(stderr, "chromosome %s not found. known chromosome names:\n", chr_name) ;
			genome->print_desc(stderr) ;
			fclose(fd) ;
			return -1 ;
		}
		for (unsigned int i=0; i<strlen(str[0]); i++)
		{
			if (str[0][i]=='/')
			{
				strcpy(str[1], &str[0][i+1]) ;
				str[0][i]=0 ;
				break ;
			}
		}
		
		for (int s=0; s<2; s++)
		{
			int spos=position ;
			if (s==1 && strcmp(str[0],str[s])==0)
				continue ;
			if (strcmp(str[s],"*")==0)
				continue ;
			
			assert(str[s][0]=='-' || str[s][0]=='+') ;

			char flank[3]="NN" ;
		
			if (str[s][0]=='+') 
			{
				strcpy(ref_str, "") ;
				if (strcmp(str[s], "*")!=0)
					strcpy(variant_str, str[s]+1) ;
				//flank[0]=genome->chromosome(chr_idx)[spos-1] ;
				//flank[1]=genome->chromosome(chr_idx)[spos-1+1] ;
			}
			else
			{
				strcpy(variant_str, "") ;
				if (strcmp(str[s], "*")!=0)
					strcpy(ref_str, str[s]+1) ;
				spos++ ;
			}

			int ref_len = strlen(ref_str) ;
			int variant_len = strlen(variant_str) ;
			
			// validate variants on genome sequence
			if (ref_len>0 && validate_variants && false)
			{
				//fprintf(stdout, "pos=%i\tref_len=%i\tvariant_len=%i\n", pos, ref_len, variant_len) ;
				for (int i=0; i<ref_len; i++)
				{
					if (genome->chromosome(chr_idx)[spos+i-1]!=ref_str[i] && 
						ref_str[i]!='N' && ref_str[i]!='Y' && ref_str[i]!='W' && ref_str[i]!='K' && ref_str[i]!='S' && ref_str[i]!='M' && ref_str[i]!='R' && ref_str[i]!='D')
					{
						if (exit_on_validation_error)
						{
							fprintf(stderr, "ERROR: variant map disagrees with genome: %i\t%i\t%c\t%c\t%s\n%s\n", i, spos+i, genome->chromosome(chr_idx)[spos+i-1], ref_str[i], ref_str, buf) ;
							exit(-1) ;
						}
						else
							fprintf(stdout, "WARNING: variant map disagrees with genome: %i\t%i\t%c\t%c\t%s\t%s\n", i, spos+i, genome->chromosome(chr_idx)[spos+i-1], ref_str[i], ref_str, buf) ;
					}
				}
				variant_lines_checked++ ;
			}
			
			insert_variant(chr_idx, spos-1, ref_len, variant_len, ref_str, variant_str, 1, 0, 0, 0, "indel", -1, -1, flank) ;
			
			variant_lines++ ;
			
			if (variant_lines%10000==0)
			{
				long pos = ftell(fd) ;
				fprintf(stdout, "num_variants=%i\t\t%ld Mb read\r", variant_lines, pos/1024/1024) ; 
			}
		}		
	}
	
	fclose(fd) ;

	fprintf(stdout, "read %i indels (checked %i)\n", variant_lines, variant_lines_checked) ;

	return 0 ;

}


int VariantMap::report_to_file(const std::string &filename) const
{
	unsigned int pos_ext=filename.rfind('.');

	if (pos_ext != std::string::npos)  
	{
		std::string extension(filename.substr(pos_ext+1));
		
		if (extension.compare("sdi")==0 ||extension.compare("SDI")==0)
			return report_to_sdi(filename) ;
		if (extension.compare("bin")==0 ||extension.compare("BIN")==0)
			return report_to_bin(filename) ; 
	}
	fprintf(stderr, "ERROR: unknown extension of variant file %s\n", filename.c_str()) ;
	exit(-1) ;
	
	return -1 ;
}

int VariantMap::report_to_sdi(const std::string &sdi_fname)  const
{
	int nb_variants=0;
	
	fprintf(stdout, "report genome variants in SDI file %s\n", sdi_fname.c_str()) ;
	
	FILE * fd=Util::openFile(sdi_fname.c_str(), "w") ;
	if (!fd)
		return -1 ;	

	const int max_len=10 ;
	int num_del=0 ;
	int num_ins=0 ;
	int num_snp=0 ;
	int num_sub=0 ;
	std::vector<int> num_del_dist(max_len+2, 0) ;
	std::vector<int> num_ins_dist(max_len+2, 0) ;

	fprintf(fd,"#chromosome\tposition\tlen_diff\tref_seq\tvariant_seq\tconf_count\tnon_conf_count\tused_count\tnon_used_count\tsource\tread position/len\tstatus\n");

	for (unsigned int i=0; i<genome->nrChromosomes(); i++)
	{
		const char * chr= genome->get_desc(i);
		std::deque<Variant>::iterator it;
		
		for (it=variantlist[i].begin(); it!=variantlist[i].end(); it++)
		{			
			if ((*it).type==pt_SNP)
				num_snp++ ;
			if ((*it).type==pt_deletion)
			{
				num_del++ ;
				if ((*it).ref_len>max_len)
					num_del_dist[max_len+1]++ ;
				else
					num_del_dist[(*it).ref_len]++ ;
			}
			if ((*it).type==pt_insertion)
			{
				num_ins++ ;
				if ((*it).variant_len>max_len)
					num_ins_dist[max_len+1]++ ;
				else
					num_ins_dist[(*it).variant_len]++ ;
			}
			if ((*it).type==pt_substitution)
				num_sub++ ;
			
			std::string ref_str = (*it).ref_str ;
			if (ref_str.size()==0)
				ref_str+='-' ;
			std::string variant_str = (*it).variant_str ;
			if (variant_str.size()==0)
				variant_str+='-' ;
			
			fprintf(fd,"%s\t%i\t%i\t%s\t%s\t%i\t%i\t%i\t%i\t%s\t%i/%i\t%s\n",
					chr, (*it).position+1, (*it).variant_len-(*it).ref_len, ref_str.c_str(), variant_str.c_str(), (*it).conf_count, (*it).non_conf_count, (*it).used_count,(*it).non_used_count, 
					(*it).read_id.c_str(),(*it).read_pos+1, (*it).read_len, (*it).id<=known_variants_limit?"known":"discovered");
			nb_variants++;
		}
	}
	fclose(fd) ;

	fprintf(stdout, "reported %i variants\n", nb_variants) ;	
	fprintf(stdout, "* %i SNPs\n", num_snp) ;
	fprintf(stdout, "* %i substitutions\n", num_sub) ;
	fprintf(stdout, "* %i deletions\n\t", num_del) ;
	for (int i=0; i<max_len+1; i++)
		fprintf(stdout, "%i:%i ", i, num_del_dist[i]) ;
	fprintf(stdout, " >%i:%i\n", max_len, num_del_dist[max_len+1]) ;
	fprintf(stdout, "* %i insertions\n\t", num_ins) ;
	for (int i=0; i<max_len+1; i++)
		fprintf(stdout, "%i:%i ", i, num_ins_dist[i]) ;
	fprintf(stdout, " >%i:%i\n", max_len, num_ins_dist[max_len+1]) ;
	
	return 0;
	
}

int VariantMap::init_from_files(std::string &fnames)
{

	int previousfound=0;
	int found=fnames.find(",");
	std::string filename;
	bool has_maf_file = false ;
	
	while (true)
	{
		if (found >=0)
			filename = fnames.substr(previousfound, found-previousfound);
		else
			filename=fnames.substr(previousfound);

		unsigned int pos_ext=filename.rfind('.');
		//fprintf(stdout, "fn=%s,%i\n", filename.c_str(), pos_ext) ;
		VariantInputEnum ext=unknown;
		if (pos_ext != std::string::npos)
		{
			std::string extension(filename.substr(pos_ext+1));
			//fprintf(stdout, "ext=%s\n", extension.c_str()) ;
			
			if (extension.compare("sdi")==0 ||extension.compare("SDI")==0)
				ext=sdi;
			if (extension.compare("maf")==0 ||extension.compare("MAF")==0)
				ext=maf;
			if (extension.compare("indel")==0 ||extension.compare("samtools")==0)
				ext=samtools;
			if (extension.compare("snp")==0 ||extension.compare("SNP")==0)
				ext=snp;
			if (extension.compare("bin")==0 ||extension.compare("BIN")==0)
				ext=bingz;
		}
		if (ext==unknown)
		{
			fprintf(stderr,	"ERROR: Variant input %s has an unknown format\n", (char*)filename.c_str());
			exit(1) ;
		}
		if ( ext==maf && _config.MAF_REF_NAME.length()<=0)
		{
			fprintf(stderr,	"ERROR: Need reference genome name (with -maf-ref option) to initialize variants from the maf file %s\n", (char*)filename.c_str());
			exit(1) ;
		}
		
		int ret = 0 ;
		if (ext == sdi)
			init_from_sdi(filename);
		if (ext == bingz)
			init_from_bin(filename);
		if (ext == samtools)
			init_from_samtools(filename);
		if (ext == snp)
			init_from_snp(filename, _config.VARIANT_SNP_TAKE_LINES);
		if (ext == maf)
		{
			init_from_maf(filename, _config.MAF_REF_NAME);
			has_maf_file = true ;
		}

		if (ret!=0)
			return ret;
		check_variant_order() ;

		if (found<=0)
			break ;
		
		previousfound=found+1;
		found=fnames.find(",", found+1);
	}
	 
	known_variants_limit=next_variant_id-1;
	
	if (!has_maf_file && _config.MAF_REF_NAME.length()>0)
		fprintf(stdout, "WARNING: maf reference given, but no maf file as input\n") ;
	
	return 0;
}


int VariantMap::insert_variants_from_multiple_alignments(std::string &ref_align,int ref_len, std::vector<std::string> &variant_align, std::vector<std::string> &variant_name, 
														 int start_position, int chr_len, int chr_idx, char strand)
{

	int start_variant=-1;
	enum polytype type_variant = pt_unknown ;
	std::string diff_seq;
	int num_variants=0;
	
	//Compare reference sequence with alignments from different genomes 
	for (unsigned int i=0;i<variant_align.size();i++){

		std::string variant_seq=variant_align[i];
		std::string source_name=variant_name[i] ;
		
		if (variant_seq.length() != ref_align.length())
			continue;
		int ref_position=-1;

		//fprintf(stdout, "%s\n%s\n",(char*)ref_align.c_str(), (char*)variant_seq.c_str()) ;
		bool ref_started = false ;
		bool variant_started = false ;
		char flank_variant='N' ;
		
		for (unsigned j=0;j<ref_align.length();j++){

			//Gap on both sequences
			if (ref_align[j]=='-' && variant_seq[j]=='-')
				continue;
			if (ref_align[j]!='-' && !ref_started)
			{
				ref_started=true ;
				start_variant=-1 ;
				flank_variant=ref_align[j] ;
				diff_seq.clear() ;
			}
			if (variant_seq[j]!='-' && !variant_started)
			{
				variant_started=true ;
				start_variant=-1 ;
				diff_seq.clear() ;
			}
			
			//Match
			if (ref_align[j]==variant_seq[j])
			{
				ref_position++;
				//Next position after an insertion or deletion
				if (start_variant!=-1)
				{
					if (type_variant == pt_insertion)
					{
						char flank[3]="NN" ;
						
						if (strand =='+')
						{
							//flank[0]=ref_align[start_variant] ;
							flank[1]=ref_align[j] ; 
							
							insert_variant(chr_idx, start_position + ref_position, 0, diff_seq.length(), "", diff_seq, 1, 0, 0, 0,source_name, -1, -1, flank);
						}
						else
						{
							flank[0]=ref_align[j] ;
							std::string rflank=QPalma::complement(flank) ;
							
							insert_variant(chr_idx, chr_len - start_position - ref_position , 0,  diff_seq.length(), "", QPalma::reverse(QPalma::complement(diff_seq)), 
										   1, 0, 0, 0,source_name, -1, -1, rflank.c_str());
						}
					}
					if (type_variant == pt_deletion){
						if (strand =='+')
							insert_variant(chr_idx, start_position +start_variant, diff_seq.length(), 0,  diff_seq, "", 1, 0, 0, 0,source_name,-1,-1);
						else
							insert_variant(chr_idx, chr_len - start_position - ref_position , diff_seq.length(), 0, QPalma::reverse(QPalma::complement(diff_seq)), "", 1, 0, 0, 0,
										   source_name, -1, -1);
					}
					start_variant=-1;
					flank_variant='N' ;
					diff_seq.clear();
					num_variants++;
				}
				continue;
			}
			
			//Insertion on ref
			if (ref_align[j]=='-')
			{
				//This position corresponds to the end of a deletion
				if (start_variant!=-1 && type_variant == pt_deletion)
				{
					if (strand =='+')
						insert_variant(chr_idx, start_position +start_variant,  diff_seq.length(), 0,  diff_seq, "", 1, 0, 0, 0,source_name, -1, -1);
					else
						insert_variant(chr_idx, chr_len - start_position - ref_position -1, diff_seq.length(), 0, QPalma::reverse(QPalma::complement(diff_seq)), 
									   "", 1, 0, 0, 0,source_name, -1, -1);
					diff_seq.clear();
					start_variant=ref_position;
					type_variant=pt_insertion;
					diff_seq+=variant_seq[j];
					num_variants++;
				}
				else{
					//New insertion
					if (start_variant==-1 && ref_started && variant_started)
					{
						start_variant=ref_position;
						type_variant=pt_insertion;
						assert(j>=0) ;
						int k=j ;
						while (ref_align[k-1]=='-' && k>1)
							k-- ;
						flank_variant=ref_align[k-1] ;
					}
					diff_seq+=variant_seq[j];
				}
				
				continue;
			}
			
			//Deletion on ref
			if (variant_seq[j]=='-')
			{
				ref_position++;
				//This position corresponds to the end of an insertion
				if (start_variant!=-1 && type_variant == pt_insertion)
				{
					char flank[3]="NN" ;
					
					if (strand =='+')
					{
						// TODO: check flanks
						insert_variant(chr_idx, start_position + ref_position, 0, diff_seq.length(), "", diff_seq, 1, 0, 0, 0,source_name, -1, -1, flank);
					}
					else
					{
						// TODO: check flanks
						insert_variant(chr_idx, chr_len - start_position - ref_position - 1, 0,  diff_seq.length(), "", QPalma::reverse(QPalma::complement(diff_seq)), 
									   1, 0, 0, 0,source_name,-1,-1, flank);
					}
					
					diff_seq.clear();
					start_variant=ref_position;
					type_variant=pt_deletion;
					diff_seq+=ref_align[j];
					num_variants++;
				}
				else{
					//New deletion
					if (start_variant==-1 && variant_started && ref_started)
					{
						start_variant=ref_position;
						type_variant=pt_deletion;
						assert(j>=0) ;
						int k=j ;
						while (ref_align[j-1]=='-' && k>1)
							k-- ;
						flank_variant=ref_align[k-1] ;
					}
					diff_seq+=ref_align[j];
				}
				
				continue;
			}


			//Mismatch => SNP
			ref_position++;
			//Next position after an insertion or deletion
			if (start_variant!=-1)
			{
				if (type_variant == pt_insertion)
				{
					char flank[3]="NN" ;
					
					if (strand =='+')
					{
						//TODO: check flanks
						insert_variant(chr_idx, start_position + ref_position, 0, diff_seq.length(), "", diff_seq, 1, 0, 0, 0,source_name,-1,-1, flank);
					}
					else
					{
						//TODO: check flanks
						insert_variant(chr_idx, chr_len - start_position - ref_position - 1, 0,  diff_seq.length(), "", QPalma::reverse(QPalma::complement(diff_seq)), 
									   1, 0, 0, 0,source_name,-1,-1, flank);
					}
				}
				if (type_variant == pt_deletion){
					if (strand =='+')
						insert_variant(chr_idx, start_position +start_variant,  diff_seq.length(), 0,  diff_seq, "", 1, 0, 0, 0,source_name,-1,-1);
					else
						insert_variant(chr_idx, chr_len - start_position - ref_position, diff_seq.length(), 0, QPalma::reverse(QPalma::complement(diff_seq)), 
									   "", 1, 0, 0, 0,source_name,-1,-1);
				}
				start_variant=-1;
				diff_seq.clear();
				num_variants++;
			}

			//Insert SNP
			if (strand =='+')
			{
				std::string ref;
				std::string var;
				char r=ref_align[j];
				char v=variant_seq[j];
				ref+=r;
				var+=v;
				
				insert_variant(chr_idx, start_position + ref_position, 1, 1,  ref, var, 1, 0, 0, 0,source_name,-1,-1);
			}
			else{
				std::string ref;
				std::string var;
				char r=QPalma::complement(ref_align[j]);
				char v=QPalma::complement(variant_seq[j]);
				ref+=r;
				var+=v;
				insert_variant(chr_idx, chr_len - start_position - ref_position -1, 1,1, ref, var, 1, 0, 0, 0,source_name,-1,-1);
			}
			
			num_variants++;
		}

		// Next position after an insertion or deletion
		// disabled: do not insert indels at the trailing ends
		if (false && start_variant!=-1)
		{
			if (type_variant == pt_insertion)
			{
				char flank[3]="NN" ;
				
				if (strand =='+')
				{
					//TODO: check flanks
					insert_variant(chr_idx, start_position + ref_position, 0, diff_seq.length(), "", diff_seq, 1, 0, 0, 0,source_name,-1,-1);
				}
				else
				{
					//TODO: check flanks
					insert_variant(chr_idx, chr_len - start_position - ref_position - 1, 0,  diff_seq.length(), "", QPalma::reverse(QPalma::complement(diff_seq)), 
								   1, 0, 0, 0,source_name,-1,-1, flank);
				}
			}
			if (type_variant == pt_deletion)
			{
				if (strand =='+')
					insert_variant(chr_idx, start_position +start_variant,  diff_seq.length(), 0,  diff_seq, "", 1, 0, 0, 0,source_name,-1,-1);
				else
					insert_variant(chr_idx, chr_len - start_position - ref_position -1, diff_seq.length(), 0, QPalma::reverse(QPalma::complement(diff_seq)), 
								   "", 1, 0, 0, 0,source_name,-1,-1);
			}
			start_variant=-1;
			diff_seq.clear();
			num_variants++;
		}//End comparison reference against one other genome
		

	}//End loop over non reference genomes
	
	return num_variants;
	
}


bool compare_variants(const Variant &a, const Variant &b)
{
	return (VariantMap::variant_cmp(a,b)<0) ;
}

int VariantMap::init_from_maf(const std::string &maf_fname, const std::string &ref_genome)
{

	fprintf(stdout, "initializing genome variant list with MAF file %s\n", maf_fname.c_str()) ;
	insert_unsorted = true ;
	
	FILE * fd=Util::openFile(maf_fname.c_str(), "r") ;
	if (!fd)
		return -1 ;
	int variant_lines = 0;
	
	std::string ref_align="";
	std::vector<std::string> variant_align;
	std::vector<std::string> variant_name;
	char strand='+';
	int start_position=-1;
	int chr_idx=-1;
	int ref_len=-1;
	int ref_chr_len=-1 ;
	int num_blocks = 0 ;
	
	while (!feof(fd))
	{
		
		Util::skip_comment_lines(fd) ;
		
		/**************************************/
		/* Parse MAF line                     */
		/* Initialize multiple alignments     */
		/**************************************/

		//Scan MAF line type
		char name_maf[1000], alignment_maf[100000];
		int position_maf, len_maf, len_chr ;
		char strand_maf;
		char type_maf;

		int num = fscanf(fd, "%c\t%1000s\t%i\t%i\t%c\t%i\t%100000s\n", &type_maf, name_maf, &position_maf, &len_maf, &strand_maf, &len_chr, alignment_maf) ;  
		//fprintf(stdout, "%c\t%s\t%i\t%i\t%c\t%i\t%s\n", type_maf, name_maf, position_maf, len_maf, strand_maf, len_chr, alignment_maf) ;  
		if (num<1)
		{
			if (feof(fd))
				break ;
			fprintf(stdout, "maf line only contained %i columns (7 expected), aborting\n", num) ;
			continue;
			
		}

		//New alignment block: get variants from the previous block 
		if (type_maf == 'a')
		{
			//Block with reference name
			if(!ref_align.empty())
				variant_lines+=insert_variants_from_multiple_alignments(ref_align, ref_len, variant_align, variant_name, start_position, ref_chr_len, chr_idx,strand);
			ref_align.clear();
			variant_align.clear();
			variant_name.clear();
			num_blocks++ ;
			
			if (num_blocks%10000==0)
			{
				long pos = ftell(fd) ;
				fprintf(stdout, "num_blocks=%i\t\t%ld Mb read\r", num_blocks, pos/1024/1024) ; 
			}
			//if (num_blocks>1000)
			//	break ;
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
			ref_chr_len = len_chr ;
			ref_align.assign(alignment_maf);
			std::transform(ref_align.begin(), ref_align.end(),ref_align.begin(), ::toupper);			
		}
		else{
			std::string var_seq(alignment_maf);
			std::transform(var_seq.begin(), var_seq.end(),var_seq.begin(), ::toupper);			
			variant_align.push_back(var_seq);
			variant_name.push_back(genome_name) ;
			var_seq.clear();
		}
		

	}		

	//Last block if exists
	if (!ref_align.empty())
	{
		variant_lines+=insert_variants_from_multiple_alignments(ref_align,ref_len,variant_align,variant_name, start_position, ref_chr_len, chr_idx,strand);
		ref_align.clear();
		variant_align.clear();
		variant_name.clear();
	}
	
	fclose(fd) ;

	//lock() ;
	fprintf(stdout, "Inserted %i variants\n", variant_lines) ;
	fflush(stdout) ;
	
	fprintf(stdout, "Sorting...") ;
	for (int i=0; i<(int)genome->nrChromosomes(); i++)
	{
		sort(variantlist[i].begin(), variantlist[i].end(), compare_variants);
		fprintf(stdout, ".") ;
		fflush(stdout) ;
	}
	fprintf(stdout, "Done.\n") ;
	insert_unsorted = false ;
	//unlock() ;

	fflush(stdout) ;

	return 0 ;

}

void VariantMap::transcribe_gff(const std::string & gff_input, const std::string & fasta_output)
{
    
	fprintf(stdout, "\ninitializing trancripts from GFF file %s\n", gff_input.c_str()) ;

	FILE * fd=Util::openFile(gff_input.c_str(), "r") ;
	FILE * fdo=Util::openFile(fasta_output.c_str(), "w") ;
    std::string last_parent;
    std::string exon_seq;
    std::string transcript_seq;

    char * ret;

	if (!fd || !fdo)
		return ;
    
    char curr_strand = '*';
    char last_strand = '*';
    int start = 0, end = 0;

	while (!feof(fd))
	{
		Util::skip_comment_lines(fd) ;

        char line[10000]="", chr_name[1000]="", source[1000]="", type[1000]="", properties[1000]="", strand='+' ;
        ret = fgets(line, sizeof(line), fd);
        if (!ret)
            break;
        char * sl = strtok(line, "\t");

        int idx = 0;
        while (sl != NULL) {
            if (idx == 0) {
                strcpy(chr_name, sl) ; 
            } else if (idx == 1) {
                strcpy(source, sl) ; 
            } else if (idx == 2) {
                strcpy(type, sl) ; 
            } else if (idx == 3) {
                start = atoi(sl) - 1;
            } else if (idx == 4) {
                end = atoi(sl) - 1;
            } else if (idx == 6) {
                strand = sl[0]; 
            } else if (idx == 8) {
                strcpy(properties, sl) ; 
            }
            sl = strtok(NULL, "\t");    
            idx++;
        }

		if (idx<9)
		{
			if (feof(fd))
				break ;
			fprintf(stdout, "gff line in file %s only contained %i columns, aborting\n", gff_input.c_str(), idx) ;
            continue ;
			//exit(-1) ;
		}


        int chr_idx = genome->find_desc(chr_name) ;
        if (chr_idx==-1)
        {
            fprintf(stderr, "chromosome %s not found. known chromosome names:\n", chr_name) ;
            genome->print_desc(stderr) ;
            fclose(fd) ;
            return ;
        }

		if (strcmp(type, "exon")==0)
		{
			//Get Parent name
			std::string parent;
			std::string tmp(properties);
			int pos_parent=tmp.find("Parent=");
			if (pos_parent <0){
				fprintf(stderr, "No parent information for exon in gff3 file\n") ;
				fclose(fd) ;
				exit(-1) ;
			}
			else
				pos_parent += 7;
		    if (tmp.find("Transcript:") > 0) 
                pos_parent += 11;

			int pos_parent_end = tmp.find(";", pos_parent);
			if (pos_parent_end < 0)
				parent=tmp.substr(pos_parent);
			else
				parent=tmp.substr(pos_parent,pos_parent_end-pos_parent);
            
            if (parent[parent.size() - 1] == '\n')
                parent = parent.substr(0, parent.size() - 1);

            curr_strand = strand;
            if (last_strand == '*')
                last_strand = curr_strand;

            int e_idx = start;
            // collect all variants overlapping to this exon
            std::vector<std::deque<Variant>::iterator> variant_list;
            std::deque<Variant>::iterator it = my_lower_bound(this->variantlist[chr_idx].begin(), this->variantlist[chr_idx].end(), e_idx-100) ;
            for (; it != this->variantlist[chr_idx].end(); it++) {
                if (it->position + it->ref_len >= start && it->position <= end) {
                    // possibility to filter variants based on confirmation
                    //if (it->conf_count > 0) 
                        variant_list.push_back(it);
                }
                if (it->position > end) 
                    break;
            }

            // define set of combatible variants in the variant list

            // incorporate compatible variants into exon
            size_t v_idx = 0;
            while (e_idx <= end && v_idx < variant_list.size()) 
			{
                it = variant_list.at(v_idx);
                if (it->position <= e_idx) 
				{
                    int offset = e_idx - it->position;
                    switch (it->type) 
					{
					case pt_substitution: 
					case pt_SNP: 
					{ 
						e_idx += it->ref_len - offset;
                        if (it->ref_len - offset >= it->variant_len)
                            exon_seq = exon_seq.append(it->variant_str);
                        else if (offset > 0)
                            exon_seq = exon_seq.append(it->variant_str.substr(it->variant_len + offset - it->ref_len, it->ref_len - offset));
                        else
                            exon_seq = exon_seq.append(it->variant_str);
                        v_idx++;
					}
					break; 
					case pt_deletion: e_idx = e_idx + it->ref_len - offset <=end?e_idx + it->ref_len - offset:end; v_idx++; break;
					case pt_insertion: 
					{
						if (offset == 0) {
							exon_seq = exon_seq.append(it->variant_str);
                            v_idx++;
                        }
					}
					break;
					default: assert(0) ;
                    } 
                } else {
                    exon_seq += genome->chromosome(chr_idx)[e_idx]; 
                    e_idx++;
                }
            }
            while (e_idx <=end) {
                exon_seq += genome->chromosome(chr_idx)[e_idx]; 
                e_idx++;
            }

            if (last_parent.size() == 0)
                last_parent = parent;

            if (!parent.compare(last_parent)) {
                transcript_seq = transcript_seq.append(exon_seq);
            } else {
                fprintf(fdo, ">%s\n", last_parent.c_str());
                if (transcript_seq.size() > 0) {
                    if (last_strand == '-')
                        transcript_seq = QPalma::reverse(QPalma::complement(transcript_seq));
                    for (size_t i = 0; i < transcript_seq.size(); i+= 80)
                    {
                        fprintf(fdo,"%s\n", transcript_seq.substr(i, 80).c_str());
                    }
                    transcript_seq = exon_seq;
                    last_strand = curr_strand;
                }
                last_parent = parent;
            }
            exon_seq.clear();
        } 
    }
    if (transcript_seq.size() > 0) {
        fprintf(fdo, ">%s\n", last_parent.c_str());
        if (last_strand == '-')
            transcript_seq = QPalma::reverse(QPalma::complement(transcript_seq));
        for (size_t i = 0; i < transcript_seq.size(); i+= 80)
        {
            fprintf(fdo,"%s\n", transcript_seq.substr(i, 80).c_str());
        }
        transcript_seq.clear();
    }
}

