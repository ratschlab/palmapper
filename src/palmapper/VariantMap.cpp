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
#include <vector>

VariantMap::VariantMap(Genome const &genome_, bool p_merge_variant_source_ids)
{
	genome = &genome_ ;
	unsigned int nbchr = genome->nrChromosomes();
	
	variantlist = new std::vector<Variant>[nbchr];
	next_variant_id =0;
	known_variants_limit =-1;
	
	int ret = pthread_mutex_init(&variant_mutex, NULL) ;
	assert(ret==0) ;
    
	validate_variants=true ;
	exit_on_validation_error=false ;
	insert_unsorted=false ;
    
#ifndef PMINDEX
	max_variant_len = _config.FILTER_VARIANT_VLEN ;
#endif
	merge_variant_source_ids=p_merge_variant_source_ids ;
}

VariantMap::~VariantMap()
{
	for (unsigned int i=0; i<genome->nrChromosomes(); i++)
		variantlist[i].clear() ;
	
	delete[] variantlist;
}

void VariantMap::filter_variants_junctions(JunctionMap & junctionmap)
{
	fprintf(stdout, "Filtering/penalizing variants in intronic regions...\n") ;
	int nbchr = genome->nrChromosomes();
	int N=0, T=0, P=0 ;
    
	const int intron_win=10 ;
	const int grey_win=2 ;
	const int inner_intron_penalty=10 ;
	
    for (int i=0; i<nbchr; i++)
    {
		//fprintf(stdout, ".") ;
		std::vector<short int> map(genome->chromosome(i).length(), 0) ;
		std::deque<Junction> & junctions = junctionmap.junctionlist[i] ;
		std::deque<Junction>::iterator it=junctions.begin() ;
		for (;it!=junctions.end(); it++)
		{
			for (int k=(*it).start-grey_win; k<(*it).end && k<(*it).start+intron_win; k++)
			{
				if (k<0)
					continue ;
				map[k]=1 ; // intron
			}
			for (int k=(*it).start+intron_win; k<(*it).end-intron_win; k++)
			{
				if (map[k]==0)
					map[k] = 2 ; // inner intron
			}
			for (int k=(*it).end+grey_win-1; k>(*it).start && k>(*it).end-intron_win; k--)
			{
				if (k>=(int)map.size())
					continue ;
				map[k]=1 ; // intron
			}
		}
		
		int n=0, t=0, p=0;
		std::vector<Variant> filtered ;
		for (unsigned int j=0; j<variantlist[i].size(); j++)
		{
			n++ ;
			bool take=false ;
			for (int l=variantlist[i][j].position; l<=variantlist[i][j].end_position && !take ; l++)
				if (l>=0 && l<(int)map.size() && (map[l]!=1))
					take = true ;
			bool apply_penalty=true ;
			for (int l=variantlist[i][j].position; l<=variantlist[i][j].end_position && apply_penalty ; l++)
				if (l>=0 && l<(int)map.size() && (map[l]==0))
					apply_penalty = false ;
			if (take)
			{
				t++ ;
				if (apply_penalty)
				{
					p++ ;
					variantlist[i][j].non_conf_count+=inner_intron_penalty ;
				}
				filtered.push_back(variantlist[i][j]) ;
			}
		}
		fprintf(stdout, "%s: analyzed %i variants, accepted %i variants (%i with penalty)\n", genome->chromosome(i).desc(), n, t, p) ;
		variantlist[i].clear() ;
		variantlist[i]=filtered ;
		N+=n ;
		T+=t ;
		P+=p ;
	}
	fprintf(stdout, "All: analyzed %i variants, accepted %i variants (%i with penalty)\n", N, T, P) ;
}

void VariantMap::unique_variant_source_ids()
{
	fprintf(stdout, "Making variant source ids unique\n") ;
    
	int nbchr = genome->nrChromosomes();
    
	for (int i=0; i<nbchr; i++)
	{
		int num_before=0 ;
		int num_after=0 ;
		
		for (unsigned int j=0; j<variantlist[i].size(); j++)
		{
			//fprintf(stdout, "before: %s\n", variantlist[i][j].read_id.c_str()) ;
            
			std::vector<std::string> sources ;
			std::string source="" ;
			for (unsigned int p=0; p<variantlist[i][j].read_id.length(); p++)
				if (variantlist[i][j].read_id[p]==',')
				{
					sources.push_back(source) ;
					source="" ;
				}
				else
					source+=variantlist[i][j].read_id[p] ;
			if (variantlist[i][j].read_id.length()>=0)
				sources.push_back(source) ;
            
			num_before+=sources.size() ;
			
			std::vector<std::string>::iterator it;
			sort(sources.begin(), sources.end()) ;
			it = unique(sources.begin(), sources.end()) ;
			sources.resize( it - sources.begin() );
            
			num_after+=sources.size() ;
            
			std::string read_id="" ;
			for (unsigned int p=0; p<sources.size(); p++)
			{
				if (p>0)
					read_id+="," ;
				read_id+=sources[p] ;
			}
			variantlist[i][j].read_id=read_id ;
			//fprintf(stdout, "after: %s\n", read_id.c_str()) ;
		}
		
		fprintf(stdout, "%s: analyzed %i variants with %i sources and %i unique sources\n", genome->chromosome(i).desc(), (int)variantlist[i].size(), num_before, num_after) ;
	}
}

void VariantMap::convert_substitutions()
{
	fprintf(stdout, "Converting substitutions to insertion+deletion\n") ;
	
	int nbchr = genome->nrChromosomes();
	int N=0, T=0;
    
	for (int i=0; i<nbchr; i++)
	{
		int n=0,t=0 ;
		
		std::vector<Variant> filtered ;
		for (unsigned int j=0; j<variantlist[i].size(); j++)
		{
			if (variantlist[i][j].type==pt_substitution)
			{
				Variant ins=variantlist[i][j], del=variantlist[i][j] ;
				ins.type=pt_insertion ;
				ins.ref_str="" ;
				ins.ref_len=0 ;
				ins.end_position=ins.position ;
				
				del.type=pt_deletion ;
				del.variant_str="" ;
				del.variant_len=0 ;
                
				filtered.push_back(ins) ;
				filtered.push_back(del) ;
				t++ ;
			}
			else
				filtered.push_back(variantlist[i][j]) ;
			n++ ; t++ ;
		}
		fprintf(stdout, "%s: analyzed %i variants, accepted %i variants\n", genome->chromosome(i).desc(), n, t) ;
		variantlist[i].clear() ;
		variantlist[i]=filtered ;
		N+=n ;
		T+=t ;
	}
	
}


void VariantMap::filter_variants(int min_source_count, int min_conf_count, double max_nonconf_ratio,
								 int min_use_count,
								 std::vector<std::string> & accept_sources, std::vector<std::string> & required_sources,
								 int max_len, int filter_by_map, const GenomeMaps & genomemaps)
{
	fprintf(stdout, "Filtering variants, requiring\n* %i as minimum confirmation count\n* %1.2f as the ratio of confirmed vs. non-confirmed\n* Accepting %ld specific sources (independent of conditions above)\n",
			min_conf_count, max_nonconf_ratio, accept_sources.size()) ;
	if (required_sources.size()>0)
		fprintf(stdout, "* requiring variation to have at least %i specific sources (e.g., reads/alleles)\n", (int)required_sources.size()) ;
	if (min_source_count>0)
		fprintf(stdout, "* requiring variation to have at least %i sources (e.g., reads/alleles)\n", min_source_count) ;
	if (min_use_count>0)
		fprintf(stdout, "* requiring variation to have at least %i uses/non-uses\n", min_use_count) ;
	if (max_len>0)
		fprintf(stdout, "* requiring variation to be shorter than %i bp\n", max_len) ;
	if (filter_by_map>=0)
		fprintf(stdout, "* requiring variation next to mapped read or annotated exon with distance at most %i bp\n", filter_by_map) ;
    
	int nbchr = genome->nrChromosomes();
	int N=0, T=0;
    
	for (int i=0; i<nbchr; i++)
	{
		int n=0,t=0 ;
		
		std::vector<Variant> filtered ;
		for (unsigned int j=0; j<variantlist[i].size(); j++)
		{
			n++ ;
			bool take=false ;
			if (variantlist[i][j].conf_count>=min_conf_count)
				take = true ;
			if (min_source_count>0 || required_sources.size()>0)
			{
				std::vector<bool> req_source_found(required_sources.size(), false) ;
				/*for (unsigned int r=0; r<required_sources.size(); r++)
                 req_source_found.push_back(false) ;*/
                
				int num_sources=0 ;
				if (variantlist[i][j].read_id.length()>0)
					num_sources++ ;
				std::string source ;
				for (unsigned int p=0; p<variantlist[i][j].read_id.length(); p++)
					if (variantlist[i][j].read_id[p]==',')
					{
						num_sources++ ;
						for (unsigned int r=0; r<required_sources.size(); r++)
							if (source==required_sources[r])
								req_source_found[r]=true ;
						source.assign("") ;
					}
					else
						source+=variantlist[i][j].read_id[p] ;
				if (source.size()>0)
					for (unsigned int r=0; r<required_sources.size(); r++)
						if (source==required_sources[r])
							req_source_found[r]=true ;
                
				if (num_sources<min_source_count)
					take = false ;
				for (unsigned int r=0; r<required_sources.size(); r++)
					if (!req_source_found[r])
						take=false ;
			}
			if (variantlist[i][j].conf_count>0)
			{
				if (((double)variantlist[i][j].non_conf_count/(double)variantlist[i][j].conf_count)>max_nonconf_ratio)
					take = false ;
			}
			else
			{
				if (variantlist[i][j].non_conf_count>0)
					take = false ;
			}
			if (take && variantlist[i][j].used_count<min_use_count && variantlist[i][j].non_used_count<min_use_count)
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
			if (take && max_len>=0)
				if (variantlist[i][j].variant_len>max_len || variantlist[i][j].ref_len>max_len)
					take=false ;
			
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

void VariantMap::report_variant(int rank, int total, Variant & j, int chr, const char* flank, bool update_only, bool ignore_variant_str_in_cmp)
{
	if (rank!=0)
		return ;
	
	insert_variant(j, chr, flank, update_only, ignore_variant_str_in_cmp) ;
}

void VariantMap::report_non_variant(int rank, int total, const Chromosome * chr, std::vector<int> & aligned_positions, std::vector<int> & exons, int no_gap_end)
{
	if (rank!=0)
		return ;
    
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
    
	lock() ;
	std::vector<Variant>::iterator it = my_lower_bound(variantlist[chr->nr()].begin(), variantlist[chr->nr()].end(), start_pos) ;
	
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
	unlock() ;
	
}

void VariantMap::insert_variant(int chr, int pos, int ref_len, int variant_len, const std::string & ref_str, const std::string & variant_str, int conf_count, int non_conf_count, int used_count, int non_used_count, const std::string & read_id, int read_pos, int read_len, const char* flank, bool update_only, bool ignore_variant_str_in_cmp)
{
	enum polytype pt = pt_unknown ;
	if (ref_len==1 && variant_len==1)
		pt = pt_SNP ;
	else if (ref_len==0 && variant_len>0)
	{
		pt = pt_insertion ;
		//fprintf(stdout, "insert: insertion\n") ;
	}
	else if (variant_len==0 && ref_len>0)
	{
		pt = pt_deletion ;
		//fprintf(stdout, "insert: deletion\n") ;
	}
	else if (variant_len>0 && ref_len>0)
	{
		pt = pt_substitution ;
		//fprintf(stdout, "insert: substitution\n") ;
	}
    
    
	/*if (variant_len>=32768 || ref_len>=32768)
     {
     fprintf(stderr, "Drop variant with invalid lengths %i/%i\n", variant_len, ref_len) ;
     return ;
     }*/
	if (ref_str==variant_str && ref_str!="N" && ref_str!="*")
    {
        fprintf(stderr, "Drop variant with equal variant and reference string %i:%i:'%s'/'%s'\n", chr, pos, variant_str.c_str(), ref_str.c_str()) ;
		return ;
    }
	
    
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
	
	insert_variant(j, chr, flank, update_only, ignore_variant_str_in_cmp) ;
}

inline int min(int a, int b)
{
	if (a<b)
		return a ;
	return b ;
}

int VariantMap::update_variant(int rank, int total, int index, int chr, const Variant &v,const char *flank)
{
	if (rank!=0)
		return -1 ;
	
	if (validate_variants)
		if (!validate_variant(v, chr, flank))
			return 0;
	
	if (v.variant_len>max_variant_len)
		return 0;
    
	lock() ;
    
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
			unlock() ;
			return 1;
		}
		
		if (v.position > variantlist[chr][i].position)
			break;
	}
    
	unlock() ;
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
			if (i>=0 && i<(int)genome->chromosome(chr).length())
				genome_str+=genome->chromosome(chr)[j.position+i] ;
			else
				genome_str+='_' ;
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
			if (j.position+i>=0 && j.position+i<(int)genome->chromosome(chr).length())
				genome_str+=genome->chromosome(chr)[j.position+i] ;
			else
				genome_str+='-' ;
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


void VariantMap::insert_variant(Variant & j, int chr, const char* flank, bool update_only, bool ignore_variant_str_in_cmp)
{
	if (validate_variants)
		if (!validate_variant(j, chr, flank))
		{
			return ;
		}
    
	if (j.variant_len>max_variant_len)
	{
		//fprintf(stdout, ">max_variant_len=%i:%i:%i\n",max_variant_len,j.variant_len,j.ref_len) ;
		return ;
	}
	if (j.variant_len<0 || j.ref_len<0)
	{
		fprintf(stderr, "Drop variant with invalid lengths %i/%i\n", j.variant_len, j.ref_len) ;
		return ;
	}
	
	lock() ;
    
	if (insert_unsorted || variantlist[chr].empty())
	{
		if (!update_only)
		{
			j.id=next_variant_id;
			
			next_variant_id++;
			variantlist[chr].push_back(j);
		}
		else
		{
			fprintf(stdout, "droping variant (1) at %i:%i:%s:%s\n", chr, j.position, j.ref_str.c_str(), j.variant_str.c_str()) ;
		}
		
		unlock() ;
		return;
	}
    
	std::vector<Variant>::iterator it = my_lower_bound(variantlist[chr].begin(), variantlist[chr].end(), j.position);
	if (it !=variantlist[chr].begin())
		it--;
	
	for (; it!=variantlist[chr].end(); it++)
	{
		if (variant_cmp(j, *it, ignore_variant_str_in_cmp)<0)
		{
			if (!update_only)
			{
				j.id=next_variant_id;
				next_variant_id++;
				variantlist[chr].insert(it, j);
			}
			else
			{
				fprintf(stdout, "droping variant (2) at %i:%i:%s:%s  >  %i:%i:[%i]:[%i]\n", chr, j.position, j.ref_str.c_str(), j.variant_str.c_str(),
						chr, (*it).position, (int)strlen((*it).ref_str.c_str()), (int)strlen((*it).variant_str.c_str())) ;
			}
			unlock() ;
			return;
		}
		if (variant_cmp(j, *it, ignore_variant_str_in_cmp)==0)
		{
			(*it).used_count += j.used_count ;
			(*it).non_used_count += j.non_used_count ;
			(*it).conf_count += j.conf_count ;
			(*it).non_conf_count += j.non_conf_count ;
			int old_dist = min((*it).read_pos, (*it).read_len-(*it).read_pos) ;
			int new_dist = min(j.read_pos, j.read_len-j.read_pos) ;
			if	((new_dist>old_dist && (*it).read_len>0) || (*it).read_id.size()==0 || j.read_len<=0 )
			{
				if (!merge_variant_source_ids)
					(*it).read_id = j.read_id;
				(*it).read_pos = j.read_pos;
				(*it).read_len = j.read_len;
			}
			if (merge_variant_source_ids)
			{
				if ((*it).read_id.size()>0)
					(*it).read_id = (*it).read_id + "," + j.read_id ;
				else
					(*it).read_id = j.read_id ;
			}
			
			unlock() ;
			return;
		}
		continue;
	}
    
	if (!update_only)
	{
		j.id=next_variant_id;
		
		next_variant_id++;
		variantlist[chr].push_back(j);
	}
	else
	{
		fprintf(stdout, "droping variant (3) at %i:%i:%s:%s\n", chr, j.position, j.ref_str.c_str(), j.variant_str.c_str()) ;
	}
	
	unlock() ;
    
	return ;
}

int VariantMap::init_from_vcf(const std::string &vcf_fname)
{
    fprintf(stdout, "initializing genome variant list with VCF file %s\n", vcf_fname.c_str()) ;
    
    FILE * fd=Util::openFile(vcf_fname.c_str(), "r") ;
    if (!fd)
        return -1 ;
    int variant_lines = 0, variant_lines_checked = 0 ;
    const int max_buf_len = 10000000 ;
	const int max_field_len = 500000 ;
    
    char * buf=(char*)malloc(max_buf_len+1) ;
    strcpy(buf, "") ;
    
    while (!feof(fd))
    {
        //variant object requirements for palmapper
        char chr_name[1001]="", ref_str[max_field_len+1]="",
        variant_str[max_field_len+1]="", source_id[1001]="",
        varElems[1001]="" ;
        char * varPch ;
		int position=0, lendiff=0, read_pos=-1, read_len=-1, conf_count=0,
        non_conf_count=0, used_count=0, non_used_count=0, chr_idx=0,
        variant_len=0, ref_len=0;
        std::vector<std::string> strainRefVec, variantVec, strainVec ;
        
        if (fgets(buf, max_buf_len, fd)==NULL)
			break ;
        
        //skip lines that contain VCF metadata
        if (buf[0] == '#' && buf[1] == '#')
        {
            continue ;
        }
        
        //store the strain names
        else if (buf[0] == '#' && isalpha(buf[1]))
        {
            char * headerPch ;
            int headCntInt = 0 ;
            headerPch = strtok (buf, "\t") ;
            
            while (headerPch != NULL)
            {
                if (headCntInt > 8)   //strains begin in column 9
                {
                    std::string strainNmStr = headerPch ;
                    std::string whitespaces (" \t\f\v\n\r") ;
                    size_t found ;
                    found = strainNmStr.find_last_not_of(whitespaces) ;
                    if (found!=std::string::npos)
                        strainNmStr.erase(found+1) ;
                    else
                        strainNmStr.clear() ;
                    strainRefVec.push_back(strainNmStr) ;
                }
                headCntInt++ ;
                
                headerPch = strtok (NULL, "\t") ;
            }
        } else if (buf[0] != '#')
        {
            char * elemPch ;
            int elemCntInt = 0 ;
            
            elemPch = strtok (buf, "\t") ;
            
            while (elemPch != NULL)
            {
                if (elemCntInt == 0) {
                    chr_idx = atoi(elemPch) ;
                    //strcpy(chr_name, elemPch) ;
                } else if (elemCntInt == 1)
                {
                    //position is 1 based; palmapper is 0 based
                    position = atoi(elemPch) - 1 ;
                } else if (elemCntInt == 3)
                {
                    strcpy(ref_str, elemPch) ;
                } else if (elemCntInt == 4)     //determine number of variants
                {
                    strcpy(varElems, elemPch) ;
                }
                else if (elemCntInt > 8)
                {
                    strainVec.push_back(elemPch) ; //store strain information in strainVec
                } else {
                }
                
                elemCntInt++;
                
                elemPch = strtok (NULL, "\t") ;
            }   // process variant elements in first columns
            
            // Tokenize different variants
            varPch = strtok(varElems, ",") ;
            while (varPch != NULL)
            {
                variantVec.push_back(varPch);
                varPch = strtok (NULL, ",") ;
            }
            //update variants
            
            ref_len = strlen(ref_str) ;
            if (ref_len < 0)
            {
                continue ;
            } else
            {
                for (unsigned int variantCntInt=0;
                	 variantCntInt < variantVec.size() ;
                     variantCntInt++)
                {
					fprintf(stdout, "variant_lines=%i\n", variant_lines) ;
					
                    //find length difference
                    variant_len = variantVec[variantCntInt].size() ;
                    std::string srcIdStr ;
                    lendiff = (variant_len - ref_len) ;
                    
                    //match strains with variants
                    for (unsigned int strainCntInt = 0 ;
                         strainCntInt != strainVec.size() ;
                         strainCntInt++)
                    {
                        //skip blank strains
                        if (strainVec[strainCntInt].size() >= 3)
                        {
							int strainCharInt=0, strainCharInt2=0, strainQualInt=0 ;
							size_t num = sscanf(strainVec[strainCntInt].c_str(), "%i/%i:%i", &strainCharInt, &strainCharInt2, &strainQualInt) ;
							assert(num>=2) ;
							if (strainCharInt!=strainCharInt2)
								fprintf(stdout, "Warning: heterozygous polymorphism\n") ;
							
                        	if (strainCharInt == variantCntInt+1 || strainCharInt2 == variantCntInt+1)
                        	{
                                //fprintf(stdout, "what i'm looking at: %i\n",
                                //       strainRefVec[strainCntInt]) ;

                                assert(variantCntInt<strainRefVec.size()) ;
								
                                //update string for source_id
                                if (srcIdStr.size() == 0)
                        		{
                        			srcIdStr.append(strainRefVec[strainCntInt]) ;
                        		} else
                        		{
                        			srcIdStr.append(",");
                        			srcIdStr.append(strainRefVec[strainCntInt]) ;
                        		}
                        	}
                        } //completed all strains for a variant
                    }
                    strcpy(variant_str, variantVec[variantCntInt].c_str()) ;
                    strcpy(source_id, srcIdStr.c_str()) ;
                    
                    insert_variant(chr_idx, position-1, ref_len, variant_len, ref_str, variant_str, 0, 0, 0,0, "", -2, -1);
                    variant_lines++ ;
                } //completed all variants in a line
            }
        }
    }
    fclose(fd) ;
    /*
    if (variant_lines%10000==0)
    {
        long pos = ftell(fd) ;
        fprintf(stdout, "num_variants=%i\t\t%ld Mb read\r",
                variant_lines, pos/1024/1024) ;
    }
    fprintf(stdout, "read %i variants (checked %i)\n", variant_lines,
            variant_lines_checked) ;
    */
    free(buf) ;
    return 0;
}

int VariantMap::init_from_sdi(const std::string &sdi_fname)
{
    
    fprintf(stdout, "initializing genome variant list with SDI file %s\n", sdi_fname.c_str()) ;
    
    FILE * fd=Util::openFile(sdi_fname.c_str(), "r") ;
    if (!fd)
        return -1 ;
    int variant_lines = 0, variant_lines_checked = 0 ;
    const int max_buf_len = 10000000 ;
    const int max_field_len = 500000 ;
    
    char * buf=(char*)malloc(max_buf_len+1) ;
    strcpy(buf, "") ;
    
    while (!feof(fd))
    {
        char chr_name[1001]="", ref_str[max_field_len+1]="", variant_str[max_field_len+1]="", prop[1001]="", source_id[1001]=""  ;
        int position, lendiff, read_pos=-1, read_len=-1, conf_count=0, non_conf_count=0, used_count=0, non_used_count=0 ;
        
        Util::skip_comment_lines(fd) ;
        
        if (fgets(buf, max_buf_len, fd)==NULL)
            break ;
        //fprintf(stdout, "%s", buf) ;
        
        //Scan sdi line
        //int num = sscanf(buf, "%1000s\t%i\t%i\t%100000s\t%100000s\t%1000s\t%1000s\n", chr_name, &position, &lendiff, ref_str, variant_str, tmp, tmp) ;
        
        int num1 = sscanf(buf,"%1000s\t%i\t%i\t%500000s\t%500000s\t%i\t%i\t%i\t%i\t%1000s\t%i/%i\t%1000s\n",
                          chr_name, &position, &lendiff, ref_str, variant_str, &conf_count, &non_conf_count, &used_count,&non_used_count, source_id, &read_pos, &read_len, prop);
        // compatibility with old format
        int num = sscanf(buf,"%1000s\t%i\t%i\t%500000s\t%500000s\t%i\t%i\t%i\t%1000s\t%i/%i\t%1000s\n",
                         chr_name, &position, &lendiff, ref_str, variant_str, &conf_count, &non_conf_count, &used_count, source_id, &read_pos, &read_len, prop);
        if (num1>num)
        {
            num = sscanf(buf,"%1000s\t%i\t%i\t%500000s\t%500000s\t%i\t%i\t%i\t%i\t%1000s\t%i/%i\t%1000s\n",
                         chr_name, &position, &lendiff, ref_str, variant_str, &conf_count, &non_conf_count, &used_count,&non_used_count, source_id, &read_pos, &read_len, prop);
        }
        
        if (num<5)
        {
            if (feof(fd))
                break ;
            fprintf(stdout, "sdi line only contained %i columns (5 expected), aborting (%s)\nftell=%ld\n%s\n", num, chr_name, ftell(fd), buf) ;
            exit(1) ;
        }
        
        int chr_idx = genome->find_desc(chr_name) ;
        if (chr_idx==-1)
        {
            fprintf(stderr, "chromosome %s not found. known chromosome names:\n", chr_name) ;
            genome->print_desc(stderr) ;
            fclose(fd) ;
            free(buf) ;
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
    
    free(buf) ;
    return 0 ;
    
}

int VariantMap::init_from_info(const std::string &info_fname)
{
    
    fprintf(stdout, "initializing genome variant list with SDI file %s\n", info_fname.c_str()) ;
    
    FILE * fd=Util::openFile(info_fname.c_str(), "r") ;
    if (!fd)
        return -1 ;
    int variant_lines = 0, variant_lines_checked = 0 ;
    const int max_var_len=500000 ;
    
    while (!feof(fd))
    {
        char chr_name[1001]="", ref_str[max_var_len+1]="", variant_str[max_var_len+1]="", buf[2*max_var_len]="", source_id[1001]="", diff_str[2*max_var_len+1]=""  ;
        int position, lendiff ;
        int non_ref_counts=0, read_pair_support=0, bp_range1=0, bp_range2=0, four_gamete_left=0, four_gamete_right=0, len=0;
        char anc_status='N' ;
        
        Util::skip_comment_lines(fd) ;
        
        if (fgets(buf, max_var_len*2, fd)==NULL)
            break ;
        
        int num = sscanf(buf,"%1000s\t%i\t%i\t%500000s\t%i\t%c\t%i\t%i\t%i\t%i\t%i\t%1000s\n",
                         chr_name, &position, &len, diff_str, &non_ref_counts, &anc_status, &read_pair_support, &bp_range1, &bp_range2, &four_gamete_left, &four_gamete_right, source_id);
        if (num<5)
        {
            if (feof(fd))
                break ;
            fprintf(stdout, "info line only contained %i columns (5 expected), aborting (%s)\n", num, chr_name) ;
            continue ;
        }
        if (len>max_var_len)
        {
            fprintf(stdout, "dropping line (variant too long):\n%s\n", buf) ;
            continue ;
        }
        
        /*if (strcmp(diff_str, "TD")==0)
         continue ;*/
        /*if (strcmp(diff_str, "Inv")==0)
         continue ;
         if (strcmp(diff_str, "INV")==0)
         continue ;*/
        /*if (strcmp(diff_str, "Ins")==0)
         continue ;
         if (strcmp(diff_str, "INS")==0)
         continue ;*/
        /*if (strcmp(diff_str, "Del")==0)
         continue ;
         if (strcmp(diff_str, "DEL")==0)
         continue ;*/
        /*if (strcmp(diff_str, "CTX")==0)
         continue ;
         if (strcmp(diff_str, "ITX")==0)
         continue ;*/
        
        int chr_idx = genome->find_desc(chr_name) ;
        if (chr_idx==-1)
            chr_idx = genome->find_desc((std::string("Chr")+std::string(chr_name)).c_str()) ;
        if (chr_idx==-1)
        {
            fprintf(stderr, "chromosome %s not found. known chromosome names:\n", chr_name) ;
            genome->print_desc(stderr) ;
            fclose(fd) ;
            return -1 ;
        }
        
        std::string diff_str2;
        for (int ii=0; ii<(int)strlen(diff_str); ii++)
        {
            assert(diff_str[ii]!='>') ;
            if (diff_str[ii]!='<')
                diff_str2+=diff_str[ii] ;
            else
            {
                int numN=atoi(&diff_str[ii+1]) ;
                for (int kk=0; kk<numN; kk++)
                    diff_str2+='N' ;
                
                while (diff_str[ii]!='>' && ii<(int)strlen(diff_str))
                    ii++ ;
            }
        }
        strcpy(diff_str, diff_str2.c_str()) ;
        if (diff_str[0]=='-')
        {
            strcpy(ref_str, &diff_str[1]) ;
            lendiff=-len ;
            
            for (int k=0; k<len; k++)
            {
                if (toupper(ref_str[k])!='A' && toupper(ref_str[k])!='C' && toupper(ref_str[k])!='G' && toupper(ref_str[k])!='T')
                    ref_str[k]=genome->chromosome(chr_idx)[position+k] ;
                else
                    if (!(ref_str[k]==genome->chromosome(chr_idx)[position+k]))
                    {
                        fprintf(stderr, "error: mismatch in line\n%s\n", buf) ;
                        exit(-1) ;
                    }
            }
            ref_str[len]=0 ;
        }
        else if (diff_str[0]=='+')
        {
            strcpy(variant_str, &diff_str[1]) ;
            lendiff=len ;
        }
        else if (strcmp(diff_str,"Del")==0 || strcmp(diff_str,"DEL")==0)
        {
            for (int k=0; k<len; k++)
                ref_str[k]=genome->chromosome(chr_idx)[position+k] ;
            ref_str[len]=0 ;
            lendiff=-len ;
        }
        else if (strcmp(diff_str,"Ins")==0 || strcmp(diff_str,"INS")==0)
        {
            lendiff=len ;
            for (int k=0; k<len; k++)
                variant_str[k]='N' ;
            variant_str[len]=0 ;
        }
        else if (strcmp(diff_str,"CTX")==0 || strcmp(diff_str,"ITX")==0)
        {
            lendiff=len ;
            for (int k=0; k<len; k++)
                variant_str[k]='N' ;
            variant_str[len]=0 ;
        }
        else if (strcmp(diff_str,"TD")==0)
        {
            lendiff=len ;
            for (int k=0; k<len; k++)
                variant_str[k]=genome->chromosome(chr_idx)[bp_range1+k] ;
            variant_str[len]=0 ;
        }
        else if (strcmp(diff_str,"Inv")==0 || strcmp(diff_str,"INV")==0)
        {
            lendiff=0 ;
            for (int k=0; k<len; k++)
                ref_str[k]=genome->chromosome(chr_idx)[position+k] ;
            ref_str[len]=0 ;
            strcpy(variant_str, QPalma::reverse(QPalma::complement(ref_str)).c_str()) ;
            assert(strlen(variant_str)==strlen(ref_str)) ;
        }
        else
        {
            fprintf(stderr, "Error:\nline=%s\ndiff_str[0]=%s\n", buf, diff_str) ;
            assert(0) ;
        }
        assert(len=(int)strlen(variant_str) || len==(int)strlen(ref_str)) ;
        
        
        int ref_len = strlen(ref_str) ;
        int variant_len = strlen(variant_str) ;
        
        if (!(lendiff==variant_len-ref_len))
        {
            fprintf(stderr, "Error:\nline=%s\ndiff_str[0]=%s\n", buf, diff_str) ;
            assert(0) ;
        }
        
        
        // validate variants on genome sequence
        if (ref_len>0 && validate_variants)
        {
            //fprintf(stdout, "pos=%i\tref_len=%i\tvariant_len=%i\n", pos, ref_len, variant_len) ;
            for (int i=0; i<ref_len; i++)
            {
                if (genome->chromosome(chr_idx)[position+i]!=ref_str[i] &&
                    ref_str[i]!='N' && ref_str[i]!='Y' && ref_str[i]!='W' && ref_str[i]!='K' && ref_str[i]!='S' && ref_str[i]!='M' && ref_str[i]!='R' && ref_str[i]!='D')
                {
                    if (exit_on_validation_error)
                    {
                        fprintf(stderr, "ERROR: variant map disagrees with genome: %i\t%i\t%c\t%c\t%s\n%s\n", i, position+i, genome->chromosome(chr_idx)[position+i], ref_str[i], ref_str, buf) ;
                        exit(-1) ;
                    }
                    else
                        fprintf(stdout, "WARNING: variant map disagrees with genome: %i\t%i\t%c\t%c\t%s\t%s\n", i, position+i, genome->chromosome(chr_idx)[position+i], ref_str[i], ref_str, buf) ;
                }
            }
            variant_lines_checked++ ;
        }
        
        insert_variant(chr_idx, position, ref_len, variant_len, ref_str, variant_str, non_ref_counts+read_pair_support, 0, 0,0, "", -2, -1);
        
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

int VariantMap::init_from_snp_info(const std::string &info_fname)
{
    
    fprintf(stdout, "initializing genome variant list with snp_info file %s\n", info_fname.c_str()) ;
    
    FILE * fd=Util::openFile(info_fname.c_str(), "r") ;
    if (!fd)
        return -1 ;
    int variant_lines = 0, variant_lines_checked = 0 ;
    const int max_buf_len = 10000000 ;
    const int max_field_len = 500000 ;
    
    char *buf=(char*)malloc(max_buf_len+1) ;
    strcpy(buf, "") ;
    char *rest=(char*)malloc(max_buf_len+1) ;
    strcpy(rest, "") ;
    while (!feof(fd))
    {
        char chr_name[1001]="", ref_str[max_field_len+1]="", variant_str[max_field_len+1]=""; ;
        int position=0, non_ref_counts=0 ;
        float q_value ;
        
        Util::skip_comment_lines(fd) ;
        
        if (fgets(buf, max_buf_len, fd)==NULL)
            break ;
        
        int num = sscanf(buf,"%1000s\t%i\t%f\t%500000s\t%i\t%c\t%10000000s\n", chr_name, &position, &q_value, ref_str, &non_ref_counts, variant_str, rest);
        if (num<5)
        {
            if (feof(fd))
                break ;
            fprintf(stdout, "info line only contained %i columns (5 expected), aborting (%s)\n%s\n", num, chr_name, buf) ;
            continue ;
        }
        
        int chr_idx = genome->find_desc(chr_name) ;
        if (chr_idx==-1)
            chr_idx = genome->find_desc((std::string("Chr")+std::string(chr_name)).c_str()) ;
        if (chr_idx==-1)
        {
            fprintf(stderr, "chromosome %s not found. known chromosome names:\n", chr_name) ;
            genome->print_desc(stderr) ;
            fclose(fd) ;
            free(buf) ;
            free(rest) ;
            return -1 ;
        }
        
        int ref_len = strlen(ref_str) ;
        int variant_len = strlen(variant_str) ;
        
        insert_variant(chr_idx, position-1, ref_len, variant_len, ref_str, variant_str, 0, 0, 0, 0, "", -2, -1);
        
        variant_lines++ ;
        
        if (variant_lines%10000==0)
        {
            long pos = ftell(fd) ;
            fprintf(stdout, "num_variants=%i\t\t%ld Mb read\r", variant_lines, pos/1024/1024) ;
        }
    }
    
    fclose(fd) ;
    free(buf) ;
    free(rest) ;
    
    fprintf(stdout, "read %i variants (checked %i)\n", variant_lines, variant_lines_checked) ;
    
    return 0 ;
    
}


int VariantMap::init_from_csv(const std::string &sdi_fname, const std::vector<std::string> & take_lines, VariantInputEnum ext)
{
    std::string ref_line ="" ;
    int ref_line_idx=-1 ;
    
    std::vector<std::string> take_lines_names ;
    std::vector<int> take_lines_idx ;
    
    const int max_buf_len = 10000000 ;
    const int max_field_len = 500000 ;
    
    fprintf(stdout, "initializing genome variant list with %s file %s\n", (ext==snpcsv)?"SNPCSV":"SVCSV", sdi_fname.c_str()) ;
    
    FILE * fd=Util::openFile(sdi_fname.c_str(), "r") ;
    if (!fd)
        return -1 ;
    
    char *buf=(char*)malloc(max_buf_len+1) ;
    strcpy(buf, "") ;
    
    int variant_lines = 0 ;
    
    {
        char list[max_field_len+1];
        if (fgets(buf, max_buf_len, fd)==NULL)
        {
            fprintf(stdout, "Error reading file\n") ;
            exit(-1) ;
        }
        int num = sscanf(buf,"Chromosome,Position,%1000000s", list);
        if (num<1)
            num = sscanf(buf,"Chromosome,Positions,%1000000s", list);
        if (num<1)
            num = sscanf(buf,"chr,loc,%1000000s", list);
        if (num<1)
            num = sscanf(buf,"chr,pos,%1000000s", list);
        if (num<1)
        {
            fprintf(stderr, "%s header only contained %i columns (>=3 expected), aborting \n%s\n", (ext==snpcsv)?"SNPCSV":"SVCSV", num, buf) ;
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
        char chr_name[1001]="", ref_str[max_field_len+1]="", variant_str[max_field_len+1]="" ;
        int position=0 ;
        char list[max_field_len+1] ;
        
        Util::skip_comment_lines(fd) ;
        
        if (fgets(buf, max_buf_len, fd)==NULL)
            break ;
        
        int num = sscanf(buf,"%1000[^,],%i,%500000s", chr_name, &position, list);
        
        if (num<3)
        {
            if (feof(fd))
                break ;
            fprintf(stderr, "%s line only contained %i columns (>=3 expected), aborting (%s)\n%s\n", (ext==snpcsv)?"SNPCSV":"SVCSV", num, chr_name, buf) ;
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
                free(buf) ;
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
            
            int ref_len ;
            int variant_len ;
            std::string source = "" ;
            
            if (ext==snpcsv)
            {
                strcpy(ref_str, "N") ;
                if (ref!='N')
                    ref_str[0]=ref ;
                else
                    ref_str[0]=genome->chromosome(chr_idx)[spos-1] ;
                strcpy(variant_str, "N") ;
                variant_str[0]=SNP[s] ;
                
                if (ref_str[0]==variant_str[0])
                    continue ;
                ref_len = strlen(ref_str) ; assert(ref_len==1) ;
                variant_len = strlen(variant_str) ; assert(variant_len==1) ;
                
                source=(variant_str[0]=='N')?"MISS:":"SNP:"	;
            }
            else
            {
                assert(SNP[s]=='1' || SNP[s]=='0') ;
                if (SNP[s]=='0')
                    continue ;
                
                strcpy(variant_str, "*") ;
                strcpy(ref_str, "*") ;
                variant_len=0 ;
                ref_len=0 ;
                spos++ ;
                source="SV:" ;
                
            }
            source+= take_lines_names[s] ;
            char flank[3]="NN" ;
            
            insert_variant(chr_idx, spos-1, ref_len, variant_len, ref_str, variant_str, 1, 0, 0, 0, source.c_str(), -1, -1, flank, true, true) ;
            
            variant_lines++ ;
            
            if (variant_lines%10000==0)
            {
                long pos = ftell(fd) ;
                fprintf(stdout, "num_variants=%i\t\t%ld Mb read\r", variant_lines, pos/1024/1024) ;
            }
        }
    }
    
    fclose(fd) ;
    
    fprintf(stdout, "read %i %s liness\n", variant_lines, (ext==snpcsv)?"SNPCSV":"SVCSV") ;
    
    free(buf) ;
    
    return 0 ;
    
}


void parse_pileupstr(int coverage, char * str, char variant, int & conf, int & non_conf)
{
    int len=strlen(str) ;
    //assert(coverage<=len) ;
    conf=0 ;
    non_conf=0 ;
    for (int i=0; i<len; i++)
    {
        if (str[i]=='$')
            continue ;
        if (str[i]=='^')
        {
            i++ ;
            continue ;
        }
        if (str[i]=='-' || str[i]=='+')
        {
            std::string ns="" ;
            for (int j=i+1; j<len; j++)
                if (str[j]>='0' && str[j]<='9')
                    ns+=str[j] ;
                else
                    break ;
            int n=atoi(ns.c_str()) ;
            i+=strlen(ns.c_str()) ;
            assert(n>0) ;
            //fprintf(stdout, "n=%i, %i\n", n, (int)strlen(ns.c_str())) ;
            for (int j=0; j<n; j++)
            {
                i++ ;
                assert(i<len) ;
            }
            continue ;
        }
        if (toupper(str[i])==toupper(variant))
            conf++ ;
        else
            non_conf++ ;
    }
    if (conf+non_conf!=coverage)
    {
        fprintf(stdout, "warning: pileup string inconsistency: %s\t%i\t%i\t%i\n", str, conf, non_conf, coverage) ;
    }
}


int VariantMap::init_from_samtools(const std::string &sdi_fname)
{
    fprintf(stdout, "initializing genome variant list with indel/samtools file %s\n", sdi_fname.c_str()) ;
    
    FILE * fd=Util::openFile(sdi_fname.c_str(), "r") ;
    if (!fd)
        return -1 ;
    int variant_lines = 0, variant_lines_checked = 0 ;
    int line_num = 0 ;
    const int max_buf_len = 10000000 ;
    const int max_field_len = 500000 ;
    
    char *buf=(char*)malloc(max_buf_len+1) ;
    strcpy(buf, "") ;
    char *rest=(char*)malloc(max_buf_len+1) ;
    strcpy(rest, "") ;
    
    while (!feof(fd))
    {
        char chr_name[1001]="", ref_str[max_field_len+1]="", variant_str[max_field_len+1]="" ;
        int position=0 ;
        char str[2][max_field_len+1]={"", ""} ;
        char allele_str[2][max_field_len+1]={"", ""} ;
        int allele_reads[2]={0,0} ;
        int xxx ;
        
        
        line_num+=Util::skip_comment_lines(fd) ;
        
        if (fgets(buf, max_buf_len, fd)==NULL)
            break ;
        line_num++ ;
        //fprintf(stdout, "%s", buf) ;
        
        int num = sscanf(buf,"%1000s\t%i\t%500000s\t%500000s\t%10000000s", chr_name, &position, ref_str, str[0], rest);
        bool is_snp=false ;
        
        if (num>=3)
        {
            if (strcmp(ref_str, "*")==0)
            {
                // indel
                num = sscanf(buf,"%1000s\t%i\t%500000s\t%500000s\t%i\t%i\t%i\t%i\t%500000s\t%500000s\t%i\t%i\t%10000000s", chr_name, &position, ref_str, str[0], &xxx, &xxx, &xxx, &xxx, allele_str[0], allele_str[1], &allele_reads[0], &allele_reads[1], rest);
                //fprintf(stdout, "NUM=%i\n", num) ;
                
                is_snp=false ;
                if (num<10)
                    strcpy(allele_str[0], "") ;
                if (num<11)
                    strcpy(allele_str[1], "") ;
                if (num<12)
                    allele_reads[0]=0 ;
                if (num<13)
                    allele_reads[1]=0 ;
            }
            else
            {
                // snp
                char tmp[max_field_len+1] ;
                allele_reads[1]=0 ;
                num = sscanf(buf,"%1000s\t%i\t%500000s\t%500000s\t%i\t%i\t%i\t%i\t%500000s\t%500000s", chr_name, &position, ref_str, variant_str, &xxx, &xxx, &xxx, &allele_reads[0], tmp, rest);
                if (num<8)
                {
                    fprintf(stdout, "samtools/indel line %i only contained %i columns (>=8 expected), aborting (%s)\n%s\n", line_num, num, chr_name, buf) ;
                    continue ;
                }
                
                assert(strlen(variant_str)==1) ;
                parse_pileupstr(allele_reads[0], tmp, variant_str[0], allele_reads[0], allele_reads[1]) ;
                
                is_snp=true ;
            }
        }
        
        if (num<3)
        {
            if (feof(fd))
                break ;
            fprintf(stdout, "samtools/indel line %i only contained %i columns (>=4 expected), aborting (%s)\n%s\n", line_num, num, chr_name, buf) ;
            continue ;
        }
        
        
        int chr_idx = genome->find_desc(chr_name) ;
        if (chr_idx==-1)
        {
            fprintf(stderr, "chromosome %s not found. known chromosome names:\n", chr_name) ;
            genome->print_desc(stderr) ;
            fclose(fd) ;
            
            free(rest) ;
            free(buf) ;
            return -1 ;
        }
        /*strcpy(str[1], str[0]) ;
         for (unsigned int i=0; i<strlen(str[0]); i++)
         {
         if (str[0][i]=='/')
         {
         strcpy(str[1], &str[0][i+1]) ;
         str[0][i]=0 ;
         break ;
         }
         }*/
        
        for (int s=0; s<2; s++)
        {
            //fprintf(stdout, "num=%i\tis_snp=%i\ts=%i\tstr=%s\n", num, is_snp, s, allele_str[s]) ;
            
            int spos=position ;
            if (s==1 && is_snp)
                continue ;
            
            if (!is_snp && strcmp(allele_str[s], "*")==0)
                continue ;
            
            char flank[3]="NN" ;
            if (allele_str[s][0]=='-' || allele_str[s][0]=='+')
            {
                assert(!is_snp) ;
                
                if (allele_str[s][0]=='+')
                {
                    strcpy(ref_str, "") ;
                    assert(strcmp(allele_str[s], "*")!=0) ;
                    strcpy(variant_str, allele_str[s]+1) ;
                    //flank[0]=genome->chromosome(chr_idx)[spos-1] ;
                    //flank[1]=genome->chromosome(chr_idx)[spos-1+1] ;
                }
                else
                {
                    strcpy(variant_str, "") ;
                    assert(strcmp(allele_str[s], "*")!=0) ;
                    strcpy(ref_str, allele_str[s]+1) ;
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
                std::string info_str = "indel" ;
#ifndef PMINDEX
                if (_config.VARIANT_SNP_TAKE_LINES.size()>0 && _config.VARIANT_SNP_TAKE_LINES[0]!="*")
                    info_str=_config.VARIANT_SNP_TAKE_LINES[0] ;
#endif
                int conf=allele_reads[s] ;
                int non_conf=allele_reads[1-s] ;
                insert_variant(chr_idx, spos-1, ref_len, variant_len, ref_str, variant_str, conf, non_conf, 0, 0, info_str.c_str(), -1, -1, flank) ;
            }
            else
            {
                assert(is_snp) ;
                //assert(strcmp(allele_str[0], allele_str[1])==0) ;
                strcpy(variant_str, str[0]) ;
                
                int ref_len = strlen(ref_str) ;
                int variant_len = strlen(variant_str) ;
                
                std::string info_str = "snp" ;
#ifndef PMINDEX
                if (_config.VARIANT_SNP_TAKE_LINES.size()>0 && _config.VARIANT_SNP_TAKE_LINES[0]!="*")
                    info_str=_config.VARIANT_SNP_TAKE_LINES[0] ;
#endif
                int conf=allele_reads[s] ;
                int non_conf=allele_reads[1-s] ;
                insert_variant(chr_idx, spos-1, ref_len, variant_len, ref_str, variant_str, conf, non_conf, 0, 0, info_str.c_str(), -1, -1, flank) ;
            }
            
            variant_lines++ ;
            
            if (variant_lines%1000==0)
            {
                long pos = ftell(fd) ;
                fprintf(stdout, "num_variants=%i\t\t%ld Mb read\r", variant_lines, pos/1024/1024) ;
            }
        }
    }
    
    fclose(fd) ;
    
    fprintf(stdout, "read %i indels (checked %i)\n", variant_lines, variant_lines_checked) ;
    
    free(buf) ;
    free(rest) ;
    
    return 0 ;
    
}


int VariantMap::report_to_file(const std::string &filename) const
{
    size_t pos_ext=filename.rfind('.');
    
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
    
    fprintf(fd,"#chromosome\tposition\tlen_diff\tref_seq\tvariant_seq\tconf_count\tnon_conf_count\tused_count\tnon_used_count\tsource\tread position/len\tstatus\n");
    
    for (unsigned int i=0; i<genome->nrChromosomes(); i++)
    {
        const char * chr= genome->get_desc(i);
        std::vector<Variant>::iterator it;
        
        for (it=variantlist[i].begin(); it!=variantlist[i].end(); it++)
        {
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
    
    return 0;
    
}

int VariantMap::stats_to_file(const std::string &stats, int max_len)  const
{
    int nb_variants=0;
    
    fprintf(stdout, "report genome variants stats to file %s\n", stats.c_str()) ;
    
    FILE * fd=stdout ;
    if (!stats.empty())
        fd=Util::openFile(stats.c_str(), "w") ;
    if (!fd)
        return -1 ;
    
    int num_del=0 ;
    int num_ins=0 ;
    int num_snp=0 ;
    int num_sub=0 ;
    std::vector<int> num_del_dist(max_len+2, 0) ;
    std::vector<int> num_ins_dist(max_len+2, 0) ;
    
    for (unsigned int i=0; i<genome->nrChromosomes(); i++)
    {
        //const char * chr= genome->get_desc(i);
        std::vector<Variant>::iterator it;
        
        for (it=variantlist[i].begin(); it!=variantlist[i].end(); it++)
        {
            if ((*it).ref_len<0 || (*it).variant_len<0)
            {
                fprintf(stderr, "dropping variant with invalid sizes %i/%i\n", (*it).ref_len, (*it).variant_len) ;
                continue ;
            }
            
            if ((*it).type==pt_SNP)
                num_snp++ ;
            if ((*it).type==pt_deletion)
            {
                num_del++ ;
                assert((*it).ref_len>=0) ;
                if ((*it).ref_len>max_len)
                    num_del_dist[max_len+1]++ ;
                else
                    num_del_dist[(*it).ref_len]++ ;
            }
            if ((*it).type==pt_insertion)
            {
                num_ins++ ;
                assert((*it).variant_len>=0) ;
                if ((*it).variant_len>max_len)
                    num_ins_dist[max_len+1]++ ;
                else
                    num_ins_dist[(*it).variant_len]++ ;
            }
            if ((*it).type==pt_substitution)
                num_sub++ ;
            
            nb_variants++;
        }
    }
    
    fprintf(fd, "reported %i variants\n", nb_variants) ;
    fprintf(fd, "* %i SNPs (%1.2f%%)\n", num_snp, (100.0*num_snp)/nb_variants) ;
    fprintf(fd, "* %i substitutions (%1.2f%%)\n", num_sub, (100.0*num_sub)/nb_variants) ;
    fprintf(fd, "* %i deletions (%1.2f%%)\n", num_del, (100.0*num_del)/nb_variants) ;
    for (int i=1; i<max_len+1; i++)
        if (num_del_dist[i]>0)
            fprintf(fd, "D\t%i\t%i\t%1.2f%%\n", i, num_del_dist[i], (100.0*num_del_dist[i])/num_del) ;
    fprintf(fd, "D\t>%i\t%i\t%1.2f%%\n", max_len, num_del_dist[max_len+1], (100.0*num_del_dist[max_len+1])/num_del) ;
    fprintf(fd, "* %i insertions (%1.2f%%)\n", num_ins, (100.0*num_ins)/nb_variants) ;
    for (int i=1; i<max_len+1; i++)
        if (num_ins_dist[i]>0)
            fprintf(fd, "I\t%i\t%i\t%1.2f%%\n", i, num_ins_dist[i], (100.0*num_ins_dist[i])/num_ins) ;
    fprintf(fd, "I\t>%i\t%i\t%1.2f%%\n", max_len, num_ins_dist[max_len+1], (100.0*num_ins_dist[max_len+1])/num_ins) ;
    
    if (fd!=stdout)
        fclose(fd) ;
    
    return 0;
    
}

int VariantMap::init_from_files(std::string &fnames)
{
    
    int previousfound=0;
    int found=fnames.find(",");
    std::string filename;
#ifndef PMINDEX
    bool has_maf_file = false ;
#endif
    
    while (true)
    {
        if (found >=0)
            filename = fnames.substr(previousfound, found-previousfound);
        else
            filename=fnames.substr(previousfound);
        
        size_t pos_ext=filename.rfind('.');
        //fprintf(stdout, "fn=%s,%i\n", filename.c_str(), pos_ext) ;
        VariantInputEnum ext=unknown;
        if (pos_ext != std::string::npos)
        {
            std::string extension(filename.substr(pos_ext+1));
            //fprintf(stdout, "ext=%s\n", extension.c_str()) ;
            
            if (extension.compare("sdi")==0 ||extension.compare("SDI")==0)
                ext=sdi;
            if (extension.compare("vcf")==0 || extension.compare("VCF")==0)
                ext=vcf;
#ifndef PMINDEX
            if (extension.compare("maf")==0 ||extension.compare("MAF")==0)
                ext=maf;
#endif
            if (extension.compare("indel")==0 || extension.compare("samtools")==0 || extension.compare("var")==0)
                ext=samtools;
            if (extension.compare("snp")==0 ||extension.compare("SNP")==0)
                ext=snpcsv;
            if (extension.compare("snpcsv")==0 ||extension.compare("SNPCSV")==0)
                ext=snpcsv;
            if (extension.compare("svcsv")==0 ||extension.compare("SVCSV")==0)
                ext=svcsv;
            if (extension.compare("bin")==0 ||extension.compare("BIN")==0)
                ext=bingz;
            if (extension.compare("snpinfo")==0 ||extension.compare("SNPINFO")==0)
                ext=snpinfo;
            else
                if (extension.compare("info")==0 ||extension.compare("INFO")==0)
                    ext=info;
        }
        if (ext==unknown)
        {
            fprintf(stderr,	"ERROR: Variant input %s has an unknown format\n", (char*)filename.c_str());
            exit(1) ;
        }
#ifndef PMINDEX
        if ( ext==maf && _config.MAF_REF_NAME.length()<=0 )
        {
            fprintf(stderr,	"ERROR: Need reference genome name (with -maf-ref option) to initialize variants from the maf file %s\n", (char*)filename.c_str());
            exit(1) ;
        }
#endif
        
        int ret = 0 ;
        if (ext == sdi)
            init_from_sdi(filename);
        if (ext == snpinfo)
            init_from_snp_info(filename);
        if (ext == info)
            init_from_info(filename);
        if (ext == bingz)
            init_from_bin(filename);
        if (ext == samtools)
            init_from_samtools(filename);
        if (ext == vcf)
            init_from_vcf(filename);
#ifndef PMINDEX
        if (ext == snpcsv)
            init_from_csv(filename, _config.VARIANT_SNP_TAKE_LINES, ext);
        if (ext == svcsv)
            init_from_csv(filename, _config.VARIANT_SNP_TAKE_LINES, ext);
        if (ext == maf)
        {
            init_from_maf(filename, _config.MAF_REF_NAME);
            has_maf_file = true ;
        }
#endif
        
        if (ret!=0)
            return ret;
        check_variant_order() ;
        
        if (found<=0)
            break ;
        
        previousfound=found+1;
        found=fnames.find(",", found+1);
    }
    
    known_variants_limit=next_variant_id-1;
    
#ifndef PMINDEX
    if (!has_maf_file && _config.MAF_REF_NAME.length()>0)
        fprintf(stdout, "WARNING: maf reference given, but no maf file as input\n") ;
#endif
    
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
        size_t position_tmp = name_tmp.find_last_of('.');
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
            std::vector<std::vector<Variant>::iterator> variant_list;
            std::vector<Variant>::iterator it = my_lower_bound(this->variantlist[chr_idx].begin(), this->variantlist[chr_idx].end(), e_idx-100) ;
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
