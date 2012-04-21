#pragma once

#include <assert.h>
#include <palmapper/Genome.h>
#include <palmapper/GenomeMaps.h>
#include <deque>
#include <map>
#include <string>
#include <stdlib.h> 
#include <iostream>
#include <gzstream/gzstream.h>
#include "JunctionMap.h"

enum polytype
{
	pt_unknown,
	pt_SNP,
	pt_deletion,
	pt_insertion,
	pt_substitution
}  ;

struct variant_str 
{
	int id ;
	int position ;
	int end_position ;
	int ref_len ;
	int variant_len ;
	enum polytype type ;
	short int read_pos;
	short int read_len;
	short unsigned int used_count ;
	short unsigned int non_used_count ;
	short unsigned int conf_count ;
	short unsigned int non_conf_count ;
	std::string ref_str, variant_str ;
	std::string read_id ;
#ifdef PMINDEX
	std::vector<int> read_id_num ;
#endif
} ;

typedef struct variant_str Variant;
	
inline igzstream& operator>>(igzstream & os, struct variant_str & a)
{
	os.read((char*)&a.id, sizeof(a.id)) ;
	os.read((char*)&a.position, sizeof(a.position)) ;
	os.read((char*)&a.end_position, sizeof(a.end_position)) ;
	os.read((char*)&a.ref_len, sizeof(a.ref_len)) ;
	os.read((char*)&a.variant_len, sizeof(a.variant_len)) ;
	os.read((char*)&a.type, sizeof(a.type)) ;
	os.read((char*)&a.read_pos, sizeof(a.read_pos)) ;
	os.read((char*)&a.read_len, sizeof(a.read_len)) ;
	os.read((char*)&a.used_count, sizeof(a.used_count)) ;
	os.read((char*)&a.non_used_count, sizeof(a.non_used_count)) ;
	os.read((char*)&a.conf_count, sizeof(a.conf_count)) ;
	os.read((char*)&a.non_conf_count, sizeof(a.non_conf_count)) ;

	int size=0 ;
	os.read((char*)&size, sizeof(size)) ;
	{
		char buf[size+1] ;
		os.read(buf, size) ;
		buf[size]=0 ;
		a.ref_str.assign(buf) ;
	}
	
	size=0 ;
	os.read((char*)&size, sizeof(size)) ;
	{
		char buf[size+1] ;
		os.read(buf, size) ;
		buf[size]=0 ;
		a.variant_str.assign(buf) ;
	}

	size=0 ;
	os.read((char*)&size, sizeof(size)) ;
	{
		char buf[size+1] ;
		os.read(buf, size) ;
		buf[size]=0 ;
		a.read_id.assign(buf, size) ;
	}

	return os ;
}

inline ogzstream& operator<<(ogzstream & os, const struct variant_str & a)
{
	os.write((char*)&a.id, sizeof(a.id)) ;
	os.write((char*)&a.position, sizeof(a.position)) ;
	os.write((char*)&a.end_position, sizeof(a.end_position)) ;
	os.write((char*)&a.ref_len, sizeof(a.ref_len)) ;
	os.write((char*)&a.variant_len, sizeof(a.variant_len)) ;
	os.write((char*)&a.type, sizeof(a.type)) ;
	os.write((char*)&a.read_pos, sizeof(a.read_pos)) ;
	os.write((char*)&a.read_len, sizeof(a.read_len)) ;
	os.write((char*)&a.used_count, sizeof(a.used_count)) ;
	os.write((char*)&a.non_used_count, sizeof(a.non_used_count)) ;
	os.write((char*)&a.conf_count, sizeof(a.conf_count)) ;
	os.write((char*)&a.non_conf_count, sizeof(a.non_conf_count)) ;

	int size=a.ref_str.size() ;
	os.write((char*)&size, sizeof(size)) ;
	os.write(a.ref_str.c_str(), size) ;

	size=a.variant_str.size() ;
	os.write((char*)&size, sizeof(size)) ;
	os.write(a.variant_str.c_str(), size) ;

	size=a.read_id.size() ;
	os.write((char*)&size, sizeof(size)) ;
	os.write(a.read_id.c_str(), size) ;
	
	return os ;
}

inline igzstream& operator>>(igzstream & os, std::vector<Variant> & list)
{
	unsigned int size=0 ;
	//os.read((char*)&next_variant_id, sizeof(next_variant_id)) ;
	os.read((char*)&size, sizeof(size)) ;

	list.clear() ;
	for (unsigned int i=0; i<size; i++)
	{
		Variant v ;
		os >> v ;
		list.push_back(v) ;
	}
	return os ;
}

inline ogzstream& operator<<(ogzstream & os, const std::vector<Variant> & list)
{
	unsigned int size=list.size() ;
	//os.write((char*)&next_variant_id, sizeof(next_variant_id)) ;
	os.write((char*)&size, sizeof(size)) ;
	//fprintf(stdout, "size=%i\n", size) ;
	
	for (unsigned int i=0; i<size; i++)
		os << list[i] ;
	return os ;
}

struct found_variant_str {
	int start_pos ;
	int end_pos ;
	int id ;
};
typedef struct found_variant_str FoundVariant;

struct variant_cache_str
{
	std::vector<int>  end_positions ;
	std::vector<int>  id_dels ;
	std::vector<int>  id_snps ;	
	std::vector<char> snps ;
	int insertion;
} ;
typedef struct variant_cache_str variant_cache_t ;


enum VariantInputEnum {
	sdi,
	maf,
	samtools,
	snpcsv,
	svcsv,
	bingz,
	info,
	snpinfo,
	unknown,
};

class VariantMap
{

private:
	int next_variant_id;
	int known_variants_limit;
	bool validate_variants ;
	bool exit_on_validation_error ;
	bool insert_unsorted ;
	int max_variant_len ;
	int merge_variant_source_ids ;
	
public:
	std::vector<Variant> * variantlist ;

	VariantMap(Genome const &genome_, bool merge_variant_source_ids=false) ;
	~VariantMap() ;
	int insert_variants_from_multiple_alignments(std::string & ref_align,int ref_len, std::vector<std::string> & variant_align, std::vector<std::string> & variant_name, 
												 int start_position, int ref_chr_len, int chr_idx, char strand);
	void insert_variant(int chr, int pos, int ref_len, int variant_len, const std::string & ref_str, const std::string & variant_str, int conf_count, int non_conf_count, int used_count,int non_used_count, 
						const std::string & read_id, int read_pos, int read_len, const char* flank="NN", 
						bool update_only=false, bool ignore_variant_str_in_cmp=false);
	
	bool validate_variant(const Variant & j, int chr, const char *flank="NN") const ;
	void insert_variant(Variant & j, int chr, const char* flank="NN", bool update_only=false, bool ignore_variant_str_in_cmp=false) ;
	int init_from_files(std::string &sdi_fname);

	int stats_to_file(const std::string &stats_fname, int max_len) const ;
	int report_to_file(const std::string &sdi_fname) const ;
	int report_to_sdi(const std::string &sdi_fname) const ;
	int report_to_bin(const std::string &sdi_fname) const 
	{
		fprintf(stdout, "report genome variants in BIN file %s ...", sdi_fname.c_str()) ;

		ogzstream s(sdi_fname.c_str()) ;
		for (int i=0; i<(int)genome->nrChromosomes(); i++)
		{
			fprintf(stdout, ".") ;
			s << variantlist[i] ;
		}
		//fprintf(stdout, "next_variant_id=%i\n", next_variant_id) ;
		std::string end_marker="sdi_bin_end_marker" ;
		s << end_marker ;
		
		int num_variants =0;
		for (int i=0; i<(int)genome->nrChromosomes(); i++)
			for (unsigned int j=0; j<variantlist[i].size(); j++, num_variants++)  ;

		fprintf(stdout, " Done.\nWrote %i variants\n", num_variants) ;

		return 0 ;
	}
	int init_from_bin(const std::string &sdi_fname)
	{
		fprintf(stdout, "init genome variants in BIN file %s ...", sdi_fname.c_str()) ;

		igzstream s(sdi_fname.c_str()) ;
		next_variant_id=0 ;
		
		int num_variants =0;
		for (int i=0; i<(int)genome->nrChromosomes(); i++)
		{
			s >> variantlist[i] ;
			fprintf(stdout, ".") ;
			for (unsigned int j=0; j<variantlist[i].size(); j++, num_variants++) 
				if (next_variant_id<variantlist[i][j].id+1)
					next_variant_id = variantlist[i][j].id+1 ;
		}
		std::string end_marker="" ;
		s >> end_marker ;
		if (end_marker!=std::string("sdi_bin_end_marker"))
			fprintf(stdout, "\nWarning: EOF marker missing\n") ;

		fprintf(stdout, " Done.\nRead %i variants\n", num_variants) ;
		
		return 0 ;
	}

#ifdef PMINDEX
	int get_read_id_num()
	{
		lock() ;
		int counter=0 ;
		std::map<std::string,int> map ;
		for (int i=0; i<(int)genome->nrChromosomes(); i++)
		{
			for (unsigned int j=0; j<variantlist[i].size(); j++)
			{
				std::string source_ids=variantlist[i][j].read_id ;
				int found=source_ids.find(",");
				int previousfound=0;

				while (true)
				{
					std::string source_id ;
					if (found >=0)
						source_id = source_ids.substr(previousfound, found-previousfound);
					else
						source_id = source_ids.substr(previousfound);
					
					if (map.count(source_id)==0)
					{
						/*if (counter<30)
						  fprintf(stdout, "%s  %i\n", source_id.c_str(), counter) ;*/
						map[source_id] = counter ;
						variantlist[i][j].read_id_num.push_back(counter) ;
						counter++ ;
					}
					else
						variantlist[i][j].read_id_num.push_back(map[source_id]) ;

					if (found<=0)
						break ;
					
					previousfound=found+1;
					found=source_ids.find(",", found+1);
					
				}
				
			}
		}
		unlock() ;
		
		return counter ;
	}
#endif

	void filter_variants(int min_source_count, int min_conf_count, double max_nonconf_ratio, 
						 int min_use_count, std::vector<std::string> & accept_sources, std::vector<std::string> & required_sources, 
						 int max_len, int filter_by_map, const GenomeMaps & genomemaps) ;
	void filter_variants_junctions(JunctionMap & junctions) ;
	void unique_variant_source_ids() ;
	
	void lock() 
	{ 
		pthread_mutex_lock( &variant_mutex) ; 
	}

	void unlock() 
	{ 
		pthread_mutex_unlock( &variant_mutex) ; 
	}
	
	bool variant_identical(const Variant &a, const Variant &b)
	{
		if (a.position!=b.position || a.end_position!=b.end_position || a.ref_len!=b.ref_len || a.variant_len!=b.variant_len)
			return false ;
		if (a.ref_str!=b.ref_str || a.variant_str!=b.variant_str)
			return false ;
		return true ;
	}

	static int variant_cmp(const Variant &a, const Variant &b, bool ignore_variant_str=false) 
	{
		if (a.position<b.position)
			return -1  ;
		if (a.position>b.position)
			return 1  ;

		/* special case of updating an existing indel */
		if (a.ref_str=="*" && a.variant_str=="*")
			return 0 ;
		if (b.ref_str=="*" && b.variant_str=="*")
			return 0 ;

		if (a.end_position<b.end_position)
			return -1 ;
		if (a.end_position>b.end_position)
			return 1 ;
		
		if ((a.ref_str.size()!=1 || b.ref_str.size()!=1) || 
			((toupper(a.ref_str[0])=='A' || toupper(a.ref_str[0])=='C' || toupper(a.ref_str[0])=='G' || toupper(a.ref_str[0])=='T') &&
			 (toupper(b.ref_str[0])=='A' || toupper(b.ref_str[0])=='C' || toupper(b.ref_str[0])=='G' || toupper(b.ref_str[0])=='T')))
		{
			if (a.ref_str<b.ref_str)
				return -1 ;
			if (a.ref_str>b.ref_str)
				return 1 ;
		}
		if (!ignore_variant_str && a.variant_str<b.variant_str)
			return -1 ;
		if (!ignore_variant_str && a.variant_str>b.variant_str)
			return 1 ;
		return 0 ;
	}

	
	void report_SNP_variant(std::vector<Variant> & variants, const Chromosome & chr, int dna_pos, char ref, char variant, const std::string & read_id, int read_pos, int read_len) const
	{
		Variant v ;
		v.type = pt_SNP ;
		v.position = dna_pos ;
		v.end_position = dna_pos+1 ;
		v.ref_len=1 ;
		v.variant_len=1;
		v.ref_str+=ref ;
		v.variant_str+=variant ;
		v.conf_count = 1 ;
		v.used_count = 0 ;
		v.non_used_count = 0 ;
		v.read_id=read_id ;
		v.non_conf_count = 0 ;
		v.read_pos=read_pos;
		v.read_len=read_len;
		
		if (validate_variants)
			if (!validate_variant(v, chr.nr()))
				return ;

		variants.push_back(v) ;
	}
	
	void report_sub_variant(std::vector<Variant> & variants, const Chromosome & chr, int dna_pos, int len_ref, std::string & ref_str, int len_var, std::string & variant_str, const std::string & read_id, int read_pos, int read_len) const
	{
		Variant v ;
		v.type = pt_deletion ;
		v.position = dna_pos ;
		v.end_position = dna_pos+ref_str.size() ;
		v.ref_len=ref_str.size() ;
		v.variant_len=variant_str.size();
		v.ref_str=ref_str ;
		v.variant_str=variant_str ;
		v.conf_count = 1 ;
		v.used_count = 0 ;
		v.non_used_count = 0 ;
		v.read_id=read_id ;
		v.non_conf_count = 0 ;
		v.read_pos=read_pos;
		v.read_len=read_len;

		if (validate_variants)
			if (!validate_variant(v, chr.nr()))
				return ;

		variants.push_back(v) ;
	}

	void report_del_variant(std::vector<Variant> & variants, const Chromosome & chr, int dna_pos, int len, std::string & ref_str, const std::string & read_id, int read_pos, int read_len) const
	{
		Variant v ;
		v.type = pt_deletion ;
		v.position = dna_pos ;
		v.end_position = dna_pos+ref_str.size() ;
		v.ref_len=ref_str.size() ;
		v.variant_len=0;
		v.ref_str=ref_str ;
		v.variant_str="" ;
		v.conf_count = 1 ;
		v.used_count = 0 ;
		v.non_used_count = 0 ;
		v.read_id=read_id ;
		v.non_conf_count = 0 ;
		v.read_pos=read_pos;
		v.read_len=read_len;

		if (validate_variants)
			if (!validate_variant(v, chr.nr()))
				return ;

		variants.push_back(v) ;
	}

	void report_ins_variant(std::vector<Variant> & variants, const Chromosome & chr, int dna_pos, int len, std::string & variant_str, const std::string & read_id, 
							int read_pos, int read_len, const char* flank) const
	{
		Variant v ;
		v.type = pt_insertion ;
		v.position = dna_pos ;
		v.end_position = dna_pos ;
		v.ref_len=0 ;
		v.variant_len=variant_str.size();
		v.ref_str="" ;
		v.variant_str=variant_str ;
		v.conf_count = 1 ;
		v.used_count = 0 ;
		v.non_used_count = 0 ;
		v.read_id=read_id ;
		v.non_conf_count = 0 ;
		v.read_pos=read_pos;
		v.read_len=read_len;

		if (validate_variants)
			if (!validate_variant(v, chr.nr(), flank))
				return ;

		variants.push_back(v) ;
	}
	
	void report_non_variant(int rank, int total, const Chromosome * chr, std::vector<int> & aligned_positions, std::vector<int> & exons, int no_gap_end) ;
	void report_variant(int rank, int total, Variant & j, int chr, const char* flank="NN", bool update_only=false, bool ignore_variant_str_in_cmp=false) ;
	int update_variant(int rank, int total, int index, int chr, const Variant &v,const char *flank="NN");
	
	void check_variant_order() const
	{
		for (int i=0; i<(int)genome->nrChromosomes(); i++)
			for (int j=0; j<(int)(variantlist[i].size()-1); j++)
			{
				if (variant_cmp(variantlist[i][j], variantlist[i][j+1])>0) 
				{
					fprintf(stdout, "ERROR: wrong order %i/%i-%i\n", i, j, j+1) ;
				}
			}		
	}

	void transcribe_gff(const std::string & gff_input, const std::string & fasta_output) ;
	

protected:
	
	int init_from_samtools(const  std::string &sam_fname);
	int init_from_csv(const  std::string &snp_fname, const std::vector<std::string> & take_lines, VariantInputEnum ext);
	int init_from_sdi(const  std::string &sdi_fname);
	int init_from_info(const  std::string &info_fname);
	int init_from_snp_info(const  std::string &info_fname);
	int init_from_maf(const  std::string &gff_fname, const std::string &ref_genome);
	pthread_mutex_t variant_mutex;

public:
	Genome const *genome;

};

inline std::vector<Variant>::iterator  my_lower_bound ( std::vector<Variant>::iterator first, std::vector<Variant>::iterator  last, const int& value )
{
	std::vector<Variant>::iterator it;
	long int count, step;
	count = distance(first,last);

	while (count>0)
	{
		it = first;
		step = count/2;
		advance (it,step);
		
		if ( (*it).position < value) 
		{
			first=it; //++it;
			count-=step+1;
		}
		else count=step;
	}
	return first;
}
