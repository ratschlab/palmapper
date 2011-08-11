#pragma once

#include <palmapper/Genome.h>
#include <deque>
#include <string>
#include <stdlib.h> 
#include <assert.h>

enum polytype
{
	pt_unknown,
	pt_SNP,
	pt_deletion,
	pt_insertion,
	pt_substitution
}  ;

struct variant_str {
	int position ;
	int end_position ;
	int ref_len ;
	int variant_len ;
	int id ;
	enum polytype type ;
	std::string ref_str, variant_str ;
	std::string read_id ;
	int read_pos;
	int used_count ;
	int conf_count ;
	int non_conf_count ;
};

typedef struct variant_str Variant;

struct found_variant_str {
	int read_position ;
	int id ;
	enum polytype type ;
};
typedef struct found_variant_str FoundVariant;

struct super_variant_str {
	int variant_id ;
	int position ;
	int end_position ;
	char SNP[2] ;
	enum polytype type ;
};
typedef struct super_variant_str SuperVariant ;


class VariantMap
{

private:
	int next_variant_id;
	int limit_known_variants;

public:
	std::deque<Variant> * variantlist ;
	VariantMap(Genome const &genome_) ;
	~VariantMap() ;

	void insert_variant(int chr, int pos, int ref_len, int variant_len, const std::string & ref_str, const std::string & variant_str, int conf_count, int used_count, const std::string & read_id,int read_pos);
	void insert_variant(Variant & j, int chr) ;
	int init_from_sdis(std::string &sdi_fname);
	int report_to_sdi(std::string &sdi_fname);



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

	int variant_cmp(const Variant &a, const Variant &b)
	{
		if (a.position<b.position)
			return -1  ;
		if (a.position>b.position)
			return 1  ;
		if (a.end_position<b.end_position)
			return -1 ;
		if (a.end_position>b.end_position)
			return 1 ;
		if (a.ref_str<b.ref_str)
			return -1 ;
		if (a.ref_str>b.ref_str)
			return 1 ;
		if (a.variant_str<b.variant_str)
			return -1 ;
		if (a.variant_str>b.variant_str)
			return 1 ;
		return 0 ;
	}

	
	static void report_SNP_variant(std::vector<Variant> & variants, const Chromosome & chr, int dna_pos, char ref, char variant, const std::string & read_id, int read_pos) 
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
		v.read_id=read_id ;
		v.non_conf_count = 0 ;
		v.read_pos=read_pos;
		variants.push_back(v) ;
	}
	
	static void report_sub_variant(std::vector<Variant> & variants, const Chromosome & chr, int dna_pos, int len_ref, std::string & ref_str, int len_var, std::string & variant_str, const std::string & read_id, int read_pos)
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
		v.read_id=read_id ;
		v.non_conf_count = 0 ;
		v.read_pos=read_pos;
		variants.push_back(v) ;
	}

	static void report_del_variant(std::vector<Variant> & variants, const Chromosome & chr, int dna_pos, int len, std::string & ref_str, const std::string & read_id, int read_pos)
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
		v.read_id=read_id ;
		v.non_conf_count = 0 ;
		v.read_pos=read_pos;
		variants.push_back(v) ;
	}

	static void report_ins_variant(std::vector<Variant> & variants, const Chromosome & chr, int dna_pos, int len, std::string & variant_str, const std::string & read_id, int read_pos)
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
		v.read_id=read_id ;
		v.non_conf_count = 0 ;
		v.read_pos=read_pos;
		variants.push_back(v) ;
	}
	
	void report_non_variant(const Chromosome * chr, std::vector<int> & aligned_positions, std::vector<int> & exons) ;

	void check_variant_order()
	{
		for (unsigned int i=0; i<genome->nrChromosomes(); i++)
			for (unsigned j=0; j<variantlist[i].size()-1; j++)
				assert(variant_cmp(variantlist[i][j], variantlist[i][j+1])<0) ;
	}

protected:
	
	int init_from_sdi(std::string &gff_fname);

	Genome const *genome;

	pthread_mutex_t variant_mutex;
};

inline std::deque<Variant>::iterator  my_lower_bound ( std::deque<Variant>::iterator first, std::deque<Variant>::iterator  last, const int& value )
{
	std::deque<Variant>::iterator it;
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
