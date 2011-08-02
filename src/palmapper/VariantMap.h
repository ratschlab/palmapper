#pragma once

#include <palmapper/Genome.h>
#include <deque>
#include <string>
#include <stdlib.h> 

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

inline void report_SNP_variant(std::vector<Variant> & variants, int chr, int dna_pos, char ref, char variant) 
{
	Variant v ;
	v.type = pt_SNP ;
	v.position = dna_pos ;
	v.ref_len=1 ;
	v.variant_len=1;
	v.ref_str+=ref ;
	v.variant_str+=variant ;
}

inline void report_del_variant(std::vector<Variant> & variants, int chr, int dna_pos, int len, std::string & ref_str) 
{
}

inline void report_ins_variant(std::vector<Variant> & variants, int chr, int dna_pos, int len, std::string & variant_str) 
{
}

class VariantMap
{

public:
	VariantMap(Genome const &genome_) ;
	~VariantMap() ;

	void insert_variant(int chr, int pos, int ref_len, int variant_len, const std::string & ref_str, const std::string & variant_str);
	int init_from_sdis(std::string &sdi_fname);
	int report_to_sdi(std::string &sdi_fname);

	std::deque<Variant> * variantlist ;

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
