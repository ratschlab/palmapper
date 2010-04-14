#pragma once

#include <genomemapper/Config.h>
#include <genomemapper/Statistics.h>
#include <genomemapper/Util.h>

#include <string>

class Read {
public:
	Read();
	Read(const Read & read); // copy constructor

	unsigned int length() const {
		return READ_LENGTH;
	}
	unsigned int max_length() const {
		return ASSUMED_READ_LENGTH;
	}

	char *data() {
		return READ;
	}

	char const *data() const {
		return READ;
	}

	char const *id() const {
		return READ_ID;
	}

	char const * const *quality() const {
		return READ_QUALITY;
	}

	int pe_flag() const {
		return READ_PE_FLAG;
	}

	char format() const {
		return READ_FORMAT;
	}

	void cutOffLast() {
		--READ_LENGTH;
		READ[READ_LENGTH] = '\0';
		READ_QUALITY[0][READ_LENGTH] = '\0';
	}

	void cutOffFirst() {
		for (unsigned int i=0; i<READ_LENGTH; i++)
		{
			READ[i] = READ[i+1];
			READ_QUALITY[0][i] = READ_QUALITY[0][i+1];
		}
		--READ_LENGTH;
	}

	void trim_read_start(Read * read, int trim_start) 
	{
		for (unsigned int i=trim_start; i<read->length(); i++)
		{
			READ[i-trim_start] = read->READ[i];
			READ_QUALITY[0][i-trim_start] = read->READ_QUALITY[0][i];
		}
		READ_LENGTH = read->length()-trim_start ;

		READ[READ_LENGTH] = '\0';
		READ_QUALITY[0][READ_LENGTH] = '\0';
	}

	void trim_read_end(Read * read, int trim_end) 
	{
		for (unsigned int i=0; i<read->length()-trim_end; i++)
		{
			READ[i] = read->READ[i];
			READ_QUALITY[0][i] = read->READ_QUALITY[0][i];
		}
		READ_LENGTH = read->length()-trim_end ;

		READ[READ_LENGTH] = '\0';
		READ_QUALITY[0][READ_LENGTH] = '\0';
	}


	void find_poly(unsigned int &poly_length_start, unsigned int &poly_length_end, float frac=0.8)
	{
		int num_a_start = 0 ;
		int num_t_start = 0 ;
		int num_a_end = 0 ;
		int num_t_end = 0 ;
		poly_length_start = 0 ;
		poly_length_end = 0 ;
		
		for (unsigned int i=0; i<READ_LENGTH; i++)
		{
			if (READ[i]=='T' || READ[i]=='t')
				num_t_start++ ;
			if (READ[i]=='A' || READ[i]=='a')
				num_a_start++ ;
			if (READ[READ_LENGTH-i]=='A' || READ[READ_LENGTH-i]=='a')
				num_a_end++ ;
			if (READ[READ_LENGTH-i]=='T' || READ[READ_LENGTH-i]=='t')
				num_t_end++ ;
			if (((float)num_t_start)/i >= frac)
				poly_length_start = i ;
			if (((float)num_a_start)/i >= frac)
				poly_length_start = i ;
			if (((float)num_a_end)/i >= frac)
				poly_length_end = i ;
			if (((float)num_t_end)/i >= frac)
				poly_length_end = i ;
		}
	}

	bool is_full_poly(float frac=0.9)
		{
			int num_t=0, num_a=0 ;
			
			for (unsigned int i=0; i<READ_LENGTH; i++)
			{
				if (READ[i]=='T' || READ[i]=='t')
					num_t++ ;
				if (READ[i]=='A' || READ[i]=='a')
					num_a++ ;
			}
			if (((float)num_a/READ_LENGTH)>=frac)
				return true ;
			if (((float)num_t/READ_LENGTH)>=frac)
				return true ;
			return false ;
		}
	
	int read_short_read(FILE *QUERY_FP);

	int determine_read_length(const std::string & query_fname) ;

	Read* get_orig() { return orig_read ; } ;
	void set_orig(Read* orig) { orig_read=orig ; } ;

private:
//	Statistics _stats;
	unsigned long int linenr;
	unsigned int READ_LENGTH;
//	Config &_config;

	char *READ_QUALITY[3];
	char READ[Config::MAX_READ_LENGTH + 1];
	char READ_FORMAT;	// 0: fq, 1: fa, 2: flat
	char *READ_ID;
	int READ_PE_FLAG;
	
	Read* orig_read ;

	unsigned int ASSUMED_READ_LENGTH ;
};
