#pragma once

#include <genomemapper/Config.h>
#include <genomemapper/Statistics.h>
#include <genomemapper/Structs.h>
#include <genomemapper/Util.h>

class Read {
public:
	Read();
	Read(const Read & read); // copy constructor

	unsigned int length() {
		return READ_LENGTH;
	}

	char const *data() {
		return READ;
	}

	char const *id() {
		return READ_ID;
	}

	char const * const *quality() {
		return READ_QUALITY;
	}

	int pe_flag() {
		return READ_PE_FLAG;
	}

	char format() {
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


	void find_poly(int &poly_length_start, int &poly_length_end, float frac=0.8)
	{
		int num_t_start = 0 ;
		int num_a_end = 0 ;
		poly_length_start = 0 ;
		poly_length_end = 0 ;
		
		for (unsigned int i=0; i<READ_LENGTH; i++)
		{
			if (READ[i]=='T' || READ[i]=='t')
				num_t_start++ ;
			if (READ[READ_LENGTH-i]=='A' || READ[READ_LENGTH-i]=='a')
				num_a_end++ ;
			if (((float)num_t_start)/i >= frac)
				poly_length_start = i ;
			if (((float)num_a_end)/i >= frac)
				poly_length_end = i ;
		}
	}
	
	int read_short_read(FILE *QUERY_FP);
	int map_fast(int & firstslot, int & firstpos);
	int map_short_read(unsigned int num, int first_slot, int first_pos);
	int get_slot(int pos);


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

/*
#ifndef BinaryStream_MAP
	STORAGE_ENTRY *INDEX_REV_MMAP;
	STORAGE_ENTRY *INDEX_FWD_MMAP;
#else
	CBinaryStream<STORAGE_ENTRY>* INDEX_REV_MMAP;
	CBinaryStream<STORAGE_ENTRY>* INDEX_FWD_MMAP;
#endif
	INDEX_ENTRY *INDEX;
	INDEX_ENTRY *INDEX_REV;

	bool HAS_SLOT;
	int SLOT;
*/
};
