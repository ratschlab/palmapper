#pragma once

#include <genomemapper/Config.h>
#include <genomemapper/Statistics.h>
#include <genomemapper/Structs.h>
#include <genomemapper/Util.h>

class Read {
public:
	Read();

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
