#pragma once

#include <palmapper/Hits.h>
#include <palmapper/QPalma.h>
#include <palmapper/Read.h>

class Mapper {

	friend class Hits;

public:

	class Result {
	public:
		Read const &_read;
		QPalma::Result _result;
		Hits _readMappings;
	};

	Mapper(Genome &genome, GenomeMaps &genomemaps, QueryFile &queryFile);
	~Mapper();

	int map_reads(Genome &genome, GenomeMaps &genomeMaps, QPalma* qpalma) ;
	int REDUNDANT;

protected:
	int init_constants()  ;
	int init_statistic_vars() ;
	int init_operators() ;
	int init_alignment_structures(Config * config);
	void map_read(Read &read, int count_reads, QPalma::Result &qpalmaResult);
	CHROMOSOME_ENTRY **GENOME; // doppelt

	Genome &_genome;
	GenomeMaps &_genomeMaps ;

	unsigned int LONGEST_HIT;
	static unsigned int SLOTS[2];

public:
	static unsigned int MAX_USED_SLOTS;
private:
	QueryFile &_queryFile;
    FILE *LEFTOVER_FP;
	FILE *ADAPTERTRIM_LOG_FP;
	unsigned int MAXHITS;
	int c_map_fast;
	int c_map_short_read;
};
