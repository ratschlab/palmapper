#pragma once

#include <palmapper/Hits.h>
#include <palmapper/QPalma.h>
#include <palmapper/Read.h>
#include <palmapper/TopAlignments.h>

class Mapper {

	friend class Hits;

public:

	class Result {
	public:
		Result(int nr, Read &read, Mapper &mapper)
		: 	_read(read),
		  	_readMappings(mapper._genome, mapper._genomeMaps, mapper, read),
		  	_qpalma(read, mapper._qpalma),
		  	_topAlignments(&mapper._genomeMaps)
		{
			_nr = nr;
		}

		int _nr;
		Read &_read;
		Hits _readMappings;
		QPalma::Result _qpalma;
		TopAlignments _topAlignments;
	};

	Mapper(Genome &genome, GenomeMaps &genomemaps, QueryFile &queryFile, QPalma &qpalma);
	~Mapper();

	int map_reads(Genome &genome, GenomeMaps &genomeMaps, QPalma* qpalma) ;
	int REDUNDANT;

protected:
	int init_constants()  ;
	int init_statistic_vars() ;
	int init_operators() ;
	int init_alignment_structures(Config * config);
	void map_read(Result &result);
	CHROMOSOME_ENTRY **GENOME; // doppelt

	Genome &_genome;
	GenomeMaps &_genomeMaps ;

	unsigned int LONGEST_HIT;
	static unsigned int SLOTS[2];

public:
	static unsigned int MAX_USED_SLOTS;
private:
	QueryFile &_queryFile;
	QPalma &_qpalma;
    FILE *LEFTOVER_FP;
	FILE *ADAPTERTRIM_LOG_FP;
	unsigned int MAXHITS;
	int c_map_fast;
	int c_map_short_read;
};
