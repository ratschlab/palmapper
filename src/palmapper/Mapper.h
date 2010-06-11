#pragma once

#include <palmapper/Genome.h>
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
		  	_qpalma(read, mapper._qpalma)
		{
			_nr = nr;
		}

		int _nr;
		Read &_read;
		Hits _readMappings;
		QPalma::Result _qpalma;
	};

	class Reporter {
	public:

	};

	Mapper(Genome &genome, GenomeMaps &genomemaps, QueryFile &queryFile, QPalma &qpalma);
	~Mapper();

	int map_reads(FILE *OUT_FP, FILE *SP_OUT_FP) ;
	int REDUNDANT;

protected:
	void map_reads_timing(int count_reads, float this_read=-1);
	int init_constants()  ;
	int init_statistic_vars() ;
	int init_operators() ;
	int init_alignment_structures(Config * config);
	void map_read(Result &result, clock_t start_time);
	CHROMOSOME_ENTRY **GENOME; // doppelt

	Genome &_genome;
	GenomeMaps &_genomeMaps ;

	static unsigned int SLOTS[2];

public:
	static unsigned int MAX_USED_SLOTS;
private:
	QueryFile &_queryFile;
	QPalma &_qpalma;

	FILE *_OUT_FP;
	FILE *_SP_OUT_FP;
    FILE *_LEFTOVER_FP;
	FILE *_ADAPTERTRIM_LOG_FP;
	FILE *_TRIGGERED_LOG_FP;

	clock_t time1, time2a, time2b, time2c, time3;
	clock_t last_timing_report;

	unsigned int MAXHITS;
	int c_map_fast;
	int c_map_short_read;
};
