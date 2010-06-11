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

		int _rtrim_cut;
		int _polytrim_cut_start;
		int _polytrim_cut_end;
	};

	class Reporter {
	public:
		virtual void report(Result &result) = 0;
	};

	Mapper(Genome const &genome, GenomeMaps &genomemaps, QueryFile &queryFile, QPalma &qpalma, Reporter &reporter);
	~Mapper();
	void setProgressChar(char c) {
		_progressChar = c;
	}

	int map_reads() ;
	int REDUNDANT;

protected:
	void map_reads_timing(int count_reads, float this_read=-1);
	int init_constants()  ;
	int init_statistic_vars() ;
	int init_operators() ;
	int init_alignment_structures(Config * config);
	void map_read(Result &result, clock_t start_time);
	CHROMOSOME_ENTRY **GENOME;

	Genome const &_genome;
	GenomeMaps &_genomeMaps ;

private:
	QueryFile &_queryFile;
	QPalma const &_qpalma;
	Reporter &_reporter;

	int num_spliced_alignments_triggered;
	char _progressChar;

    FILE *_LEFTOVER_FP;
	FILE *_ADAPTERTRIM_LOG_FP;
	FILE *_TRIGGERED_LOG_FP;

	clock_t time1, time2a, time2b, time2c, time3;
	clock_t last_timing_report;

	unsigned int MAX_USED_SLOTS;
	unsigned int MAXHITS;
	int c_map_fast;
	int c_map_short_read;
	std::vector<bool> seed_covered ;
};
