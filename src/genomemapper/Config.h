#pragma once

#include <stdint.h>

#include <string>
#include <vector>

//#define MAX_READ_LENGTH 1000
//#define MAX_READ_ID_LENGTH 100
#define MAX_INDEX_DEPTH 13
#define ASSUMED_READ_LENGTH 76 // TODO Refactoring needed. This should be a parameter.
#define VERSION "0.4"


class Config {
public:
	Config();
	int parseCommandLine(int argc, char *argv[]);

	static int const REPORT_REPETITIVE_SEED_COUNT = 2 ;
	static int const REPORT_REPETITIVE_SEED_COUNT_MANY1 = 10 ;
	static int const REPORT_REPETITIVE_SEED_COUNT_MANY2 = 100 ;
	static int const REPORT_MAPPED_REGIONS_MIN_LENGTH = 25 ;

	static int const QPALMA_USE_MAP_WINDOW = 10 ; //Window used to go all away the extended long region to find mapped positions

	static int const MAX_EDIT_OPS = 10;
	static size_t const INDEX_SIZE_12 = 16777216; //4^12
	static size_t const INDEX_SIZE_13 = 67108864; //4^13
	static size_t const MAX_READ_LENGTH = 1000;
	static int const MAX_READ_ID_LENGTH = 100;
	static unsigned int const NUM_TOP_ALIGNMENTS = 10 ;

	int NUM_THREADS;
	char ALL_HIT_STRATEGY;
	char SUMMARY_HIT_STRATEGY;

	unsigned int RTRIM_STRATEGY;
	unsigned int RTRIM_STRATEGY_MIN_LEN;
	unsigned int POLYTRIM_STRATEGY;
	unsigned int POLYTRIM_STRATEGY_MIN_LEN;

	//int SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS[2] ;
	//int SUMMARY_HIT_STRATEGY_HIT_FOUND[2] ;
	int HITLEN_LIMIT;
	char VERBOSE;
	char MAP_REVERSE;
	char REPEATMAP;
	char STRINGENT_GAPLIMIT;
	int PRINT_SEQ;
	int INDEX_DEPTH;
	int INDEX_DEPTH_EXTRA;

	int INDEX_DEPTH_EXTRA_THRESHOLD;
	int SEED_HIT_CANCEL_THRESHOLD;
	bool NOT_MAXIMAL_HITS;
	std::string CHR_INDEX_FILE_NAME;
	std::string INDEX_FWD_FILE_NAME;
	std::string INDEX_REV_FILE_NAME;
	std::string META_INDEX_FILE_NAME;
	std::string QUERY_FILE_NAME;
	std::string OUT_FILE_NAME;
	std::string SPLICED_OUT_FILE_NAME;
	std::string GENOME_FILE_NAME;
	std::string LEFTOVER_FILE_NAME;
	std::string TRIGGERED_LOG_FILE;  // #A#
	char OUTPUT_FORMAT;

	char * REPORT_FILE;
	int REPORT_FILE_READONLY;
	int REPORT_REPETITIVE_SEEDS;
	int REPORT_MAPPED_REGIONS;
	int REPORT_MAPPED_READS;
	int REPORT_SPLICED_READS;

	int REPORT_RESET;

	int QPALMA_USE_MAP;
	int QPALMA_USE_MAP_MAX_SIZE;
	int QPALMA_USE_SPLICE_SITES;
	float QPALMA_USE_SPLICE_SITES_THRESH_DON;
	float QPALMA_USE_SPLICE_SITES_THRESH_ACC;
	float QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC;

	int READ_COUNT_LIMIT; // limits the number of reads for alignment


	bool LOG_TRIGGERED;  // #A#
	unsigned int FILTER_BY_MAX_MISMATCHES;
	unsigned int FILTER_BY_MAX_GAPS;
	unsigned int FILTER_BY_SPLICE_SITES;
	int FILTER_BY_SPLICE_SITES_REGION;
	unsigned int FILTER_BY_SPLICE_SITES_EDIT_MIN;
	float FILTER_BY_SPLICE_SITES_THRESH_ACC;
	float FILTER_BY_SPLICE_SITES_THRESH_DON;
	float FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC;

	std::string QPALMA_FILE;
	std::string ACC_FILES;
	std::string DON_FILES;
	int NO_SPLICE_PREDICTIONS;

	int INDEX_PRECACHE;
	unsigned int FLANKING;
	int NUM_EDIT_OPS;
	int NUM_MISMATCHES;
	int NUM_GAPS;
	double MM_SCORE;
	double M_SCORE;
	double GAP_SCORE;
	char GAPS_MOST_RIGHT;
	char OVERHANG_ALIGNMENT;
	char SCORES_OUT;

	bool SPLICED_HITS;
	int SPLICED_HIT_MIN_LENGTH_SHORT;
	int SPLICED_HIT_MIN_LENGTH_COMB;
	int SPLICED_HIT_MIN_LENGTH_LONG;
	int SPLICED_LONGEST_INTRON_LENGTH;
	int SPLICED_MAX_NUM_ALIGNMENTS;
	int SPLICED_CLUSTER_TOLERANCE;
	int SPLICED_MAX_INTRONS;

	char STATISTICS;
	unsigned int CHROM_CONTAINER_SIZE;

private:
	int usage();
	int getInt(int &i, char *argv[]) const;
	int getString(int &i, char *argv[]) const;
};
