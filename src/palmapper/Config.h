#pragma once

#include <stdint.h>
#include <string.h>
#include <string>
#include <vector>

#define MAX_INDEX_DEPTH 16
#define VERSION "0.6"

enum OutputFormatEnum
{
	OUTPUT_FORMAT_DEFAULT=-1,
	OUTPUT_FORMAT_SHORE=0,
	OUTPUT_FORMAT_BED=1,
	OUTPUT_FORMAT_BEDX=2,
	OUTPUT_FORMAT_SAM=3,
	OUTPUT_FORMAT_BAM=4,
}  ;

enum OutputFormatOptionEnum
{
	OUTPUT_FORMAT_OPTION_DEFAULT=-1,
	OUTPUT_FORMAT_OPTION_SORTNAME=1,
	OUTPUT_FORMAT_OPTION_SORTPOS=2,
}  ;

enum OutputFormatFlagsEnum
{
	OUTPUT_FORMAT_FLAGS_DEFAULT=15,
	OUTPUT_FORMAT_FLAGS_READ=1,
	OUTPUT_FORMAT_FLAGS_QUALITY=2,
	OUTPUT_FORMAT_FLAGS_SAMFLAGS=4,
	OUTPUT_FORMAT_FLAGS_MORESAMFLAGS=8,
	OUTPUT_FORMAT_FLAGS_MAQQUALITY=16,
}  ;

enum OutputFilterEnum
{
	OUTPUT_FILTER_DEFAULT=-1,
	OUTPUT_FILTER_ALL=0,
	OUTPUT_FILTER_TOP=1,
	OUTPUT_FILTER_LIMIT=2,
	OUTPUT_FILTER_RANDOM=3,
	OUTPUT_FILTER_MMDROP=4
}  ;

enum ProtocolEnum
{
	PROTOCOL_UNSTRANDED=0,
	PROTOCOL_FIRST=1,
	PROTOCOL_SECOND=2
}  ;

enum Personality {
	Palmapper,
	GenomeMapper,
};

const int DEFAULT_SETTING=123456789 ;

class Genome;

class Config {
public:
	
	Config();
	int parseCommandLine(int argc, char *argv[]);
	int applyDefaults(Genome* genome) ;
	int checkConfig() ;

	static int const REPORT_REPETITIVE_SEED_COUNT = 2 ;
	static int const REPORT_REPETITIVE_SEED_COUNT_MANY1 = 10 ;
	static int const REPORT_REPETITIVE_SEED_COUNT_MANY2 = 100 ;
	static int const REPORT_MAPPED_REGIONS_MIN_LENGTH = 25 ;

	static int const QPALMA_USE_MAP_WINDOW = 10 ; //Window used to go all away the extended long region to find mapped positions

	static int const MAX_EDIT_OPS = 55;
	static size_t const INDEX_SIZE_12 = 16777216; //4^12
	static size_t const INDEX_SIZE_13 = 67108864; //4^13
	static size_t const INDEX_SIZE_14 = 67108864*4; //4^14
	static size_t const INDEX_SIZE_15 = 67108864*4*4; //4^15
	static size_t const INDEX_SIZE_16 = 67108864L*4*4*4; //4^16
	static size_t const MAX_READ_LENGTH = 20000;
	static int const MAX_READ_ID_LENGTH = 1000;
	//static unsigned int const NUM_TOP_ALIGNMENTS = 10 ;

	Personality _personality;
	unsigned int BWA_INDEX;
	int NUM_THREADS;
	OutputFilterEnum OUTPUT_FILTER ;
	unsigned int OUTPUT_FILTER_NUM_TOP ;
	int OUTPUT_FILTER_NUM_LIMIT ;
	int OUTPUT_FILTER_DELTA_MMDROP ;
	int OUTPUT_FILTER_NUM_MMDROP ;
	
	//char ALL_HIT_STRATEGY;
	//char BEST_HIT_STRATEGY;
	//char SUMMARY_HIT_STRATEGY;

	char BSSEQ;

	unsigned int FIXTRIM_STRATEGY_LEN;
	unsigned int FIXTRIMRIGHT_STRATEGY_LEN;
	unsigned int FIXTRIMLEFT_STRATEGY_LEN;
	unsigned int RTRIM_STRATEGY;
	unsigned int RTRIM_STRATEGY_MIN_LEN;
	unsigned int RTRIM_STRATEGY_STEP;
	unsigned int POLYTRIM_STRATEGY;
	unsigned int POLYTRIM_STRATEGY_MIN_LEN;
	unsigned int POLYTRIM_STRATEGY_STEP ;
	unsigned int POLYTRIM_STRATEGY_POLY_MIN_LEN;
	unsigned int ADAPTERTRIM_STRATEGY ;
	unsigned int ADAPTERTRIM_STRATEGY_MIN_LEN;
	std::string ADAPTERTRIM_STRATEGY_LOG ;
	
	bool REALIGN_READS ;

	//int SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS[2] ;
	//int SUMMARY_HIT_STRATEGY_HIT_FOUND[2] ;
	unsigned int HITLEN_LIMIT;
	char VERBOSE;
	char MAP_REVERSE;
	//char REPEATMAP;
	char STRINGENT_GAPLIMIT;
	int PRINT_SEQ;
	unsigned int INDEX_DEPTH;
	unsigned int INDEX_DEPTH_EXTRA;
	unsigned int INDEX_SEED_STEP ;
	std::string READ_ID_PREFIX ;
	unsigned int FIRST_READ_NR;
	char* FIRST_READ_ID;
	unsigned int LAST_READ_NR;
	char* LAST_READ_ID;
	std::string ONLY_READ_ID ;

	unsigned int INDEX_DEPTH_EXTRA_THRESHOLD;
	unsigned int SEED_HIT_CANCEL_THRESHOLD;
	unsigned int SEED_HIT_TRUNCATE_THRESHOLD ;

	bool NOT_MAXIMAL_HITS;
	std::string CHR_INDEX_FILE_NAME;
	std::string INDEX_FWD_FILE_NAME;
	std::string INDEX_REV_FILE_NAME;
	std::string META_INDEX_FILE_NAME;
	std::string Q_QUERY_FILE_NAMES;
	std::string Q1_QUERY_FILE_NAMES;
	std::string Q2_QUERY_FILE_NAMES;
	std::vector<std::string> QUERY_FILE_NAMES;
	std::vector<int> QUERY_FILE_STRANDS;
	std::string OUT_FILE_NAME;
	std::string SPLICED_OUT_FILE_NAME;
	std::string GENOME_FILE_NAME;
	std::string LEFTOVER_FILE_NAME;
	std::string TRIGGERED_LOG_FILE;  // #A#
	OutputFormatEnum OUTPUT_FORMAT;
	std::string SAMTOOLS_PATH_NAME;
	OutputFormatOptionEnum OUTPUT_FORMAT_OPTION;
	unsigned int OUTPUT_FORMAT_FLAGS;
	bool INCLUDE_UNMAPPED_READS_SAM;
	
	std::string LOCK_FILE_NAME ;

	int READ_LENGTH_PARAM ;

	char * REPORT_FILE;
	int REPORT_FILE_READONLY;
	int REPORT_MAPPED_REGIONS;
	int REPORT_MAPPED_READS;
	int REPORT_SPLICED_READS;
	int REPORT_GENOME_COVERAGE ;
	
	char * REPORT_GENOME_COVERAGE_FILE;
	std::string REPORT_GFF_FILE_NAME ;

	std::string REPORT_JUNCTIONS_FILE ;
	int REPORT_JUNCTIONS;
	
	std::string MAP_JUNCTIONS_FILE ;
	int MAP_JUNCTIONS;
	int MAP_JUNCTIONS_COVERAGE;
	int MAP_JUNCTIONS_MIN_SEGMENT_LENGTH ;
	int MAP_JUNCTIONS_PSEUDO_ANNO_COV ;
	int MAP_JUNCTIONS_MAP_WINDOW ;
	bool MAP_JUNCTIONS_ONLY ;
	
	std::string ANNOTATED_SPLICE_SITES_FILE;
	int SCORE_ANNOTATED_SPLICE_SITES ;

	int REPORT_RESET;

	int QPALMA_USE_MAP;
	int QPALMA_USE_MAP_MAX_SIZE;
	int QPALMA_USE_SPLICE_SITES;
	int QPALMA_MIN_NUM_MATCHES ;
	float QPALMA_USE_SPLICE_SITES_THRESH_DON;
	float QPALMA_USE_SPLICE_SITES_THRESH_ACC;
	float QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC;
	bool QPALMA_PRB_OFFSET_FIX ;

	unsigned int READ_COUNT_LIMIT; // limits the number of reads for alignment


	bool LOG_TRIGGERED;  // #A#
	unsigned int FILTER_BY_MAX_MISMATCHES;
	unsigned int FILTER_BY_MAX_GAPS;
	bool FILTER_BY_SPLICE_SITES;
	int FILTER_BY_SPLICE_SITES_REGION;
	unsigned int FILTER_BY_SPLICE_SITES_EDIT_MIN;
	float FILTER_BY_SPLICE_SITES_THRESH_ACC;
	float FILTER_BY_SPLICE_SITES_THRESH_DON;
	float FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC;

	std::string QPALMA_FILE;
	std::string ACC_FILES;
	std::string DON_FILES;

	int NO_SPLICE_PREDICTIONS;
	int NO_SPLICE_DISCOVERY;
	bool NO_QPALMA;
	bool NO_GENOMEMAPPER ;
	double QPALMA_INDEL_PENALTY ;
	bool NO_QPALMA_SCORE_FIX ;
	
	std::vector<const char*> ACC_CONSENSUS, DON_CONSENSUS, ACC_CONSENSUS_REV, DON_CONSENSUS_REV ;
	bool non_consensus_search ;
	int non_consensus_search_gap ;
	float non_consensus_search_discount ;
	int MIN_NUM_MATCHES_PEN;

	bool NO_LOAD_GENOME ;
	bool NO_READ_MAPPING ;

	bool USE_VARIANTS ;
	bool DISCOVER_VARIANTS ;
	bool REPORT_VARIANTS ;
	bool REPORT_USED_VARIANTS ;
	bool IUPAC_SNPS ;
	bool IUPAC_GENOME ;
	bool CONVERT_SUBSTITUTION_VARIANTS ;
	bool VALIDATE_VARIANTS ;
	bool FILTER_VARIANTS ;
	bool MERGE_VARIANT_SOURCE_IDS ;
	bool UNIQUE_VARIANT_SOURCE_IDS ;
	bool USE_VARIANTS_EDITOP_FILTER ;

	std::string USE_VARIANT_FILE_NAME ;
	std::string REPORT_USED_VARIANT_FILE_NAME;
	std::string MGF_REF_NAME;
	std::string REPORT_VARIANTS_FILE_NAME ;
	std::string REPORT_VARIANTS_STATS_FILE_NAME ;
	std::vector<std::string> FILTER_VARIANT_SOURCES ;
	std::vector<std::string> VARIANT_SNP_TAKE_LINES ;
	
	std::string MAF_SOURCE_ID ;

	int FILTER_VARIANT_MINSOURCECOUNT ;
	std::vector<std::string> FILTER_VARIANT_REQSOURCES ;
	int FILTER_VARIANT_MINCONFCOUNT ;
	double FILTER_VARIANT_MAXNONCONFRATIO ;
	int FILTER_VARIANT_MINUSECOUNT ;
	int FILTER_VARIANT_MAP_WINDOW ;
	bool FILTER_VARIANT_JUNCTIONS ;
	int DP_MAX_DELETIONS ;
	int FILTER_VARIANT_VLEN ;
	int FILTER_VARIANT_MAXLEN ;
	
	int REPORT_SNP_TERMINAL_DIST ;
	int REPORT_INDEL_TERMINAL_DIST ;
	
	int NO_GAP_END;

	bool TRANSCRIBE_GFF ;
	std::string TRANSCRIBE_GFF_FILE, TRANSCRIBE_FASTA_FILE ;
	
	int INDEX_PRECACHE;
	unsigned int FLANKING;
	int NUM_EDIT_OPS;
	int NUM_MISMATCHES;
	int NUM_GAPS;
    int STRAND;
    ProtocolEnum PROTOCOL;
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
	int SPLICED_SHORTEST_INTRON_LENGTH;
	int UNSPLICED_MAX_NUM_ALIGNMENTS;
	int JUNCTION_MAX_NUM_ALIGNMENTS;
	int SPLICED_MAX_NUM_ALIGNMENTS;
	int SPLICED_CLUSTER_TOLERANCE;
	int SPLICED_MAX_INTRONS;
	int SPLICED_MIN_SEGMENT_LENGTH ;
	int JUNCTION_REMAP_MAX_INTRONS;

	char STATISTICS;
	unsigned int CHROM_CONTAINER_SIZE;

	static void VersionHeader() ;
	static int usage();
	static Personality getPersonality();

    std::string COMMAND_LINE;

    double TAG_MULTIMAPPERS ;
    int TAG_MULTIMAPPERS_QPALMA ;

private:
	int getInt(int &i, char *argv[]) const;
	int getString(int &i, char *argv[]) const;
	int postprocess_consensus_list(std::vector<const char *> & consensus_list) ;
	int postprocess_query_filenames(const std::string filenames, int strand) ;
};


inline void split_string(std::string text, std::vector<const char *>& words, char sep)
{
  int i=0;
  char ch;
  std::string word;

  while ((ch=text[i++]))
    {
      if (ch==sep)
	{
	  if (!word.empty())
	    {
	      words.push_back(strdup(word.c_str()));
	    }
	  word = "";
	}
      else
	{
	  word += ch;
	}
    }
  if (!word.empty())
    {
      words.push_back(strdup(word.c_str())) ;
    }
}
