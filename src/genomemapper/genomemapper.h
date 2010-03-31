// Authors: Korbinian Schneeberger, Stephan Ossowski and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#ifndef GENOMEMAPPER_H_
#define GENOMEMAPPER_H_

#include "shogun/init.h"
#include "shogun/common.h"
#include "shogun/io.h"
#include "shogun/Array.h"
#include "DNAArray.h"
#include "DNAArray4.h"
#include "dyn_prog/qpalma_dp.h"

#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <string>


//#define BinaryStream_MAP

#ifdef BinaryStream_MAP
#include "BinaryStream.h"
#endif


#define THREADS 1

// ##############################################################
// ####### GLOBAL VARIABLES #####################################
// ##############################################################

#include <genomemapper/Config.h>
#include <genomemapper/Genome.h>
#include <genomemapper/Hit.h>
#include <genomemapper/Read.h>
#include <genomemapper/Statistics.h>
#include <genomemapper/Util.h>

extern Config _config;
extern Statistics _stats;
extern Read _read;
extern Genome _genome;

#if 1 // dd
extern char HAS_SLOT;
extern int SLOT;
#endif // dd
extern int BINARY_CODE[4];

extern unsigned int NUM_POS;

extern unsigned int LONGEST_HIT;

extern int REPORT_REPETITIVE_SEED_DEPTH_EXTRA ;

// ##############################################################
// ####### FILE HANDLING ########################################
// ##############################################################

//extern FILE *CHR_INDEX_FP;
//extern FILE *META_INDEX_FP;
//extern FILE *QUERY_FP;
extern FILE *OUT_FP;
extern FILE *SP_OUT_FP;
//extern FILE *LEFTOVER_FP;
extern FILE *TRIGGERED_LOG_FP; // #A#

// ##############################################################
// ####### MEMORY MANAGEMENT ####################################
// ##############################################################

extern size_t INDEX_SIZE ;

typedef struct meta_idx_file_entry {
	int slot;
	unsigned int num;
} META_INDEX_ENTRY;

typedef struct position_structure {
	unsigned int pos;
	unsigned int chr;
} POS;

extern POS *BLOCK_TABLE;
extern unsigned int BLOCK_TABLE_SIZE;


#if 1 // dd
extern INDEX_ENTRY *INDEX;
extern INDEX_ENTRY *INDEX_REV;

#ifndef BinaryStream_MAP
extern STORAGE_ENTRY *INDEX_REV_MMAP;
extern STORAGE_ENTRY *INDEX_FWD_MMAP;
#else
extern CBinaryStream<STORAGE_ENTRY>* INDEX_REV_MMAP;
extern CBinaryStream<STORAGE_ENTRY>* INDEX_FWD_MMAP;
#endif
#endif // dd

extern unsigned long int MAX_POSITIONS;

// ##############################################################
// ####### SHORT READ ###########################################
// ##############################################################

extern int REDUNDANT;
extern unsigned long int linenr;

extern char FLANK_SEQ[Config::MAX_READ_LENGTH + 200];

// ##############################################################
// ####### ALIGNMENT ############################################
// ##############################################################

extern char *chrseq;	// for debugging
extern char *ALIGNSEQ;
extern int NUM_MATCHES; // for new version of QPalma 

#define DIAGONAL 'D'
#define LEFT 'L'
#define UP 'U'
extern double WORST_SCORE;
extern double WORST_MM_SCORE;

// ##############################################################
// ####### SPLICED HITTING ######################################
// ##############################################################

typedef struct alignment_t {
  double qpalma_score;
  uint32_t num_matches;
  uint32_t num_gaps;
  char read_anno[4*Config::MAX_READ_LENGTH] ;
  std::vector<int> exons;
  Chromosome const *chromosome;
  char orientation;
  char strand ;
  char read_id[Config::MAX_READ_ID_LENGTH];
  int min_exon_len ;
  int max_intron_len ;
  bool spliced ;
  //bool rescued ;
  //int rescue_start, rescue_end ;
} ALIGNMENT;

extern std::vector<alignment_t *> top_alignments ;

struct region_t {
	int32_t start;
	int32_t end;
	bool erased ;
	bool from_map ;
	char orientation ;
	bool* read_map ;
	//int32_t chromosome ;
	//char strand ;
};
extern std::vector<std::vector<region_t *> > regions[];


// ##############################################################
// ####### HIT/EDIT #############################################
// ##############################################################

typedef struct hits_by_score_structure {
	HIT *hitpointer;
	int num;
} HITS_BY_SCORE_STRUCT;

#define SCORE_INTERVAL 1
extern int NUM_SCORE_INTERVALS;

extern HIT **HIT_LISTS_OPERATOR;
extern HIT **READSTART_BINS;
extern HITS_BY_SCORE_STRUCT *HITS_BY_SCORE;
extern unsigned int HITS_IN_SCORE_LIST;

typedef struct mapping_entry_structure {
	unsigned int readpos;
	HIT *hit;
	struct mapping_entry_structure *pred;
	struct mapping_entry_structure *succ;
} MAPPING_ENTRY;

typedef struct chromosome_entry {

	int chromosome; // @TODO? Minus values indicate reverse hits (would save strand var = 1 byte, ~15MB)
	unsigned int genome_pos;
	char strand;

	struct chromosome_entry *next;
	MAPPING_ENTRY *mapping_entries;

	// It seems to be cheaper to store the back-pointer information (pos)
	// in each of these entries rather than having a superior structure.
} CHROMOSOME_ENTRY;

extern CHROMOSOME_ENTRY **GENOME;

extern unsigned int LONGEST_CHROMOSOME;

// ##############################################################
// ####### MEMORY CONTAINER #####################################
// ##############################################################

#define CONTAINER_SIZE 100000

typedef struct mapping_entry_container_structure {
	struct mapping_entry_structure entries[CONTAINER_SIZE];
	unsigned int used;
	struct mapping_entry_container_structure *next;
} MAPPING_ENTRY_CONTAINER;

typedef struct hit_container_structure {
	HIT entries[CONTAINER_SIZE];
	unsigned int used;
	struct hit_container_structure *next;
} HIT_CONTAINER;

typedef struct chromosome_entry_container_structure {
	CHROMOSOME_ENTRY* entries;
	unsigned int used;
} CHROMOSOME_ENTRY_CONTAINER;

extern MAPPING_ENTRY_CONTAINER* MAPPING_ENTRY_OPERATOR_FIRST;
extern MAPPING_ENTRY_CONTAINER* MAPPING_ENTRY_OPERATOR;
extern HIT_CONTAINER* HIT_OPERATOR_FIRST;
extern HIT_CONTAINER* HIT_OPERATOR;
extern CHROMOSOME_ENTRY_CONTAINER* CHROMOSOME_ENTRY_OPERATOR;

extern unsigned int MAX_USED_SLOTS;
extern unsigned int NUM_MAPPING_ENTRIES;

// ##############################################################
// ####### ROUTINES #############################################
// ##############################################################

// genomemmapper.c
extern int map_reads();

// usage.c
extern int usage();

// init.c
extern int init(int argc, char *argv[]);
extern int init_from_meta_index();

// index.c
//dd extern int build_index();
#ifndef BinaryStream_MAP
extern void index_pre_buffer(STORAGE_ENTRY* index_mmap, STORAGE_ENTRY* buffer, long index, long size) ;
#endif

// alloc.c
extern int alloc_index_memory();
extern int alloc_genome_memory();
extern MAPPING_ENTRY* alloc_mapping_entry();
extern CHROMOSOME_ENTRY* alloc_chromosome_entry(unsigned int pos, Chromosome const &chr, char strand);
extern HIT* alloc_hit();
extern int alloc_hit_lists_operator();
extern int alloc_hits_by_score();
extern MAPPING_ENTRY_CONTAINER* alloc_mapping_entry_container();
extern CHROMOSOME_ENTRY_CONTAINER* alloc_chromosome_entry_container();
extern HIT_CONTAINER* alloc_hit_container();
extern int dealloc_mapping_entries();
extern int dealloc_hits();
extern int dealloc_chromosome_entries();
extern int dealloc_hit_lists_operator();
extern int dealloc_hits_by_score();

//hit.c
extern int map_reads();
extern int seed2genome(unsigned int num, unsigned int slot, unsigned int readpos, char reverse);
extern int insert_into_scorelist(HIT* hit, char d);
extern char *get_seq(unsigned int n);
extern void update_num_edit_ops(int num_edit_ops, char & all_hit_strategy, int & NUM_EDIT_OPS_);

//align.c
#if THREADS == 1
int check_mm(Chromosome const &chr, int genome_pos, int readpos, int ori);
extern int align_hit_simple(HIT* hit, int start, int end, int readpos, Chromosome const &chromosome, int orientation, unsigned char mismatches);
extern int prepare_kbound_alignment(HIT* hit, int start, int end, int readpos, Chromosome const &chromosome, char orientation, char mismatches);
#else
extern int align_hit_simple(HIT* hit);
extern int prepare_kbound_alignment (HIT* hit);
#endif

//print.c

extern int compare_editops(const void *a, const void *b) ;
extern char get_compl_base(char c);
extern void print_stats();
extern int print_alignment(HIT* hit, unsigned int num);
extern int print_hits();
extern int print_perfect_hits(unsigned int num);
extern int print_largest_hit();
extern void print_leftovers(const char *tag, FILE *LEFTOVER_FP);
extern void printhits();
extern void printhit(HIT* hit);
extern void print_alignment_matrix(int chrstart, int readstart, int length, int offset_front, int offset_end, Chromosome const &chr, char ori, int K);
extern void print_alignment_records(std::vector<alignment_t *> &hits, int num_unspliced_alignments, int num_spliced_alignments, int RTRIM_STRATEGY_CUT) ;
extern int compare_int(const void *a, const void *b) ;


// qpalma_init.cpp

struct alignment_parameter_struct {
	struct penalty_struct h, a, d, *qualityPlifs;
	int num_qualityPlifs;
	double *matchmatrix;
	int matchmatrix_dim[2];
	int quality_offset ;
} ;

extern struct alignment_parameter_struct *alignment_parameters;

extern int init_qpalma() ;
extern int clean_qpalma() ;
extern int map_splice_sites(std::string file_template, char type, float &splice_site_threshold, bool estimate_thresh, bool do_report) ;

// qpalma_filter.cpp

extern int qpalma_filter(struct alignment_t *ali, int num_N) ;
extern void qpalma_filter_stat(bool spliced) ;
extern void qpalma_filter_stat_report() ;

// qpalma_alignment.cpp

extern int get_splicesite_positions(std::string file_template, const char *type, Chromosome const &chr, char strand, int start, int end, float thresh, bool store_pos,
									std::vector<int> &positions) ;
extern int capture_hits();
extern int perform_alignment(std::string &read_string, std::string &read_quality, std::string &dna, std::vector<region_t *> &regions, std::vector<int> &positions,
							 Chromosome const &contig_id, char strand, int ori, int & num_reported,int hit_read, int hit_dna, int hit_length) ;
extern void capture_hits_timing(int read_count=-1, float this_read=-1.0) ;

extern float score_unspliced(const char * read_anno) ;


// filter.cpp

extern u_int8_t report_unspliced_hit(HIT *hit)  ;

extern void start_best_alignment_record();
extern void end_best_alignment_record(int RTRIM_STRATEGY_CUT);
extern void add_alignment_record(alignment_t *alignment, int num_alignments);


#include "report_maps.h"

#endif
