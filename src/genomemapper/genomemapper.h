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


// ##############################################################
// ####### GLOBAL VARIABLES #####################################
// ##############################################################

#include <genomemapper/Config.h>
#include <genomemapper/Genome.h>
#include <genomemapper/Hit.h>
#include <genomemapper/Read.h>
#include <genomemapper/Statistics.h>
#include <genomemapper/Util.h>
#include <genomemapper/TopAlignments.h>
#include <genomemapper/QPalma.h>
#include <genomemapper/GenomeMaps.h>

extern Config _config;
extern Statistics _stats;
extern Read _read;
extern Hits _hits;
extern Genome _genome;
extern TopAlignments _topalignments ;
extern QPalma _qpalma ;
extern GenomeMaps _genomemaps ;



#if 1 // dd
extern char HAS_SLOT;
extern int SLOT;
#endif // dd


// ##############################################################
// ####### FILE HANDLING ########################################
// ##############################################################

extern FILE *OUT_FP;
extern FILE *SP_OUT_FP;
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
// ####### HIT/EDIT #############################################
// ##############################################################


#define SCORE_INTERVAL 1
extern unsigned int NUM_SCORE_INTERVALS;



// ##############################################################
// ####### MEMORY CONTAINER #####################################
// ##############################################################


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
int check_mm(Chromosome const &chr, int genome_pos, int readpos, int ori);
extern int align_hit_simple(HIT* hit, int start, int end, int readpos, Chromosome const &chromosome, int orientation, unsigned char mismatches);
extern int prepare_kbound_alignment(HIT* hit, int start, int end, int readpos, Chromosome const &chromosome, char orientation, char mismatches);

//print.c

extern int compare_editops(const void *a, const void *b) ;
extern char get_compl_base(char c);
extern void print_stats();
extern int print_alignment(HIT* hit, unsigned int num);
extern int print_hits();
extern int print_perfect_hits(unsigned int num);
extern int print_largest_hit();
extern void print_leftovers(const char *tag, FILE *LEFTOVER_FP);
extern void print_alignment_matrix(int chrstart, int readstart, int length, int offset_front, int offset_end, Chromosome const &chr, char ori, int K);
extern void print_alignment_stats(int num_unspliced_best, int num_unspliced_suboptimal, int num_spliced_best, int num_spliced_suboptimal) ;
extern int compare_int(const void *a, const void *b) ;

#endif
