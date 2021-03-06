// Authors: Korbinian Schneeberger, Stephan Ossowski and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#if defined CONF
#else
#define CONF

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <palmapper/VariantMap.h>
#include <palmapper/GenomeMaps.h>

#define VERSION "0.3.0"

// ##############################################################
// ####### GLOBAL VARIABLES #####################################
// ##############################################################

#define MIN_INDEX_DEPTH 5
#define MAX_INDEX_DEPTH 16

extern char VERBOSE;
extern char HAS_SLOT;
extern unsigned int SLOT;

extern int NUM_THREADS ;

extern int INDEX_DEPTH;
extern long int POWER[MAX_INDEX_DEPTH];
extern long int BINARY_CODE[5];

extern unsigned int LONGEST_CHROMOSOME;

extern int debug;

// ##############################################################
// ####### FILE HANDLING ########################################
// ##############################################################

extern char GENOME_FILE_NAME[500];
extern char GENOME_VARIANTS_FILE_NAME[500];
extern char CHR_INDEX_FILE_NAME[500];
extern char MAPFWD_INDEX_FILE_NAME[500];
extern char META_INDEX_FILE_NAME[500];
extern char GENOME_OUT_FILE_NAME[500];
//char OCC_FILE_NAME[500];

extern FILE *CHR_INDEX_FP;
extern FILE *MAPFWD_INDEX_FP;
extern FILE *META_INDEX_FP;
extern FILE *GENOME_FP;
extern FILE *GENOME_OUT_FP;
//FILE *OCC_FP;

// ##############################################################
// ######## SEQUENCE STORAGE ####################################
// ##############################################################

#define CHR_DESC_LENGTH 50

extern unsigned int CHR_LENGTH;
extern char *CHR_SEQ;

extern char CHR_DESC[2000];
extern char CHR_DESC_TMP[2000];

extern unsigned long int POSITION_COUNTER;
extern int MAX_SOURCE_COMBINATIONS ;

// ##############################################################
// ####### INDEX MANAGEMENT #####################################
// ##############################################################

#define BIN_SIZE 3
#define BIN_SIZE_EXT 20
#define INDEX_SIZE_12 16777216 //4^12
#define INDEX_SIZE_13 67108864 //4^13
#define INDEX_SIZE_14 268435456L //4^14
#define INDEX_SIZE_15 (268435456L*4) //4^15
#define INDEX_SIZE_16 (268435456L*4*4) //4^16
extern unsigned long int INDEX_SIZE ;

#define BLOCK_TABLE_SIZE 16777216	// 2^24 (3 Byte)
#define BLOCK_SIZE 256	// 2^8 (1 Byte)

typedef struct id_structure {
	unsigned char id[4];
} ID;

typedef struct position_structure {
	unsigned int pos;
	unsigned int chr;
} POS;

extern POS *BLOCK_TABLE;
extern unsigned int BLOCK;
extern unsigned short int POSITION;

typedef struct bin_extension_structure {
	ID ids[BIN_SIZE_EXT];
	struct bin_extension_structure *bin_ext;
} BIN_EXT;

typedef struct bin_structure {
	unsigned int num_pos;
	ID ids[BIN_SIZE];
	BIN_EXT *bin_ext;
	BIN_EXT *last_bin_ext;
} BIN;

//BIN *INDEX;
//BIN *INDEX_REV;

//extern BIN *INDEX[INDEX_SIZE];
extern BIN **INDEX;//[INDEX_SIZE];

extern long int NUM_USED_SLOTS; //different to SLOT_COUNTER! This counts the number of different used slots, used by reverse and forward Index.
extern unsigned int *USED_SLOTS;//[INDEX_SIZE];


// ##############################################################
// ####### MEMORY MANAGEMENT ####################################
// ##############################################################

#define BIN_EXT_PER_NUGGET 100000

typedef struct nugget_str {
	BIN_EXT buffer[BIN_EXT_PER_NUGGET];
	struct nugget_str *next;
} NUGGET;

typedef struct stmg_str {
	int curr_num;
	NUGGET *nuggets;
} STORAGE;

extern STORAGE *MEM_MGR;

// ##############################################################
// ####### ROUTINES #############################################
// ##############################################################

//init.c
extern int init(int argc, char *argv[]);

//usage.c
extern int usage();

//load.c
extern int load_chromosomes(Genome & genome, VariantMap & variants, GenomeMaps & genome_mask);
extern int desc_parsing(char *c);

//indec.c
extern int index_chromosome(unsigned int chr, Genome & genome, VariantMap & variants, GenomeMaps & genome_mask, bool mask_do_alloc=true, bool mask_do_secondary=false, bool mask_do_add=true);
extern int index_chromosome_novariants(unsigned int chr, Genome & genome, GenomeMaps & genome_mask, bool mask_do_alloc=true, bool use_secondary_regions=false, bool mask_do_add=true)  ;

//alloc.c
extern int alloc_bin(unsigned int slot);
extern BIN_EXT *alloc_bin_ext() ;
extern void alloc_blocktable();
extern int dealloc_chr();

//write.c
extern int write_meta_index(unsigned int num_chr);
//int write_index();
extern int write_chr_desc(unsigned int chr);

//printindex.c
extern void printindex();

extern int seed_hit_cancel_threshold ;

extern int has_genome_mask ;
extern char GENOME_MASK_FILE_NAME[500] ;
extern char GENOME_MASK_GFF_FILE_NAME[500] ;

#define MASK_REGION_PRIMARY MASK_MAPPED_READ_BEST
#define MASK_REGION_SECONDARY MASK_MAPPED_READ
#define MASK_REGION_SECONDARY_REGION MASK_SPLICED_READ_BEST

extern bool genome_mask_use_rep_seeds; 
extern bool genome_mask_use_secondary_regions ;
extern int secondary_min_num_hits ;
extern int secondary_region_extra ;
extern int genome_mask_gff_extra ;
extern int secondary_min_hits_perc ;
#endif
