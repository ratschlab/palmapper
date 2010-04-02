#pragma once

#include <genomemapper/Config.h>
#include <genomemapper/Chromosome.h>
#include <genomemapper/Read.h>

#define CONTAINER_SIZE 100000

typedef struct edit_op_structure {
	signed int pos;
	int mm: 1;
} EDIT_OPS;

typedef struct hit_structure {
	unsigned short int readpos;
	unsigned int start;
	unsigned int end;
	Chromosome const *chromosome;
	char orientation;
	unsigned char mismatches;	// including gaps!
	unsigned char gaps;
	signed char start_offset;
	signed char end_offset;
	EDIT_OPS edit_op[Config::MAX_EDIT_OPS];
	char aligned;
	struct hit_structure *same_eo_succ;	// the list of HITS_BY_SCORE - only forward pointer for now
	struct hit_structure *next;
	struct hit_structure *last;
} HIT;



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


typedef struct mapping_entry_container_structure {
	struct mapping_entry_structure entries[CONTAINER_SIZE];
	unsigned int used;
	struct mapping_entry_container_structure *next;
} MAPPING_ENTRY_CONTAINER;


typedef struct chromosome_entry_container_structure {
	CHROMOSOME_ENTRY* entries;
	unsigned int used;
} CHROMOSOME_ENTRY_CONTAINER;

typedef struct hits_by_score_structure {
	HIT *hitpointer;
	int num;
} HITS_BY_SCORE_STRUCT;

typedef struct hit_container_structure {
	HIT entries[CONTAINER_SIZE];
	unsigned int used;
	struct hit_container_structure *next;
} HIT_CONTAINER;

class Hits {

public:
	
	Hits();
	int map_reads() ;
	int size_hit(HIT *hit, unsigned int *oldlength, char num) ;
	int seed2genome(unsigned int num, unsigned int index_slot, unsigned int readpos, char reverse) ;
	void printgenome() ;
	int alloc_genome_memory() ;
	int duplicate(HIT* hit) ;
	int insert_into_scorelist(HIT* hit, char d) ;
	int browse_hits() ;
	int map_fast(Read & read, int & firstslot, int & firstpos);
	int map_short_read(Read & read, unsigned int num, int first_slot, int first_pos);
	int get_slot(Read & read, int pos);

	void printhits() ;
	int init_from_meta_index() ;
	int init_hit_lists()  ;

protected:
	int init_constants()  ;
	int init_statistic_vars() ;
	int init_operators() ;
	int dealloc_mapping_entries() ;
	int dealloc_hits() ;
	int dealloc_chromosome_entries() ;
	int dealloc_hit_lists_operator() ;
	int dealloc_hits_by_score() ;
	int alloc_readstart_bins() ;
	int alloc_hits_by_score() ;
	int alloc_hit_lists_operator() ;
	HIT_CONTAINER* alloc_hit_container() ;
	HIT* alloc_hit() ;
	CHROMOSOME_ENTRY_CONTAINER* alloc_chromosome_entry_container() ;
	CHROMOSOME_ENTRY* alloc_chromosome_entry(unsigned int pos, Chromosome const &chr, char strand) ;
	MAPPING_ENTRY_CONTAINER* alloc_mapping_entry_container() ;
	MAPPING_ENTRY* alloc_mapping_entry() ;
	
	

	CHROMOSOME_ENTRY **GENOME;

/*
	unsigned short int readpos;
	unsigned int start;
	unsigned int end;
	Chromosome const *chromosome;
	char orientation;
	unsigned char mismatches;	// including gaps!
	unsigned char gaps;
	signed char start_offset;
	signed char end_offset;
	EDIT_OPS edit_op[Config::MAX_EDIT_OPS];
	char aligned;
	HIT *same_eo_succ;	// the list of HITS_BY_SCORE - only forward pointer for now
	HIT *next;
	HIT *last;
*/

	unsigned int LONGEST_HIT;

public:
	HIT **HIT_LISTS_OPERATOR;
	HITS_BY_SCORE_STRUCT *HITS_BY_SCORE;
	unsigned int HITS_IN_SCORE_LIST;

protected:
	HIT **READSTART_BINS;


	MAPPING_ENTRY_CONTAINER* MAPPING_ENTRY_OPERATOR_FIRST;
	MAPPING_ENTRY_CONTAINER* MAPPING_ENTRY_OPERATOR;
	HIT_CONTAINER* HIT_OPERATOR_FIRST;
	HIT_CONTAINER* HIT_OPERATOR;
	CHROMOSOME_ENTRY_CONTAINER* CHROMOSOME_ENTRY_OPERATOR;
	
	unsigned int MAX_USED_SLOTS;
	unsigned int NUM_MAPPING_ENTRIES;
};
