#pragma once

#include <assert.h>
#include <genomemapper/Config.h>
#include <genomemapper/Chromosome.h>
#include <genomemapper/Read.h>

#define CONTAINER_SIZE 100000
#define SCORE_INTERVAL 1

extern Config _config ;

extern char *get_seq(unsigned int n);


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

struct CHROMOSOME_ENTRY {
	int chromosome; // @TODO? Minus values indicate reverse hits (would save strand var = 1 byte, ~15MB)
	unsigned int genome_pos;
	char strand;

	CHROMOSOME_ENTRY *next;
	MAPPING_ENTRY *mapping_entries;

	// It seems to be cheaper to store the back-pointer information (pos)
	// in each of these entries rather than having a superior structure.
};


typedef struct mapping_entry_container_structure {
	struct mapping_entry_structure entries[CONTAINER_SIZE];
	unsigned int used;
	struct mapping_entry_container_structure *next;
} MAPPING_ENTRY_CONTAINER;


struct CHROMOSOME_ENTRY_CONTAINER {
	CHROMOSOME_ENTRY_CONTAINER(int nrEntries) {
		entries = new CHROMOSOME_ENTRY[nrEntries];
		used = 0;
	}
	CHROMOSOME_ENTRY* entries;
	unsigned int used;
};

typedef struct hits_by_score_structure {
	HIT *hitpointer;
	int num;
} HITS_BY_SCORE_STRUCT;

typedef struct hit_container_structure {
	HIT entries[CONTAINER_SIZE];
	unsigned int used;
	struct hit_container_structure *next;
} HIT_CONTAINER;

class QPalma ;
class Genome ;
class GenomeMaps ;
class TopAlignments ;


class Hits {

public:
	
	Hits(Genome &genome, GenomeMaps &genomemaps);

	int map_reads(TopAlignments* topalignments, QPalma* qpalma) ;
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

	int analyze_hits(TopAlignments* topalignments, QPalma * qpalma) ;
	int report_read_alignment(HIT* hit, int nbest)  ;

	int REDUNDANT;

	int align_hit_simple(HIT* hit, int start, int end, int readpos, Chromosome const &chromosome, int orientation, unsigned char mismatches) ;
	int prepare_kbound_alignment(HIT* hit, int start, int end, int readpos, Chromosome const &chromosome, char orientation, char mismatches) ;

	inline void reset_num_edit_ops(int num_edit_ops)
	{
		SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.resize(0) ;
	}

	inline void update_num_edit_ops(int num_edit_ops, char & all_hit_strategy, int & NUM_EDIT_OPS_)
	{
		assert(num_edit_ops<Config::MAX_EDIT_OPS) ;
		
		if (!all_hit_strategy && num_edit_ops < NUM_EDIT_OPS_)
			NUM_EDIT_OPS_ = num_edit_ops ;
		
		if (_config.OUTPUT_FILTER==OUTPUT_FILTER_TOP)
		{
			// keep a list of minimal edit operation (assuming that each hit is only reported once)
			// the list is twice as long as NUM_TOP_ALIGNMENTS, as will be filtered later according to
			// the qpalma-alignment score (this is a heuristic, which may work well enough; alternatively, one
			// would need to compute the qpalma score for each hit, which will be too expensive
			
			bool inserted = false ;
			
			std::vector<int>::iterator it = SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.begin();
			
			for (uint8_t i = 0; i < SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.size(); i++, it++)
			{
				if ( num_edit_ops < SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS[i] )
				{
					SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.insert(it, num_edit_ops);
					inserted = true;
					
					if (SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.size() > _config.OUTPUT_FILTER_NUM_TOP*2)
						SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.pop_back();
					
					break;
				}
			}
			if (!inserted && SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.size() < _config.OUTPUT_FILTER_NUM_TOP*2)
			{
				SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.push_back(num_edit_ops) ;
				inserted = true ;
			}
		}
		
		if (!all_hit_strategy && SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.size()>0)
			NUM_EDIT_OPS_ = SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS[SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.size()-1] ;
	}
	
protected:
	int init_alignment_structures(Config * config);
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
	CHROMOSOME_ENTRY* alloc_chromosome_entry(unsigned int pos, Chromosome const &chr, char strand) ;
	MAPPING_ENTRY_CONTAINER* alloc_mapping_entry_container() ;
	MAPPING_ENTRY* alloc_mapping_entry() ;
	
	CHROMOSOME_ENTRY **GENOME;
	char HAS_SLOT;
	int SLOT;

	Genome * genome ;
	GenomeMaps * genomemaps ;

	std::vector<int> SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS ;
	
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

	unsigned int NUM_SCORE_INTERVALS;
	unsigned int MAX_USED_SLOTS;

protected:
	HIT **READSTART_BINS;


	MAPPING_ENTRY_CONTAINER* MAPPING_ENTRY_OPERATOR_FIRST;
	MAPPING_ENTRY_CONTAINER* MAPPING_ENTRY_OPERATOR;
	HIT_CONTAINER* HIT_OPERATOR_FIRST;
	HIT_CONTAINER* HIT_OPERATOR;
	CHROMOSOME_ENTRY_CONTAINER CHROMOSOME_ENTRY_OPERATOR;
	
	unsigned int NUM_MAPPING_ENTRIES;
};


