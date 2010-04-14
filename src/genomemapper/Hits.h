#pragma once

#include <assert.h>
#include <vector>
#include <genomemapper/Config.h>
#include <genomemapper/Chromosome.h>
#include <genomemapper/Read.h>

#define CONTAINER_SIZE 100000
#define SCORE_INTERVAL 1

extern Config _config ;

extern char *get_seq(Read const &read, unsigned int n);


typedef struct edit_op_structure {
	signed int pos;
	int mm: 1;
} EDIT_OPS;

struct HIT {
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
};

template <bool clearNew, class Entry> class Container {
public:
	Container() {
		_current = _end = NULL;
	}
	~Container() {
		clear();
	}
	Entry *newEntry() {
		if (_current >= _end) {
			_entries.push_back(_current = new Entry[CONTAINER_SIZE]);
			_end = _current + CONTAINER_SIZE;
		}
		::memset(_current, 0, sizeof(Entry));
		return _current++;
	}
	void clear() {
		while (_entries.size() > 0) {
			delete[] _entries.back();
			_entries.pop_back();
		}
		_current = _end = NULL;
	}

private:
	std::vector<Entry*> _entries;
	Entry *_current;
	Entry *_end;
};

struct MAPPING_ENTRY {
	unsigned int readpos;
	HIT *hit;
	struct MAPPING_ENTRY *pred;
	struct MAPPING_ENTRY *succ;
};

struct CHROMOSOME_ENTRY {
	int chromosome; // @TODO? Minus values indicate reverse hits (would save strand var = 1 byte, ~15MB)
	unsigned int genome_pos;
	char strand;

	CHROMOSOME_ENTRY *next;
	MAPPING_ENTRY *mapping_entries;

	// It seems to be cheaper to store the back-pointer information (pos)
	// in each of these entries rather than having a superior structure.
};

struct CHROMOSOME_ENTRY_CONTAINER {
	CHROMOSOME_ENTRY_CONTAINER(int nrEntries) {
		entries = new CHROMOSOME_ENTRY[nrEntries];
		used = 0;
	}
	~CHROMOSOME_ENTRY_CONTAINER() {
		delete[] entries;
	}
	CHROMOSOME_ENTRY* entries;
	unsigned int used;
};

typedef struct hits_by_score_structure {
	HIT *hitpointer;
	int num;
} HITS_BY_SCORE_STRUCT;

class QPalma ;
class Genome ;
class GenomeMaps ;
class TopAlignments ;


class Hits {

	friend class ReadMappings;

public:
	
	Hits(Genome &genome, GenomeMaps &genomemaps, Read &read);
	~Hits();

	int map_reads(Genome &genome, GenomeMaps &genomeMaps, TopAlignments* topalignments, QPalma* qpalma) ;
	int REDUNDANT;

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
	int init_constants()  ;
	int init_statistic_vars() ;
	int init_operators() ;
	int init_alignment_structures(Config * config);
	CHROMOSOME_ENTRY* alloc_chromosome_entry(unsigned int pos, Chromosome const &chr, char strand) ;
	CHROMOSOME_ENTRY **GENOME; // doppelt

	Genome &_genome;
	GenomeMaps &_genomeMaps ;

	std::vector<int> SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS ;
	
	unsigned int LONGEST_HIT;

public:
	static unsigned int MAX_USED_SLOTS;
	CHROMOSOME_ENTRY_CONTAINER CHROMOSOME_ENTRY_OPERATOR;

protected:
	Read &_read;
};

class ReadMappings {
public:
	ReadMappings(Genome &genome, GenomeMaps &genomeMaps, Hits &hits, Read const &read);
	~ReadMappings();
	int dealloc_hit_lists_operator() ;
	int dealloc_hits_by_score() ;
	int align_hit_simple(HIT* hit, int start, int end, int readpos, Chromosome const &chromosome, int orientation, unsigned char mismatches) ;
	int prepare_kbound_alignment(HIT* hit, int start, int end, int readpos, Chromosome const &chromosome, char orientation, char mismatches) ;
	int size_hit(HIT *hit, unsigned int *oldlength, char num) ;
	int seed2genome(unsigned int num, unsigned int index_slot, unsigned int readpos, char reverse) ;
	void printgenome() ;
//	int alloc_genome_memory() ;
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

	HIT **HIT_LISTS_OPERATOR;
	unsigned int HITS_IN_SCORE_LIST;
	void dealloc_mapping_entries() {_mappings.clear();}
	void dealloc_hits() {_hits.clear();}
private:
	int alloc_hits_by_score() ;
	int alloc_hit_lists_operator() ;
	CHROMOSOME_ENTRY* alloc_chromosome_entry(unsigned int pos, Chromosome const &chr, char strand) {
		return _outer.alloc_chromosome_entry(pos, chr, strand);
	}

	Genome &_genome;
	GenomeMaps &_genomeMaps;
	HITS_BY_SCORE_STRUCT *HITS_BY_SCORE;
	unsigned int NUM_SCORE_INTERVALS;
	Container<false, MAPPING_ENTRY> _mappings;
	Container<true, HIT> _hits;
	CHROMOSOME_ENTRY **GENOME; // doppelt
	Hits &_outer;
	Read const &_read;
	CHROMOSOME_ENTRY_CONTAINER &CHROMOSOME_ENTRY_OPERATOR;

	static char HAS_SLOT;
	static int SLOT;
};
