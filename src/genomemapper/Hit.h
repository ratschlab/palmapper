#pragma once

#include <genomemapper/Config.h>
#include <genomemapper/Chromosome.h>

typedef struct edit_op_structure {
	signed int pos;
	int mm: 1;
} EDIT_OPS;

class HIT {
public:
	HIT();

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
