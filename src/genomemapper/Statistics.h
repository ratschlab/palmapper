#pragma once

#include <palmapper/Config.h>

class Statistics {
public:
	Statistics();

	unsigned int PERFECT_READS;
	unsigned int PERFECT_HITS;
	unsigned int PERFECT_HITS_REV;
	unsigned long int NUM_HITS;
	unsigned int HITS_LEN[Config::MAX_READ_LENGTH]; // TODO Refactoring needed. Should be a parameter.
	unsigned int HITS_MM[Config::MAX_EDIT_OPS + 1];
	unsigned int READS_MAPPED;
	unsigned long int NUM_ALIGNMENTS;
	unsigned long int NUM_WHOLE_ALIGNMENTS;
	unsigned int ENDSTART_MAPPED[2];
	unsigned int NOT_ALIGNED[2];
	unsigned int NUM_READS;
	unsigned int HITS_PER_READ;
	unsigned long int GAPS_ENCOUNTERED[3];
	unsigned long int TOO_MANY_MMS[2];
	unsigned long int BREAK_GLOBAL_ALIGNMENT[2];
	unsigned long int BREAK_TB_IN_GLOBAL_ALIGNMENT;
	unsigned long int CELLS_GLOBAL;
	unsigned long int CELLS_OVERHANG;
	unsigned long int W;
	unsigned int listcount, listocc;
};
