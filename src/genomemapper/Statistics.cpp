#include <genomemapper/Statistics.h>

Statistics::Statistics() 
{
	PERFECT_READS = 0;
	PERFECT_HITS = 0;
	PERFECT_HITS_REV = 0;
	NUM_HITS = 0;
	for (unsigned int i = 0; i < Config::MAX_READ_LENGTH; ++i)
		HITS_LEN[i] = 0;
	for (int i = 0; i != Config::MAX_EDIT_OPS + 1; ++i)
		HITS_MM[i] = 0;
	READS_MAPPED = 0;
	NUM_ALIGNMENTS = 0;
	NUM_WHOLE_ALIGNMENTS = 0;
	ENDSTART_MAPPED[0] = 0;
	ENDSTART_MAPPED[1] = 0;
	NOT_ALIGNED[0] = 0;
	NOT_ALIGNED[1] = 0;
	NUM_READS = 0;
	HITS_PER_READ = 0;
	GAPS_ENCOUNTERED[0] = 0;
	GAPS_ENCOUNTERED[1] = 0;
	GAPS_ENCOUNTERED[2] = 0;
	TOO_MANY_MMS[0] = 0;
	TOO_MANY_MMS[1] = 0;
	BREAK_GLOBAL_ALIGNMENT[0] = 0;
	BREAK_GLOBAL_ALIGNMENT[1] = 0;
	BREAK_TB_IN_GLOBAL_ALIGNMENT = 0;
	CELLS_GLOBAL = 0;
	CELLS_OVERHANG = 0;
	W = 0;
	listcount = 0;
	listocc = 0;
}
