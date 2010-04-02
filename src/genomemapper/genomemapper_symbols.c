#if 1 // dd
char HAS_SLOT;
int SLOT;
#endif // dd
size_t INDEX_SIZE = Config::INDEX_SIZE_13 ;

// Genome
int BINARY_CODE[4];
unsigned int NUM_POS;
unsigned int LONGEST_HIT;

//config
int REPORT_REPETITIVE_SEED_DEPTH_EXTRA = 31 - MAX_INDEX_DEPTH ;

FILE *OUT_FP;
FILE *SP_OUT_FP;
FILE *TRIGGERED_LOG_FP; // #A#

POS *BLOCK_TABLE;
unsigned int BLOCK_TABLE_SIZE;
INDEX_ENTRY *INDEX;
INDEX_ENTRY *INDEX_REV;
#if 1 // dd
#ifndef BinaryStream_MAP
STORAGE_ENTRY *INDEX_REV_MMAP;
STORAGE_ENTRY *INDEX_FWD_MMAP;
#else
CBinaryStream<STORAGE_ENTRY>* INDEX_REV_MMAP;
CBinaryStream<STORAGE_ENTRY>* INDEX_FWD_MMAP;
#endif
#endif // dd
unsigned long int MAX_POSITIONS;
int REDUNDANT;
char FLANK_SEQ[Config::MAX_READ_LENGTH + 200];
char *chrseq;	// for debugging
char *ALIGNSEQ;
int NUM_MATCHES; // For new version of QPalma
double WORST_SCORE;
double WORST_MM_SCORE;

unsigned int NUM_SCORE_INTERVALS;
HIT **HIT_LISTS_OPERATOR;
HIT **READSTART_BINS;
HITS_BY_SCORE_STRUCT *HITS_BY_SCORE; ///< Collection of bins, each of which contains all hits with the same score for the current read.
unsigned int HITS_IN_SCORE_LIST;
CHROMOSOME_ENTRY **GENOME;
unsigned int LONGEST_CHROMOSOME;
// dd unsigned int CHROM_CONTAINER_SIZE;

//Hit resultat
MAPPING_ENTRY_CONTAINER* MAPPING_ENTRY_OPERATOR_FIRST;
MAPPING_ENTRY_CONTAINER* MAPPING_ENTRY_OPERATOR;
HIT_CONTAINER* HIT_OPERATOR_FIRST;
HIT_CONTAINER* HIT_OPERATOR;

CHROMOSOME_ENTRY_CONTAINER* CHROMOSOME_ENTRY_OPERATOR;
unsigned int MAX_USED_SLOTS;
unsigned int NUM_MAPPING_ENTRIES;

Read _read;

//QPalma
struct alignment_parameter_struct *alignment_parameters = NULL ;

