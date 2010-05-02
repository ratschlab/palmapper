#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <palmapper/Config.h>
#include <palmapper/Read.h>
#include <palmapper/Genome.h>

Config::Config() {
	NUM_THREADS = 4;
	OUTPUT_FILTER = OUTPUT_FILTER_DEFAULT ;
	OUTPUT_FILTER_NUM_TOP = 10 ;
	OUTPUT_FILTER_NUM_LIMIT = 0 ; // all
	RTRIM_STRATEGY=0 ;
	RTRIM_STRATEGY_MIN_LEN = 25 ;
	RTRIM_STRATEGY_STEP = DEFAULT_SETTING ;
	POLYTRIM_STRATEGY=0 ;
	POLYTRIM_STRATEGY_MIN_LEN=25 ;
	POLYTRIM_STRATEGY_STEP = DEFAULT_SETTING ;
	POLYTRIM_STRATEGY_POLY_MIN_LEN = 10 ;
	ADAPTERTRIM_STRATEGY = 0 ;
	ADAPTERTRIM_STRATEGY_MIN_LEN = 40 ;

	//int SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS[2] ;
	HITLEN_LIMIT = 0;
	VERBOSE = 0;
	MAP_REVERSE = 1 ;
	STRINGENT_GAPLIMIT = 1 ;
	PRINT_SEQ = 0;
	INDEX_DEPTH = 0;
	INDEX_DEPTH_EXTRA = 3 ;
	INDEX_DEPTH_EXTRA_THRESHOLD = 100000000 ;
	SEED_HIT_CANCEL_THRESHOLD = 100000000 ;
	OUTPUT_FORMAT = OUTPUT_FORMAT_DEFAULT ;
	REPORT_FILE = NULL;
	REPORT_FILE_READONLY = 0 ;
	REPORT_REPETITIVE_SEEDS = 0 ;
	REPORT_MAPPED_REGIONS = 0 ;
	REPORT_MAPPED_READS = 0 ;
	REPORT_SPLICED_READS = 0 ;
	REPORT_RESET = 0 ;
	QPALMA_USE_MAP = 1 ;
	QPALMA_USE_MAP_MAX_SIZE = 10000 ;
	QPALMA_USE_SPLICE_SITES = 0 ;
	QPALMA_USE_SPLICE_SITES_THRESH_DON = 0.0 ;
	QPALMA_USE_SPLICE_SITES_THRESH_ACC = 0.0 ;
	QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC = 0.0 ;
	QPALMA_MIN_NUM_MATCHES = 3 ;
	READ_COUNT_LIMIT = 0 ; // limits the number of reads for alignment
	LOG_TRIGGERED = false;  // #A#
	FILTER_BY_MAX_MISMATCHES = 1 ;
	FILTER_BY_MAX_GAPS = 0 ;
	FILTER_BY_SPLICE_SITES = 5 ;
	FILTER_BY_SPLICE_SITES_REGION = 5 ;
	FILTER_BY_SPLICE_SITES_EDIT_MIN = 1 ;
	FILTER_BY_SPLICE_SITES_THRESH_ACC=0.8 ;
	FILTER_BY_SPLICE_SITES_THRESH_DON=0.8 ;
	FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC = 0.0;
	NO_SPLICE_PREDICTIONS=0 ;
	INDEX_PRECACHE = 0 ;
	FLANKING = 0;
	NUM_EDIT_OPS = DEFAULT_SETTING ;
	NUM_MISMATCHES = DEFAULT_SETTING ;
	NUM_GAPS = DEFAULT_SETTING ;
	MM_SCORE = 4;
	M_SCORE = 0;
	GAP_SCORE = 5;
	GAPS_MOST_RIGHT = 0;
	OVERHANG_ALIGNMENT = 0;
	SCORES_OUT = 1;
	SPLICED_HITS = 0 ;
	SPLICED_HIT_MIN_LENGTH_SHORT = DEFAULT_SETTING ;
	SPLICED_HIT_MIN_LENGTH_COMB = DEFAULT_SETTING ;
	SPLICED_HIT_MIN_LENGTH_LONG = DEFAULT_SETTING ;
	SPLICED_LONGEST_INTRON_LENGTH = DEFAULT_SETTING ;
	SPLICED_MAX_NUM_ALIGNMENTS = 10 ;
	SPLICED_CLUSTER_TOLERANCE = 10 ;
	SPLICED_MAX_INTRONS = DEFAULT_SETTING ;
	STATISTICS = 0;
	CHROM_CONTAINER_SIZE = 15000000 ;

	ALL_HIT_STRATEGY = 0 ;
	SUMMARY_HIT_STRATEGY= 0 ;
	HITLEN_LIMIT = 0 ;

};

int Config::applyDefaults(Genome * genome)
{
	{
		Read read ;
		int read_length = read.determine_read_length(QUERY_FILE_NAME) ;
		if ((SPLICED_HITS && (SPLICED_HIT_MIN_LENGTH_SHORT == DEFAULT_SETTING || SPLICED_HIT_MIN_LENGTH_LONG == DEFAULT_SETTING || SPLICED_HIT_MIN_LENGTH_COMB == DEFAULT_SETTING || SPLICED_MAX_INTRONS == DEFAULT_SETTING)) || 
			NUM_EDIT_OPS == DEFAULT_SETTING || NUM_MISMATCHES == DEFAULT_SETTING || NUM_GAPS == DEFAULT_SETTING)
			fprintf(stdout, "Automatically determining alignment parameters based on read length (%int):", read_length) ;
		if (SPLICED_HITS && SPLICED_HIT_MIN_LENGTH_SHORT == DEFAULT_SETTING)
		{
			SPLICED_HIT_MIN_LENGTH_SHORT = 15 ;
			fprintf(stdout, " -K %i", SPLICED_HIT_MIN_LENGTH_SHORT) ;
		}
		if (SPLICED_HITS && SPLICED_HIT_MIN_LENGTH_LONG == DEFAULT_SETTING)
		{
			SPLICED_HIT_MIN_LENGTH_LONG = (read_length/4<20) ? 20 : read_length/4 ;
			fprintf(stdout, " -L %i", SPLICED_HIT_MIN_LENGTH_LONG) ;
		}
		if (SPLICED_HITS && SPLICED_HIT_MIN_LENGTH_COMB == DEFAULT_SETTING)
		{
			SPLICED_HIT_MIN_LENGTH_COMB = (read_length/2<30) ? 30 : read_length/2 ;
			fprintf(stdout, " -C %i", SPLICED_HIT_MIN_LENGTH_COMB) ;
		}
		if (SPLICED_HITS && SPLICED_MAX_INTRONS == DEFAULT_SETTING)
		{
			SPLICED_MAX_INTRONS = (read_length>50) ? (read_length>=100 ? 3 : 2) : 1 ;
			if (read_length>200) SPLICED_MAX_INTRONS = 6 ;
			fprintf(stdout, " -NI %i", SPLICED_MAX_INTRONS) ;
		}
		if (NUM_EDIT_OPS == DEFAULT_SETTING)
		{
			NUM_EDIT_OPS = read_length*0.08 ;
			fprintf(stdout, " -E %i", NUM_EDIT_OPS) ;
		}
		if (NUM_MISMATCHES == DEFAULT_SETTING)
		{
			NUM_MISMATCHES = read_length*0.08 ;
			fprintf(stdout, " -M %i", NUM_MISMATCHES) ;
		}
		if (NUM_GAPS == DEFAULT_SETTING)
		{
			NUM_GAPS = read_length*0.03 ;
			fprintf(stdout, " -G %i", NUM_GAPS) ;
		}
		fprintf(stdout, "\n") ;
	}
	
	if (SPLICED_HITS && SPLICED_LONGEST_INTRON_LENGTH == DEFAULT_SETTING)
	{
		unsigned long int genome_size = 0 ;
		for (unsigned int i=0; i<genome->nrChromosomes(); i++)
			genome_size+=genome->chromosome(i).length() ;

		if (genome_size<900000000)
			SPLICED_LONGEST_INTRON_LENGTH = 10000 ;
		else if (genome_size<500000000)
			SPLICED_LONGEST_INTRON_LENGTH = 50000 ;
		else
			SPLICED_LONGEST_INTRON_LENGTH = 200000 ;
		fprintf(stdout, "Automatically determined maximal intron size based on genome size (%ikb)\n", SPLICED_LONGEST_INTRON_LENGTH/1000) ;
	}
	
	// determine default output format
	if (OUTPUT_FORMAT==OUTPUT_FORMAT_DEFAULT)
		if (SPLICED_HITS)
		{
			OUTPUT_FORMAT=OUTPUT_FORMAT_BEDX ;
			fprintf(stdout, "Selecting BEDX output format\n") ;
		}
		else
		{
			OUTPUT_FORMAT=OUTPUT_FORMAT_SHORE ;
			fprintf(stdout, "Selecting SHORE output format\n") ;
		}

	if (POLYTRIM_STRATEGY)
	{
		if ((int)POLYTRIM_STRATEGY_STEP == DEFAULT_SETTING)
		{
			POLYTRIM_STRATEGY_STEP = (NUM_EDIT_OPS>=1) ? NUM_EDIT_OPS : 1 ;
			fprintf(stdout, "Automatically selecting polytrim step size: %int\n", POLYTRIM_STRATEGY_STEP) ;
		}
	}

	if (RTRIM_STRATEGY)
	{
		if ((int)RTRIM_STRATEGY_STEP == DEFAULT_SETTING)
		{
			RTRIM_STRATEGY_STEP = (NUM_EDIT_OPS>=2) ? NUM_EDIT_OPS/2 : 1 ;
			fprintf(stdout, "Automatically selecting rtrim step size: %int\n", RTRIM_STRATEGY_STEP) ;
		}
	}

	return 0 ;
}

int Config::checkConfig()
{
	if (SPLICED_OUT_FILE_NAME.length()>0 && !SPLICED_HITS)
	{
		fprintf(stderr, "ERROR: output files for spliced hits provided, but no spliced alignment is performed\n");
		exit(1);
	}

	if (SPLICED_HITS && !(NO_SPLICE_PREDICTIONS || (ACC_FILES.length()>0 && DON_FILES.length()>0)))
	{
		fprintf(stderr, "ERROR: for spliced alignments either -acc and -don or -no-ss-pred need to be given as argument\n");
		exit(1);
	}

	if (NO_SPLICE_PREDICTIONS && (ACC_FILES.length()>0 || DON_FILES.length()>0))
	{
		fprintf(stderr, "ERROR: the options -acc/-don and -no-ss-pred have to be used exclusively\n");
		exit(1);
	}

	if (RTRIM_STRATEGY && POLYTRIM_STRATEGY)
	{
		fprintf(stderr, "ERROR: RTRIM and POLYTRIM cannot be combined\n") ;
		exit(1) ;
	}

	if (SPLICED_HITS && (OUTPUT_FORMAT==OUTPUT_FORMAT_SHORE || OUTPUT_FORMAT==OUTPUT_FORMAT_BED))
	{
		fprintf(stderr, "ERROR: SHORE or BED format currently do not support spliced alignments (choose BEDX or SAM) %i\n", (int)OUTPUT_FORMAT) ;
		exit(1) ;
	}

/*	if (OUTPUT_FORMAT==OUTPUT_FORMAT_SAM)
	{
		fprintf(stderr, "ERROR: SAM format not implemented yet\n") ;
		exit(1) ;
	}*/
	return 0 ;
}

int Config::parseCommandLine(int argc, char *argv[]) 
{
	int i;
	char not_defined;
	char has_index = 0;
	char has_query = 0;
	char has_genome = 0;
	char output[5];

	for (i = 1; i < argc; i++) {
		not_defined = 1;

		//genome file
		if (strcmp(argv[i], "-i") == 0) {

			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -i\n") ;
				usage();
				exit(1);
			}
			i++;
			GENOME_FILE_NAME.assign(argv[i]);
			has_genome = 1;
			has_index = 1;

			CHR_INDEX_FILE_NAME.assign(argv[i]) ;
			CHR_INDEX_FILE_NAME += ".cid";

			META_INDEX_FILE_NAME.assign(argv[i]);
			META_INDEX_FILE_NAME += ".mta";

			INDEX_FWD_FILE_NAME.assign(argv[i]);
			INDEX_FWD_FILE_NAME += ".mfd";

			INDEX_REV_FILE_NAME.assign(argv[i]);
			INDEX_REV_FILE_NAME += ".mrc";

		}

		/**

		 //chr index file
		 if(strcmp(argv[i],"-x")==0){
		 not_defined = 0;
		 if(i+1 > argc - 1) 
		 { 
		 fprintf(stderr, "ERROR: Argument missing for option -i\n") ;
		 usage(); exit(1); 
		 }
		 i++;
		 strcpy(CHR_INDEX_FILE_NAME, argv[i]);
		 has_index = 1;
		 }

		 //meta index file
		 if(strcmp(argv[i],"-t")==0){
		 not_defined = 0;
		 if(i+1 > argc - 1) { usage(); exit(1); }
		 i++;
		 strcpy(META_INDEX_FILE_NAME, argv[i]);
		 has_meta_index = 1;
		 }

		 //index fwd file
		 if(strcmp(argv[i],"-z")==0){
		 not_defined = 0;
		 if(i+1 > argc - 1) { usage(); exit(1); }
		 i++;
		 strcpy(INDEX_FWD_FILE_NAME, argv[i]);
		 has_fwd_index = 1;
		 }

		 //index rev file
		 if(strcmp(argv[i],"-y")==0){
		 not_defined = 0;
		 if(i+1 > argc - 1) { usage(); exit(1); }
		 i++;
		 strcpy(INDEX_REV_FILE_NAME, argv[i]);
		 has_rev_index = 1;
		 }

		 */

		 //meta index file
		if(strcmp(argv[i],"-threads")==0)
		{
			not_defined = 0;
			if(i+1 > argc - 1) 
			{ 
				fprintf(stderr, "ERROR: Argument missing for option -threads\n") ;
				usage(); 
				exit(1); 
			}
			i++;
		    NUM_THREADS = atoi(argv[i]) ;
			if (NUM_THREADS<1)
			{
				fprintf(stderr,	"ERROR: number of threads too small\n");
				exit(1) ;
			}
		}

		//query file
		if (strcmp(argv[i], "-q") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -q\n") ;
				usage();
				exit(1);
			}
			i++;
			QUERY_FILE_NAME.assign(argv[i]);
			has_query = 1;
		}

		//output file
		if (strcmp(argv[i], "-o") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -o\n") ;
				usage();
				exit(1);
			}
			i++;
			OUT_FILE_NAME.assign(argv[i]);
		}

		//partial alignments for rtrim
		if (strcmp(argv[i], "-rtrim") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -rtrim\n") ;
				usage();
				exit(1);
			}
			i++;
			RTRIM_STRATEGY = 1 ;
			RTRIM_STRATEGY_MIN_LEN = atoi(argv[i]) ;
			if (RTRIM_STRATEGY_MIN_LEN<INDEX_DEPTH)
			{
				fprintf(stderr,	"ERROR: minimal rtrim alignment length too short\n");
				exit(1) ;
			}
		}

		//partial alignments for rtrim
		if (strcmp(argv[i], "-rtrim-step") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -rtrim-step\n") ;
				usage();
				exit(1);
			}
			i++;
			RTRIM_STRATEGY = 1 ;
			RTRIM_STRATEGY_STEP = atoi(argv[i]) ;
			if (RTRIM_STRATEGY_STEP<1)
			{
				fprintf(stderr,	"ERROR: rtrim step size too short\n");
				exit(1) ;
			}
		}

		//polyA/T trimming
		if (strcmp(argv[i], "-polytrim") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -polytrim\n") ;
				usage();
				exit(1);
			}
			i++;
			POLYTRIM_STRATEGY = 1 ;
			POLYTRIM_STRATEGY_MIN_LEN = atoi(argv[i]) ;
			if (POLYTRIM_STRATEGY_MIN_LEN<INDEX_DEPTH)
			{
				fprintf(stderr,	"ERROR: minimal polytrim alignment length too short\n");
				exit(1) ;
			}
		}

		//adapter trimming
		if (strcmp(argv[i], "-adaptertrim") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -polytrim\n") ;
				usage();
				exit(1);
			}
			i++;
			ADAPTERTRIM_STRATEGY = 1 ;
			ADAPTERTRIM_STRATEGY_MIN_LEN = atoi(argv[i]) ;
			if (ADAPTERTRIM_STRATEGY_MIN_LEN<INDEX_DEPTH)
			{
				fprintf(stderr,	"ERROR: minimal polytrim alignment length too short\n");
				exit(1) ;
			}
		}

		//report output file
		if (strcmp(argv[i], "-report") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -report\n") ;
				usage();
				exit(1);
			}
			i++;
			REPORT_FILE=strdup(argv[i]) ;
			REPORT_FILE_READONLY = 0 ;
		}

		//report output file
		if (strcmp(argv[i], "-report-ro") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -report-ro\n") ;
				usage();
				exit(1);
			}
			i++;
			REPORT_FILE=strdup(argv[i]) ;
			REPORT_FILE_READONLY = 1 ;
		}

		//report repetitive seeds
		if (strcmp(argv[i], "-report-rep-seed") == 0) {
			not_defined = 0;
			REPORT_REPETITIVE_SEEDS = 1 ;
			//assert(REPORT_SPLICE_SITES==0) ; // currently not supported
		}

		//report mapped regions
		if (strcmp(argv[i], "-report-map-region") == 0) {
			not_defined = 0;
			REPORT_MAPPED_REGIONS  = 1 ;
		}

		//report mapped regions
		if (strcmp(argv[i], "-report-map-read") == 0) {
			not_defined = 0;
			REPORT_MAPPED_READS = 1 ;
		}

		//report mapped regions
		if (strcmp(argv[i], "-report-spliced-read") == 0) {
			not_defined = 0;
			REPORT_SPLICED_READS = 1 ;
		}

		//report splice sites - confidence threshold
		if (strcmp(argv[i], "-report-splice-sites") == 0) {
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -report-splice-sites\n") ;
				usage();
				exit(1);
			}
			i++;
			QPALMA_USE_SPLICE_SITES_THRESH_ACC = atof(argv[i]);
			QPALMA_USE_SPLICE_SITES_THRESH_DON = atof(argv[i]);
			QPALMA_USE_SPLICE_SITES= 1 ;
			not_defined = 0;
			//assert(REPORT_REPETITIVE_SEEDS==0) ; // currently not supported
			//assert(REPORT_SPLICE_SITES_THRESH_TOP_PERC==0.0) ;
		}

		//report splice sites - percentile threshold
		if (strcmp(argv[i], "-report-splice-sites-top-perc") == 0) {
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -report-splice-sites-top-perc\n") ;
				usage();
				exit(1);
			}
			i++;
			QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC = atof(argv[i]);
			assert(QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC>=0 && QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC<=1.0) ;
			QPALMA_USE_SPLICE_SITES= 1 ;
			not_defined = 0;
			//assert(REPORT_REPETITIVE_SEEDS==0) ; // currently not supported
			//assert(REPORT_SPLICE_SITES_THRESH==0.0) ;
		}

		// filter by splice sites - percentile threshold
		if (strcmp(argv[i], "-filter-splice-sites-top-perc") == 0) {
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -filter-splice-sites-top-perc\n") ;
				usage();
				exit(1);
			}
			i++;
			FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC = atof(argv[i]);
			assert(FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC>=0 && FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC<=1.0) ;
			not_defined = 0;
			FILTER_BY_SPLICE_SITES_THRESH_ACC = 0.0 ;
			FILTER_BY_SPLICE_SITES_THRESH_DON = 0.0 ;
		}

		// filter by mismatches
		if (strcmp(argv[i], "-filter-max-mismatches") == 0) {
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -filter-max-mismatches\n") ;
				usage();
				exit(1);
			}
			i++;
			FILTER_BY_MAX_MISMATCHES = atoi(argv[i]);
			assert(FILTER_BY_MAX_MISMATCHES>=0) ;
			not_defined = 0;
		}

		// splice-site based filter: require at least this many edits
		if (strcmp(argv[i], "-filter-splice-min-edit") == 0) {
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -filter-splice-min-edit\n") ;
				usage();
				exit(1);
			}
			i++;
			FILTER_BY_SPLICE_SITES_EDIT_MIN = atoi(argv[i]);
			assert(FILTER_BY_SPLICE_SITES_EDIT_MIN>=0) ;
			not_defined = 0;
		}

		// splice-site based filter: consider a region of this length around the read; -1 switches off the splice site-based filter
		if (strcmp(argv[i], "-filter-splice-region") == 0)
		{
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -filter-splice-region\n") ;
				usage();
				exit(1);
			}
			i++;
			FILTER_BY_SPLICE_SITES_REGION = atoi(argv[i]);
			if (FILTER_BY_SPLICE_SITES_REGION==-1)
				FILTER_BY_SPLICE_SITES = 0 ;

			assert(FILTER_BY_SPLICE_SITES_REGION>=-1) ;
			not_defined = 0;
		}

		// filter by mismatches
		if (strcmp(argv[i], "-filter-max-gaps") == 0) {
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -filter-max-gaps\n") ;
				usage();
				exit(1);
			}
			i++;
			FILTER_BY_MAX_GAPS = atoi(argv[i]);
			assert(FILTER_BY_MAX_GAPS>=0) ;
			not_defined = 0;
		}

		// reset the report (i.e. don't load it from file, even when available)
		if (strcmp(argv[i], "-report-reset") == 0) {
			not_defined = 0;
			REPORT_RESET = 1 ;
		}

		// use regions around mapped reads for qpalma alignment
		/*if (strcmp(argv[i], "-qpalma-use-map") == 0) {
			not_defined = 0;
			QPALMA_USE_MAP = 1 ;
			}*/

		// use regions around mapped reads for qpalma alignment
		if (strcmp(argv[i], "-qpalma-use-map-max-len") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -qpalma-use-map-max-len\n") ;
				usage();
				exit(1);
			}
			i++;
			QPALMA_USE_MAP_MAX_SIZE = atoi(argv[i]);
		}

		//spliced output file
		if (strcmp(argv[i], "-report-gff-init") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -report-gff-init\n") ;
				usage();
				exit(1);
			}
			i++;
			REPORT_GFF_FILE_NAME.assign(argv[i]);
		}

		//spliced output file
		if (strcmp(argv[i], "-H") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -H\n") ;
				usage();
				exit(1);
			}
			i++;
			SPLICED_OUT_FILE_NAME.assign(argv[i]);
		}

		if (strcmp(argv[i], "-read-id-prefix") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -read-id-prefix\n") ;
				usage();
				exit(1);
			}
			i++;
			READ_ID_PREFIX.assign(argv[i]);
		}

		//spliced hits min combined length
		if (strcmp(argv[i], "-C") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -C\n") ;
				usage();
				exit(1);
			}
			i++;
			SPLICED_HIT_MIN_LENGTH_COMB = atoi(argv[i]);
		}

		//spliced hits min length
		if (strcmp(argv[i], "-K") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -K\n") ;
				usage();
				exit(1);
			}
			i++;
			SPLICED_HIT_MIN_LENGTH_SHORT = atoi(argv[i]);
		}

		//spliced hits min length of longer hit
		if (strcmp(argv[i], "-L") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -L\n") ;
				usage();
				exit(1);
			}
			i++;
			SPLICED_HIT_MIN_LENGTH_LONG = atoi(argv[i]);
		}

		// longest intron length
		if (strcmp(argv[i], "-I") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -I\n") ;
				usage();
				exit(1);
			}
			i++;
			SPLICED_LONGEST_INTRON_LENGTH = atoi(argv[i]);
		}

		// maximal number of spliced alignments to be performed per read
		if (strcmp(argv[i], "-SA") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -SA\n") ;
				usage();
				exit(1);
			}
			i++;
			int tmp = atoi(argv[i]);
			if (tmp < 0) {
				fprintf(stderr, "ERROR: Argument for option -SA too small\n") ;
				usage();
				exit(1);
			}
			SPLICED_MAX_NUM_ALIGNMENTS = tmp;
		}

		// maximal number of introns in spliced alignments
		if (strcmp(argv[i], "-NI") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -NI\n") ;
				usage();
				exit(1);
			}
			i++;
			int tmp = atoi(argv[i]);
			if (tmp < 0) {
				fprintf(stderr, "ERROR: Argument for option -NI too small\n") ;
				usage();
				exit(1);
			}
			SPLICED_MAX_INTRONS = tmp;
		}

		// How many matches are necessary for identifying a possible splice site in QPALMA recursive alignment algorithm?
		if (strcmp(argv[i], "-QMM") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -QMM\n") ;
				usage();
				exit(1);
			}
			i++;
			int tmp = atoi(argv[i]);
			if (tmp < 0) {
				fprintf(stderr, "ERROR: Argument for option -QMM too small\n") ;
				usage();
				exit(1);
			}
			QPALMA_MIN_NUM_MATCHES = tmp;
		}

		// how much distance to tolerate between a hit and an existing
		// hit cluster
		if (strcmp(argv[i], "-CT") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -CT\n") ;
				usage();
				exit(1);
			}
			i++;
			int tmp = atoi(argv[i]);
			if (tmp < 0) {
				fprintf(stderr, "ERROR: Argument for option -CT too small\n") ;
				usage();
				exit(1);
			}
			SPLICED_CLUSTER_TOLERANCE = tmp;
		}

		// limit the number of reads for alignment
		if (strcmp(argv[i], "-rlim") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -rlim\n") ;
				usage();
				exit(1);
			}
			i++;
			int tmp = atoi(argv[i]);
			if (tmp < 0) {
				fprintf(stderr, "ERROR: Argument for option -rlim too small\n") ;
				usage();
				exit(1);
			}
			READ_COUNT_LIMIT = tmp;
		}

		// limit the number of reads for alignment
		if (strcmp(argv[i], "-index-precache") == 0) {
			not_defined = 0;
			INDEX_PRECACHE = 1;
		}

		//output format
		if (strcmp(argv[i], "-f") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -f\n") ;
				usage();
				exit(1);
			}
			i++;
			//printf("argv .%s.\n", argv[i]);
			strcpy(output, argv[i]);
			if (strcmp(output, "shore") == 0 || strcmp(output, "SHORE") == 0) {
				OUTPUT_FORMAT = OUTPUT_FORMAT_SHORE ;
			} else if (strcmp(output, "bed") == 0 || strcmp(output, "BED") == 0) {
				OUTPUT_FORMAT = OUTPUT_FORMAT_BED;
			} else if (strcmp(output, "bedx") == 0 || strcmp(output, "BEDX") == 0) {
				OUTPUT_FORMAT = OUTPUT_FORMAT_BEDX;
			} else if (strcmp(output, "sam") == 0 || strcmp(output, "SAM") == 0) {
				OUTPUT_FORMAT = OUTPUT_FORMAT_SAM;
			} else {
				fprintf(stderr,
						"ERROR: Output file format must either be \"shore\" or \"bed\"\n");
				exit(1);
			}
		}

		//leftover file
		if (strcmp(argv[i], "-u") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -u\n") ;
				usage();
				exit(1);
			}
			i++;
			LEFTOVER_FILE_NAME.assign(argv[i]);
		}

		//verbose
		if (strcmp(argv[i], "-v") == 0) {
			not_defined = 0;
			VERBOSE = 1;
		}

		/*//scores out
		if (strcmp(argv[i], "-e") == 0) {
			not_defined = 0;
			SCORES_OUT = 0;
			}*/

		// do not align using reverse index
		if (strcmp(argv[i], "-r") == 0) {
			not_defined = 0;
			MAP_REVERSE = 0;
		}

		//spliced hits analysis
		if (strcmp(argv[i], "-S") == 0) {
			not_defined = 0;
			SPLICED_HITS = 1;
		}

		// extend the seed-length if too many seed-matches were found
		if (strcmp(argv[i], "-index-extend") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -index-extend\n") ;
				usage();
				exit(1);
			}
			i++;
			int tmp = atoi(argv[i]);
			if (tmp < 0) {
				fprintf(stderr, "ERROR: Argument for option -index-extend too small\n") ;
				usage();
				exit(1);
			}
			INDEX_DEPTH_EXTRA = tmp ;
		}

		// extend the seed-length if too many seed-matches were found
		if (strcmp(argv[i], "-index-extend-threshold") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -index-extend-threshold\n") ;
				usage();
				exit(1);
			}
			i++;
			int tmp = atoi(argv[i]);
			if (tmp < 0) {
				fprintf(stderr, "ERROR: Argument for option -index-extend-threshold too small\n") ;
				usage();
				exit(1);
			}
			INDEX_DEPTH_EXTRA_THRESHOLD = tmp ;
		}

		// cancel seed processing if more seed matches exist than a threshold
		if (strcmp(argv[i], "-seed-hit-cancel-threshold") == 0)
		{
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -seed-hit-cancel-threshold\n") ;
				usage();
				exit(1);
			}
			i++;
			int tmp = atoi(argv[i]);
			if (tmp < 0) {
				fprintf(stderr, "ERROR: Argument for option -seed-hit-cancel-threshold too small\n") ;
				usage();
				exit(1);
			}
			SEED_HIT_CANCEL_THRESHOLD = tmp ;
		}

		//report all hits-strategy
		if (strcmp(argv[i], "-a") == 0) {
			not_defined = 0;
			if (OUTPUT_FILTER!=OUTPUT_FILTER_DEFAULT)
			{
				fprintf(stderr, "ERROR: output filter already defined\n") ;
				exit(1) ;
			}
			OUTPUT_FILTER=OUTPUT_FILTER_ALL ;
			ALL_HIT_STRATEGY = 1 ;
		}

		if (strcmp(argv[i], "-z") == 0) {
			not_defined = 0;
			if (OUTPUT_FILTER!=OUTPUT_FILTER_DEFAULT)
			{
				fprintf(stderr, "ERROR: output filter already defined\n") ;
				exit(1) ;
			}
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -z\n") ;
				usage();
				exit(1);
			}
			i++;
			
			if ((OUTPUT_FILTER_NUM_TOP = atoi(argv[i])) == 0) {
				if (argv[i][0] != '0' || OUTPUT_FILTER_NUM_TOP<0) 
				{
					fprintf(stderr,
							"ERROR: Number of alignments must be a positive integer value!\n");
					exit(1);
				}
			}
			OUTPUT_FILTER=OUTPUT_FILTER_TOP ;
			SUMMARY_HIT_STRATEGY=1 ;
		}

		//max nr of alignments per read
		if (strcmp(argv[i], "-ar") == 0) {
			not_defined = 0;
			if (OUTPUT_FILTER!=OUTPUT_FILTER_DEFAULT)
			{
				fprintf(stderr, "ERROR: output filter already defined\n") ;
				exit(1) ;
			}
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -ar\n") ;
				usage();
				exit(1);
			}
			i++;
			if (OUTPUT_FILTER_NUM_LIMIT != 0) {
				fprintf(stderr,
						"ERROR: options -a and -ar exclude themselves!\n");
				exit(1);
			}
			if ((OUTPUT_FILTER_NUM_LIMIT = atoi(argv[i])) == 0) {
				if (argv[i][0] != '0') {
					fprintf(stderr,
							"ERROR: Number of alignments must be an integer value!\n");
					exit(1);
				}
			}
			OUTPUT_FILTER = OUTPUT_FILTER_LIMIT ;
		}

		//max nr of alignments per read, randomly chosen!
		if (strcmp(argv[i], "-n") == 0) {
			not_defined = 0;
			if (OUTPUT_FILTER!=OUTPUT_FILTER_DEFAULT)
			{
				fprintf(stderr, "ERROR: output filter already defined\n") ;
				exit(1) ;
			}
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -n\n") ;
				usage();
				exit(1);
			}
			i++;
			if (OUTPUT_FILTER_NUM_LIMIT != 0) {
				fprintf(stderr,
						"ERROR: options -a and -ar exclude themselves!\n");
				exit(1);
			}
			if ((OUTPUT_FILTER_NUM_LIMIT = atoi(argv[i])) == 0) {
				if (argv[i][0] != '0') {
					fprintf(stderr,
							"ERROR: Number of alignments must be an integer value!\n");
					exit(1);
				}
			}
			OUTPUT_FILTER_NUM_LIMIT = -OUTPUT_FILTER_NUM_LIMIT ;
			OUTPUT_FILTER = OUTPUT_FILTER_RANDOM ;
		}

		//max number of allowed edit operations
		if (strcmp(argv[i], "-E") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -E\n") ;
				usage();
				exit(1);
			}
			i++;
			if ((NUM_EDIT_OPS = atoi(argv[i])) == 0) {
				if (argv[i][0] != '0') {
					fprintf(stderr,
							"ERROR: Number of edit operations must be an integer value!\n");
					exit(1);
				}
			}
			if (NUM_EDIT_OPS > MAX_EDIT_OPS) {
				fprintf(
						stderr,
						"ERROR: Number of allowed mismatches exceeds maximal number of edit operations (=%d)! Please restart with a lower value!\n",
						MAX_EDIT_OPS);
				exit(1);
			}
		}

		//max number of allowed mismatches
		if (strcmp(argv[i], "-M") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -M\n") ;
				usage();
				exit(1);
			}
			i++;
			if ((NUM_MISMATCHES = atoi(argv[i])) == 0) {
				if (argv[i][0] != '0') {
					fprintf(stderr,
							"ERROR: Number of mismatches must be an integer value!\n");
					exit(1);
				}
			}
			if (NUM_MISMATCHES > MAX_EDIT_OPS) {
				fprintf(
						stderr,
						"ERROR: Number of allowed mismatches exceeds maximal number of edit operations (=%d)! Please restart with a lower value!\n",
						MAX_EDIT_OPS);
				exit(1);
			}
		}

		//max number of allowed gaps
		if (strcmp(argv[i], "-G") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -G\n") ;
				usage();
				exit(1);
			}
			i++;
			if ((NUM_GAPS = atoi(argv[i])) == 0) {
				if (argv[i][0] != '0') {
					fprintf(stderr,
							"ERROR: Number of gaps must be an integer value!\n");
					exit(1);
				}
			}
			if (NUM_GAPS > MAX_EDIT_OPS) {
				fprintf(
						stderr,
						"ERROR: Number of allowed gaps exceeds maximal number of edit operations (=%d)! Please restart with a lower value!\n",
						MAX_EDIT_OPS);
				exit(1);
			}
		}

		//match score
		if (strcmp(argv[i], "-match_score") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -match_score\n") ;
				usage();
				exit(1);
			}
			i++;
			M_SCORE = atof(argv[i]);
		}

		//mismatch score
		if (strcmp(argv[i], "-m") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -m\n") ;
				usage();
				exit(1);
			}
			i++;
			MM_SCORE = atof(argv[i]);
			if (MM_SCORE < 0) {
				fprintf(stderr, "ERROR: Mismatch score must be positive!\n");
				exit(1);
			}
			if (MM_SCORE == 0)
				fprintf(stderr,
						"\n!!! WARNING: Mismatch score is 0! This could lead to bad alignments!\n\n");
		}

		//gap score
		if (strcmp(argv[i], "-g") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -g\n") ;
				usage();
				exit(1);
			}
			i++;
			GAP_SCORE = atof(argv[i]);
			if (GAP_SCORE < 0) {
				fprintf(stderr, "ERROR: Gap score must be positive!\n");
				exit(1);
			}
			if (GAP_SCORE == 0)
				fprintf(stderr,
						"\n!!! WARNING: Gap score is 0! This could lead to bad alignments!\n\n");
		}

		//hitlength limit
		if (strcmp(argv[i], "-l") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -l\n") ;
				usage();
				exit(1);
			}
			i++;
			if ((HITLEN_LIMIT = atoi(argv[i])) == 0) {
				fprintf(stderr,
						"ERROR: Hitlength limit must be an integer value unequal to 0!\n");
				exit(1);
			}
		}

		//chromosome container size
		if (strcmp(argv[i], "-c") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -c\n") ;
				usage();
				exit(1);
			}
			i++;
			if ((CHROM_CONTAINER_SIZE += atoi(argv[i])) == 0) {
				fprintf(stderr,
						"ERROR: Chromosome Container Size must be an integer value unequal to 0!\n");
				exit(1);
			}
		}

		//nr of unallowed chars in reads
		/*if(strcmp(argv[i],"-n")==0){
		 not_defined = 0;
		 if (i+1 > argc - 1){ usage(); exit(1); }
		 i++;
		 if ((NUM_UNALLOWED_CHARS = atoi(argv[i])) == 0 && argv[i][0] != '0') {
		 fprintf(stderr, "ERROR: Number of non-base symbols in read must be an integer value!\n");
		 exit(1);
		 }
		 }*/

		//gaps most right
		if (strcmp(argv[i], "-d") == 0) {
			not_defined = 0;
			GAPS_MOST_RIGHT = 1;
		}

		//overhang alignment
		if (strcmp(argv[i], "-h") == 0) {
			not_defined = 0;
			OVERHANG_ALIGNMENT = 1;
		}

		//overhang alignment
		if (strcmp(argv[i], "-w") == 0) {
			not_defined = 0;
			STRINGENT_GAPLIMIT = 0;
		}

		//statistics
		if (strcmp(argv[i], "-s") == 0) {
			not_defined = 0;
			STATISTICS = 1;
		}

		//print out gene, too (for every hit) used for WMD2
		if (strcmp(argv[i], "-target_info") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr,
						"-target_info needs an integer value (1: length only, 2: length+sequence)\n");
				exit(1);
			}
			i++;
			if (((PRINT_SEQ = atoi(argv[i])) == 0) || PRINT_SEQ < 1
					|| PRINT_SEQ > 2) {
				fprintf(stderr, "-target_info value must be either 1 or 2!\n");
				exit(1);
			}
		}

		//flanking region of a hit
		if (strcmp(argv[i], "-flanking") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(
						stderr,
						"-flanking needs an integer value, the length of one of the symmetric flanking regions of a hit\n");
				exit(1);
			}
			i++;
			if (((FLANKING = atoi(argv[i])) == 0) || FLANKING < 0 || FLANKING
					> 100) {
				fprintf(stderr,
						"-flanking value must be a positive integer and must not exceed 100!\n");
				exit(1);
			}
		}

		//print out gene, too (for every hit) used for WMD2
		if (strcmp(argv[i], "-qpalma") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr,
						"-qpalma needs an argument\n");
				exit(1);
			}
			i++;
			QPALMA_FILE.assign(argv[i]);
		}

		if (strcmp(argv[i], "-acc") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr,
						"-acc needs an argument\n");
				exit(1);
			}
			i++;
			ACC_FILES.assign(argv[i]);
		}

		if (strcmp(argv[i], "-don") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr,
						"-don needs an argument\n");
				exit(1);
			}
			i++;
			DON_FILES.assign(argv[i]);
		}

		if (strcmp(argv[i], "-no-ss-pred") == 0) {
			not_defined = 0;
			NO_SPLICE_PREDICTIONS = 1 ;
		}

        // qpalma-triggered reads file                     //  #A#
		if (strcmp(argv[i], "-log-triggered-reads") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -log-triggered-reads\n") ;
				usage();
				exit(1);
			}
			i++;
            LOG_TRIGGERED = 1;
			TRIGGERED_LOG_FILE.assign(argv[i]);
		}                                                   // #A#

		if (not_defined == 1) {
			fprintf(stderr, "ERROR: unknown option %s\n", argv[i]) ;
			usage();
			exit(1);
		}
	}

	if (has_index == 0 || has_query == 0 || has_genome == 0) {
		usage();
		exit(1);
	}

	NOT_MAXIMAL_HITS = SEED_HIT_CANCEL_THRESHOLD || INDEX_DEPTH_EXTRA_THRESHOLD;

	return 0;
}

void Config::VersionHeader()
{
	printf("\nPALMapper version %s   (PALMapper is a fusion of GenomeMapper & QPALMA)\n", VERSION);
	printf("written by Korbinian Schneeberger, Joerg Hagmann, Gunnar Raetsch, Geraldine Jean, Fabio De Bona, Stephan Ossowski, and others\n");
	printf("Max Planck Institute for Developmental Biology and Friedrich Miescher Laboratory, Tuebingen, Germany, 2008-2010\n\n");
}

int Config::usage() 
{
	//VersionHeader() ;

	printf("USAGE: palmapper [options]\n");
	printf("\n");
	printf("mandatory:\n");
	printf(" -i STRING      reference sequence (fasta file and prefix to index files)\n");
	printf(" -q STRING      query filename (fasta, fastq, SHORE flat file)\n");
	printf("\n\n");
	printf("optional:\n");
	printf(" -S             report spliced alignments (detailed options below)\n\n");
	printf(" -f STRING      output format (\"shore\", \"bed\", \"bedx\", or \"sam\")\n");
	printf(" -o STRING      output filename (stdout)\n");
	printf(" -H STRING      output filename for spliced hits (stdout)\n");
	printf(" -u STRING      unmapped reads filename\n\n");

	printf(" -a             report all alignments (best alignments only)\n");
	printf(" -ar INT        report a limited number of alignments (best alignments only)\n");
	printf(" -z INT         report a number of top alignments\n\n");

	printf(" -r             disable reverse alignment\n");
	printf(" -h             perform alignment of flanking regions of hits\n");
	printf(" -d             align gaps most right (most left) (ignored for spliced alignments)\n");
	printf(" -w             allow more gaps for best hit (ignored for spliced alignments)\n");
	//printf(" -e         report edit operations (alignment scores)\n");
	printf(" -l INT         seed length (index size)\n");
	printf(" -n INT         max number of best alignments (all)\n");
	printf(" -c INT         seed container size (15.000.000)\n\n");

	printf(" -threads INT                       maximal number of threads (4) \n");
	printf(" -seed-hit-cancel-threshold INT     number of hits of a seed that lead to its ignoration\n");
	printf(" -index-extend-threshold INT        number of hits of a seed that lead to an seed-length extension\n");
	printf(" -index-extend INT                  length of seed-length extension\n");
	printf(" -index-precache                    linearly read index file to fill caches\n\n");

	printf(" -rtrim INT                         shortens the read until a hit is found or the minimal length is reached (INT)\n");
	printf(" -polytrim INT                      trims polyA or polyT ends until a hit is found or the minimal length is reached (INT)\n");
	printf(" -adaptertrim INT                   trims away known adapter sequences, read is dropped, when shorter than parameter (INT)\n\n");

	printf(" -rlim INT      limit the number of reads for alignment\n\n");

	printf(" -M INT         max number of mismatches (3)\n");
	printf(" -G INT         max number of gaps (1)\n");
	printf(" -E INT         max edit operations(3)\n");
	printf(" -m DOUBLE      mismatch penalty (4)\n");
	printf(" -g DOUBLE      gap penalty (5)\n");

	printf(" -v             verbose (silent)\n\n");

	printf("spliced hits definitions: (-S required)\n");
	printf(" -qpalma STRING                        file name with qpalma parameters (essential)\n");
	printf(" -qpalma-use-map-max-len INT           limit the map extension up- and downstream to the given length (10.000)\n\n");

	printf(" -acc STRING                           path name to acceptor splice site predictions (essential)\n");
	printf(" -don STRING                           path name to donor splice site predictions (essential)\n");
	printf(" -no-ss-pred                           indicates that no splice site predictions should be used\n\n");

	printf(" -filter-splice-sites-top-perc FLOAT   trigger spliced alignments, if read covers top percentile splice site (between 0 and 1)\n");
	printf(" -filter-max-mismatches INT            trigger spliced alignment, if unspliced alignment has at least this many mismatches\n");
	printf(" -filter-max-gaps INT                  trigger spliced alignment, if unspliced alignment has at least this many mismatches\n\n");

	printf(" -C INT                                min combined length (auto)\n");
	printf(" -L INT                                min length of long hit (auto)\n");
	printf(" -K INT                                min length of short hit (auto)\n");
	printf(" -I INT                                longest intron length  (auto)\n");
	printf(" -SA INT                               maximum number of spliced alignments per read (10)\n");
	printf(" -NI INT                               maximum number of introns in spliced alignments (auto)\n");
	printf(" -CT INT                               distance to tolerate between hit and existing hit cluster (10)\n");
	printf(" -QMM INT                              number of matches required for identifying a splice site (3)\n\n");

	printf(" -report STRING                        file for map reporting\n");
	printf(" -report-ro STRING                     file for map reporting (read only)\n");
	printf(" -report-rep-seed                      switch on reporting of repetitive seeds\n");
	printf(" -report-map-region                    switch on reporting of mapped regions\n");
	printf(" -report-map-read                      switch on reporting of mapped reads\n");
	printf(" -report-spliced-read                  switch on reporting of spliced reads\n");
	printf(" -report-splice-sites FLOAT            report splice sites with confidence not less that threshold\n");
	printf(" -report-splice-sites-top-perc FLOAT   report splice sites with confidence in top percentile (between 0 and 1)\n");
	printf(" -report-gff-init STRING               initialize map with exons from GFF file\n");
	//printf(" -qpalma-use-map                       use map for qpalma alignments\n");

	return 0;
}