// Authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include "genomemapper.h"
#include "align.h"
#include <assert.h>

int init_defaults(Config * config);
int init_output_file(Config * config);
int init_spliced_output_file(Config * config);

int valid_char[256] ;
char compl_char[256] ;
char upper_char[256] ;

int init(int argc, char *argv[], Config * config)
{
	init_defaults(config);
	config->parseCommandLine(argc, argv);
	init_output_file(config);
	init_spliced_output_file(config);

	//_qpalma.init_qpalma() ;

	return 0;
}

int init_defaults(Config * config) {

	config->NUM_THREADS = 1;
	_stats.NUM_HITS = 0;
	_stats.NUM_ALIGNMENTS = 0;

	config->ALL_HIT_STRATEGY = 0; // default: best hit strategy
	config->SUMMARY_HIT_STRATEGY = 0; //            "
	config->HITLEN_LIMIT = 0; // only temporarily 0

	config->VERBOSE = 0;

	config->MAP_REVERSE = 1;

	_stats.READS_MAPPED = 0;
	config->NUM_EDIT_OPS = 4;
	config->NUM_GAPS = 1;
	config->NUM_MISMATCHES = 3;
	config->M_SCORE = 0;
	config->MM_SCORE = 4;
	config->GAP_SCORE = 5;
	config->GAPS_MOST_RIGHT = 0;
	config->OVERHANG_ALIGNMENT = 1;
	config->STRINGENT_GAPLIMIT = 1;
	config->PRINT_SEQ = 0;
	config->FLANKING = 0;
    config->LOG_TRIGGERED = 0;   // #A#

	config->OUTPUT_FORMAT = OUTPUT_FORMAT_SHORE;
	config->SCORES_OUT = 1;
	config->OUTPUT_FILTER_NUM_RANDOM = 0;
	config->OUTPUT_FILTER_NUM_TOP=1 ;
	config->CHROM_CONTAINER_SIZE = 15000000;
	config->STATISTICS = 0;

	_stats.listcount = 0;
	_stats.listocc = 0;

	return (0);
}

int init_output_file(Config * config) {
	if (config->OUT_FILE_NAME.length() > 0) {
		if ((OUT_FP = fopen(config->OUT_FILE_NAME.c_str(), "w")) == NULL) {
			fprintf(stderr, "ERROR : Couldn't open output file %s\n",
					config->OUT_FILE_NAME.c_str());
			exit(1);
		}
	} else {
		OUT_FP = stdout;
	}

	return (0);
}

int init_spliced_output_file(Config * config) {
	if (config->SPLICED_HITS && config->SPLICED_OUT_FILE_NAME.length() > 0) {
		if ((SP_OUT_FP = fopen(config->SPLICED_OUT_FILE_NAME.c_str(), "w")) == NULL) {
			fprintf(stderr,
					"ERROR : Couldn't open output file for spliced hits %s\n",
					config->SPLICED_OUT_FILE_NAME.c_str());
			exit(1);
		}
	} else {
		SP_OUT_FP = stdout;
	}

	return (0);
}
