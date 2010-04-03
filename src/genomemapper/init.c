// Authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include "genomemapper.h"
#include "align.h"
#include <assert.h>

int init_defaults(Config * config);
int init_output_file(Config * config);
int init_spliced_output_file(Config * config);
int init_alignment_structures(Config * config, Hits * hits);

int valid_char[256] ;
char compl_char[256] ;
char upper_char[256] ;

int init(int argc, char *argv[], Config * config, Hits* hits) 
{
	init_defaults(config);
	config->parseCommandLine(argc, argv);
	init_output_file(config);
	init_spliced_output_file(config);
	init_alignment_structures(config, hits);
	hits->init_hit_lists();

	//_qpalma.init_qpalma() ;

	return 0;
}

int init_alignment_structures(Config * config, Hits* hits) {
	hits->REDUNDANT = 0;

	/////////////////////////////
	////// if you change config->NUM_EDIT_OPS here, then put alloc_hits_by_editops in init_hit_lists_operator() after this function!!!
	/////////////////////////////

	if (config->NUM_MISMATCHES > config->NUM_EDIT_OPS) {
		fprintf(
				stderr,
				"\n!!! WARNING: Nr. of allowed mismatches is set to nr. of allowed edit operations!\n\n");
		config->NUM_MISMATCHES = config->NUM_EDIT_OPS;
	}
	if (config->NUM_GAPS > config->NUM_EDIT_OPS) {
		fprintf(
				stderr,
				"\n!!! WARNING: Nr. of allowed gaps is set to nr. of allowed edit operations!\n\n");
		config->NUM_GAPS = config->NUM_EDIT_OPS;
	}

	if (config->M_SCORE == config->MM_SCORE) {
		fprintf(
				stderr,
				"ERROR: Sorry, until now, program cannot handle equal mismatch- and match-scores! Please restart with different values!\n");
		exit(1);
	}

	if (config->VERBOSE)
		printf(
				"allowed edit operations: %d / allowed mismatches: %d / allowed gaps: %d\n%s strategy%s\t%s overhang alignment\n------------------\n",
				config->NUM_EDIT_OPS, config->NUM_MISMATCHES, config->NUM_GAPS,
				config->ALL_HIT_STRATEGY ? "all hit" : "best hit",
				config->SUMMARY_HIT_STRATEGY ? " (summary)" : "",
				config->OVERHANG_ALIGNMENT ? "with" : "without");

	if (config->GAP_SCORE > config->MM_SCORE)
	WORST_SCORE = config->NUM_GAPS * config->GAP_SCORE + (config->NUM_EDIT_OPS - config->NUM_GAPS)
				* config->MM_SCORE;
	else
		WORST_SCORE = config->NUM_MISMATCHES * config->MM_SCORE + (config->NUM_EDIT_OPS
				- config->NUM_MISMATCHES) * config->GAP_SCORE;
	WORST_MM_SCORE = config->NUM_MISMATCHES * config->MM_SCORE;
  

	//only for debugging purposes:
	//chrseq = (char *) malloc (config->MAX_READ_LENGTH * sizeof(char));
	//chrseq[0] = '\0';

	return (0);
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

	// lt: Previously, there was no minimum distance for spliced hits.
	// In effect, this meant that a pair of spliced hits could overlap.
	// In order to preserve the original behaviour if the
	// user does not explicitly define a minimum distance for spliced hits,
	// we define the default minimum distance to be -ASSUMED_READ_LENGTH.
	// The effect is that overlapping pairs of spliced hits will still pass
	// the test in print.c:comp_hits_4_splicing:
	// if (hit1->start < hit2->start) {
	//	  dist = (hit2->start - hit1->end + 1);
	// } else {
	// 	  dist = (hit1->start - hit2->end + 1);
	// }
	// if (dist > SPLICED_HIT_MAX_DIST || dist < SPLICED_HIT_MIN_DIST) {
	// 	  return count;
	// }

	config->OUTPUT_FORMAT = OUTPUT_FORMAT_SHORE;

	config->SCORES_OUT = 1;

	config->REPEATMAP = 0;

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
