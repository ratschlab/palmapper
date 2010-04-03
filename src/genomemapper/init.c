// Authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include "genomemapper.h"
#include "align.h"
#include <assert.h>

int init_defaults();
int init_output_file();
int init_spliced_output_file();
int init_alignment_structures();

int valid_char[256] ;
char compl_char[256] ;
char upper_char[256] ;

int init(int argc, char *argv[]) {
	init_defaults();
	_config.parseCommandLine(argc, argv);
	init_output_file();
	init_spliced_output_file();
	init_alignment_structures();
	_hits.init_hit_lists();

	//_qpalma.init_qpalma() ;

	return 0;
}

int init_alignment_structures() {
	_hits.REDUNDANT = 0;

	/////////////////////////////
	////// if you change _config.NUM_EDIT_OPS here, then put alloc_hits_by_editops in init_hit_lists_operator() after this function!!!
	/////////////////////////////

	if (_config.NUM_MISMATCHES > _config.NUM_EDIT_OPS) {
		fprintf(
				stderr,
				"\n!!! WARNING: Nr. of allowed mismatches is set to nr. of allowed edit operations!\n\n");
		_config.NUM_MISMATCHES = _config.NUM_EDIT_OPS;
	}
	if (_config.NUM_GAPS > _config.NUM_EDIT_OPS) {
		fprintf(
				stderr,
				"\n!!! WARNING: Nr. of allowed gaps is set to nr. of allowed edit operations!\n\n");
		_config.NUM_GAPS = _config.NUM_EDIT_OPS;
	}

	if (_config.M_SCORE == _config.MM_SCORE) {
		fprintf(
				stderr,
				"ERROR: Sorry, until now, program cannot handle equal mismatch- and match-scores! Please restart with different values!\n");
		exit(1);
	}

	if (_config.VERBOSE)
		printf(
				"allowed edit operations: %d / allowed mismatches: %d / allowed gaps: %d\n%s strategy%s\t%s overhang alignment\n------------------\n",
				_config.NUM_EDIT_OPS, _config.NUM_MISMATCHES, _config.NUM_GAPS,
				_config.ALL_HIT_STRATEGY ? "all hit" : "best hit",
				_config.SUMMARY_HIT_STRATEGY ? " (summary)" : "",
				_config.OVERHANG_ALIGNMENT ? "with" : "without");

	if (_config.GAP_SCORE > _config.MM_SCORE)
	WORST_SCORE = _config.NUM_GAPS * _config.GAP_SCORE + (_config.NUM_EDIT_OPS - _config.NUM_GAPS)
				* _config.MM_SCORE;
	else
		WORST_SCORE = _config.NUM_MISMATCHES * _config.MM_SCORE + (_config.NUM_EDIT_OPS
				- _config.NUM_MISMATCHES) * _config.GAP_SCORE;
	WORST_MM_SCORE = _config.NUM_MISMATCHES * _config.MM_SCORE;
  

	//only for debugging purposes:
	//chrseq = (char *) malloc (_config.MAX_READ_LENGTH * sizeof(char));
	//chrseq[0] = '\0';

	return (0);
}


int init_defaults() {

	_config.NUM_THREADS = 1;
	_stats.NUM_HITS = 0;
	_stats.NUM_ALIGNMENTS = 0;

	_config.ALL_HIT_STRATEGY = 0; // default: best hit strategy
	_config.SUMMARY_HIT_STRATEGY = 0; //            "
	_config.HITLEN_LIMIT = 0; // only temporarily 0

	_config.VERBOSE = 0;

	_config.MAP_REVERSE = 1;

	_stats.READS_MAPPED = 0;
	_config.NUM_EDIT_OPS = 4;
	_config.NUM_GAPS = 1;
	_config.NUM_MISMATCHES = 3;
	_config.M_SCORE = 0;
	_config.MM_SCORE = 4;
	_config.GAP_SCORE = 5;
	_config.GAPS_MOST_RIGHT = 0;
	_config.OVERHANG_ALIGNMENT = 1;
	_config.STRINGENT_GAPLIMIT = 1;
	_config.PRINT_SEQ = 0;
	_config.FLANKING = 0;
    _config.LOG_TRIGGERED = 0;   // #A#

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

	_config.OUTPUT_FORMAT = OUTPUT_FORMAT_SHORE;

	_config.SCORES_OUT = 1;

	_config.REPEATMAP = 0;

	_config.CHROM_CONTAINER_SIZE = 15000000;

	_config.STATISTICS = 0;

	_stats.listcount = 0;
	_stats.listocc = 0;

	return (0);
}


int init_output_file() {
	if (_config.OUT_FILE_NAME.length() > 0) {
		if ((OUT_FP = fopen(_config.OUT_FILE_NAME.c_str(), "w")) == NULL) {
			fprintf(stderr, "ERROR : Couldn't open output file %s\n",
					_config.OUT_FILE_NAME.c_str());
			exit(1);
		}
	} else {
		OUT_FP = stdout;
	}

	return (0);
}

int init_spliced_output_file() {
	if (_config.SPLICED_HITS && _config.SPLICED_OUT_FILE_NAME.length() > 0) {
		if ((SP_OUT_FP = fopen(_config.SPLICED_OUT_FILE_NAME.c_str(), "w")) == NULL) {
			fprintf(stderr,
					"ERROR : Couldn't open output file for spliced hits %s\n",
					_config.SPLICED_OUT_FILE_NAME.c_str());
			exit(1);
		}
	} else {
		SP_OUT_FP = stdout;
	}

	return (0);
}
