// Authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include "genomemapper.h"
#include <assert.h>

int init_defaults();
int init_constants();
int init_output_file();
int init_spliced_output_file();
int init_hit_lists();
int init_operators();
int init_statistic_vars();
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
	init_hit_lists();

	//_qpalma.init_qpalma() ;

	return 0;
}

int init_alignment_structures() {
	REDUNDANT = 0;

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

	ALIGNSEQ = (char *) malloc((_config.MAX_READ_LENGTH + 3 * Config::MAX_EDIT_OPS)
			* sizeof(char));
	if (ALIGNSEQ == NULL) {
		fprintf(stderr, "[init_alignment_structures] Could not allocate memory\n");
		exit(1);
	}

	//only for debugging purposes:
	//chrseq = (char *) malloc (_config.MAX_READ_LENGTH * sizeof(char));
	//chrseq[0] = '\0';

	return (0);
}

int init_from_meta_index() {
	init_constants(); // updated

	alloc_genome_memory(); // updated

	init_operators(); //updated

	if (_config.STATISTICS) {
		init_statistic_vars(); //updated
	}
	if (!_config.HITLEN_LIMIT)
		_config.HITLEN_LIMIT = _config.INDEX_DEPTH;

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
	NUM_MATCHES = 3;
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

int init_constants() {
	if (_config.INDEX_DEPTH == 5) {
		BINARY_CODE[0] = 0; //binary number: 0000 0000 0000 0000 0000 0000
		BINARY_CODE[1] = 256; //binary number: 0000 0000 0000 0001 0000 0000
		BINARY_CODE[2] = 512; //binary number: 0000 0000 0000 0010 0000 0000
		BINARY_CODE[3] = 768; //binary number: 0000 0000 0000 0011 0000 0000
	}
	if (_config.INDEX_DEPTH == 6) {
		BINARY_CODE[0] = 0; //binary number: 0000 0000 0000 0000 0000 0000
		BINARY_CODE[1] = 1024; //binary number: 0000 0000 0000 0100 0000 0000
		BINARY_CODE[2] = 2048; //binary number: 0000 0000 0000 1000 0000 0000
		BINARY_CODE[3] = 3072; //binary number: 0000 0000 0000 1100 0000 0000
	}
	if (_config.INDEX_DEPTH == 7) {
		BINARY_CODE[0] = 0; //binary number: 0000 0000 0000 0000 0000 0000
		BINARY_CODE[1] = 4096; //binary number: 0000 0000 0001 0000 0000 0000
		BINARY_CODE[2] = 8192; //binary number: 0000 0000 0010 0000 0000 0000
		BINARY_CODE[3] = 12288; //binary number: 0000 0000 0011 0000 0000 0000
	}
	if (_config.INDEX_DEPTH == 8) {
		BINARY_CODE[0] = 0; //binary number: 0000 0000 0000 0000 0000 0000
		BINARY_CODE[1] = 16384; //binary number: 0000 0000 0100 0000 0000 0000
		BINARY_CODE[2] = 32768; //binary number: 0000 0000 1000 0000 0000 0000
		BINARY_CODE[3] = 49152; //binary number: 0000 0000 1100 0000 0000 0000
	}
	if (_config.INDEX_DEPTH == 9) {
		BINARY_CODE[0] = 0; //binary number: 0000 0000 0000 0000 0000 0000
		BINARY_CODE[1] = 65536; //binary number: 0000 0001 0000 0000 0000 0000
		BINARY_CODE[2] = 131072; //binary number: 0000 0010 0000 0000 0000 0000
		BINARY_CODE[3] = 196608; //binary number: 0000 0011 0000 0000 0000 0000
	}
	if (_config.INDEX_DEPTH == 10) {
		BINARY_CODE[0] = 0; //binary number: 0000 0000 0000 0000 0000 0000
		BINARY_CODE[1] = 262144; //binary number: 0000 0100 0000 0000 0000 0000
		BINARY_CODE[2] = 524288; //binary number: 0000 1000 0000 0000 0000 0000
		BINARY_CODE[3] = 786432; //binary number: 0000 1100 0000 0000 0000 0000
	}
	if (_config.INDEX_DEPTH == 11) {
		BINARY_CODE[0] = 0; //binary number: 0000 0000 0000 0000 0000 0000
		BINARY_CODE[1] = 1048576; //binary number: 0001 0000 0000 0000 0000 0000
		BINARY_CODE[2] = 2097152; //binary number: 0010 0000 0000 0000 0000 0000
		BINARY_CODE[3] = 3145728; //binary number: 0011 0000 0000 0000 0000 0000
	}
	if (_config.INDEX_DEPTH == 12) {
		BINARY_CODE[0] = 0; //binary number: 0000 0000 0000 0000 0000 0000
		BINARY_CODE[1] = 4194304; //binary number: 0100 0000 0000 0000 0000 0000
		BINARY_CODE[2] = 8388608; //binary number: 1000 0000 0000 0000 0000 0000
		BINARY_CODE[3] = 12582912; //binary number: 1100 0000 0000 0000 0000 0000
	}
	if (_config.INDEX_DEPTH == 13) {
		BINARY_CODE[0] = 0; //binary number: 0000 0000 0000 0000 0000 0000 0000
		BINARY_CODE[1] = 16777216; //binary number: 0001 0000 0000 0000 0000 0000 0000
		BINARY_CODE[2] = 33554432; //binary number: 0010 0000 0000 0000 0000 0000 0000
		BINARY_CODE[3] = 50331648; //binary number: 0011 0000 0000 0000 0000 0000 0000
	}
	if (_config.INDEX_DEPTH == 14) {
		BINARY_CODE[0] = 0; //binary number: 0000 0000 0000 0000 0000 0000 0000
		BINARY_CODE[1] = 67108864; //binary number: 0001 0000 0000 0000 0000 0000 0000
		BINARY_CODE[2] = 134217728; //binary number: 0010 0000 0000 0000 0000 0000 0000
		BINARY_CODE[3] = 201326592; //binary number: 0011 0000 0000 0000 0000 0000 0000
	}
	if (_config.INDEX_DEPTH>14 || _config.INDEX_DEPTH<5)
	  {
	    fprintf(stderr, "ERROR: _config.INDEX_DEPTH out of range\n") ;
	    exit(1) ;
	  }

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

int init_hit_lists() {
	alloc_hit_lists_operator();
	alloc_hits_by_score(); // A T T E N T I O N ! ! !   needs correct _config.NUM_EDIT_OPS which can be changed in init_alignment_structures() !!!

	return (0);
}

int init_operators() {
	MAPPING_ENTRY_OPERATOR = alloc_mapping_entry_container();
	assert(MAPPING_ENTRY_OPERATOR!=NULL) ;

	MAPPING_ENTRY_OPERATOR_FIRST = MAPPING_ENTRY_OPERATOR;

	HIT_OPERATOR = alloc_hit_container();
	assert(HIT_OPERATOR!=NULL) ;
	HIT_OPERATOR_FIRST = HIT_OPERATOR;

	alloc_chromosome_entry_container();

	return (0);
}

int init_statistic_vars() {
	new (&_stats) Statistics();
	return (0);
}

