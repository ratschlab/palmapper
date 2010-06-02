// Authors: Korbinian Schneeberger, Stephan Ossowski, Joerg Hagmann, Gunnar Raetsch
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#ifndef GENOMEMAPPER_H_
#define GENOMEMAPPER_H_

#include "shogun/init.h"
#include "shogun/common.h"
#include "shogun/io.h"
#include "shogun/Array.h"
#include "DNAArray.h"
#include "DNAArray4.h"

#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <string>


//#define BinaryStream_MAP

#ifdef BinaryStream_MAP
#include "BinaryStream.h"
#endif


// ##############################################################
// ####### GLOBAL VARIABLES #####################################
// ##############################################################

#include <palmapper/Config.h>
#include <palmapper/Genome.h>
#include <palmapper/Hits.h>
#include <palmapper/Read.h>
#include <palmapper/Statistics.h>
#include <palmapper/Util.h>
#include <palmapper/TopAlignments.h>
#include <palmapper/QPalma.h>
#include <palmapper/GenomeMaps.h>

extern Config _config;
extern Statistics _stats;
//extern Read _read;

// ##############################################################
// ####### FILE HANDLING ########################################
// ##############################################################

extern FILE *OUT_FP;
extern FILE *SP_OUT_FP;
extern FILE *TRIGGERED_LOG_FP; // #A#

// init.c
extern int init_output_file(Config * config);
extern int init_spliced_output_file(Config * config);

#endif
