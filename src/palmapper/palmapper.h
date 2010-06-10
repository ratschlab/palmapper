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

//#define BinaryStream_MAP

#ifdef BinaryStream_MAP
#include "BinaryStream.h"
#endif


// ##############################################################
// ####### GLOBAL VARIABLES #####################################
// ##############################################################

#include <palmapper/Config.h>
#include <palmapper/Statistics.h>

extern Config _config;
extern Statistics _stats;

// ##############################################################
// ####### FILE HANDLING ########################################
// ##############################################################

//extern FILE *OUT_FP;
//extern FILE *SP_OUT_FP;

// init.c
extern FILE *init_output_file(Config &config);
extern FILE *init_spliced_output_file(Config &config, FILE *OUT_FP);

#endif
