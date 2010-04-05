// Authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include "genomemapper.h"
#include "align.h"
#include <assert.h>

int init_output_file(Config * config);
int init_spliced_output_file(Config * config);

int valid_char[256] ;
char compl_char[256] ;
char upper_char[256] ;


int init_output_file(Config * config) 
{
	if (config->OUT_FILE_NAME.length() > 0) 
	{
		if ((OUT_FP = fopen(config->OUT_FILE_NAME.c_str(), "w")) == NULL) 
		{
			fprintf(stderr, "ERROR : Couldn't open output file %s\n",
					config->OUT_FILE_NAME.c_str());
			exit(1);
		}
	} else 
	{
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
