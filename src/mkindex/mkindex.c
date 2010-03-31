// Authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mkindex.h"
#include "mkindex_symbols.c"

int load_chromosomes();
char* get_seq(unsigned int n);

int main(int argc, char *argv[])
{
	clock_t start, end;
	double elapsed;

	start = clock();
	
	if (VERBOSE) { printf("Start initializing\n"); }

	init(argc, argv);

	if (VERBOSE) { printf("Start loading\n"); }

	load_chromosomes();

	if (VERBOSE) printf("\nTotal number of seed occurrences: %lu\n\n", POSITION_COUNTER);

	if (VERBOSE) { printf("Finish.\n"); }
	
	end = clock();
	elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
	if (VERBOSE) printf ("Time needed: %g s\n",elapsed);
	
	return EXIT_SUCCESS;
}


char *get_seq(unsigned int n)
{
	char *seq = (char *) malloc (INDEX_DEPTH*sizeof(char));
	int i, c;
	for (i=INDEX_DEPTH-1; i>=0; --i) {
		c = (int) (n / POWER[i]);
		switch (c)
		{
			case 0: seq[i] = 'A';
					break;
			case 1: seq[i] = 'C';
					break;
			case 2: seq[i] = 'G';
					break;
			case 3: seq[i] = 'T';
					break;
		}
		n -= (int) (c * POWER[i]);
	}
	return seq;
}

