// Authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "pmindex.h"
#include "pmindex_symbols.c" 

int seed_hit_cancel_threshold = -1 ;

int load_chromosomes();
char* get_seq(unsigned int n);

int main(int argc, char *argv[])
{
	clock_t start, end;
	double elapsed;

	start = clock();
	
	if (VERBOSE) { printf("Start initializing\n"); }

	init(argc, argv);

	INDEX = (BIN**)malloc(INDEX_SIZE*sizeof(BIN*)) ;
	if (INDEX==NULL)
	{
		fprintf(stderr, "Failed to allocate %ld bytes\n", INDEX_SIZE*sizeof(BIN*)) ;
		exit(-1) ;
	}
	memset(INDEX, 0, INDEX_SIZE*sizeof(BIN*)) ;

	USED_SLOTS =(unsigned int*)malloc(INDEX_SIZE*sizeof(unsigned int)) ;
	if (USED_SLOTS==NULL)
	{
		fprintf(stderr, "Failed to allocate %ld bytes\n", INDEX_SIZE*sizeof(unsigned int)) ;
		exit(-1) ;
	}
	memset(INDEX, 0, INDEX_SIZE*sizeof(unsigned int)) ;
	
	if (VERBOSE) { printf("Start loading\n"); }

	Genome genome(0) ; 
	if (strlen(GENOME_VARIANTS_FILE_NAME)>0 || strlen(GENOME_MASK_FILE_NAME)>0 || strlen(GENOME_MASK_GFF_FILE_NAME)>0)
	{
		fprintf(stdout, "loading complete genome\n") ;
		genome.load_genome() ;
	}
	VariantMap variants(genome, false, false) ;
	if (strlen(GENOME_VARIANTS_FILE_NAME)>0)
	{
		std::string fnames = std::string(GENOME_VARIANTS_FILE_NAME) ;
		variants.init_from_files(fnames) ;
		int cnt=variants.get_read_id_num() ;
		fprintf(stdout, "found %i source ids\n", cnt) ;
	}
	assert(variants.genome!=NULL) ;

	GenomeMaps genome_mask(genome) ;
	if (strlen(GENOME_MASK_FILE_NAME)>0)
	  genome_mask.read_reporting(GENOME_MASK_FILE_NAME) ;
	if (strlen(GENOME_MASK_GFF_FILE_NAME)>0)
	  {
	    std::string fname=std::string(GENOME_MASK_GFF_FILE_NAME) ;
	    genome_mask.init_with_gff(fname) ;
	  }

	load_chromosomes(genome, variants, genome_mask);

	if (strlen(GENOME_MASK_FILE_NAME)>0)
	  genome_mask.write_reporting(GENOME_MASK_FILE_NAME) ;

	if (VERBOSE) printf("\nTotal number of seed occurrences: %lu\n\n", POSITION_COUNTER);

	if (VERBOSE) { printf("Finish.\n"); }
	
	end = clock();
	elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
	if (VERBOSE) printf ("Time needed: %g s\n",elapsed);

	if (GENOME_OUT_FP) 
	  fclose(GENOME_OUT_FP) ;
	if (GENOME_FP)
	  fclose(GENOME_FP) ;
	if (CHR_INDEX_FP)
	  fclose(CHR_INDEX_FP) ;
	if (META_INDEX_FP)
	  fclose(META_INDEX_FP) ;
	if (MAPFWD_INDEX_FP)
	  fclose(MAPFWD_INDEX_FP) ;

	return EXIT_SUCCESS;
}


char *get_seq(unsigned int n)
{
	char *seq = (char *) malloc (INDEX_DEPTH*sizeof(char));
	memset(seq, 0, INDEX_DEPTH*sizeof(char)) ;
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

