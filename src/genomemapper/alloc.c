// Authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include "genomemapper.h"

int alloc_index_memory()
{

        //if ( (MEM_MGR.first_entry = (STORAGE_ENTRY *) malloc (NUM_POS * (1+_config.MAP_REVERSE) * sizeof(STORAGE_ENTRY)) ) == NULL) {
        //        fprintf(stderr, "ERROR : not enough memory for mem_master (1)\n");
        //        exit(1);
        //}

	//MEM_MGR.num_bins = 0;
	//MEM_MGR.next_unused_entry = MEM_MGR.first_entry;

	INDEX_SIZE=_config.INDEX_SIZE_13 ;
	if (_config.INDEX_DEPTH==12)
		INDEX_SIZE = _config.INDEX_SIZE_12 ;
	
	fprintf(stdout, "_config.INDEX_DEPTH=%i, INDEX_SIZE=%ld, sizeof(INDEX_ENTRY)=%ld, index size=%ld\n",  (int)_config.INDEX_DEPTH, (long int)INDEX_SIZE, sizeof(INDEX_ENTRY), sizeof(INDEX_ENTRY)*INDEX_SIZE) ;

	if ( (INDEX = (INDEX_ENTRY *) calloc (INDEX_SIZE, sizeof(INDEX_ENTRY)) ) == NULL) {
		fprintf(stderr, "ERROR : not enough memory for mem_master (2)\n");
		exit(1);
	}

	if ( _config.MAP_REVERSE && ((INDEX_REV = (INDEX_ENTRY *) calloc (INDEX_SIZE, sizeof(INDEX_ENTRY))) == NULL) ) {
		fprintf(stderr, "ERROR : not enough memory for mem_master (3)\n");
		exit(1);
	}

	return(0);
}


