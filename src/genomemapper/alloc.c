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

int alloc_genome_memory()
{
	if ((GENOME = (CHROMOSOME_ENTRY **) calloc (LONGEST_CHROMOSOME, sizeof(CHROMOSOME_ENTRY**))) == NULL) {
		fprintf(stderr, "ERROR : not enough memory for genome memory\n");
		exit(1);
	}
		
	//@TODO is this really necessary? why isn't it already NULL
	unsigned int i;
	for (i=0; i!=LONGEST_CHROMOSOME; ++i) {
		*(GENOME+i) = NULL;
	}

	return(0);
}

MAPPING_ENTRY* alloc_mapping_entry()
{
	MAPPING_ENTRY_CONTAINER *container;
	MAPPING_ENTRY *entry;

	if (MAPPING_ENTRY_OPERATOR->used >= CONTAINER_SIZE - 1) {
		container = alloc_mapping_entry_container();
		if (!container)
			return NULL ;
		MAPPING_ENTRY_OPERATOR->next = container;
		MAPPING_ENTRY_OPERATOR = container;
	}
	
	entry = &(MAPPING_ENTRY_OPERATOR->entries[MAPPING_ENTRY_OPERATOR->used]);
	MAPPING_ENTRY_OPERATOR->used++;

	entry->hit = NULL;
	entry->pred = NULL;
	entry->succ = NULL;
	entry->readpos = -1;
	
	return(entry);
}

MAPPING_ENTRY_CONTAINER* alloc_mapping_entry_container()
{
	MAPPING_ENTRY_CONTAINER *container;

	if ((container =  (MAPPING_ENTRY_CONTAINER *) calloc(1, sizeof(MAPPING_ENTRY_CONTAINER))) == NULL) {
		//fprintf(stderr, "ERROR : not enough memory for mapping entry container memory (alloc_mapping_entry_container)\n");
		fprintf(stderr, "WARNING : not enough memory for mapping entry container memory (alloc_mapping_entry_container)\n");
		return NULL ;
	}

	container->used = 0;
	container->next = NULL;

	return(container);
}


CHROMOSOME_ENTRY* alloc_chromosome_entry(unsigned int pos, Chromosome const &chr, char strand)
{
	if (CHROMOSOME_ENTRY_OPERATOR->used > _config.CHROM_CONTAINER_SIZE - 1) {
		fprintf(stderr, "\n!!! WARNING: Chromosome container size of %d is too small! Hits for read %s cannot be reported any more!\n\n", _config.CHROM_CONTAINER_SIZE, _read.id());
		return(NULL);
	}
	
	CHROMOSOME_ENTRY *entry;

	entry = &(CHROMOSOME_ENTRY_OPERATOR->entries[CHROMOSOME_ENTRY_OPERATOR->used]);
	entry->chromosome = chr.nr();
	entry->genome_pos = pos;
	entry->strand = strand;
	entry->next = NULL;
    entry->mapping_entries = NULL;
    
    CHROMOSOME_ENTRY_OPERATOR->used++;

    if (_config.STATISTICS && CHROMOSOME_ENTRY_OPERATOR->used > MAX_USED_SLOTS)
    	MAX_USED_SLOTS = CHROMOSOME_ENTRY_OPERATOR->used;

	return(entry);
}

CHROMOSOME_ENTRY_CONTAINER* alloc_chromosome_entry_container()
{

    if ((CHROMOSOME_ENTRY_OPERATOR = (CHROMOSOME_ENTRY_CONTAINER *) calloc (1, sizeof(CHROMOSOME_ENTRY_CONTAINER))) == NULL) {
    	fprintf(stderr, "ERROR : not enough memory for chromosome entry container memory\n");
		exit(1);
	}

	if ((CHROMOSOME_ENTRY_OPERATOR->entries = (CHROMOSOME_ENTRY *) calloc (_config.CHROM_CONTAINER_SIZE, sizeof(CHROMOSOME_ENTRY))) == NULL) {
		fprintf(stderr, "ERROR : not enough memory for chromosome entry container memory\n");
		exit(1);
	}
	
	CHROMOSOME_ENTRY_OPERATOR->used = 0;

	return(CHROMOSOME_ENTRY_OPERATOR);
}

HIT* alloc_hit()
{
	HIT_CONTAINER *container;
	HIT *hit;

	if (HIT_OPERATOR->used >= CONTAINER_SIZE - 1) {
		container = alloc_hit_container();
		if (!container)
			return NULL ;
		HIT_OPERATOR->next = container;
		HIT_OPERATOR = container;
	}
	
	hit = &(HIT_OPERATOR->entries[HIT_OPERATOR->used]);
	HIT_OPERATOR->used++;

	new (hit) HIT();
	return(hit);
}

HIT_CONTAINER* alloc_hit_container()
{
	HIT_CONTAINER *container;

	if ((container =  (HIT_CONTAINER *) calloc (1, sizeof(HIT_CONTAINER))) == NULL) {
		fprintf(stderr, "WARNING: not enough memory for mapping entry container memory (alloc_hit_container)\n");
		return NULL ;
	}

	container->used = 0;
	container->next = NULL;

	return(container);
}

int alloc_hit_lists_operator()
{
	if ((HIT_LISTS_OPERATOR = (HIT **) calloc (_config.MAX_READ_LENGTH, sizeof(HIT*))) == NULL) {
		fprintf(stderr, "ERROR : not enough memory for hitlist (alloc_hit_lists_operator)\n");
		exit(1);
	}

	return(0);
}

int alloc_hits_by_score()
{
	double max_score = _config.NUM_GAPS * _config.GAP_SCORE + (_config.NUM_EDIT_OPS - _config.NUM_GAPS) * _config.MM_SCORE;
	NUM_SCORE_INTERVALS = max_score / SCORE_INTERVAL;
	if (NUM_SCORE_INTERVALS * SCORE_INTERVAL != max_score) ++NUM_SCORE_INTERVALS;
	NUM_SCORE_INTERVALS++;
	
	if ((HITS_BY_SCORE = (HITS_BY_SCORE_STRUCT *) calloc (NUM_SCORE_INTERVALS, sizeof(HITS_BY_SCORE_STRUCT))) == NULL) {
		fprintf(stderr, "ERROR : not enough memory for hitlist by score (alloc_hits_by_score)\n");
		exit(1);
	}

	int i;
	for (i=0; i!=NUM_SCORE_INTERVALS; ++i) {
		HITS_BY_SCORE[i].hitpointer = NULL;
		HITS_BY_SCORE[i].num = 0;
	}

	return(0);
}

int alloc_readstart_bins()
{
	if ((READSTART_BINS = (HIT **) calloc ((LONGEST_CHROMOSOME / ((_config.NUM_GAPS==0)? 1: _config.NUM_GAPS)), sizeof(HIT*) )) == NULL) {
		fprintf(stderr, "ERROR : not enough memory for readstart bin structure (alloc_readstart_bins)\n");
		exit(1);
	}

	return(0);
} 

int dealloc_mapping_entries()
{
	NUM_MAPPING_ENTRIES = 0;
	
	MAPPING_ENTRY_CONTAINER *container, *next;
	container = MAPPING_ENTRY_OPERATOR_FIRST->next;

	while (container != NULL) {
		next = container->next;
		NUM_MAPPING_ENTRIES += container->used;
		free(container);
		container = next;
	}

	NUM_MAPPING_ENTRIES += MAPPING_ENTRY_OPERATOR_FIRST->used;
	
	MAPPING_ENTRY_OPERATOR = MAPPING_ENTRY_OPERATOR_FIRST;
	MAPPING_ENTRY_OPERATOR->used = 0;
	MAPPING_ENTRY_OPERATOR->next = NULL;

	return(0);
}

int dealloc_hits()
{
	HIT_CONTAINER *container, *next;
	container = HIT_OPERATOR_FIRST->next;

	while (container != NULL) {
		next = container->next;
		free(container);
		container = next;
	}


	HIT_OPERATOR = HIT_OPERATOR_FIRST;
	HIT_OPERATOR->used = 0;
	HIT_OPERATOR->next = NULL;

	return(0);
}

int dealloc_chromosome_entries()
{
	unsigned int i;
	
	for (i=0; i < _config.CHROM_CONTAINER_SIZE; i++) {
		free(CHROMOSOME_ENTRY_OPERATOR->entries[i].mapping_entries);
		free(CHROMOSOME_ENTRY_OPERATOR->entries+i);
	}
	//free(CHROMOSOME_ENTRY_OPERATOR);

	return(0);
}

int dealloc_hit_lists_operator()
{
	unsigned int i;

	//for (i = _config.INDEX_DEPTH; i <= LONGEST_HIT; i++) {
	for (i = 0; i <= LONGEST_HIT; i++) 
	{
		*(HIT_LISTS_OPERATOR+i) = NULL;
	}
	
	return(0);
}

int dealloc_hits_by_score()
{
	unsigned int i;

	for (i = 0; i != NUM_SCORE_INTERVALS; i++) {
		HITS_BY_SCORE[i].hitpointer = NULL;
		HITS_BY_SCORE[i].num = 0;
	}

	return(0);
}

