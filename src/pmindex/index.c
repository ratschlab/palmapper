// Authors: Gunnar Raetsch, Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include <palmapper/VariantMap.h>
#include "pmindex.h"
#include <iterator>


int pos2bin(unsigned int slot, unsigned int chr);
int pos2bin_rev(unsigned int slot, unsigned int chr);
unsigned int get_slot(const char *seq, int pos);

const bool perform_extra_checks=true ;


inline int insert_variants(std::vector<Variant>::iterator start, std::vector<Variant>::iterator end, char * seq, unsigned int pos, int chr, 
						   int counter, std::map<int,bool> & source_ids, std::map<unsigned int,bool> & slots)  
{
	if (counter < 0)
		return 0 ;

	std::vector<Variant>::iterator it ;
	for (it=start; it<end; it++)
	{
		if ((*it).type == pt_SNP)
		{
			int r1=0, r2=0 ;
			bool found = false, used=false ;

			for (unsigned int j=0; j<(*it).read_id_num.size(); j++)
			{
				found = (source_ids.count((*it).read_id_num[j])>0) ;
				if (found)
				{
					used = source_ids[(*it).read_id_num[j]] ;
					break ;
				}
			}
			
			if (found)
			{
				if (!used)
					r2+=insert_variants( it+1, end, seq, pos, chr, counter, source_ids, slots) ;
				else
				{
					char c = (*it).variant_str[0] ;
					if (c=='A' || c=='C' || c=='G' || c=='T')
					{
						if (perform_extra_checks)
							assert((*it).position-(int)pos>=0 && (*it).position-(int)pos<INDEX_DEPTH) ;
						
						char o=seq[(*it).position-pos] ;
						seq[(*it).position-pos] = c ;
						
						r1+=insert_variants( it+1, end, seq, pos, chr, counter, source_ids, slots) ;
						
						// cleanup
						seq[(*it).position-pos] = o ;
					}
				}
			}
			else
			{
				if (counter<1)
					continue ; // spare preparing the sequence

				for (unsigned int j=0; j<(*it).read_id_num.size(); j++)
				{
					source_ids[(*it).read_id_num[j]]=false ;
					r2+=insert_variants( it+1, end, seq, pos, chr, counter, source_ids, slots) ;
					
					char c = (*it).variant_str[0] ;
					if ((c=='A' || c=='C' || c=='G' || c=='T'))
					{
						source_ids[(*it).read_id_num[j]]=true ;
						if (perform_extra_checks)
							assert((*it).position-(int)pos>=0 && (*it).position-(int)pos<INDEX_DEPTH) ;
						
						char o=seq[(*it).position-pos] ;
						seq[(*it).position-pos] = c ;
						
						r1+=insert_variants( it+1, end, seq, pos, chr, counter-1, source_ids, slots) ;
						
						// cleanup
						seq[(*it).position-pos] = o ;
					}
					
					source_ids.erase((*it).read_id_num[j]) ;
				}
			}
			
			return r1+r2 ;
		}
	}
	
	if (it==end)
	{
		if (false && slots.size()>=2)
		{
			fprintf(stdout, "variant at pos %i\nsources\n", pos) ;
			for (std::map<int,bool>::iterator it2=source_ids.begin(); it2!=source_ids.end(); it2++)
				fprintf(stdout, "%i\n", it2->first) ;
		}
		
		HAS_SLOT = 0;
		
		unsigned int slot = get_slot(seq, 0);
		slots[slot]=true ;

		HAS_SLOT = 0 ;
		
		return 1 ;
	}
	
	return 0 ;
}

int index_chromosome_novariants(unsigned int chr, Genome & genome, GenomeMaps & genome_mask, bool mask_do_alloc, bool mask_do_secondary, bool mask_do_add) 
{
	if (VERBOSE) { printf("\tBuilding index ..."); fflush(stdout); }

	unsigned int pos = 0;
	int spacer = 0;
	unsigned int slot = 0;
	POS p;
	int num_seeds_total=0 ;
	int num_positions_total=0 ;
	HAS_SLOT = 0;

	while (spacer < (int)CHR_LENGTH) 
	{ 
		if (spacer < (int)pos + INDEX_DEPTH - 1) 
		{
			if (CHR_SEQ[spacer]=='A' || CHR_SEQ[spacer]=='T' || CHR_SEQ[spacer]=='C' || CHR_SEQ[spacer]=='G') {
				spacer++;
			}
			else {
				spacer++;
				POSITION += (unsigned short int)(spacer - pos);
				pos = spacer;
				HAS_SLOT = 0;
			}
		}
		else {
			if (CHR_SEQ[spacer]=='A' || CHR_SEQ[spacer]=='T' || CHR_SEQ[spacer]=='C' || CHR_SEQ[spacer]=='G') 
			{
				slot = get_slot(CHR_SEQ, pos);
				
				if (!has_genome_mask)
				  {
				    if(INDEX[slot] == NULL) 
					alloc_bin(slot);
				    assert(INDEX[slot]!=NULL) ;
				    pos2bin(slot, chr);	// 0-initialized
				    
				    num_seeds_total++ ;
				    num_positions_total++ ;
				  }
				else
				  {
				    char elem=genome_mask.CHR_MAP(genome.chromosome(chr), pos) ;
				    if (mask_do_alloc && ((elem & MASK_REGION_PRIMARY)>0))
				      {
					if(INDEX[slot] == NULL) 
					    alloc_bin(slot);
					pos2bin(slot, chr);	// 0-initialized

					num_seeds_total++ ;
					num_positions_total++ ;
				      }
				    if (mask_do_secondary && INDEX[slot] != NULL)
				      {
					if ((elem & MASK_REGION_SECONDARY) == 0)
					  {
					    elem+=MASK_REGION_SECONDARY ;
					    genome_mask.CHR_MAP_set(genome.chromosome(chr), pos, elem) ;
					  }
					if (!mask_do_alloc)
					  {
					    num_seeds_total++ ;
					    num_positions_total++ ;
					  }
				      }
				    if (mask_do_add && 
					( /* ((elem & MASK_REGION_PRIMARY)>0) ||  */
					  ((elem & MASK_REGION_SECONDARY)>0) || ((elem & MASK_REGION_SECONDARY_REGION)>0)))
				      {
					if (INDEX[slot] == NULL) 
					  alloc_bin(slot);
					pos2bin(slot, chr);	// 0-initialized
					
					if (!mask_do_alloc && !mask_do_secondary)
					  {
					    num_seeds_total++ ;
					    num_positions_total++ ;
					  }
				      }
				  }
				POSITION++;
				HAS_SLOT = 1;
				spacer++;
				pos++;
			}
			else {
				spacer++;
				POSITION += spacer - pos;
				pos = spacer;
				HAS_SLOT = 0;
			}
		}
		
		if (POSITION >= BLOCK_SIZE) 
		{
			if (BLOCK == BLOCK_TABLE_SIZE - 1) 
			{
				fprintf(stderr, "ERROR: Too large chromosomes/contigs or too many chromosomes/contigs! Split input file into many smaller ones!\n");
				exit(0);
			}
			BLOCK++;
			POSITION %= BLOCK_SIZE;
			
			p.chr = chr;
			p.pos = pos - POSITION;
			BLOCK_TABLE[BLOCK] = p;
		}
	}
	if (VERBOSE)
	  fprintf(stdout, "found %i seeds at %i positions (chr_len=%i; without variants)", num_seeds_total,  num_positions_total, CHR_LENGTH) ;
	if (VERBOSE) 
	  printf("... done\n");

	return 0;
}

struct thread_data
{
	unsigned int chr ;
	int from ;
	int to ;
	VariantMap * variants ;
	std::map<int, std::vector<unsigned int> > chr_slots  ;
	int thread_id ;
} ;

std::map<int, std::vector<unsigned int> > get_slots_from_chromosome(unsigned int chr, int from, int to, VariantMap & variants, int thread_id)  ;

void * get_slots_from_chromosome_wrapper(void * data)
{
	struct thread_data *pdata=(struct thread_data *) data ;
	
	pdata->chr_slots = get_slots_from_chromosome(pdata->chr, pdata->from, pdata->to, *(pdata->variants), pdata->thread_id) ;

	return data ;
}


int index_chromosome(unsigned int chr, Genome & genome, VariantMap & variants, GenomeMaps & genome_mask, bool mask_do_alloc, bool mask_do_secondary, bool mask_do_add) 
{
	if (VERBOSE) { printf("\tBuilding index ..."); fflush(stdout); }
	unsigned int pos = 0;
	int spacer = 0;
	unsigned int slot = 0;
	POS p;
	
	HAS_SLOT = 0;
	assert(variants.genome!=NULL) ;

	// verify genome sequence
	assert(variants.genome->chromosome(chr).length()==CHR_LENGTH) ;
	for (unsigned int i=0; i<CHR_LENGTH; i++)
		assert(variants.genome->chromosome(chr)[i]==CHR_SEQ[i]) ;

	// compute slots (in parallel)
	std::map<int, std::vector<unsigned int> > chr_slots ;
#ifdef NILL
	{
		const int num_threads=NUM_THREADS ;
		struct thread_data thread_data[num_threads] ;
		pthread_t thread_desc[num_threads] ;
		int step = (CHR_LENGTH/num_threads)+1 ;
		
		if (num_threads>1)
			fprintf(stdout, "starting %i threads\n", num_threads) ;
		int from=0;
		for (int i=0; i<num_threads; i++)
		{
			assert(from<=(int)CHR_LENGTH) ;
			int to=from+step+INDEX_DEPTH ;
			if (to>(int)CHR_LENGTH)
				to=CHR_LENGTH ;
			
			thread_data[i].chr=chr ;
			thread_data[i].from=from ;
			thread_data[i].to=to ;
			thread_data[i].variants=&variants ;
			thread_data[i].thread_id = i ;
			//fprintf(stdout, "starting thread %i: %i-%i\n", i, from, to) ;
			
			if (num_threads==1)
				get_slots_from_chromosome_wrapper(&thread_data[i]) ;
			else
				pthread_create(&thread_desc[i], NULL, get_slots_from_chromosome_wrapper, &thread_data[i]);

			from+=step ;
		}
		if (num_threads>1)
		{
			for (int i=0; i<num_threads; i++)
			{
				void *value_ptr ;
				pthread_join(thread_desc[i], &value_ptr);
			}
			fprintf(stdout, "%i threads finished\n", num_threads) ;
		}

		for (int i=0; i<num_threads; i++)
		{
			std::map<int, std::vector<unsigned int> >::iterator it=thread_data[i].chr_slots.begin() ;
			for (; it!=thread_data[i].chr_slots.end(); it++)
				chr_slots[it->first]=it->second ;
			thread_data[i].chr_slots.clear() ;
		}
	}
#endif

	int num_seeds_total = 0 ;
	int num_positions_total = 0 ;

	int step = 1000000 ;
	struct thread_data slot_data ;
		
	int from = 0 ;
	int to=0 ;
	while (spacer < (int)CHR_LENGTH) 
	{ 

		if (spacer < (int)pos + INDEX_DEPTH - 1) 
		{
			if (CHR_SEQ[spacer]=='A' || CHR_SEQ[spacer]=='T' || CHR_SEQ[spacer]=='C' || CHR_SEQ[spacer]=='G') {
				spacer++;
			}
			else {
				spacer++;
				POSITION += spacer - pos;
				pos = spacer;
			}
		}
		else {
			if (CHR_SEQ[spacer]=='A' || CHR_SEQ[spacer]=='T' || CHR_SEQ[spacer]=='C' || CHR_SEQ[spacer]=='G') 
			{
			  if ((signed)pos>=to-INDEX_DEPTH)
			    {
			      from=pos ;
			      to=from+step+INDEX_DEPTH ;
			      fprintf(stdout, "producing slots for %i-%i: ", from, to) ;
			      if (to>(int)CHR_LENGTH)
				to=CHR_LENGTH ;
			      slot_data.chr=chr ;
			      slot_data.from=from ;
			      slot_data.to=to ;
			      slot_data.variants=&variants ;
			      slot_data.thread_id = 0 ;
			      
			      get_slots_from_chromosome_wrapper(&slot_data) ;
			      
			      std::map<int, std::vector<unsigned int> >::iterator it=slot_data.chr_slots.begin() ;
			      for (; it!=slot_data.chr_slots.end(); it++)
				chr_slots[it->first]=it->second ;
			      slot_data.chr_slots.clear() ;
			    }
			  
			  assert(chr_slots.count(pos)>0) ;
			  std::vector<unsigned int> & slots = chr_slots[pos] ;
			  bool used_pos=false ;

				if (!has_genome_mask)
				  {
				    for (std::vector<unsigned int>::iterator it=slots.begin(); it != slots.end(); it++)
				      {
					slot = (*it) ;
					if(INDEX[slot] == NULL) 
					  alloc_bin(slot);
					pos2bin(slot, chr);	// 0-initialized
					num_seeds_total++ ;
					used_pos=true ;
				      }
				  }
				else
				  {
				    char elem=genome_mask.CHR_MAP(genome.chromosome(chr), pos) ;
				    for (std::vector<unsigned int>::iterator it=slots.begin(); it != slots.end(); it++)
				      {
					slot = (*it) ;
					if (mask_do_alloc && ((elem & MASK_REGION_PRIMARY)>0))
					  {
					    if(INDEX[slot] == NULL)
						alloc_bin(slot);
					    pos2bin(slot, chr);   // 0-initialized
					    num_seeds_total++ ;
					    used_pos=true ;
					  }
					if (mask_do_secondary && INDEX[slot] != NULL)
					  {
					    if ((elem & MASK_REGION_SECONDARY) == 0)
					      {
						elem+=MASK_REGION_SECONDARY ;
						genome_mask.CHR_MAP_set(genome.chromosome(chr), pos, elem) ;
					      }
					    if (!mask_do_alloc)
					      num_seeds_total++ ;
					  }
					if (mask_do_add && 
					    ( /* ((elem & MASK_REGION_PRIMARY)>0) || */ 
					     ((elem & MASK_REGION_SECONDARY)>0) || ((elem & MASK_REGION_SECONDARY_REGION)>0)) )
					  {
					    if(INDEX[slot] == NULL) 
					      alloc_bin(slot);
					    pos2bin(slot, chr);   // 0-initialized                                                                                                                                                                                                                                                                   
					    used_pos=true ;

					    if (!mask_do_alloc && !mask_do_secondary)
						num_seeds_total++ ;
					  }
				      }
				  }
				POSITION++;
				spacer++;
				pos++;
				if (used_pos)
				  num_positions_total++ ;
			}
			else {
				spacer++;
				POSITION += spacer - pos;
				pos = spacer;
			}
		}
		assert( POSITION%BLOCK_SIZE == pos%BLOCK_SIZE ) ;
		
		if (POSITION >= BLOCK_SIZE) 
		{
			if (BLOCK == BLOCK_TABLE_SIZE - 1) 
			{
				fprintf(stderr, "ERROR: Too large chromosomes/contigs or too many chromosomes/contigs! Split input file into many smaller ones!\n");
				exit(0);
			}
			BLOCK++;
			POSITION %= BLOCK_SIZE;
			
			p.chr = chr;
			p.pos = pos - POSITION;
			BLOCK_TABLE[BLOCK] = p;
		}
	}
	fprintf(stdout, "found %i seeds at %i positions (chr_len=%i; with variants)", num_seeds_total,  num_positions_total, CHR_LENGTH) ;

	if (VERBOSE) printf("... done\n");

	return 0;
}


std::map<int, std::vector<unsigned int> > get_slots_from_chromosome(unsigned int chr, int from, int to, VariantMap & variants, int thread_id) 
{
	unsigned int pos = from ;
	int spacer = from ;
	unsigned int slot = 0;

	std::map<int, std::vector<unsigned int> > chr_slots ;
	
	HAS_SLOT = 0;
	assert(variants.genome!=NULL) ;

	assert(variants.genome->chromosome(chr).length()==CHR_LENGTH) ;
	for (unsigned int i=0; i<CHR_LENGTH; i++)
		assert(variants.genome->chromosome(chr)[i]==CHR_SEQ[i]) ;

	int num_seeds_total = 0 ;
	int num_positions_total = 0 ;

	std::vector<Variant>::iterator curr_start = variants.variantlist[chr].begin() ;
	std::vector<Variant>::iterator curr_stop = curr_start ;

	int num_considered_variants = 0 ; 
	while (spacer < to) 
	{ 
		
		while (curr_start != variants.variantlist[chr].end() && (int)(*curr_start).position < (int)pos) 
			advance(curr_start, 1) ;
		curr_stop = curr_start ;
		while (curr_stop != variants.variantlist[chr].end() && (*curr_stop).position < (int)pos+INDEX_DEPTH)
			advance(curr_stop, 1) ;

		if (spacer < (int)pos + INDEX_DEPTH - 1) 
		{
			if (CHR_SEQ[spacer]=='A' || CHR_SEQ[spacer]=='T' || CHR_SEQ[spacer]=='C' || CHR_SEQ[spacer]=='G') {
				spacer++;
			}
			else {
				spacer++;
				pos = spacer;
			}
		}
		else {
			if (CHR_SEQ[spacer]=='A' || CHR_SEQ[spacer]=='T' || CHR_SEQ[spacer]=='C' || CHR_SEQ[spacer]=='G') 
			{
				if (distance(curr_start, curr_stop)>0)
				{
					char seq[INDEX_DEPTH] ;
					strncpy(seq, &CHR_SEQ[pos], INDEX_DEPTH) ;
					
					std::map<int,bool> source_ids ;
					std::map<unsigned int,bool> slots ;
					
					num_considered_variants+=curr_stop-curr_start ;
					int num_seeds=insert_variants(curr_start, curr_stop, seq, pos, chr, MAX_SOURCE_COMBINATIONS, source_ids, slots) ;
					int num_slots=slots.size() ;
					
					if (perform_extra_checks)
						for (unsigned int i=pos; i<pos+INDEX_DEPTH; i++)
							assert(seq[i-pos]==CHR_SEQ[i]) ;
					
					if (slots.size()>100)
						fprintf(stdout, "num_seeds=%i, num_slots=%i, #variants=%ld\n", num_seeds, num_slots, distance(curr_start, curr_stop)) ;
					num_seeds_total += num_slots ;
					
					for (std::map<unsigned int,bool>::iterator it=slots.begin(); it != slots.end(); it++)
					{
						slot = it->first ;
						chr_slots[pos].push_back(slot) ;
					}
				}
				else
				{
					//HAS_SLOT=0 ;
					
					slot = get_slot(CHR_SEQ, pos);
					chr_slots[pos].push_back(slot) ;
					
					HAS_SLOT = 1;
					num_seeds_total++ ;
				}
				num_positions_total++ ;
				spacer++;
				pos++;
			}
			else {
				spacer++;
				pos = spacer;
				HAS_SLOT = 0;
			}
		}
	}
	/* fprintf(stdout, "distance to start to end: %ld\n", distance(variants.variantlist[chr].begin(), variants.variantlist[chr].end())) ;
	fprintf(stdout, "distance to curr_start to curr_end: %ld\n", distance(curr_start, curr_stop)) ;
	fprintf(stdout, "distance to start: %ld\n", distance(variants.variantlist[chr].begin(), curr_start)) ;
	fprintf(stdout, "distance to end: %ld\n", distance(curr_stop, variants.variantlist[chr].end())) ; */
	fprintf(stdout, "[%i: %2.1f%%, %i, %i, %i: %i]\n ", thread_id, 100.0*(double)(pos-from)/(double)(to-from), spacer-from, num_positions_total, num_seeds_total, num_considered_variants) ;

	return chr_slots;
}


int pos2bin(unsigned int slot, unsigned int chr) 
{
	BIN_EXT **bin_ext;
	BIN *bin;
	unsigned int num;

	bin = INDEX[slot];
	num = bin->num_pos;

	POSITION_COUNTER++;

	if (num == 0) {
		USED_SLOTS[NUM_USED_SLOTS] = slot;
		NUM_USED_SLOTS++;
	}

	if (num < BIN_SIZE) {
		memcpy(&(bin->ids[num].id[0]), &BLOCK, 3 * sizeof(char));
		memcpy(&(bin->ids[num].id[3]), &POSITION, sizeof(unsigned char));
		
		bin->num_pos++;
	}
	else {
		bin_ext = &(bin->bin_ext);
		if (*bin_ext == 0) {
			*bin_ext = alloc_bin_ext();
			memcpy(&(*bin_ext)->ids[0].id, &BLOCK, 3 * sizeof(char));
			memcpy(&(*bin_ext)->ids[0].id[3], &POSITION, sizeof(unsigned char));
			bin->num_pos++;
			bin->last_bin_ext = *bin_ext;
			return 0;
		}
		else {
			bin_ext = &(bin->last_bin_ext);
			if ((num % BIN_SIZE_EXT) != BIN_SIZE) {
				memcpy(&(*bin_ext)->ids[(num-BIN_SIZE) % BIN_SIZE_EXT].id, &BLOCK, 3 * sizeof(char));
				memcpy(&(*bin_ext)->ids[(num-BIN_SIZE) % BIN_SIZE_EXT].id[3], &POSITION, sizeof(unsigned char));
				bin->num_pos++;
				return 0;
			}
			else  {
				bin_ext = &((*bin_ext)->bin_ext);
				*bin_ext = alloc_bin_ext();
				memcpy(&(*bin_ext)->ids[0].id, &BLOCK, 3 * sizeof(char));
				memcpy(&(*bin_ext)->ids[0].id[3], &POSITION, sizeof(unsigned char));
				bin->num_pos++;
				bin->last_bin_ext = *bin_ext;
			}
		}
	}

	return(0);
}


unsigned int get_slot(const char *seq, int pos)
{
	unsigned int slot = 0;
	unsigned int i;
	int c = 0;
	
	if (HAS_SLOT == 0) 
	{
		for (i=0; (int)i<INDEX_DEPTH; i++) 
		{
			if (seq[pos+i] == 'A') {
				c = 0;
			}
			else {
				if (seq[pos+i] == 'C') {
					c = 1;
				}
				else {
					if (seq[pos+i] == 'G') {
						c = 2;
					}
					else {
						if (seq[pos+i] == 'T') {
							c = 3;
						}
						else {
							fprintf(stderr, "ERROR: Reached unallowed character in calculateSlot %c at position %d \n", seq[pos+i], (pos+i));
							exit(1);
						}
					}
				}
			}
			slot = slot + (unsigned int)(POWER[i] * c);
		}
	}
	else 
	{

		slot = SLOT;
		slot >>= 2;

		if (seq[pos+INDEX_DEPTH-1] == 'A') {
			slot = slot | (unsigned int)(BINARY_CODE[0]);
		}
		else if (seq[pos+INDEX_DEPTH-1] == 'C') {
			slot = slot | (unsigned int)(BINARY_CODE[1]);
		}
		else if (seq[pos+INDEX_DEPTH-1] == 'G') {
			slot = slot | (unsigned int)(BINARY_CODE[2]);
		}
		else { //if (seq[pos+INDEX_DEPTH-1] == 'T') {
			slot = slot | (unsigned int)(BINARY_CODE[3]);
		}

	}

	SLOT = slot;

	return(slot);
}
