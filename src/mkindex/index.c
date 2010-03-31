// Authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include "mkindex.h"

int pos2bin(unsigned int slot, unsigned int chr);
int pos2bin_rev(unsigned int slot, unsigned int chr);
int get_slot(char *seq, int pos);

int index_chromosome(unsigned int chr) 
{
	if (VERBOSE) { printf("\tBuilding index ..."); fflush(stdout); }

	unsigned int pos = 0;
	int spacer = 0;
	int slot = 0;
	POS p;
	
	HAS_SLOT = 0;

	while (spacer < CHR_LENGTH) {
		
		if (spacer < pos + INDEX_DEPTH - 1) {
			if (CHR_SEQ[spacer]=='A' || CHR_SEQ[spacer]=='T' || CHR_SEQ[spacer]=='C' || CHR_SEQ[spacer]=='G') {
				spacer++;
			}
			else {
				spacer++;
				POSITION += spacer - pos;
				pos = spacer;
				HAS_SLOT = 0;
			}
		}
		else {
			if (CHR_SEQ[spacer]=='A' || CHR_SEQ[spacer]=='T' || CHR_SEQ[spacer]=='C' || CHR_SEQ[spacer]=='G') {

				slot = get_slot(CHR_SEQ, pos); // yields also SLOT_REV if -r option wasn't set

				if(INDEX[slot] == NULL) {
					alloc_bin(slot);
					if (BUILD_REVERSE_INDEX) alloc_bin_rev(SLOT_REV);
				}
				pos2bin(slot, chr);	// 0-initialized
				
				if (BUILD_REVERSE_INDEX) {
					pos2bin_rev(SLOT_REV, chr); // 0-initialized
				}
				
				spacer++;
				pos++;
				POSITION++;
				HAS_SLOT = 1;
			}
			else {
				spacer++;
				POSITION += spacer - pos;				
				pos = spacer;
				HAS_SLOT = 0;
			}
		}
		
		if (POSITION >= BLOCK_SIZE) {
			if (BLOCK == BLOCK_TABLE_SIZE - 1) {
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

	if (VERBOSE) printf("... done\n");

	return 0;
}

int pos2bin(unsigned int slot, unsigned int chr) 
{
	BIN_EXT **bin_ext;
	BIN *bin, *binrev;
	unsigned int num, numrev;

	bin = INDEX[slot];
	//bin->num_all++;
	num = bin->num_pos;

	numrev = 0;
	if (BUILD_REVERSE_INDEX && INDEX_REV[slot] != NULL) {
		binrev = INDEX_REV[slot];
		numrev = binrev->num_pos;
	}

	POSITION_COUNTER++;

	if (num == 0 && numrev == 0) {
		USED_SLOTS[NUM_USED_SLOTS] = slot;
		NUM_USED_SLOTS++;
	}

	if (num == 0) {
		SLOT_COUNTER++;
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


int pos2bin_rev(unsigned int slot, unsigned int chr) 
{
	BIN_EXT **bin_ext;
	BIN *bin, *binfwd;
	unsigned int num, numfwd;
	
	bin = INDEX_REV[slot];
	//bin->num_all++;
	num = bin->num_pos;

	numfwd = 0;
	if (INDEX[slot] != NULL) {
		binfwd = INDEX[slot];
		numfwd = binfwd->num_pos;
	}

	if (num == 0 && numfwd == 0) {
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

int get_slot(char *seq, int pos)
{
	unsigned int slot = 0;
	unsigned int i;
	int c = 0;
	
	if (HAS_SLOT == 0) { 
		SLOT_REV = 0;
		for (i=0; i<INDEX_DEPTH; i++) {
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
			slot = slot + POWER[i] * c;
			if (BUILD_REVERSE_INDEX) SLOT_REV += POWER[INDEX_DEPTH - i - 1] * (c ^ 3);
		}
	}
	else {

		slot = SLOT;
		slot >>= 2;

		if (BUILD_REVERSE_INDEX) {
			SLOT_REV <<= 34 - INDEX_DEPTH * 2;
			SLOT_REV >>= 32 - INDEX_DEPTH * 2;
		}

		if (seq[pos+INDEX_DEPTH-1] == 'A') {
			slot = slot | BINARY_CODE[0];
			if (BUILD_REVERSE_INDEX) SLOT_REV |= 3;
		}
		else if (seq[pos+INDEX_DEPTH-1] == 'C') {
			slot = slot | BINARY_CODE[1];
			if (BUILD_REVERSE_INDEX) SLOT_REV |= 2;
		}
		else if (seq[pos+INDEX_DEPTH-1] == 'G') {
			slot = slot | BINARY_CODE[2];
			if (BUILD_REVERSE_INDEX) SLOT_REV |= 1;
		}
		else { //if (seq[pos+INDEX_DEPTH-1] == 'T') {
			slot = slot | BINARY_CODE[3];
		}

	}

	SLOT = slot;

	return(slot);
}
