// Authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include <genomemapper/Util.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <iostream>

#include <genomemapper/Hit.h>
#include <genomemapper/Read.h>
#include <genomemapper/genomemapper.h>

using namespace std;

extern clock_t time2a_seek ;
extern clock_t time2a_seed2genome  ;

int get_slot(int pos);

Read::Read() {
	READ_ID = new char[_config.MAX_READ_ID_LENGTH];
	for (int i = 0; i < 3; ++i) {
		READ_QUALITY[i] = (char*) calloc(_config.MAX_READ_LENGTH + 1, sizeof(char));
		if (READ_QUALITY[i] == NULL) {
			fprintf(stderr, "[init_defaults] Could not allocate memory\n");
			exit(1);
		}
	}
	linenr = 0;
	READ_LENGTH = 0;
}

Read::Read(const Read& read) 
{
	READ_ID = new char[_config.MAX_READ_ID_LENGTH];
	for (int i = 0; i < 3; ++i) {
		READ_QUALITY[i] = (char*) calloc(_config.MAX_READ_LENGTH + 1, sizeof(char));
		if (READ_QUALITY[i] == NULL) {
			fprintf(stderr, "[init_defaults] Could not allocate memory\n");
			exit(1);
		}
	}
	linenr = read.linenr ;
	READ_LENGTH = read.READ_LENGTH ;
	strncpy(READ_ID, read.READ_ID, _config.MAX_READ_ID_LENGTH) ;
	strncpy(READ, read.READ, _config.MAX_READ_LENGTH+1) ;
	for (int i = 0; i < 3; ++i) {
		strncpy(READ_QUALITY[i], read.READ_QUALITY[i], _config.MAX_READ_LENGTH + 1) ;
	}
	READ_FORMAT=read.READ_FORMAT ;
	READ_PE_FLAG = read.READ_PE_FLAG ;
	
}

/** Parses one line of read descriptions in either FASTA, FASTQ or FLATFILE
 *  format and sets a number of global variables describing this read.
 *
 *	The rest of the GenomeMapper code always operates on the current read.
 *	Information on this read is available through a number of global variables
 *	that are set in this routine.
 *
 *	\return A status code; 0 on success
 */
int Read::read_short_read(FILE *QUERY_FP)
{
	char line[10000];
	char *tmp;
	int linelen;

	if (fgets(line, 10000, QUERY_FP) == NULL) {
		if (READ_LENGTH == 0)
			cerr << "\n!!! WARNING: Input read file '" << /*_config.*/_config.QUERY_FILE_NAME << "' is empty!\n\n";
		return 1;
	}
	++linenr;

	if (strcspn(line, " \n\t") == 0) {
		do {
			if (fgets(line, 10000, QUERY_FP) == NULL) {
				if (READ_LENGTH == 0)
					cerr << "\n!!! WARNING: Input read file '" << /*_config.*/_config.QUERY_FILE_NAME << "' is empty!\n\n";
				return 1;
			}
			++linenr;
		} while (strcspn(line, " \n\t") == 0);
	}

	linelen = strlen(line);
	if (linelen < 3) {
		cerr << "ERROR: Unknown read input format! Do all the reads have an identifier?\n";
		exit(0);
	}

	if (line[0] == '@') {
		/////// FastQ input ///////

		// R E A D _ I D
		memset(READ, 0, /*_config.*/_config.MAX_READ_LENGTH) ;
		memset(READ_QUALITY[0], 0, /*_config.*/_config.MAX_READ_LENGTH) ;
		memset(READ_ID, 0, /*_config.*/_config.MAX_READ_ID_LENGTH) ;
		strncpy(READ_ID, line+1, strcspn(line, " \t\n")-1);

		do {
			if (fgets(line, 10000, QUERY_FP) == NULL) {
				fprintf(stderr, "ERROR: Read '%s' in line %lu is not complete in input query file '%s'! Missing read sequence and quality!\n", READ_ID, linenr, /*_config.*/_config.QUERY_FILE_NAME.c_str());
				exit(0);
			}
			++linenr;
		} while (strcspn(line, " \t\n") == 0);

		// R E A D
		strncpy(READ, line, strcspn(line, " \t\n"));
		//READ[36]=0 ;
		if (strlen(READ) > /*_config.*/_config.MAX_READ_LENGTH) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu is longer than the max read length (=%zu)! It will be omitted!\n\n", READ_ID, linenr, /*_config.*/_config.MAX_READ_LENGTH);
			return -1;
		}
		else if (strlen(READ) == 0) {
			fprintf(stderr, "ERROR: Cannot find read sequence of read '%s' in line %lu in input query file '%s'!\n", READ_ID, linenr, /*_config.*/_config.QUERY_FILE_NAME.c_str());
			exit(0);
		}
		if (strcspn(READ, "aAcCgGtTnNrRyYmMkKwWsSbBdDhHvV") != 0) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu contains non-IUPAC characters! It will be omitted!\n\n", READ_ID, linenr);
			return -1;
		}

		if ((int)strlen(READ) < /*_config.*/_config.INDEX_DEPTH) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu is shorter than the specified seedlength! It will be omitted!\n\n", READ_ID, linenr);
			return -1;
		}

		do {
			if (fgets(line, 10000, QUERY_FP) == NULL) {
				fprintf(stderr, "ERROR: Read '%s' in line %lu is not complete in input query file '%s'! Missing quality!\n", READ_ID, linenr, /*_config.*/_config.QUERY_FILE_NAME.c_str());
				exit(0);
			}
			++linenr;
		} while (strcspn(line, " \t\n") == 0);

		// +
		if (strlen(line) < 1 || line[0] != '+') {
			fprintf(stderr, "ERROR: Read '%s' in line %lu is not in fastq format!\n", READ_ID, linenr);
			exit(0);
		}

		do {
			if (fgets(line, 10000, QUERY_FP) == NULL) {
				fprintf(stderr, "ERROR: Read '%s' in line %lu is not complete in input query file '%s'! Missing quality!\n", READ_ID, linenr, /*_config.*/_config.QUERY_FILE_NAME.c_str());
				exit(0);
			}
			++linenr;
		} while (strcspn(line, " \t\n") == 0);

		// Q U A L I T Y
		if (strlen(line) > 0)
			strncpy(READ_QUALITY[0], line, strcspn(line, " \t\n"));
		else {
			fprintf(stderr, "ERROR: Cannot find read quality of read '%s' in line %lu in input query file '%s'!\n", READ_ID, linenr, /*_config.*/_config.QUERY_FILE_NAME.c_str());
			exit(0);
		}

		/*if (strlen(READ_QUALITY[0]) != READ_LENGTH) {
			fprintf(stderr, "ERROR: Read quality 1 of read '%s' in line %lu hasn't length of read!\n", READ_ID, linenr);
			exit(0);
		}*/

		READ_LENGTH = strlen(READ);

		// O T H E R
		READ_PE_FLAG = 0;
		READ_QUALITY[1]=(char*)"" ;
		READ_QUALITY[2]=(char*)"";

		READ_FORMAT = 0;

	}
	else if (line[0] == '>') {
		/////// Fasta input ///////
		memset(READ, 0, /*_config.*/_config.MAX_READ_LENGTH) ;
		memset(READ_ID, 0, /*_config.*/_config.MAX_READ_LENGTH) ;

		strncpy(READ_ID, line+1, strcspn(line, " \t\n")-1);

		do {
			if (fgets(line, 10000, QUERY_FP) == NULL) {
				fprintf(stderr, "ERROR: Read '%s' in line %lu is not complete in input query file '%s'! Missing read sequence!\n", READ_ID, linenr, /*_config.*/_config.QUERY_FILE_NAME.c_str());
				exit(0);
			}
			++linenr;
		} while (strcspn(line, " \t\n") == 0);

		// R E A D
		strncpy(READ, line, strcspn(line, " \t\n"));
		if (strlen(READ) > /*_config.*/_config.MAX_READ_LENGTH) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu is longer than the max read length (=%zu)! It will be omitted!\n\n", READ_ID, linenr, /*_config.*/_config.MAX_READ_LENGTH);
			return -1;
		}
		else if (strlen(READ) == 0) {
			fprintf(stderr, "ERROR: Cannot find read sequence of read '%s' in line %lu in input query file '%s'!\n", READ_ID, linenr, /*_config.*/_config.QUERY_FILE_NAME.c_str());
			exit(0);
		}
		for (int i=0; i<(int)strlen(READ); i++)
			if (READ[i]=='-')
			{
				//fprintf(stderr, "replaced '-' with 'N'\n")  ;
				READ[i]='N' ;
			}

		if (strcspn(READ, "aAcCgGtTnNrRyYmMkKwWsSbBdDhHvV") != 0) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu contains non-IUPAC characters! It will be omitted!\n\n", READ_ID, linenr);
			return -1;
		}

		if ((int)strlen(READ) < /*_config.*/_config.INDEX_DEPTH) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu is shorter than the specified seedlength! It will be omitted!\n\n", READ_ID, linenr);
			return -1;
		}

		READ_LENGTH = strlen(READ);

		READ_PE_FLAG = 0;
		strcpy(READ_QUALITY[0], READ) ;
		memset(READ_QUALITY[0], 'h', strlen(READ)) ;

		//READ_QUALITY[0] = (char*)"";
		READ_QUALITY[1] = (char*)"";
		READ_QUALITY[2] = (char*)"";

		READ_FORMAT = 1;

	}
	else {
		/////// Flatfile input ///////
		READ_ID = strtok(line, "\t");

		if ((int)strlen(READ_ID) == linelen) {
			fprintf(stderr, "ERROR: wrong read input data format, line %lu!\n", linenr);
			exit(0);
		}
		if (READ_ID == NULL) {
			fprintf(stderr, "ERROR: Read ID is empty, line %lu!\n", linenr);
			exit(0);
		}

		char *tok = strtok(NULL, "\t");
		if (tok == NULL) {
			fprintf(stderr, "ERROR: Read sequence is empty, line %lu!\n", linenr);
			exit(0);
		}

		if (strlen(tok) > /*_config.*/_config.MAX_READ_LENGTH) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu is longer than the max read length (=%zu)! It will be omitted!\n\n", READ_ID, linenr, /*_config.*/_config.MAX_READ_LENGTH);
			return -1;
		}
		//printf("%s sp: %d\n",READ,(int) strcspn(READ, "A"));
		if (strcspn(tok, "aAcCgGtTnNrRyYmMkKwWsSbBdDhHvV") != 0) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu contains non-IUPAC characters! It will be omitted!\n\n", READ_ID, linenr);
			return -1;
		}

		strcpy(READ, tok);
		READ_LENGTH = strlen(tok);
		if ((int)READ_LENGTH < /*_config.*/_config.INDEX_DEPTH) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu is shorter than the specified seedlength! It will be omitted!\n\n", READ_ID, linenr);
			return -1;
		}

		tmp = strtok(NULL, "\t");
		if (tmp == NULL) {
			fprintf(stderr, "ERROR: Paired-end flag is empty, line %lu!\n", linenr);
			exit(0);
		}
		READ_PE_FLAG = atoi(tmp);

		READ_QUALITY[0] = strtok('\0', "\t");

		//fprintf(stderr, "hack!!!\n") ;
		//for (int i=0; i<strlen(READ_QUALITY[0]); i++)
		//	READ_QUALITY[0][i]='h' ;

		if (READ_QUALITY[0] == NULL) {
			fprintf(stderr, "ERROR: Read Quality 1 is empty, line %lu!\n", linenr);
			exit(0);
		}
		/*if (strlen(READ_QUALITY[0]) != READ_LENGTH) {
			fprintf(stderr, "ERROR: Read quality 1 hasn't length of read, line %lu!\n", linenr);
			exit(0);
		}*/

		// TODO it is very questionable to reassign the READ_QUALITY pointer since other branches of
		// this method function are memset()ting on it...
		READ_QUALITY[1] = strtok(NULL, "\t");
		if (READ_QUALITY[1] == NULL) {
			fprintf(stderr, "ERROR: Read Quality 2 is empty, line %lu!\n", linenr);
			exit(0);
		}
		/*if (strlen(READ_QUALITY[1]) != READ_LENGTH) {
			fprintf(stderr, "ERROR: Read quality 2 hasn't length of read, line %lu!\n", linenr);
			exit(0);
		}*/

		READ_QUALITY[2] = strtok(NULL, "\n");
		if (READ_QUALITY[2] == NULL) {
			fprintf(stderr, "ERROR: Read Quality 3 is empty, line %lu!\n", linenr);
			exit(0);
		}
		/*if (strlen(READ_QUALITY[2]) != READ_LENGTH) {
			fprintf(stderr, "ERROR: Read quality 3 hasn't length of read, line %lu!\n", linenr);
			exit(0);
		}*/

		READ_FORMAT = 2;
	}

	return 0;
}

int Read::map_fast(int & firstslot, int & firstpos)
{
#ifndef BinaryStream_MAP
	STORAGE_ENTRY *index_mmap ;
#else
	CBinaryStream<STORAGE_ENTRY>* index_mmap;
#endif
	INDEX_ENTRY index_entry;

	unsigned int pos, i, j, p, chars, block, chrom_overlap = 0, hits_reported = 0;
	int slot, run, nr_mms, readpos, chrpos, readstart, chrstart;
	firstslot = -1 ; firstpos=0 ;

	unsigned char position;
	char nr_runs, rev, nr_seeds = READ_LENGTH / _config.INDEX_DEPTH, mm, perfect = 0, cancel=0;

	if (/*_config.*/_config.NUM_GAPS != 0)
		nr_runs = (nr_seeds > 1)? 2: 1;
	else {
		if (nr_seeds > /*_config.*/_config.NUM_MISMATCHES) nr_runs = /*_config.*/_config.NUM_MISMATCHES + 1;
			else nr_runs = nr_seeds;
	}

	bool seed_already_inspected_fwd[nr_runs+1];
	bool seed_already_inspected_rev[nr_runs+1];

	int max_mms = nr_runs - 1;
	int mmpos[max_mms];
	//fprintf(stdout, "nr_runs=%i\n", nr_runs) ;

	for (run=1; run<=nr_runs; ++run) {

		if (nr_runs == 1) nr_runs = 0;	// a bit fishy, but nr_runs and run only have to be different, thats why this assignment is due
		if (run == nr_runs) readstart = ((int)READ_LENGTH) - /*_config.*/_config.INDEX_DEPTH;
			else		 	readstart = (run-1) * /*_config.*/_config.INDEX_DEPTH;

		slot = get_slot(readstart) ;

		//fprintf(stdout, "map_fast: slot=%i HAS_SLOT=%i read_start=%i seq=%s\n", slot, HAS_SLOT, readstart, get_seq(slot)) ;

		if (run==1)
		{
			if (slot<0)
			{
				for (int ii=1; ii<((int)READ_LENGTH)-/*_config.*/_config.INDEX_DEPTH && firstslot<0; ii++)
				{
					firstpos = ii ;
					firstslot = get_slot(firstpos) ;
				}
			}
			else
			{
				firstslot = slot ;
				//firstpos = 0 ;
			}
		}

		// does not work when one of the first symbols are non-ACGT.
		//if (run == 1)
		//	firstslot = slot;


		if (slot >= 0)
		{	// tests if slot has an unallowed char!

			for (rev=0; rev <= /*_config.*/_config.MAP_REVERSE; ++rev) {
				if (!rev) {
					index_entry = INDEX[slot];
					index_mmap = INDEX_FWD_MMAP;
				}
				else  {
					index_entry = INDEX_REV[slot];
					index_mmap = INDEX_REV_MMAP;
				}

				// for each mapping position
				if (index_entry.num) {

					STORAGE_ENTRY* se_buffer=NULL;

					try
					{
						se_buffer=new STORAGE_ENTRY[index_entry.num];
					}
					catch (std::bad_alloc)
					{
						return -1;
					}

					int index_entry_num=index_entry.num ;
					if ((int)index_entry.num > _config.SEED_HIT_CANCEL_THRESHOLD) { // && !REPORT_REPETITIVE_SEEDS)
						index_entry_num=0 ;
						if (rev) seed_already_inspected_rev[run] = false;
						else seed_already_inspected_fwd[run] = false;
					}
					else {
						if (rev) seed_already_inspected_rev[run] = true;
						else seed_already_inspected_fwd[run] = true;
					}

					time_t start_time = clock() ;
#ifndef BinaryStream_MAP
					index_pre_buffer(index_mmap, se_buffer, index_entry.offset-index_entry_num, index_entry_num);
#else
					index_mmap->pre_buffer(se_buffer, index_entry.offset-index_entry_num, index_entry_num);
#endif
					time2a_seek += clock()-start_time ;

					for (i=0; (int)i<index_entry_num; i++)
					{

						block = 0; position = 0;
						unsigned char* p_block = (unsigned char*) &block;
						STORAGE_ENTRY se;
						unsigned char* p_id;

						if (!rev) {
							//memcpy(&block, &((index_mmap+(index_entry.offset-(i+1)))->id[0]), 3 * sizeof(char));
							//memcpy(&position, &((index_mmap+(index_entry.offset-(i+1)))->id[3]), sizeof(unsigned char));
							se = se_buffer[index_entry.num-(i+1)];
						}
						else {
							//memcpy(&block, &((index_mmap+(index_entry.offset-(index_entry.num-i)))->id[0]), 3 * sizeof(char));
							//memcpy(&position, &((index_mmap+(index_entry.offset-(index_entry.num-i)))->id[3]), sizeof(unsigned char));
							se = se_buffer[i];//index_entry.offset-(index_entry.num-i)];
						}

						p_id=se.id;
						p_block[0]=p_id[0];
						p_block[1]=p_id[1];
						p_block[2]=p_id[2];
						position = p_id[3];
						pos = (unsigned int) position + BLOCK_TABLE[block].pos;	// 0-initialized
						Chromosome &chr = _genome.chromosome(BLOCK_TABLE[block].chr);

						//if (_config.REPORT_REPETITIVE_SEEDS)
						//	report_repetitive_seed(chr, pos, index_entry.num)  ;

						//if (index_entry.num>_config.SEED_HIT_CANCEL_THRESHOLD)
						//	continue ;

						if (!rev) {
							chrstart = pos - (run!=nr_runs) * (run-1) * /*_config.*/_config.INDEX_DEPTH - (run==nr_runs) * (((int)READ_LENGTH) - /*_config.*/_config.INDEX_DEPTH);
						}	// 0-initialized
						else {
							chrstart = pos + (run!=nr_runs) *   run   * /*_config.*/_config.INDEX_DEPTH + (run==nr_runs) * ((int)READ_LENGTH) - 1;
						}

						// check if read can map on position in genome:
						if ( (!rev && chrstart < 0) ||
								(rev && chrstart < ((int)READ_LENGTH) - 1) ) {
							if (/*_config.*/_config.STATISTICS) chrom_overlap++;
						}
						else if ( (!rev && chrstart + ((int)READ_LENGTH) > (int)chr.length()) ||
								  ( rev && chrstart > (int)chr.length() - 1) ) {
							if (/*_config.*/_config.STATISTICS) chrom_overlap++;
						}
						else {

							nr_mms = 0;
							chars = 0;
							cancel = 0;

							readpos = 0;
							chrpos = chrstart;

							for (j=1; (int)j!=run; ++j) {

								mm = 0;

								for (p=0; (int)p!=/*_config.*/_config.INDEX_DEPTH; ++p) {

									READ[readpos+p] = mytoupper(READ[readpos+p]);
									if ( ( rev && get_compl_base(chr[chrpos-p]) != READ[readpos+p]) ||
											(!rev &&                chr[chrpos+p]  != READ[readpos+p]) ||
											!(unique_base(READ[readpos+p])) ) {

										if (nr_mms < max_mms) {
											mmpos[nr_mms] = readpos + p + 1;
											++nr_mms;
										}
										else {
											cancel = 1;
											break;
										}
										++mm;

									}

									++chars;

								}

								if (!mm && ((rev && seed_already_inspected_rev[j]) || (!rev && seed_already_inspected_fwd[j]))) {
									cancel = 1;
									break;
								}

								if (cancel) break;

								chrpos += /*_config.*/_config.INDEX_DEPTH * (rev? -1: 1);
								readpos += /*_config.*/_config.INDEX_DEPTH;

							}


							chrpos  = chrpos  + (run!=nr_runs) * /*_config.*/_config.INDEX_DEPTH * (rev? -1: 1);
							readpos = readpos + (run!=nr_runs) * /*_config.*/_config.INDEX_DEPTH;
							while (!cancel && (int)chars != ((int)READ_LENGTH) - /*_config.*/_config.INDEX_DEPTH) {

								READ[readpos] = mytoupper(READ[readpos]);
								if ( ( rev && get_compl_base(chr[chrpos]) != READ[readpos]) ||
										(!rev && 				 chr[chrpos]  != READ[readpos]) ||
										!(unique_base(READ[readpos])) )
								{
									if (nr_mms == max_mms) {
										cancel = 1;
										break;
									}
									mmpos[nr_mms] = readpos + 1;
									++nr_mms;
								}

								readpos++;
								chrpos += rev? -1: 1;
								chars++;

							}


							if ( !cancel && nr_mms <= max_mms ) {
								// create hit
								HIT* hit = new HIT();
								if (!hit)
								{
									delete[] se_buffer;
									return -1 ;
								}

								hit->chromosome = &chr;
								hit->readpos = 1;

								if (!rev) {
									hit->orientation = '+';
									hit->start = chrstart + 1;					// 1-initialized
									hit->end = chrstart + ((int)READ_LENGTH);
								}
								else {
									hit->orientation = '-';
									hit->start = chrstart - ((int)READ_LENGTH) + 2;	// 1-initialized
									hit->end = chrstart + 1;
								}

								mm = 0;
								// create possible mismatches
								for (j=0; (int)j!=nr_mms; ++j) {
									hit->edit_op[j].mm = 1;
									if (hit->orientation == '+') hit->edit_op[j].pos = mmpos[j];
									else					 hit->edit_op[j].pos = ((int)READ_LENGTH) - mmpos[j] + 1;
									assert(hit->edit_op[j].pos >= -((int)READ_LENGTH) && hit->edit_op[j].pos<=(int)READ_LENGTH) ;
									hit->mismatches++;
									mm = 1;
								}

								//if (!_config.ALL_HIT_STRATEGY && nr_mms < max_mms)
								//	max_mms = nr_mms;

								update_num_edit_ops(nr_mms, _config.ALL_HIT_STRATEGY, max_mms) ;

								// perfect matching read
								if (/*_config.*/_config.STATISTICS) {
									if (!mm) {
										if (rev) /*_stats.*/_stats.PERFECT_HITS_REV++;
										else /*_stats.*/_stats.PERFECT_HITS++;
										if (!perfect) /*_stats.*/_stats.PERFECT_READS++;

										perfect = 1;
									}
									else /*_stats.*/_stats.NOT_ALIGNED[1]++;
								}

								int ret = insert_into_scorelist(hit, 0);
								assert(ret>=0) ;

								if (/*_config.*/_config.STATISTICS)	_stats.HITS_LEN[READ_LENGTH]++;

								hits_reported++;

							} // end of create hit

						} // end of no hit-overlap with chrom border

					} // end of for each mapping pos

					delete[] se_buffer;
				} // end of index entry num

			} // end of forward/reverse	rev

		} // end of slot != -1

	} // end of runs = different slots

	if (!/*_config.*/_config.ALL_HIT_STRATEGY && !hits_reported)
	{	//if best hit strategy, but no mappings found -> prepare for complete mapping!
		/*_config.*/_config.ALL_HIT_STRATEGY = -1;
	}
	else
	{
		if (/*_config.*/_config.STATISTICS)
		{
			/*_stats.*/_stats.NUM_HITS += hits_reported;
			/*_stats.*/_stats.HITS_PER_READ += hits_reported;
			/*_stats.*/_stats.ENDSTART_MAPPED[0] += chrom_overlap;
		}
	}
	//fprintf(stdout, "map_fast: firstslot=%i, firstpos=%i\n", firstslot, firstpos) ;

	return 0 ;
}


int Read::map_short_read(unsigned int num, int first_slot, int first_pos)
{
	char reverse;
	unsigned int readpos = first_pos ;
	unsigned int spacer = first_pos;
	unsigned int slot;

	//fprintf(stdout, "first_slot=%i\n", first_slot) ;

	if (first_slot < 0) {	// first slot has an unallowed char!
		return -1 ;
		/*assert(0) ; //obsoleted code
		spacer = -first_slot+1 +1 ;
		readpos = -first_slot+1 + 1 ;
		HAS_SLOT = 0;*/
	}
	else if (/*_config.*/_config.ALL_HIT_STRATEGY && !/*_config.*/_config.SUMMARY_HIT_STRATEGY) {	// first slot hasn't been computed yet
		spacer = 0;
		slot = 0;
		HAS_SLOT = 0;
	}
	else {	// first slot has already been computed in map_fast and doesn't contain unallowed chars
		slot = first_slot;

		reverse = 0;
		if (INDEX[slot].num != 0)
			reverse = 1;
		if (/*_config.*/_config.MAP_REVERSE && INDEX_REV[slot].num != 0)
			reverse = (reverse + reverse) + 2;

		if (reverse > 0)
		{
			clock_t seed2genome_start=clock() ;
			//fprintf(stdout, "1: slot=%i\n", slot) ;
			int ret = seed2genome(num, slot, readpos + 1, reverse) ;
			time2a_seed2genome += clock() - seed2genome_start ;

			if (ret<0)
			{
				fprintf(stderr, "seed2genome<0 (1)\n")  ;
				return ret ; // add one to readpos for 1-initialization
			}
		}

		readpos++;
		spacer = /*_config.*/_config.INDEX_DEPTH;
		SLOT = first_slot;
		HAS_SLOT = 1;
	}

	while (spacer < READ_LENGTH)
	{
		READ[spacer] = mytoupper(READ[spacer]);
		if (spacer < readpos + /*_config.*/_config.INDEX_DEPTH - 1)
		{
			if (READ[spacer]=='A' || READ[spacer]=='T' || READ[spacer]=='C' || READ[spacer]=='G')
			{
				spacer++;
			}
			else
			{
				spacer++;
				readpos = spacer;
				HAS_SLOT = 0;
			}
		}
		else {
			if (READ[spacer]=='A' || READ[spacer]=='T' || READ[spacer]=='C' || READ[spacer]=='G')
			{
				slot = get_slot(readpos);
				//fprintf(stdout, "map_short_read: slot=%i HAS_SLOT=%i readpos=%i seq=%s\n", slot, HAS_SLOT, readpos, get_seq(slot)) ;

				if (slot<0)
				{
					spacer++;
					readpos++;
					continue ;
				}

				// reverse: 0: slot doesnt match in either index or index_rev, 1: only index, 2: only index_rev, 4: both
				reverse = 0;
				//if (INDEX[slot].last_entry != NULL) {
				if (INDEX[slot].num != 0) {
					reverse = 1;
				}
				if (/*_config.*/_config.MAP_REVERSE && INDEX_REV[slot].num != 0) {
					reverse = (reverse + reverse) + 2;
				}

				if (reverse > 0) {
					clock_t seed2genome_start=clock() ;
					int ret=seed2genome(num, slot, readpos + 1, reverse) ;
					time2a_seed2genome += clock() - seed2genome_start ;
					if (ret<0)
					{
						fprintf(stderr, "seed2genome<0 (2)\n") ;
						return ret ; // add one to readpos for 1-initialization
					}
				}

				spacer++;
				readpos++;
				HAS_SLOT = 1;
			}
			else {
				spacer++;
				readpos=spacer;
				HAS_SLOT = 0;
			}
		}
	}

	HAS_SLOT = 0;

	return 1;
}


/**
 *
 *  \param pos Position in the current read for which to get a slot
 *  \return A slot
 */

int Read::get_slot(int pos)
{
    unsigned int slot = 0;
	unsigned int i;
	int c = 0;

	if (HAS_SLOT == 0)
	{

		for (i = 0; (int)i < /*_config.*/_config.INDEX_DEPTH; i++)
		{
			READ[pos + i] = mytoupper(READ[pos + i]);

			switch (READ[pos + i])
			{
			case 'A':
				c = 0;
				break;
			case 'C':
				c = 1;
				break;
			case 'G':
				c = 2;
				break;
			case 'T':
				c = 3;
				break;
			default:
				return -pos - i - 1 ; // subtract -1 to make sure the output is indeed negative (pos=i=0..?)
			}

			slot = slot + Util::POWER[i] * c;
		}

	} else {

		slot = SLOT;
		slot >>= 2;

		READ[pos + /*_config.*/_config.INDEX_DEPTH - 1] = mytoupper(READ[pos + /*_config.*/_config.INDEX_DEPTH - 1]);

		switch (READ[pos + /*_config.*/_config.INDEX_DEPTH - 1]) {
		case 'A':
			slot = slot | BINARY_CODE[0];
			break;
		case 'C':
			slot = slot | BINARY_CODE[1];
			break;
		case 'G':
			slot = slot | BINARY_CODE[2];
			break;
		case 'T':
			slot = slot | BINARY_CODE[3];
			break;
		default:
			return -1;
		}
	}

	SLOT = slot;

	return (slot);
}
