// Authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include "genomemapper.h"
#include <assert.h>
#include <vector>

static const int MAX_EXON_LEN=100 ;

int num_unspliced_best=0 ;
int num_unspliced_suboptimal=0 ;
int num_spliced_best=0 ;
int num_spliced_suboptimal=0 ;
int num_unmapped=0 ;
clock_t last_spliced_report = 0 ;

void printhit(HIT* hit) ;
int compare_int(const void *a, const void *b) ;
void construct_aligned_string(HIT *hit, int *num_gaps_p, int *num_mismatches_p, int*num_matches_p) ;
int report_read_alignment(HIT* hit, int num)  ;

void print_stats() {
	printf("\n########### _config.STATISTICS ###########\n");
	printf("Mapped Reads: %i of all %d reads\n", _stats.READS_MAPPED, _stats.NUM_READS);
	int i;
	for (i = 0; i != _config.NUM_MISMATCHES + 1; ++i)
		printf(" Reads with %d mismatches: %u\n", i, _stats.HITS_MM[i]);
	printf(
			"  Perfect Plus-Hits: %i\n  Perfect Minus-Hits:\t%i\n  Total:\t\t%i\n",
			_stats.PERFECT_HITS, _stats.PERFECT_HITS_REV, _stats.PERFECT_HITS + _stats.PERFECT_HITS_REV);
	printf("   Perfect matching reads (+ or -): %i\n", _stats.PERFECT_READS);

	for (i = _config.INDEX_DEPTH; i != ASSUMED_READ_LENGTH; ++i)
		printf("    Total Hits of length %d: %u - %.4f%%\n", i, _stats.HITS_LEN[i],
				((double) (100 * _stats.HITS_LEN[i]) / _stats.NUM_HITS));
	printf("\n    Total number of hits:                  %lu\n", _stats.NUM_HITS);
	printf("      Number of overhang alignments:       %lu\n", _stats.NUM_ALIGNMENTS);
	printf("      Successful overhang alignments:      %lu\n", _stats.NUM_ALIGNMENTS
			- _stats.GAPS_ENCOUNTERED[0] - _stats.GAPS_ENCOUNTERED[1] - _stats.GAPS_ENCOUNTERED[2]
			- _stats.TOO_MANY_MMS[0] - _stats.TOO_MANY_MMS[1]);
	printf("      Number of global alignments:         %lu\n",
			_stats.NUM_WHOLE_ALIGNMENTS);
	printf("      Successful global alignments:        %lu\n",
			_stats.NUM_WHOLE_ALIGNMENTS - _stats.BREAK_GLOBAL_ALIGNMENT[0]
					- _stats.BREAK_TB_IN_GLOBAL_ALIGNMENT);
	printf("      Breaks in global alignments(no hit): %lu\n",
			_stats.BREAK_GLOBAL_ALIGNMENT[0]);
	printf("      Breaks in global alignments(hit):    %lu\n",
			_stats.BREAK_GLOBAL_ALIGNMENT[1]);
	printf("      Breaks after global alignment:       %lu\n\n",
			_stats.BREAK_TB_IN_GLOBAL_ALIGNMENT);
	printf("      Hits of len %d, not aligned:         %d\n", ASSUMED_READ_LENGTH - 1, _stats.NOT_ALIGNED[0]);
	printf("      Mappings of fast mapping:            %d\n", _stats.NOT_ALIGNED[1]);
	//printf("      hits of len 36, not aligned (with 1MM at start/beg of read): %d\n",_stats.NOT_ALIGNED[1]);
	//printf("        Number of unmapped hits at beg/end of chrom: %d and in align-step: %d\n", _stats.ENDSTART_MAPPED[0], _stats.ENDSTART_MAPPED[1]);
	printf("      Overlapping hits at beg/end of chr:  %d\n",
			_stats.ENDSTART_MAPPED[0]);

	printf("\nGaps encountered on diagonal during overhang-alignment:   %lu\n",
			_stats.GAPS_ENCOUNTERED[0]);
	printf("Gaps encountered due to score during overhang-alignment:  %lu\n",
			_stats.GAPS_ENCOUNTERED[1]);
	printf("Gaps encountered after overhang-alignment:                %lu\n",
			_stats.GAPS_ENCOUNTERED[2]);
	printf("Too many MMs during overhang-alignment (discarded hits):  %lu\n",
			_stats.TOO_MANY_MMS[0]);
	printf("Too many MMs after overhang-alignment (discarded hits):   %lu\n\n",
			_stats.TOO_MANY_MMS[1]);

	printf("cells processed in overhang: %lu\n", _stats.CELLS_OVERHANG);
	printf("cells processed in global:   %lu\n", _stats.CELLS_GLOBAL);
	printf("total cells:                 %lu\n\n", _stats.CELLS_OVERHANG
			+ _stats.CELLS_GLOBAL);

	printf("gap limit overcome: %lu\n\n", _stats.W);

	//printf("reads filtered out due to too much non-base chars: %i\n\n",READS_FILTERED);

	/*printf("Print in stat... ");
	 unsigned int sum=0;
	 FILE* f;
	 f = fopen("stat","w");
	 for (i=0; i!= _config.MAX_READ_LENGTH-_config.INDEX_DEPTH+1; ++i) {
	 for (j=0; j!=_config.MAX_READ_LENGTH - _config.INDEX_DEPTH - i + 1; ++j) {
	 if (*(*(_stats.HITS_READPOS+i)+j) != 0) {
	 fprintf(f, "\t(%d)%u ", j+1, *(*(_stats.HITS_READPOS+i)+j));
	 sum += *(*(_stats.HITS_READPOS+i)+j);
	 }
	 }
	 fprintf(f, "\nLENGTH %d: %d", i+_config.INDEX_DEPTH, sum);
	 if (i!=_config.MAX_READ_LENGTH-_config.INDEX_DEPTH && sum != 0 && i<36) fprintf(f, "\tfirst: %.4f%%, last: %.4f%%, other: %.4f%%\n", ((double) (100 * *(*(_stats.HITS_READPOS+i))) / sum), ((double) (100 * (*(*(_stats.HITS_READPOS+i)+(36 - _config.INDEX_DEPTH - i))))/sum), 100 - ((double) (100 * *(*(_stats.HITS_READPOS+i))) / sum) - ((double) (100 * (*(*(_stats.HITS_READPOS+i)+(36 - _config.INDEX_DEPTH - i))))/sum));
	 else fprintf(f, "\n");
	 sum = 0;
	 }
	 printf("done\n");*/
}

/*int different_hit(HIT* hit1, HIT* hit2)
{
	//fprintf(stdout, "diff %i %i\n", hit1->start, hit2->start) ;
	
	if (abs(hit1->start-hit2->start)<(int)READ_LENGTH)
		return 0 ;
	return 1 ;
}*/

// prints out all hits which have been inserted into HITS_BY_EDITOPS
// called once for each read (?)
int print_hits() {
	int i, printed = 0, nr;
	HIT *hit;

	//HIT *best_hit = NULL;
	//u_int32_t num_best_hits = 0;
	//vector<HIT *> found_hits ;
	//u_int32_t num_2nd_best_hits = 0;
	//u_int32_t num_other_hits = 0;
	int reported_reads = 0 ;

	for (i = 0; i != NUM_SCORE_INTERVALS; ++i) {

		if (printed && !_config.ALL_HIT_STRATEGY && !_config.SUMMARY_HIT_STRATEGY)
			break; // best hit strategy

		if (HITS_BY_SCORE[i].hitpointer != NULL) 
		{
			// only _config.REPEATMAP numbers of alignment will be chosen randomly:
			if (!_config.ALL_HIT_STRATEGY && !_config.SUMMARY_HIT_STRATEGY && _config.REPEATMAP < 0 && HITS_BY_SCORE[i].num > -_config.REPEATMAP) 
			{
				srand((unsigned) time(NULL));
				
				int j, k, n;
				int hits[-_config.REPEATMAP];
				for (j = 0; j != -_config.REPEATMAP; ++j) {
					n = 1;
					while (n != 0) {
						n = 0;
						hits[j] = rand() % HITS_BY_SCORE[i].num;
						for (k = 0; k != j; ++k) {
							if (hits[j] == hits[k])
								++n;
						}
					}
				}

				qsort(hits, -_config.REPEATMAP, sizeof(int), compare_int);

				hit = HITS_BY_SCORE[i].hitpointer;

				nr = 0;
				for (j = 0; j != HITS_BY_SCORE[i].num; ++j) {

					if (hits[nr] == j) {
						printed += print_alignment(hit, -_config.REPEATMAP);
						nr++;
					}

					if (nr == -_config.REPEATMAP)
						break;

					hit = hit->same_eo_succ;
				}

			} else if (_config.SUMMARY_HIT_STRATEGY) {

				hit = HITS_BY_SCORE[i].hitpointer;

				// iterate over all hits for this score and collect data to generate a summary
				// This somewhat counter-intuitive code works because of the way the pre-existing
				// code was set up: we see the hits in the order of their score here - better hits first.
				while (hit != NULL) 
				{
					_topalignments.report_unspliced_hit(hit) ;
					hit = hit->same_eo_succ;
				}

			} else { // no random selection of output alignments:

				hit = HITS_BY_SCORE[i].hitpointer;

				while (hit != NULL) {

					if (!_config.ALL_HIT_STRATEGY)
						nr = HITS_BY_SCORE[i].num;
					else
						nr = HITS_IN_SCORE_LIST;

					if (_config.REPEATMAP == 0) { // no max nr of hits per read was specified, print all
						printed += print_alignment(hit, nr);
					} else if (_config.REPEATMAP > 0 && printed < _config.REPEATMAP) {
						printed += print_alignment(hit, (nr < _config.REPEATMAP) ? nr : _config.REPEATMAP);
					} else if (_config.REPEATMAP == printed) { // repeatmap many alignments already printed out -> stop printing -> next read
						return 1;
					}

					hit = hit->same_eo_succ;
				}
			}

			if (_config.REPORT_MAPPED_READS)
			{
				hit = HITS_BY_SCORE[i].hitpointer;

				while (hit != NULL) {

					nr = HITS_IN_SCORE_LIST;

					printed += report_read_alignment(hit, reported_reads);
					reported_reads++ ;
					hit = hit->same_eo_succ;
				}
			}
		}
	}

	if (printed != 0)
		return 1; // read could have been mapped
	else
		return 0; // read couldn't be mapped
}

int print_alignment(HIT* hit, unsigned int num) {

	if (_config.STATISTICS)
		_stats.HITS_MM[hit->mismatches]++;

	int j, fstart, flen;

	int hitlength = hit->end - hit->start + 1;
	unsigned int readstart;
	if (hit->orientation == '+') {
		readstart = hit->start - hit->readpos + hit->start_offset; // start pos of read in genome	0-initialized
	} else {
		readstart = hit->start - (((int)_read.length()) - hit->readpos - hitlength + 2)
				+ hit->start_offset; // 0-initialized
	}

	// PERFECT HITS:
	if (hit->mismatches == 0) {

		if (hit->orientation == '+')
			strcpy(ALIGNSEQ, _read.data());
		else
		{
			for (size_t i=0; i<_read.length(); i++)
				ALIGNSEQ[i] = (*hit->chromosome)[readstart+i];
			//strncpy(ALIGNSEQ, CHR_SEQ[hit->chromosome] + readstart, ((int)_read.lenght()));
		}
		ALIGNSEQ[((int)_read.length())] = '\0';

		// print in file:
		if (_config.OUTPUT_FORMAT == 0) {
			/////// SHORE file ///////

			fprintf(OUT_FP, "%s\t%d\t%s\t%s\t%c", hit->chromosome->desc(),
					readstart + 1, // 1-initialized
					ALIGNSEQ,//Alignment String
					_read.id(), ((hit->orientation == '+') ? 'D' : 'P'));

			if (_config.SCORES_OUT)
				fprintf(OUT_FP, "\t%.1f", (double) (-((int)_read.length()) * _config.M_SCORE));
			else
				fprintf(OUT_FP, "\t%d", ((int)_read.length()) - hit->mismatches);
			// lt - changed this to report num matches instead of num mismatches

			fprintf(OUT_FP, "\t%d\t%d\t%d", num, // Number of hits
					((int)_read.length()), 0);

			if (_read.format() == 2)
				fprintf(OUT_FP, "\t%d", _read.pe_flag());
			if (strlen(_read.quality()[0]) != 0)
				fprintf(OUT_FP, "\t%s", _read.quality()[0]);
			if (strlen(_read.quality()[1]) != 0)
				fprintf(OUT_FP, "\t%s", _read.quality()[1]);
			if (strlen(_read.quality()[2]) != 0)
				fprintf(OUT_FP, "\t%s", _read.quality()[2]);

			if (_config.FLANKING != 0) {
				fstart = (readstart < _config.FLANKING) ? 0 : readstart - _config.FLANKING;
				flen
						= (readstart + ((int)_read.length()) + _config.FLANKING
								> hit->chromosome->length()) ? hit->chromosome->length()
								- fstart
					: readstart - fstart + ((int)_read.length()) + _config.FLANKING;
				{
					for (size_t i=0; i<flen; i++)
						FLANK_SEQ[i]=(*hit->chromosome)[fstart + i];
					//strncpy(FLANK_SEQ, CHR_SEQ[hit->chromosome] + fstart, flen);
				}
				FLANK_SEQ[flen] = '\0';
				fprintf(OUT_FP, "\t%d\t%s", hit->chromosome->length(),
						FLANK_SEQ);
			} else if (_config.PRINT_SEQ > 0)
				fprintf(OUT_FP, "\t%d", hit->chromosome->length());
			//if (_config.PRINT_SEQ == 2)
			//	fprintf(OUT_FP, "\t%s", CHR_SEQ[hit->chromosome]);

			fprintf(OUT_FP, "\n");
		} else {
			/////// BED file ///////

			if (_config.SCORES_OUT) {
				fprintf(OUT_FP, "%s\t%d\t%d\t%s\t%.1f\t%c\n",
						hit->chromosome->desc(), readstart, readstart + 1
								+ ((int)_read.length()), _read.id(), (double) (-((int)_read.length())
								* _config.M_SCORE), hit->orientation);
			} else {
				fprintf(OUT_FP, "%s\t%d\t%d\t%s\t%d\t%c\t0\t0\t0,0,0\t0\t0\t0\n",
						hit->chromosome->desc(), readstart, readstart + 1
								+ ((int)_read.length()), _read.id(), ((int)_read.length()) - hit->mismatches,
						hit->orientation);
				// lt  - changed this to report num matches instead of num mismatches
			}

		}

	}
	// HITS WITH MISMATCHES:
	else {

		int count_char = 0;
		char gap_offset = 0;
		char gap_in_read = 0;
		char gap_in_chr = 0;

		// sort mismatches in ascending order according to their abs positions and 'gap before mm'-strategy if equal
		qsort(hit->edit_op, hit->mismatches, sizeof(EDIT_OPS), compare_editops);

		ALIGNSEQ[0] = '\0';

		for (j = 0; j != hit->mismatches; ++j) {

			if (hit->edit_op[j].pos < 0) {
				hit->edit_op[j].pos = -hit->edit_op[j].pos;
				gap_in_chr = 1;
			}

			if (j == 0) 
			{
				for (size_t i=0; i< hit->edit_op[0].pos - 1; i++)
					ALIGNSEQ[i]=(*hit->chromosome)[readstart + i];
				//strncpy(ALIGNSEQ, CHR_SEQ[hit->chromosome] + (readstart),
				//		hit->edit_op[0].pos - 1);
				count_char += hit->edit_op[0].pos - 1;
			} else if (hit->edit_op[j].pos - hit->edit_op[j - 1].pos != 0) 
			{
				for (size_t i=0; i<hit->edit_op[j].pos - hit->edit_op[j - 1].pos - 1 + gap_in_read; i++)
					ALIGNSEQ[count_char+i] = (*hit->chromosome)[(readstart + hit->edit_op[j - 1].pos + gap_offset - gap_in_read) + i];
				/*strncpy(ALIGNSEQ + count_char, CHR_SEQ[hit->chromosome]
						+ (readstart + hit->edit_op[j - 1].pos + gap_offset
								- gap_in_read), hit->edit_op[j].pos
						- hit->edit_op[j - 1].pos - 1 + gap_in_read); // -1???
				*/
				count_char += hit->edit_op[j].pos - hit->edit_op[j - 1].pos - 1
						+ gap_in_read;
			} // else: edit_op[j-1] must have been a gap!

			gap_in_read = 0;

			if (hit->edit_op[j].mm) {
				if (hit->orientation == '+')
					sprintf(ALIGNSEQ + count_char, "[%c%c]",
							(*hit->chromosome)[readstart+ hit->edit_op[j].pos - 1 + gap_offset],
							_read.data()[hit->edit_op[j].pos - 1]);
				else
					sprintf(ALIGNSEQ + count_char, "[%c%c]",
							(*hit->chromosome)[readstart + hit->edit_op[j].pos - 1 + gap_offset],
							get_compl_base(_read.data()[((int)_read.length())
									- hit->edit_op[j].pos]));
			} else if (gap_in_chr) {
				if (hit->orientation == '+')
					sprintf(ALIGNSEQ + count_char, "[-%c]",
							_read.data()[hit->edit_op[j].pos - 1]);
				else
					sprintf(ALIGNSEQ + count_char, "[-%c]", get_compl_base(
							_read.data()[((int)_read.length()) - hit->edit_op[j].pos]));

				gap_offset--;
				gap_in_chr = 0;
			} else {
				sprintf(ALIGNSEQ + count_char, "[%c-]",
						(*hit->chromosome)[readstart + hit->edit_op[j].pos - 1 + gap_offset]);

				gap_offset++;
				gap_in_read = 1;
			}

			count_char += 4;

		}

		// from last mismatch to end of read:
		{
			for (size_t i=0; i<((int)_read.length()) - hit->edit_op[j - 1].pos + gap_in_read; i++)
				ALIGNSEQ[count_char+i] = (*hit->chromosome)[(readstart + hit->edit_op[j - 1].pos + gap_offset - gap_in_read) + i];
			/*strncpy(ALIGNSEQ + count_char, CHR_SEQ[hit->chromosome] + (readstart
				+ hit->edit_op[j - 1].pos + gap_offset - gap_in_read),
					  ((int)_read.lenght()) - hit->edit_op[j - 1].pos + gap_in_read);*/
		}
		count_char += ((int)_read.length()) - hit->edit_op[j - 1].pos + gap_in_read;

		ALIGNSEQ[count_char] = '\0';

		// print in file:
		if (_config.OUTPUT_FORMAT == 0) {

			/////// SHORE file ///////

			fprintf(OUT_FP, "%s\t%d\t%s\t%s\t%c", hit->chromosome->desc(),
					readstart + 1, // 1-initialized
					ALIGNSEQ,//Alignment String
					_read.id(), ((hit->orientation == '+') ? 'D' : 'P'));

			if (_config.SCORES_OUT)
				fprintf(OUT_FP, "\t%.1f", (double) (hit->gaps * _config.GAP_SCORE
						+ (hit->mismatches - hit->gaps) * _config.MM_SCORE
						- (((int)_read.length()) - hit->mismatches) * _config.M_SCORE));
			else
				fprintf(OUT_FP, "\t%d", ((int)_read.length()) - hit->mismatches);
			// lt - changed this to report num matches instead of num mismatches

			fprintf(OUT_FP, "\t%d\t%d\t%d", num, // Number of hits
					((int)_read.length()), // length of hit on genome
					0);

			if (_read.format() == 2)
				fprintf(OUT_FP, "\t%d", _read.pe_flag());
			if (strlen(_read.quality()[0]) != 0)
				fprintf(OUT_FP, "\t%s", _read.quality()[0]);
			if (strlen(_read.quality()[1]) != 0)
				fprintf(OUT_FP, "\t%s", _read.quality()[1]);
			if (strlen(_read.quality()[2]) != 0)
				fprintf(OUT_FP, "\t%s", _read.quality()[2]);

			if (_config.FLANKING != 0) {
				fstart = (readstart < _config.FLANKING) ? 0 : readstart - _config.FLANKING;
				flen
						= (readstart + ((int)_read.length()) + _config.FLANKING
								> hit->chromosome->length()) ? hit->chromosome->length()
								- fstart
					: readstart - fstart + ((int)_read.length()) + _config.FLANKING;
				{
					for (size_t i=0; i<flen; i++)
						FLANK_SEQ[i]= (*hit->chromosome)[fstart + i];
					//strncpy(FLANK_SEQ, CHR_SEQ[hit->chromosome] + fstart, flen);
				}
				FLANK_SEQ[flen] = '\0';
				fprintf(OUT_FP, "\t%d\t%s", hit->chromosome->length(),
						FLANK_SEQ);
			} else if (_config.PRINT_SEQ > 0)
				fprintf(OUT_FP, "\t%d", hit->chromosome->length());
			//if (_config.PRINT_SEQ == 2)
			//	fprintf(OUT_FP, "\t%s", CHR_SEQ[hit->chromosome]);

			fprintf(OUT_FP, "\n");
		} else {

			/////// BED file ///////

			if (_config.SCORES_OUT) {
				fprintf(OUT_FP, "%s\t%d\t%d\t%s\t%.1f\t%c\n",
						hit->chromosome->desc(), readstart, readstart + 1
								+ ((int)_read.length()) + gap_offset, _read.id(),
						(double) (hit->gaps * _config.GAP_SCORE + (hit->mismatches
								- hit->gaps) * _config.MM_SCORE - (((int)_read.length())
								- hit->mismatches) * _config.M_SCORE), hit->orientation);
			} else {
				fprintf(OUT_FP, "%s\t%d\t%d\t%s\t%d\t%c\t0\t0\t0,0,0\t0\t0\t0\n",
						hit->chromosome->desc(), readstart, readstart + 1
								+ ((int)_read.length()) + gap_offset, _read.id(),
						((int)_read.length()) - hit->mismatches, hit->orientation);
				// lt - changed this to report num matches instead of num mismatches
			}
		}
	}

	return 1;
}


int report_read_alignment(HIT* hit, int nbest) 
{
	int hitlength = hit->end - hit->start + 1;
	unsigned int readstart;
	if (hit->orientation == '+') 
	{
		readstart = hit->start - hit->readpos + hit->start_offset; // start pos of read in genome	0-initialized
	} else 
	{
		readstart = hit->start - (((int)_read.length()) - hit->readpos - hitlength + 2) + hit->start_offset; // 0-initialized
	}

	// PERFECT HITS:
	if (hit->mismatches == 0) 
	{
		report_mapped_read(*hit->chromosome, readstart, readstart+1+((int)_read.length()), ((int)_read.length()) - hit->mismatches, nbest) ;
	}
	// HITS WITH MISMATCHES:
	else 
	{
		char gap_offset = 0;
		char gap_in_read = 0;
		char gap_in_chr = 0;
		
		// sort mismatches in ascending order according to their abs positions and 'gap before mm'-strategy if equal
		qsort(hit->edit_op, hit->mismatches, sizeof(EDIT_OPS), compare_editops);
		
		int j ;
		for (j = 0; j != hit->mismatches; ++j) 
		{
			if (hit->edit_op[j].pos < 0) 
			{
				gap_in_chr = 1;
			}
			gap_in_read = 0;
			if (!hit->edit_op[j].mm) 
			{
				if (gap_in_chr) 
				{
					gap_offset--;
					gap_in_chr = 0;
				}
				else 
				{
					gap_offset++;
					gap_in_read = 1;
				}
			}
		}
		
		report_mapped_read(*hit->chromosome, readstart, readstart+1+((int)_read.length())+gap_offset, ((int)_read.length()) - hit->mismatches, nbest) ;
	}

	return 1;
}


// sorts by position, if it's identical, then gaps have to be before mismatches
int compare_editops(const void *a, const void *b) {
	EDIT_OPS *ia = (EDIT_OPS *) a;
	EDIT_OPS *ib = (EDIT_OPS *) b;
	if (ia->pos == ib->pos)
		return (int) (ib->mm - ia->mm);

	return (int) (abs(ia->pos) - abs(ib->pos));
}

int compare_int(const void *a, const void *b) 
{
	int *ia = (int *) a;
	int *ib = (int *) b;

	return (int) (*ia - *ib);
}

void print_leftovers(const char * tag, FILE *LEFTOVER_FP)
{
	num_unmapped++ ;
	
	if (_read.format() == 0)
		fprintf(LEFTOVER_FP, "@%s%s\n%s\n+\n%s\n", _read.id(), tag, _read.data(), _read.quality()[0]);
	else if (_read.format() == 1)
		fprintf(LEFTOVER_FP, ">%s%s\n%s\n", _read.id(), tag, _read.data());
	else
		fprintf(LEFTOVER_FP, "%s%s\t%s\t%d\t%s\t%s\t%s\n", _read.id(), tag, _read.data(),
				_read.pe_flag(), _read.quality()[0], _read.quality()[1], _read.quality()[2]);
}

void print_alignment_matrix(int chrstart, int readstart, int length,
		int offset_front, int offset_end, Chromosome const &chr, char ori, int K) {
	int i, j;

	printf(" k=%d |\t-\t", K);
	if (ori == '+')
		for (i = 0; i != length; i++)
			printf("%c\t", chr[chrstart + i]);
	else
		for (i = 0; i != length; i++)
			printf("%c\t", get_compl_base(chr[chrstart - i]));
	printf(
			"\n----------------------------------------------------------------------------------------------");
	printf(
			"----------------------------------------------------------------------------------------------\n");

	for (j = 0; j != length + 1; ++j) {
		if (j == 0)
			printf("  -  |\t");
		else {
			if ((readstart == 0 && (j <= offset_front || (j > length
					- offset_end && K != 1))) || (readstart != 0 && j > length
					- offset_end))
				printf("  X  |\t");
			else
				printf("  %c  |\t", _read.data()[readstart + j - (readstart == 0)
						* offset_front - 1]);
		}

		if (j > K)
			for (i = 0; i != j - K; ++i)
				printf("\t");

		for (i = 0; i != 2* K + 1; ++i) {
			if (i - j > K)
				break;
			if (i + j > length + K)
				break;
			if (i > length || j > length)
				break;
			if ((readstart != 0 && j > length - offset_end) || (readstart == 0
					&& j < offset_front))
				break;
			//printf("{i%d,j%d}",i,j);
			//printf("%.1f(%c)\t", M[i][j], T[i][j]);
			if (i == 6 && j == 1)
				printf("---");
		}

		printf(
				"\n----------------------------------------------------------------------------------------------");
		printf(
				"----------------------------------------------------------------------------------------------\n");
	}

}

void printhits() {
	//printf("list:\n");
	//printf("print hitlist with readlength %i, read %s, last[rl-2]=%c\n",((int)_read.lenght()), READ, READ[((int)_read.lenght())-2]);
	int i;
	int c = 0;
	HIT* hit;
	for (i = _config.INDEX_DEPTH; i != ((int)_read.length()) + 1; ++i) {

		if (*(HIT_LISTS_OPERATOR + i) != NULL) {
			printf("%i: ", i);
			hit = *(HIT_LISTS_OPERATOR + i);
			do {
				//if (hit->orientation == '-') {
				printf("[%i: %i-%i/%c/chr%i/rp%i/%imm", hit->end - hit->start
						+ 1, hit->start, hit->end, hit->orientation,
						hit->chromosome->nr() + 1, hit->readpos, hit->mismatches);
				if (hit->mismatches != 0) {
					printf(":");
					int j;
					for (j = 0; j != hit->mismatches; ++j)
						printf(" %i%c", hit->edit_op[j].pos,
								(hit->edit_op[j].mm) ? 'M' : 'G');
				}
				printf("/%s]", _read.id());
				++c;
				//}
				hit = hit->next;
				if (c > 4)
					hit = NULL;
			} while (hit != NULL);
			printf("\n");
		}
		c = 0;
	}
	//printf("done\n");
}

void printhit(HIT* hit) {
	printf("[%i: %i-%i/%c/chr%i/rp%i/%imm", hit->end - hit->start + 1,
			hit->start, hit->end, hit->orientation, hit->chromosome->nr() + 1,
			hit->readpos, hit->mismatches);
	if (hit->mismatches != 0) {
		printf(":");
		int i;
		for (i = 0; i != hit->mismatches; ++i) 
			printf(" %i%c", hit->edit_op[i].pos, (hit->edit_op[i].mm) ? 'M' : 'G');
	}
	printf("/ID: %s]", _read.id());
}


void print_alignment_records(std::vector<alignment_t *> & hits, int num_unspliced_alignments, int num_spliced_alignments, int RTRIM_STRATEGY_CUT)
{
	if (hits.size()==0)
		return ;
	alignment_t * best = hits[0] ;
	assert(best->exons.size() >= 2);

	FILE* MY_OUT_FP = OUT_FP ;
	if (best->spliced)
	{
		assert(_config.SPLICED_HITS) ;
		MY_OUT_FP = SP_OUT_FP ;
		num_spliced_best++ ;
	}
	else
		num_unspliced_best++ ;

	for (int i=1; i<hits.size(); i++)
		if (hits[i]->spliced)
			num_spliced_suboptimal+=1 ;
		else
			num_unspliced_suboptimal+= 1 ;

	if ((clock()-last_spliced_report)/CLOCKS_PER_SEC>=10)
	{
		last_spliced_report = clock() ;
		fprintf(stdout, "\n# %i (%i) unspliced, %i (%i) spliced alignments, %i unmapped (spliced %2.1f%%, unmapped %2.1f%%)\n", 
				num_unspliced_best, num_unspliced_suboptimal, num_spliced_best, num_spliced_suboptimal, num_unmapped,
				100.0*num_spliced_best/(num_spliced_best+num_unspliced_best), 100.0*num_unmapped/(num_spliced_best+num_unspliced_best+num_unmapped)) ;
	}
	
	// Print spliced alignment to open file MY_OUT_FP in BED format
	fprintf(MY_OUT_FP, "%s\t%d\t%d\t%s\t%d\t%c\t%i\t%i\t0,0,0\t%d\t",
			best->chromosome->desc(),
			best->exons[0],
			best->exons[best->exons.size()-1],
			best->read_id,
			best->num_matches,
			best->strand,
			num_unspliced_alignments,
			num_spliced_alignments,
			int(best->exons.size() / 2));

	fprintf(MY_OUT_FP, "%d", best->exons[1] - best->exons[0]);
	assert(best->exons[1] - best->exons[0] >  0) ;
	assert(best->exons[1] - best->exons[0] <= MAX_EXON_LEN) ;
	
	for (uint32_t i = 2; i < best->exons.size(); i+=2)
	{
		fprintf(MY_OUT_FP, ",%d", best->exons[i+1] - best->exons[i]);
		assert(best->exons[i+1] - best->exons[i] >  0) ;
		assert(best->exons[i+1] - best->exons[i] <= MAX_EXON_LEN) ;
	}

	fprintf(MY_OUT_FP, "\t0");
	for (uint32_t i = 2; i < best->exons.size(); i+=2)
	{
		fprintf(MY_OUT_FP, ",%d", best->exons[i] - best->exons[0]);
		assert(best->exons[i] - best->exons[0] > 0) ;
	}
	
	double qpalma_score = best->qpalma_score ;

	if (_config.RTRIM_STRATEGY)
	{
		fprintf(MY_OUT_FP, "\ttrimmed=%i", RTRIM_STRATEGY_CUT) ;
		fprintf(MY_OUT_FP, "\n");
		return ;
	} 
	else
	{
		if (best->orientation=='+')
			fprintf(MY_OUT_FP, "\tqpalmaScore=%1.3f;numMatches=%i;numGaps=%i;minExonLen=%i;maxIntronLen=%i;readOrientation=%c;read=%s;quality=%s", 
					qpalma_score, best->num_matches, best->num_gaps, best->min_exon_len, best->max_intron_len, best->orientation, best->read_anno, _read.quality()[0]) ;
		else
		{
			// reverse order of quality 
			char qual[500] ;
			for (int i=0; i<((int)_read.length()); i++)
				qual[i]=_read.quality()[0][((int)_read.length())-i-1] ;
			qual[((int)_read.length())]=0 ;
			
			fprintf(MY_OUT_FP, "\tqpalmaScore=%1.3f;numMatches=%i;numGaps=%i;minExonLen=%i;maxIntronLen=%i;readOrientation=%c;read=%s;quality=%s", 
					qpalma_score, best->num_matches, best->num_gaps, best->min_exon_len, best->max_intron_len, best->orientation, best->read_anno, qual) ;
		}
	}
	
	for (int j=1; j<hits.size(); j++)
	{
		alignment_t * second  = hits[j] ;
		
		assert(second->exons.size() >= 2);

		fprintf(MY_OUT_FP, "\t%s\t%d\t%d\t%d\t%c\t%d\t",
				second->chromosome->desc(),
				second->exons[0],
				second->exons[second->exons.size() - 1],
				second->num_matches,
				second->strand,
				int(second->exons.size() / 2));

		fprintf(MY_OUT_FP, "%d", second->exons[1] - second->exons[0]);
		assert(second->exons[1] - second->exons[0] >  0) ;
		assert(second->exons[1] - second->exons[0] <= MAX_EXON_LEN) ;
		
		for (uint32_t i = 2; i < second->exons.size(); i+=2)
		{
			fprintf(MY_OUT_FP, ",%d", second->exons[i+1] - second->exons[i]);
			assert(second->exons[i+1] - second->exons[i] >  0) ;
			assert(second->exons[i+1] - second->exons[i] <= MAX_EXON_LEN) ;
		}

		fprintf(MY_OUT_FP, "\t0");
		for (uint32_t i = 2; i < second->exons.size(); i+=2)
		{
			fprintf(MY_OUT_FP, ",%d", second->exons[i] - second->exons[0]);
			assert(second->exons[i] - second->exons[0] >  0 ) ;
			//assert(second->exons[i] - second->exons[0] <= MAX_EXON_LEN ) ;
		}

		fprintf(MY_OUT_FP, "\tqpalmaScore=%1.3f;numMatches=%i;numGaps=%i;minExonLen=%i;maxIntronLen=%i;readOrientation=%c;read=%s", 
				second->qpalma_score, second->num_matches, second->num_gaps, second->min_exon_len, best->max_intron_len, second->orientation, second->read_anno) ;
	}

	fprintf(MY_OUT_FP, "\n");

}
