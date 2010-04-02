#include <assert.h>
#include <string>
#include <sys/time.h>
#include <time.h>
#include <pthread.h>

#include "genomemapper.h"
#include "dyn_prog/qpalma_dp.h"

#include "TopAlignments.h"

const int TopAlignments::MAX_EXON_LEN = 100 ;

TopAlignments::TopAlignments() : top_alignments(), num_spliced_alignments(0),
		num_unspliced_alignments(0), verbosity(0)
{
	int ret = pthread_mutex_init(&top_mutex, NULL) ;// = PTHREAD_MUTEX_INITIALIZER ;
	assert(ret==0) ;

	num_unspliced_best=0 ;
	num_unspliced_suboptimal=0 ;
	num_spliced_best=0 ;
	num_spliced_suboptimal=0 ;
}

u_int8_t TopAlignments::report_unspliced_hit(HIT *hit) 
{
	alignment_t *algn_hit = gen_alignment_from_hit(hit) ;
	add_alignment_record(algn_hit, 1) ;
	return 1;
}

int TopAlignments::construct_aligned_string(HIT *hit, int *num_gaps_p, int *num_mismatches_p, int *num_matches_p)
{
	int alignment_length = _read.length() ;
	
	int j;
	int count_char = 0;
	char gap_offset = 0;
	char gap_in_read = 0;
	char gap_in_chr = 0;

	int num_gaps = 0 ;
	int num_mismatches = 0 ;
	int num_matches = _read.length() ;

	int hitlength = hit->end - hit->start + 1;
	unsigned int readstart;
	
	if (hit->orientation == '+') {
		readstart = hit->start - hit->readpos + hit->start_offset; // start pos of read in genome	0-initialized
	} else {
		readstart = hit->start - (((int)_read.length()) - hit->readpos - hitlength + 2)
				+ hit->start_offset; // 0-initialized
	}


	if (hit->mismatches == 0) {

		if (hit->orientation == '+')
			strcpy(ALIGNSEQ, _read.data());
		else
		{
			for (size_t i=0; i<_read.length(); i++)
				ALIGNSEQ[i] = hit->chromosome->operator [](readstart+i) ;
			//strncpy(ALIGNSEQ, CHR_SEQ[hit->chromosome] + readstart, ((int)_read.lenght()));
		}
		ALIGNSEQ[((int)_read.length())] = '\0';
		
		return alignment_length ;
		
	}

	//fprintf(stderr, "mismatches=%i\n", hit->mismatches) ;
	
	
	for (int i=0; i<hit->mismatches; i++)
		assert(hit->edit_op[i].pos>=-((int) _read.length()) && hit->edit_op[i].pos<=((int)_read.length())) ;
	
	// sort mismatches in ascending order according to their abs positions and 'gap before mm'-strategy if equal
	qsort(hit->edit_op, hit->mismatches, sizeof(EDIT_OPS), compare_editops);

	for (int i=0; i<hit->mismatches; i++)
		assert(hit->edit_op[i].pos>=-((int) _read.length()) && hit->edit_op[i].pos<=((int)_read.length())) ;

	ALIGNSEQ[0] = '\0';
	
	for (j = 0; j != hit->mismatches; ++j) 
	{
		//fprintf(stderr, "j=%i, edit_op[j].pos=%i\n count_char=%i", j, hit->edit_op[j].pos, count_char) ;
		assert(hit->edit_op[j].pos>= -((int)_read.length()) && hit->edit_op[j].pos<=((int)_read.length())) ;
		
		if (hit->edit_op[j].pos < 0) 
		{
			hit->edit_op[j].pos = -hit->edit_op[j].pos;
			gap_in_chr = 1;
		}

		if (j == 0) 
		{
			if (hit->edit_op[0].pos - 1>=0)
			{ // don't know what to do here ... hit->edit_op[0].pos - 1 can be negative
				for (size_t i=0; (int)i<hit->edit_op[0].pos - 1; i++)
					ALIGNSEQ[i] = (*hit->chromosome)[readstart+i] ;
				//strncpy(ALIGNSEQ, CHR_SEQ[hit->chromosome] + (readstart), hit->edit_op[0].pos - 1);
				count_char += hit->edit_op[0].pos - 1;
			}
		} 
		else if (hit->edit_op[j].pos - hit->edit_op[j - 1].pos != 0) 
		{
			for (size_t i=0; (int)i<hit->edit_op[j].pos - hit->edit_op[j - 1].pos - 1 + gap_in_read; i++)
				ALIGNSEQ[count_char+i] = (*hit->chromosome)[(readstart + hit->edit_op[j - 1].pos + gap_offset - gap_in_read)+i] ;
			//strncpy(ALIGNSEQ + count_char, CHR_SEQ[hit->chromosome] + (readstart + hit->edit_op[j - 1].pos + gap_offset - gap_in_read), 
			//hit->edit_op[j].pos - hit->edit_op[j - 1].pos - 1 + gap_in_read); // -1???
			count_char += hit->edit_op[j].pos - hit->edit_op[j - 1].pos - 1 + gap_in_read;
		} // else: edit_op[j-1] must have been a gap!

		gap_in_read = 0;

		if (hit->edit_op[j].mm) 
		{
			num_mismatches++ ;
			num_matches-- ;
			
			if (hit->orientation == '+')
			{
				sprintf(ALIGNSEQ + count_char, "[%c%c]", (*hit->chromosome)[readstart + hit->edit_op[j].pos - 1 + gap_offset], _read.data()[hit->edit_op[j].pos - 1]);
			}
			else
			{
				sprintf(ALIGNSEQ + count_char, "[%c%c]", (*hit->chromosome)[readstart + hit->edit_op[j].pos - 1 + gap_offset], get_compl_base(_read.data()[((int)_read.length()) - hit->edit_op[j].pos]));
			}
		} 
		else if (gap_in_chr) 
		{
			num_gaps++ ;
			num_matches-- ;
			alignment_length-- ;
			
			if (hit->orientation == '+')
			{
				sprintf(ALIGNSEQ + count_char, "[-%c]", _read.data()[hit->edit_op[j].pos - 1]);
			}
			else
			{
				sprintf(ALIGNSEQ + count_char, "[-%c]", get_compl_base(_read.data()[((int)_read.length()) - hit->edit_op[j].pos]));
			}
			
			gap_offset--;
			gap_in_chr = 0;
		} else 
		{
			num_gaps++ ;
			alignment_length++ ;
			
			sprintf(ALIGNSEQ + count_char, "[%c-]", (*hit->chromosome)[readstart + hit->edit_op[j].pos - 1 + gap_offset]);

			gap_offset++;
			gap_in_read = 1;
		}

		count_char += 4;
	}

	// from last mismatch to end of read:

	if (((int)_read.length()) - hit->edit_op[j - 1].pos + gap_in_read >= 0)
	{ // again some strange case ... don't know what else to do
		for (size_t i=0; (int)i<((int)_read.length())  - hit->edit_op[j - 1].pos + gap_in_read; i++)
			ALIGNSEQ[count_char+i] = (*hit->chromosome)[(readstart + hit->edit_op[j - 1].pos + gap_offset - gap_in_read)+i] ;
		//strncpy(ALIGNSEQ + count_char, CHR_SEQ[hit->chromosome] + (readstart + hit->edit_op[j - 1].pos + gap_offset - gap_in_read), 
		//((int)_read.lenght())	- hit->edit_op[j - 1].pos + gap_in_read);
		count_char += ((int)_read.length()) - hit->edit_op[j - 1].pos + gap_in_read;
	}
	//fprintf(stderr, "count_char=%i num_mismatches=%i num_matches=%i\n", count_char, num_mismatches, num_matches) ;
	
	ALIGNSEQ[count_char] = '\0';
	
	if (num_gaps_p)
		*num_gaps_p = num_gaps ;
	if (num_mismatches_p)
		*num_mismatches_p = num_mismatches ;
	if (num_matches_p)
		*num_matches_p = num_matches ;

	return alignment_length ;

}

alignment_t *TopAlignments::gen_alignment_from_hit(HIT *best_hit)
{
	int hitlength = best_hit->end - best_hit->start + 1;
	unsigned int readstart;
	unsigned int readend;

	int num_gaps=0, num_mismatches=0, num_matches=((int)_read.length()) ;
	int alignment_length=construct_aligned_string(best_hit, &num_gaps, &num_mismatches, &num_matches); // constructs the aligned string with mismatches and gaps in ALIGNSEQ

	if (best_hit->orientation == '+') 
	{
		readstart = best_hit->start - best_hit->readpos + best_hit->start_offset; // start pos of read in genome	0-initialized
		readend = readstart + alignment_length ;
	} else 
	{
		readstart = best_hit->start - (((int)_read.length()) - best_hit->readpos - hitlength + 2) + best_hit->start_offset; // 0-initialized
		readend = readstart + alignment_length ;
	}


	alignment_t *best = new alignment_t ;
	best->qpalma_score = 1000 ;
	best->num_matches = num_matches ;
	best->num_gaps = num_gaps ;
	strcpy(best->read_anno, ALIGNSEQ) ;
	best->exons.push_back(readstart) ;
	best->exons.push_back(readend) ;
	best->chromosome = best_hit->chromosome ;
	best->orientation = best_hit->orientation ;
	best->strand='+' ;
	strcpy(best->read_id, _read.id()) ;
	best->min_exon_len = _read.length() ;
	best->max_intron_len = 0 ;
	best->spliced = false ;

	best->qpalma_score = _qpalma.score_unspliced(ALIGNSEQ) ;

	return best ;
}

// returns > 0 if a1 scores better than a2
// returns < 0 if a1 scores worse than a2
// returns 0 if a1 and a2 score equally well
int32_t TopAlignments::compare_score(alignment_t *a1, alignment_t *a2) {

	assert(a1->qpalma_score!=1000) ;
	assert(a2->qpalma_score!=1000) ;
	
	if (a1->qpalma_score > a2->qpalma_score)
		return 1;
	else if (a1->qpalma_score < a2->qpalma_score)
		return -1;
	else
		// a1->score == a2->score
		return 0;
}

void TopAlignments::start_top_alignment_record() 
{
	pthread_mutex_lock( &top_mutex) ;
	
	// Cleaning up and starting over.
	for (uint32_t i = 0; i < top_alignments.size(); i++) 
	{
		delete top_alignments[i];
	}

	top_alignments.clear();
	num_spliced_alignments = 0;
	num_unspliced_alignments = 0;

	pthread_mutex_unlock( &top_mutex) ;
}


void TopAlignments::check_alignment(struct alignment_t * alignment)
{
	if (alignment->exons.size()>4)
		return ;
	
	for (int i=0; i<((int)alignment->exons.size())-1; i++)
	{
		assert(alignment->exons[i]<alignment->exons[i+1]) ;
	}
	for (int i=0; i<((int)alignment->exons.size()); i+=2)
	{
		assert(alignment->exons[i+1]-alignment->exons[i]<=MAX_EXON_LEN) ;
	}
}

void TopAlignments::end_top_alignment_record(int rtrim_cut, int polytrim_cut_start, int polytrim_cut_end) {

	if (top_alignments.empty())
		return;

	// Process collected hits and write extracted information

	pthread_mutex_lock( &top_mutex) ;

	for (unsigned int i=0; i<top_alignments.size(); i++)
	{
		check_alignment(top_alignments[i]) ;
	}

	if (_config.REPORT_SPLICED_READS)
	{
		for (unsigned int i=0; i<top_alignments.size(); i++)
			if (top_alignments[i]->exons.size()>2)
				_genomemaps.report_spliced_read(*top_alignments[i]->chromosome, top_alignments[i]->exons, 
												top_alignments[i]->num_matches, i) ;
	}
	if (_config.REPORT_MAPPED_READS)
	{
		for (unsigned int i=0; i<top_alignments.size(); i++)
		{
			if (top_alignments[i]->exons.size()<=2) 
				_genomemaps.report_mapped_read(*top_alignments[i]->chromosome, top_alignments[i]->exons[0], top_alignments[i]->exons[1], 
											   top_alignments[i]->num_matches, i) ;
		}
	}

	print_top_alignment_records(rtrim_cut, polytrim_cut_start, polytrim_cut_end);
	
	pthread_mutex_unlock( &top_mutex) ;

	start_top_alignment_record();
}

void TopAlignments::add_alignment_record(alignment_t *alignment, int num_alignments) {

	if (alignment == NULL || !alignment->spliced)
		num_unspliced_alignments += num_alignments;
	else
		num_spliced_alignments += num_alignments;

	if (alignment == NULL)
		return;
	assert(num_alignments>0);

	check_alignment(alignment) ;

	if (verbosity >= 2)
		printf("entering alignment with score %f\n", alignment->qpalma_score);

	pthread_mutex_lock( &top_mutex) ;

	// seems duplicated (its already done in end_top_alignment_record)
	//if (_config.REPORT_SPLICED_READS)
	//	report_spliced_read(*alignment->chromosome, alignment->exons, alignment->num_matches, -1) ;

	bool inserted = false;

	// Go through the list and see whether we can find a hit that has a worse
	// score than this one.

	// Special case: list of top hits is still empty. Kick-start it with the
	// current hit.

	if (top_alignments.empty()) 
	{
		top_alignments.push_back(alignment);
		pthread_mutex_unlock( &top_mutex) ;
		return;
	}

	std::vector<alignment_t *>::iterator it = top_alignments.begin();

	if (verbosity >= 2)

		printf(	"[add_alignment_record] About to go through list of top hits (%ld entries)\n",
				top_alignments.size());

	for (uint8_t i = 0; i < top_alignments.size(); i++, it++)
	{
		if (alignment==top_alignments[i]) // already present
		{
			pthread_mutex_unlock( &top_mutex) ;
			return ;
		}
		
		if (verbosity >= 2)
			printf("[add_alignment_record] looking at next alignment %d\n", i);

		if (compare_score(alignment, top_alignments[i]) > 0) 
		{
			if (verbosity >= 2)
				printf("[add_alignment_record] reference alignment scores better than current one in top_alignments\n");
			top_alignments.insert(it, alignment);
			inserted = true;

			if (top_alignments.size() > _config.NUM_TOP_ALIGNMENTS)
				top_alignments.pop_back();

			break;
		}
	}

	if (!inserted && top_alignments.size() < _config.NUM_TOP_ALIGNMENTS)
	{
		top_alignments.push_back(alignment);
		inserted = true;
	}

	if (!inserted)
		delete alignment;

	pthread_mutex_unlock( &top_mutex) ;
}

void TopAlignments::print_top_alignment_records(int rtrim_cut, int polytrim_cut_start, int polytrim_cut_end)
{
	if (_config.OUTPUT_FORMAT==OUTPUT_FORMAT_BED)
		print_top_alignment_records_bed(rtrim_cut, polytrim_cut_start, polytrim_cut_end) ;

	if (_config.OUTPUT_FORMAT==OUTPUT_FORMAT_SHORE)
		print_top_alignment_records_shore(rtrim_cut, polytrim_cut_start, polytrim_cut_end) ;

	if (_config.OUTPUT_FORMAT==OUTPUT_FORMAT_SAM)
		print_top_alignment_records_sam(rtrim_cut, polytrim_cut_start, polytrim_cut_end) ;

}


void TopAlignments::print_top_alignment_records_bed(int rtrim_cut, int polytrim_cut_start, int polytrim_cut_end)
{
	if (rtrim_cut!=0)
		assert(polytrim_cut_start==0 && polytrim_cut_end==0) ;
	if (polytrim_cut_start!=0)
		assert(polytrim_cut_end==0) ;

	if (top_alignments.size()==0)
		return ;
	alignment_t * best = top_alignments[0] ;
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

	for (int i=1; i<(int)top_alignments.size(); i++)
		if (top_alignments[i]->spliced)
			num_spliced_suboptimal+=1 ;
		else
			num_unspliced_suboptimal+= 1 ;

	print_alignment_stats(num_unspliced_best, num_unspliced_suboptimal, 
						  num_spliced_best, num_spliced_suboptimal) ;
	
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
		fprintf(MY_OUT_FP, "\ttrimmed=%i", rtrim_cut) ;
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
		if (_config.POLYTRIM_STRATEGY)
		{
			if (polytrim_cut_start)
				fprintf(MY_OUT_FP, ";polytrimStart=%i", polytrim_cut_start) ;
			if (polytrim_cut_start)
				fprintf(MY_OUT_FP, ";polytrimEnd=%i", polytrim_cut_end) ;
		}
	}
	
	for (unsigned int j=1; j<top_alignments.size(); j++)
	{
		alignment_t * second  = top_alignments[j] ;
		
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

		if (_config.POLYTRIM_STRATEGY)
		{
			if (polytrim_cut_start)
				fprintf(MY_OUT_FP, ";polytrimStart=%i", polytrim_cut_start) ;
			if (polytrim_cut_start)
				fprintf(MY_OUT_FP, ";polytrimEnd=%i", polytrim_cut_end) ;
		}
	}

	fprintf(MY_OUT_FP, "\n");
}

void TopAlignments::print_top_alignment_records_shore(int rtrim_cut, int polytrim_cut_start, int polytrim_cut_end)
{
}

void TopAlignments::print_top_alignment_records_sam(int rtrim_cut, int polytrim_cut_start, int polytrim_cut_end)
{
}


int report_read_alignment(HIT* hit, int nbest)  ;

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
					for (size_t i=0; i<(size_t)flen; i++)
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
				for (size_t i=0; (int)i< hit->edit_op[0].pos - 1; i++)
					ALIGNSEQ[i]=(*hit->chromosome)[readstart + i];
				//strncpy(ALIGNSEQ, CHR_SEQ[hit->chromosome] + (readstart),
				//		hit->edit_op[0].pos - 1);
				count_char += hit->edit_op[0].pos - 1;
			} else if (hit->edit_op[j].pos - hit->edit_op[j - 1].pos != 0) 
			{
				for (size_t i=0; (int)i<hit->edit_op[j].pos - hit->edit_op[j - 1].pos - 1 + gap_in_read; i++)
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
			for (size_t i=0; (int)i<((int)_read.length()) - hit->edit_op[j - 1].pos + gap_in_read; i++)
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
					for (size_t i=0; i<(size_t)flen; i++)
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

	for (i = 0; i != (int)NUM_SCORE_INTERVALS; ++i) {

		if (printed && !_config.ALL_HIT_STRATEGY && !_config.SUMMARY_HIT_STRATEGY)
			break; // best hit strategy

		if (_hits.HITS_BY_SCORE[i].hitpointer != NULL) 
		{
			// only _config.REPEATMAP numbers of alignment will be chosen randomly:
			if (!_config.ALL_HIT_STRATEGY && !_config.SUMMARY_HIT_STRATEGY && _config.REPEATMAP < 0 && _hits.HITS_BY_SCORE[i].num > -_config.REPEATMAP) 
			{
				srand((unsigned) time(NULL));
				
				int j, k, n;
				int hits[-_config.REPEATMAP];
				for (j = 0; j != -_config.REPEATMAP; ++j) {
					n = 1;
					while (n != 0) {
						n = 0;
						hits[j] = rand() % _hits.HITS_BY_SCORE[i].num;
						for (k = 0; k != j; ++k) {
							if (hits[j] == hits[k])
								++n;
						}
					}
				}

				qsort(hits, -_config.REPEATMAP, sizeof(int), compare_int);

				hit = _hits.HITS_BY_SCORE[i].hitpointer;

				nr = 0;
				for (j = 0; j != _hits.HITS_BY_SCORE[i].num; ++j) {

					if (hits[nr] == j) {
						printed += print_alignment(hit, -_config.REPEATMAP);
						nr++;
					}

					if (nr == -_config.REPEATMAP)
						break;

					hit = hit->same_eo_succ;
				}

			} else if (_config.SUMMARY_HIT_STRATEGY) {

				hit = _hits.HITS_BY_SCORE[i].hitpointer;

				// iterate over all hits for this score and collect data to generate a summary
				// This somewhat counter-intuitive code works because of the way the pre-existing
				// code was set up: we see the hits in the order of their score here - better hits first.
				while (hit != NULL) 
				{
					_topalignments.report_unspliced_hit(hit) ;
					hit = hit->same_eo_succ;
				}

			} else { // no random selection of output alignments:

				hit = _hits.HITS_BY_SCORE[i].hitpointer;

				while (hit != NULL) {

					if (!_config.ALL_HIT_STRATEGY)
						nr = _hits.HITS_BY_SCORE[i].num;
					else
						nr = _hits.HITS_IN_SCORE_LIST;

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
				hit = _hits.HITS_BY_SCORE[i].hitpointer;

				while (hit != NULL) {

					nr = _hits.HITS_IN_SCORE_LIST;

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
