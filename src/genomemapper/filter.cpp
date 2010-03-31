#include <assert.h>
#include <string>
#include <sys/time.h>
#include <time.h>
#include <pthread.h>

#include "genomemapper.h"
#include "dyn_prog/qpalma_dp.h"

static const int MAX_EXON_LEN=100 ;

const int verbosity = 0 ;

std::vector<alignment_t *> top_alignments;
int num_spliced_alignments = 0;
int num_unspliced_alignments = 0;

alignment_t *gen_alignment_from_hit(HIT *best_hit) ;

pthread_mutex_t top_mutex = PTHREAD_MUTEX_INITIALIZER;

u_int8_t report_unspliced_hit(HIT *hit) 
{
	alignment_t *algn_hit = gen_alignment_from_hit(hit) ;
	add_alignment_record(algn_hit, 1) ;
	return 1;
}

int construct_aligned_string(HIT *hit, int *num_gaps_p, int *num_mismatches_p, int *num_matches_p)
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
				for (size_t i=0; i<hit->edit_op[0].pos - 1; i++)
					ALIGNSEQ[i] = (*hit->chromosome)[readstart+i] ;
				//strncpy(ALIGNSEQ, CHR_SEQ[hit->chromosome] + (readstart), hit->edit_op[0].pos - 1);
				count_char += hit->edit_op[0].pos - 1;
			}
		} 
		else if (hit->edit_op[j].pos - hit->edit_op[j - 1].pos != 0) 
		{
			for (size_t i=0; i<hit->edit_op[j].pos - hit->edit_op[j - 1].pos - 1 + gap_in_read; i++)
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
		for (size_t i=0; i<((int)_read.length())  - hit->edit_op[j - 1].pos + gap_in_read; i++)
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

alignment_t *gen_alignment_from_hit(HIT *best_hit)
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

	best->qpalma_score =score_unspliced(ALIGNSEQ) ;

	return best ;
}

// returns > 0 if a1 scores better than a2
// returns < 0 if a1 scores worse than a2
// returns 0 if a1 and a2 score equally well
int32_t compare_score(alignment_t *a1, alignment_t *a2) {

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

void start_best_alignment_record() 
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


void check_alignment(struct alignment_t * alignment)
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

void end_best_alignment_record(int RTRIM_STRATEGY_CUT) {

	if (top_alignments.empty())
		return;

	// Process collected hits and write extracted information

	pthread_mutex_lock( &top_mutex) ;

	for (int i=0; i<top_alignments.size(); i++)
	{
		check_alignment(top_alignments[i]) ;
	}


	if (_config.REPORT_SPLICED_READS)
	{
		for (int i=0; i<top_alignments.size(); i++)
			report_spliced_read(*top_alignments[i]->chromosome, top_alignments[i]->exons, top_alignments[i]->num_matches, i) ;
	}

	print_alignment_records(top_alignments,	num_unspliced_alignments, num_spliced_alignments, RTRIM_STRATEGY_CUT);
	
	pthread_mutex_unlock( &top_mutex) ;

	start_best_alignment_record();
}

void add_alignment_record(alignment_t *alignment, int num_alignments) {

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
	if (_config.REPORT_SPLICED_READS)
		report_spliced_read(*alignment->chromosome, alignment->exons, alignment->num_matches, -1) ;

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

			if (top_alignments.size() > Config::NUM_TOP_ALIGNMENTS)
				top_alignments.pop_back();

			break;
		}
	}

	if (!inserted && top_alignments.size() < Config::NUM_TOP_ALIGNMENTS)
	{
		top_alignments.push_back(alignment);
		inserted = true;
	}

	if (!inserted)
		delete alignment;

	pthread_mutex_unlock( &top_mutex) ;
}

