#include <sys/time.h>
#include <assert.h>
#include <limits.h>
#include <pthread.h>
#include <string>
#include <time.h>

#include "palmapper.h"
#include "print.h"
#include "dyn_prog/qpalma_dp.h"

#include <palmapper/Genome.h>
#include <palmapper/Hits.h>
#include <palmapper/TopAlignments.h>

TopAlignments::TopAlignments(GenomeMaps* genomemaps_) 
  : top_alignments(), num_spliced_alignments(0),
		num_unspliced_alignments(0), verbosity(0)
{
	int ret = pthread_mutex_init(&top_mutex, NULL) ;// = PTHREAD_MUTEX_INITIALIZER ;
	assert(ret==0) ;

	ALIGNSEQ = (char *) malloc((_config.MAX_READ_LENGTH + 3 * Config::MAX_EDIT_OPS)
			* sizeof(char));
	if (ALIGNSEQ == NULL) {
		fprintf(stderr, "[init_alignment_structures] Could not allocate memory\n");
		exit(1);
	}

	genomemaps = genomemaps_ ;
	MAX_EXON_LEN = 200 ;
}

u_int8_t TopAlignments::report_unspliced_hit(Read const &read, HIT *hit, int num, QPalma const *qpalma)
{
	alignment_t *algn_hit = gen_alignment_from_hit(read, hit, qpalma) ;
	algn_hit->hit = hit ;
	algn_hit->num = num ;
	algn_hit->from_gm = 2;
	
	//if (_config.OUTPUT_FILTER==OUTPUT_FILTER_TOP)
	add_alignment_record(read, algn_hit, 1) ;
	return 1 ;

	//unsigned int printed = print_top_alignment_records(read, OUT_FP, SP_OUT_FP) ;
}

int TopAlignments::construct_aligned_string(Read const &read, HIT *hit, int *num_gaps_p, int *num_mismatches_p, int *num_matches_p)
{
	int alignment_length = read.length() ;
	
	int j;
	int count_char = 0;
	char gap_offset = 0;
	char gap_in_read = 0;
	char gap_in_chr = 0;

	int num_gaps = 0 ;
	int num_mismatches = 0 ;
	int num_matches = read.length() ;

	int hitlength = hit->end - hit->start + 1;
	unsigned int readstart;
	
	if (hit->orientation == '+') {
		readstart = hit->start - hit->readpos + hit->start_offset; // start pos of read in genome	0-initialized
	} else {
		readstart = hit->start - (((int)read.length()) - hit->readpos - hitlength + 2)
				+ hit->start_offset; // 0-initialized
	}


	if (hit->mismatches == 0) {

		if (hit->orientation == '+')
			strcpy(ALIGNSEQ, read.data());
		else
		{
			for (size_t i=0; i<read.length(); i++)
				ALIGNSEQ[i] = hit->chromosome->operator [](readstart+i) ;
			//strncpy(ALIGNSEQ, CHR_SEQ[hit->chromosome] + readstart, ((int)_read.lenght()));
		}
		ALIGNSEQ[((int)read.length())] = '\0';
		
		return alignment_length ;
		
	}

	//fprintf(stderr, "mismatches=%i\n", hit->mismatches) ;
	
	
	for (int i=0; i<hit->mismatches; i++)
		assert(hit->edit_op[i].pos>=-((int) read.length()) && hit->edit_op[i].pos<=((int)read.length())) ;
	
	// sort mismatches in ascending order according to their abs positions and 'gap before mm'-strategy if equal
	qsort(hit->edit_op, hit->mismatches, sizeof(EDIT_OPS), compare_editops);

	for (int i=0; i<hit->mismatches; i++)
		assert(hit->edit_op[i].pos>=-((int) read.length()) && hit->edit_op[i].pos<=((int)read.length())) ;

	ALIGNSEQ[0] = '\0';
	
	for (j = 0; j != hit->mismatches; ++j) 
	{
		//fprintf(stderr, "j=%i, edit_op[j].pos=%i\n count_char=%i", j, hit->edit_op[j].pos, count_char) ;
		assert(hit->edit_op[j].pos>= -((int)read.length()) && hit->edit_op[j].pos<=((int)read.length())) ;
		
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
				sprintf(ALIGNSEQ + count_char, "[%c%c]", (*hit->chromosome)[readstart + hit->edit_op[j].pos - 1 + gap_offset], read.data()[hit->edit_op[j].pos - 1]);
			}
			else
			{
				sprintf(ALIGNSEQ + count_char, "[%c%c]", (*hit->chromosome)[readstart + hit->edit_op[j].pos - 1 + gap_offset], get_compl_base(read.data()[((int)read.length()) - hit->edit_op[j].pos]));
			}
		} 
		else if (gap_in_chr) 
		{
			num_gaps++ ;
			num_matches-- ;
			alignment_length-- ;
			
			if (hit->orientation == '+')
			{
				sprintf(ALIGNSEQ + count_char, "[-%c]", read.data()[hit->edit_op[j].pos - 1]);
			}
			else
			{
				sprintf(ALIGNSEQ + count_char, "[-%c]", get_compl_base(read.data()[((int)read.length()) - hit->edit_op[j].pos]));
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

	if (((int)read.length()) - hit->edit_op[j - 1].pos + gap_in_read >= 0)
	{ // again some strange case ... don't know what else to do
		for (size_t i=0; (int)i<((int)read.length())  - hit->edit_op[j - 1].pos + gap_in_read; i++)
			ALIGNSEQ[count_char+i] = (*hit->chromosome)[(readstart + hit->edit_op[j - 1].pos + gap_offset - gap_in_read)+i] ;
		//strncpy(ALIGNSEQ + count_char, CHR_SEQ[hit->chromosome] + (readstart + hit->edit_op[j - 1].pos + gap_offset - gap_in_read), 
		//((int)_read.lenght())	- hit->edit_op[j - 1].pos + gap_in_read);
		count_char += ((int)read.length()) - hit->edit_op[j - 1].pos + gap_in_read;
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

alignment_t *TopAlignments::gen_alignment_from_hit(Read const &read, HIT *best_hit, QPalma const * qpalma)
{
	int hitlength = best_hit->end - best_hit->start + 1;
	unsigned int readstart;
	unsigned int readend;

	int num_gaps=0, num_mismatches=0, num_matches=((int)read.length()) ;
	int alignment_length=construct_aligned_string(read, best_hit, &num_gaps, &num_mismatches, &num_matches); // constructs the aligned string with mismatches and gaps in ALIGNSEQ

	if (best_hit->orientation == '+') 
	{
		readstart = best_hit->start - best_hit->readpos + best_hit->start_offset; // start pos of read in genome	0-initialized
		readend = readstart + alignment_length ;
	} else 
	{
		readstart = best_hit->start - (((int)read.length()) - best_hit->readpos - hitlength + 2) + best_hit->start_offset; // 0-initialized
		readend = readstart + alignment_length ;
	}


	alignment_t *best = new alignment_t ;
	best->qpalma_score = 1000 ;
	best->num_matches = num_matches ;
	best->num_mismatches = num_mismatches ;
	best->num_gaps = num_gaps ;
	strcpy(best->read_anno, ALIGNSEQ) ;
	best->exons.push_back(readstart) ;
	best->exons.push_back(readend) ;
	best->chromosome = best_hit->chromosome ;
	best->orientation = best_hit->orientation ;
	best->strand='+' ;
	strcpy(best->read_id, read.id()) ;
	best->min_exon_len = read.length() ;
	best->max_intron_len = 0 ;
	best->spliced = false ;

	if (qpalma)
		best->qpalma_score = qpalma->score_unspliced(read, ALIGNSEQ) ;
	best->rtrim_cut=0 ;
	best->polytrim_cut_start=0 ;
	best->polytrim_cut_end=0 ;

	best->hit = NULL ;
	best->num = 0 ;

	return best ;
}

// returns > 0 if a1 scores better than a2
// returns < 0 if a1 scores worse than a2
// returns 0 if a1 and a2 score equally well
int32_t TopAlignments::compare_score(alignment_t *a1, alignment_t *a2) {

	assert(a1->qpalma_score!=1000) ;
	assert(a2->qpalma_score!=1000) ;

	if (fabs(a1->qpalma_score-a2->qpalma_score)<1e-6)
	{
		// scores are numerically identical, use start position as sorting criterion
		if (a1->exons[0] > a2->exons[0])
			return 1;
		else if (a1->exons[0] < a2->exons[0])
			return -1;
		else
			return 0;
	}
	
	
	if (a1->qpalma_score > a2->qpalma_score)
		return 1;
	else if (a1->qpalma_score < a2->qpalma_score)
		return -1;
	else
		// a1->score == a2->score
		return 0;
}

bool TopAlignments::alignment_is_equal(alignment_t *a1, alignment_t *a2) 
{
	if (fabs(a1->qpalma_score-a2->qpalma_score)>1e-6)
		return false ;
	if (a1->num_matches!=a2->num_matches)
		return false ;
	if (a1->orientation!=a2->orientation)
		return false ;
	if (a1->strand!=a2->strand)
		return false ;
	if (a1->chromosome!=a2->chromosome)
		return false ;
	if (a1->exons!=a2->exons)
	{
		/*for (int i=0; i<a1->exons.size(); i++)
			fprintf(stderr, "a1.e[%i]=%i\n", i, a1->exons[i]) ;
		for (int i=0; i<a2->exons.size(); i++)
		fprintf(stderr, "a2.e[%i]=%i\n", i, a2->exons[i]) ;*/
		
		return false ;
	}

	return true ;
	
}

void TopAlignments::clean_top_alignment_record() 
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

void TopAlignments::start_top_alignment_record() 
{
	clean_top_alignment_record() ;
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

void TopAlignments::end_top_alignment_record(Read const &read, FILE *OUT_FP, FILE *SP_OUT_FP, int rtrim_cut, int polytrim_cut_start, int polytrim_cut_end) {

	if (top_alignments.empty())
		return;

	// Process collected hits and write extracted information

	pthread_mutex_lock( &top_mutex) ;

	for (unsigned int i=0; i<top_alignments.size(); i++)
	{
		top_alignments[i]->rtrim_cut=rtrim_cut ;
		top_alignments[i]->polytrim_cut_start=polytrim_cut_start ;
		top_alignments[i]->polytrim_cut_end=polytrim_cut_end ;
		check_alignment(top_alignments[i]) ;
	}

	if (_config.REPORT_SPLICED_READS)
	{
		for (unsigned int i=0; i<top_alignments.size(); i++)
			if (top_alignments[i]->exons.size()>2)
				genomemaps->report_spliced_read(*top_alignments[i]->chromosome, top_alignments[i]->exons, 
												top_alignments[i]->num_matches, i) ;
	}
	if (_config.REPORT_MAPPED_READS)
	{
		for (unsigned int i=0; i<top_alignments.size(); i++)
		{
			if (top_alignments[i]->exons.size()<=2) 
				genomemaps->report_mapped_read(*top_alignments[i]->chromosome, top_alignments[i]->exons[0], top_alignments[i]->exons[1], 
											   top_alignments[i]->num_matches, i) ;
		}
	}

	print_top_alignment_records(read, OUT_FP, SP_OUT_FP) ;
	
	pthread_mutex_unlock( &top_mutex) ;

	clean_top_alignment_record();
	//exit(1) ;
	
}

void TopAlignments::add_alignment_record(Read const &read, alignment_t *alignment, int num_alignments)
{
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

	if ((int)read.length()>MAX_EXON_LEN)
		MAX_EXON_LEN = read.length() ;

	// seems duplicated (its already done in end_top_alignment_record)
	//if (_config.REPORT_SPLICED_READS)
	//	report_spliced_read(*alignment->chromosome, alignment->exons, alignment->num_matches, -1) ;

	bool inserted = false;

	// Go through the list and see whether we can find a hit that has a worse
	// score than this one.

	// Special case: list of top hits is still empty. Kick-start it with the
	// current hit.

	if (top_alignments.empty() || (_config.OUTPUT_FILTER!=OUTPUT_FILTER_TOP))
	{
		top_alignments.push_back(alignment);
		//fprintf(stderr, "inserting beginning %s\n", alignment->read_id);
		pthread_mutex_unlock( &top_mutex) ;
		return;
	}

	std::vector<alignment_t *>::iterator it = top_alignments.begin();

	if (verbosity >= 2)

		printf(	"[add_alignment_record] About to go through list of top hits (%ld entries)\n",
				top_alignments.size());

	//fprintf(stderr, "trying to insert %s\n", alignment->read_id);
	
	for (uint8_t i = 0; i < top_alignments.size(); i++, it++)
	{
		if (alignment_is_equal(alignment,top_alignments[i])) // already present
		{
			pthread_mutex_unlock( &top_mutex) ;
			//fprintf(stderr, "stopped inserting %s\n", alignment->read_id);
			return ;
		}
		
		if (verbosity >= 2)
			printf("[add_alignment_record] looking at next alignment %d\n", i);

		if (compare_score(alignment, top_alignments[i]) > 0) 
		{
			if (verbosity >= 2)
				printf("[add_alignment_record] reference alignment scores better than current one in top_alignments\n");
			//fprintf(stderr, "inserting %s\n", alignment->read_id);
			top_alignments.insert(it, alignment);
			inserted = true;

			if (top_alignments.size() > _config.OUTPUT_FILTER_NUM_TOP)
				top_alignments.pop_back();

			break;
		}
	}

	if (!inserted && top_alignments.size() < _config.OUTPUT_FILTER_NUM_TOP)
	{
		top_alignments.push_back(alignment);
		//fprintf(stderr, "inserting at end %s\n", alignment->read_id);
		inserted = true;
	}

	if (!inserted)
		delete alignment;

	pthread_mutex_unlock( &top_mutex) ;
}

int TopAlignments::print_top_alignment_records(Read const &read, FILE *OUT_FP, FILE *SP_OUT_FP)
{
	if (_config.OUTPUT_FORMAT==OUTPUT_FORMAT_BEDX)
	{
		return print_top_alignment_records_bedx(read, OUT_FP, SP_OUT_FP) ;
	}

	if (_config.OUTPUT_FORMAT==OUTPUT_FORMAT_SHORE || _config.OUTPUT_FORMAT==OUTPUT_FORMAT_BED)
	{
		return print_top_alignment_records_shorebed(read, OUT_FP) ;
	}

	if (_config.OUTPUT_FORMAT==OUTPUT_FORMAT_SAM)
	{
		return print_top_alignment_records_sam(read, OUT_FP) ;
	}
	
	fprintf(stderr, "ERROR: unknow output format\n") ;
	exit(1) ;
	return 0 ;
}


int TopAlignments::print_top_alignment_records_bedx(Read const &read, FILE *OUT_FP, FILE *SP_OUT_FP)
{
	if (top_alignments.size()==0)
		return 0 ;
	alignment_t * best = top_alignments[0] ;
	assert(best->exons.size() >= 2);

	if (best->rtrim_cut!=0)
		assert(best->polytrim_cut_start==0 && best->polytrim_cut_end==0) ;
	if (best->polytrim_cut_start!=0)
		assert(best->polytrim_cut_end==0) ;

	FILE* MY_OUT_FP = OUT_FP ;

	pthread_mutex_lock( &_stats.alignment_num_mutex ) ;

	if (best->spliced)
	{
		assert(_config.SPLICED_HITS) ;
		MY_OUT_FP = SP_OUT_FP ;
		_stats.alignment_num_spliced_best++ ;
	}
	else
		_stats.alignment_num_unspliced_best++ ;
	

	for (int i=1; i<(int)top_alignments.size(); i++)
		if (top_alignments[i]->spliced)
			_stats.alignment_num_spliced_suboptimal+=1 ;
		else
			_stats.alignment_num_unspliced_suboptimal+= 1 ;

	_stats.print_alignment_stats() ;
	
	pthread_mutex_unlock( &_stats.alignment_num_mutex) ;

	// Print spliced alignment to open file MY_OUT_FP in BED format
	fprintf(MY_OUT_FP, "%s\t%s\t%d\t%d\t%s\t%d\t%c\t%i\t%i\t0,0,0\t%d\t",
			best->from_gm == 2? "gm":"qp",
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
	
	{
		int polytrim_cut_start=best->polytrim_cut_start ;
		int polytrim_cut_end=best->polytrim_cut_end ;

		if (_config.RTRIM_STRATEGY && !_config.POLYTRIM_STRATEGY)
		{
			assert(polytrim_cut_end==0) ;
			polytrim_cut_end = best->rtrim_cut ;
		}

		char *read_anno=new char[strlen(best->read_anno)+(polytrim_cut_start+polytrim_cut_end)*4+20] ;
		const char *read_qual = read.quality(0);
		int read_len = read.length() ;
		
		if ((_config.POLYTRIM_STRATEGY||_config.RTRIM_STRATEGY) && (polytrim_cut_start>0 || polytrim_cut_end>0))
		{
			char const * orig_read = read.get_orig()->data() ;
			int orig_len = read.get_orig()->length() ;
			read_len = orig_len ;
			read_qual = read.get_orig()->quality(0);

			if (best->orientation=='+')
			{
				for (int i=0; i<polytrim_cut_start; i++)
				{
					read_anno[4*i]='[' ;
					read_anno[4*i+1]='=' ;
					read_anno[4*i+2]=orig_read[i] ;
					read_anno[4*i+3]=']' ;
				}
				strcpy(&(read_anno[4*polytrim_cut_start]), best->read_anno) ;
				int len=strlen(read_anno) ;
				for (int i=0; i<polytrim_cut_end; i++)
				{
					read_anno[len+4*i]='[' ;
					read_anno[len+4*i+1]='=' ;
					read_anno[len+4*i+2]=orig_read[orig_len-i-1] ;
					read_anno[len+4*i+3]=']' ;
				}
				read_anno[len+4*polytrim_cut_end]=0 ;
			}
			else
			{ // in this case read_anno is reverse complemented
				for (int i=0; i<polytrim_cut_end; i++)
				{
					read_anno[4*i]='[' ;
					read_anno[4*i+1]='=' ;
					read_anno[4*i+2]=get_compl_base(orig_read[orig_len-i-1]) ;
					read_anno[4*i+3]=']' ;
				}
				strcpy(&read_anno[4*polytrim_cut_end], best->read_anno) ;
				int len=strlen(read_anno) ;
				for (int i=0; i<polytrim_cut_start; i++)
				{
					read_anno[len+4*i]='[' ;
					read_anno[len+4*i+1]='=' ;
					read_anno[len+4*i+2]=get_compl_base(orig_read[i]) ;
					read_anno[len+4*i+3]=']' ;
				}
				read_anno[len+4*polytrim_cut_start]=0 ;
			}
		}
		else
			strcpy(read_anno, best->read_anno) ;
		
		double qpalma_score = best->qpalma_score ;
		
		{
			if (best->orientation=='+')
				fprintf(MY_OUT_FP, "\tqpalmaScore=%1.3f;numMatches=%i;numGaps=%i;minExonLen=%i;maxIntronLen=%i;readOrientation=%c;read=%s;quality=%s", 
						qpalma_score, best->num_matches, best->num_gaps, best->min_exon_len, best->max_intron_len, best->orientation, read_anno, read_qual) ;
			else
			{
				// reverse order of quality 
				char qual[read_len+1] ;
				for (int i=0; i<read_len; i++)
					qual[i]=read_qual[((int)read_len)-i-1] ;
				qual[read_len]=0 ;
				
				fprintf(MY_OUT_FP, "\tqpalmaScore=%1.3f;numMatches=%i;numGaps=%i;minExonLen=%i;maxIntronLen=%i;readOrientation=%c;read=%s;quality=%s", 
						qpalma_score, best->num_matches, best->num_gaps, best->min_exon_len, best->max_intron_len, best->orientation, read_anno, qual) ;
			}
			if (_config.POLYTRIM_STRATEGY)
			{
				if (best->polytrim_cut_start)
					fprintf(MY_OUT_FP, ";polytrimStart=%i", best->polytrim_cut_start) ;
				if (best->polytrim_cut_end)
					fprintf(MY_OUT_FP, ";polytrimEnd=%i", best->polytrim_cut_end) ;
			}
			if (_config.RTRIM_STRATEGY)
			{
				if (best->rtrim_cut)
					fprintf(MY_OUT_FP, ";rtrimEnd=%i", best->rtrim_cut) ;
			} 
		}
		delete[] read_anno ;
	}


	for (unsigned int j=1; j<top_alignments.size(); j++)
	{
		alignment_t * second  = top_alignments[j] ;
		
		assert(second->exons.size() >= 2);

		fprintf(MY_OUT_FP, "\t%s\t%s\t%d\t%d\t%d\t%c\t%d\t",
				second->from_gm == 2? "gm":"qp",
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
			if (second->polytrim_cut_start)
				fprintf(MY_OUT_FP, ";polytrimStart=%i", second->polytrim_cut_start) ;
			if (second->polytrim_cut_start)
				fprintf(MY_OUT_FP, ";polytrimEnd=%i", second->polytrim_cut_end) ;
		}
	}

	fprintf(MY_OUT_FP, "\n");
	return top_alignments.size() ;
}

int TopAlignments::print_top_alignment_records_shorebed(Read const &read, FILE *OUT_FP)
{
	int printed = 0 ;
	for (unsigned int i=0; i<top_alignments.size(); i++)
	{
		printed+= print_alignment_shorebed(read, OUT_FP, top_alignments[i]->hit, top_alignments[i]->num)  ;
	}
	return printed ;
}

int TopAlignments::print_alignment_shorebed(Read const &read, FILE *OUT_FP, HIT* hit, unsigned int num)
{
	if (_config.STATISTICS)
		_stats.HITS_MM[hit->mismatches]++;

	int fstart, flen;

	int hitlength = hit->end - hit->start + 1;
	unsigned int readstart;
	if (hit->orientation == '+') {
		readstart = hit->start - hit->readpos + hit->start_offset; // start pos of read in genome	0-initialized
	} else {
		readstart = hit->start - (((int)read.length()) - hit->readpos - hitlength + 2)
				+ hit->start_offset; // 0-initialized
	}

	char FLANK_SEQ[Config::MAX_READ_LENGTH + 200];


	fprintf(OUT_FP, "%s\t%d\t%s\t%s\t%c",
		hit->chromosome->desc(),
		readstart + 1, 	// 1-initialized
		ALIGNSEQ,	//Alignment String
		read.id(),
		((hit->orientation == '+') ? 'D' : 'P'));

	if (_config.SCORES_OUT)
		fprintf(OUT_FP, "\t%.1f", (double) (hit->gaps * _config.GAP_SCORE
				+ (hit->mismatches - hit->gaps) * _config.MM_SCORE
				- (((int)read.length()) - hit->mismatches) * _config.M_SCORE));
	else
		fprintf(OUT_FP, "\t%d", ((int)read.length()) - hit->mismatches);
	// lt - changed this to report num matches instead of num mismatches

	fprintf(OUT_FP, "\t%d\t%d\t%d",
		num, 			// Number of hits
		((int)read.length()),	// length of hit on genome
		0);

	if (read.format() == 2)
		fprintf(OUT_FP, "\t%d", read.pe_flag());
	if (strlen(read.quality(0)) != 0)
		fprintf(OUT_FP, "\t%s", read.quality(0));
	if (strlen(read.quality(1)) != 0)
		fprintf(OUT_FP, "\t%s", read.quality(1));
	if (strlen(read.quality(2)) != 0)
		fprintf(OUT_FP, "\t%s", read.quality(2));

	if (_config.FLANKING != 0) {

		fstart = (readstart < _config.FLANKING) ? 0 : readstart - _config.FLANKING;
		flen = (readstart + 
				((int)read.length()) + _config.FLANKING > hit->chromosome->length()) ?
					hit->chromosome->length() - fstart
				: 	readstart - fstart + ((int)read.length())
			+ _config.FLANKING;
		{
			for (size_t i=0; i<(size_t)flen; i++)
				FLANK_SEQ[i]= (*hit->chromosome)[fstart + i];
			//strncpy(FLANK_SEQ, CHR_SEQ[hit->chromosome] + fstart, flen);
		}
		FLANK_SEQ[flen] = '\0';
		fprintf(OUT_FP, "\t%d\t%s",
			hit->chromosome->length(),
			FLANK_SEQ);

	} else if (_config.PRINT_SEQ > 0) {
	
		fprintf(OUT_FP, "\t%d", hit->chromosome->length());
	
	}
	//if (_config.PRINT_SEQ == 2)
	//	fprintf(OUT_FP, "\t%s", CHR_SEQ[hit->chromosome]);

	fprintf(OUT_FP, "\n");

	return 1;
}

/* In the following (commented) ALIGNSEQ is generated a second time. Since by
   creation of ALIGNSEQ the hit structure is changed (gap encoding), this won't work.
   The function above just uses the already generated ALIGNSEQ variable.
*/
	// PERFECT HITS:
/*	if (hit->mismatches == 0) {

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
				//strncpy(ALIGNSEQ + count_char, CHR_SEQ[hit->chromosome]
				//		+ (readstart + hit->edit_op[j - 1].pos + gap_offset
				//				- gap_in_read), hit->edit_op[j].pos
				//		- hit->edit_op[j - 1].pos - 1 + gap_in_read); // -1???

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
			//strncpy(ALIGNSEQ + count_char, CHR_SEQ[hit->chromosome] + (readstart
			//	+ hit->edit_op[j - 1].pos + gap_offset - gap_in_read),
			//		  ((int)_read.lenght()) - hit->edit_op[j - 1].pos + gap_in_read);
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
}*/

// SAM format
int TopAlignments::print_top_alignment_records_sam(Read const &read, FILE *OUT_FP)
{
	if (top_alignments.size()==0)
		return 0 ;
   // if (_config.RTRIM_STRATEGY)
   //     return 0;

    // pre compute H0, H1, and H2 tags
    uint32_t H0 = 0 ;
    uint32_t H1 = 0 ;
    uint32_t H2 = 0 ;
	for (unsigned int j=0; j<top_alignments.size(); j++)
	{
        if (top_alignments[j]->num_gaps + top_alignments[j]->num_mismatches == 0)
            H0 += 1 ;
        if (top_alignments[j]->num_gaps + top_alignments[j]->num_mismatches == 1)
            H1 += 1 ;
        if (top_alignments[j]->num_gaps + top_alignments[j]->num_mismatches == 2)
            H2 += 1 ;
    }

    FILE* MY_OUT_FP = OUT_FP ;
    Read const * curr_read;

	for (unsigned int j=0; j<top_alignments.size(); j++)
	{
		alignment_t * curr_align  = top_alignments[j] ;

		assert(curr_align->exons.size() >= 2);
		
		int min_exon_len = curr_align->exons[1]-curr_align->exons[0] ;
		int max_intron_len = 0 ;
		int min_intron_len = INT_MAX ;
		for (unsigned int k=2; k<curr_align->exons.size(); k+=2)
		  {
		    if (min_exon_len>curr_align->exons[k+1]-curr_align->exons[k])
		      min_exon_len = curr_align->exons[k+1]-curr_align->exons[k] ;
		    if (max_intron_len<curr_align->exons[k]-curr_align->exons[k-1])
		      max_intron_len = curr_align->exons[k]-curr_align->exons[k-1] ;
		    if (min_intron_len>curr_align->exons[k]-curr_align->exons[k-1])
		      min_intron_len = curr_align->exons[k]-curr_align->exons[k-1] ;
		  }

        if (j == 0) 
        {
            if (curr_align->rtrim_cut!=0)
                assert(curr_align->polytrim_cut_start==0 && curr_align->polytrim_cut_end==0) ;
            if (curr_align->polytrim_cut_start!=0)
                assert(curr_align->polytrim_cut_end==0) ;
			
			pthread_mutex_lock( &_stats.alignment_num_mutex) ;
			
            if (curr_align->spliced)
            {
                assert(_config.SPLICED_HITS) ;
                //MY_OUT_FP = SP_OUT_FP ;
                _stats.alignment_num_spliced_best++ ;
            }
            else
                _stats.alignment_num_unspliced_best++ ;
			
            for (int i=1; i<(int)top_alignments.size(); i++)
                if (top_alignments[i]->spliced)
                    _stats.alignment_num_spliced_suboptimal+=1 ;
                else
                    _stats.alignment_num_unspliced_suboptimal+= 1 ;
			
			
            _stats.print_alignment_stats() ;

			pthread_mutex_unlock( &_stats.alignment_num_mutex ) ;
        }

		int polytrim_cut_start=curr_align->polytrim_cut_start ;
		int polytrim_cut_end=curr_align->polytrim_cut_end ;

		if (_config.RTRIM_STRATEGY && !_config.POLYTRIM_STRATEGY)
		{
			assert(polytrim_cut_end==0) ;
			polytrim_cut_end = curr_align->rtrim_cut ;
		}

        if ((_config.POLYTRIM_STRATEGY || _config.RTRIM_STRATEGY) && (polytrim_cut_start>0 || polytrim_cut_end>0))
			if (read.get_orig())
				curr_read = read.get_orig() ;
			else
				curr_read = &read ;
        else
            curr_read = &read ;

		fprintf(MY_OUT_FP, "%s", curr_align->read_id) ;
		uint32_t flag=0 ;
		flag+=((curr_align->orientation=='-')*16) ;
        /* flag+=_config.SEQUENCING_WAS_PAIRED ;
         * flag+=(curr_read.MAPPED_AS_PAIR*2) ;
         * flag+=(IS_UNMAPPED*4) ;
         * flag+=(curr_read.MATE_IS_UNMAPPED*8) ;
         * flag+=(curr_read.STRAND_OF_MATE*32) ;
         * flag+=(curr_read.FIRST_IN_PAIR)?64:128 ;
         */
        flag+=((top_alignments.size()>1)*256) ;
		fprintf(MY_OUT_FP, "\t%d\t%s\t%d\t%i", 
				flag, 
				curr_align->chromosome->desc(),
				curr_align->exons[0] + 1,
                254 - j);

	    //	double qpalma_score = best->qpalma_score ;
        
        // determine CIGAR
        char __cigar[500] ; 
        char cigar[500] ;
        char cig_buf[255] ;

        uint32_t pos = 0;
        uint32_t idx = 0 ;
        
        // handle trimmed start as soft clips
		if (_config.POLYTRIM_STRATEGY && polytrim_cut_start>0 )
        {
            snprintf (cig_buf, (size_t) 255, "%d", polytrim_cut_start) ;
            for (uint32_t ii=0; ii < strlen(cig_buf); ii++)
                cigar[pos + ii] = cig_buf[ii] ;
            pos += strlen(cig_buf) ;
            cigar[pos++] = 'S' ;
        }

        for (uint32_t i = 0; i < strlen(curr_align->read_anno); i++)
        {
            if (curr_align->read_anno[i] != '[') 
                __cigar[idx] = 'M' ;	
            else
            {
                if (curr_align->read_anno[i+1] == '-')
                    __cigar[idx] = 'I' ;
                else if (curr_align->read_anno[i+2] == '-')
                    __cigar[idx] = 'D' ;
                else
                    __cigar[idx] = 'M' ;
                i += 3 ;
            }
            idx += 1 ;
        }
        __cigar[idx] = 0 ;

        uint32_t last = __cigar[0] ;
        uint32_t count = 1 ;
        uint32_t ii = 0;
        uint32_t indel_offset = 0 ;
        idx = 0 ;
        uint32_t cum_size = (curr_align->exons[idx + 1] - curr_align->exons[idx]) ;


        for (uint32_t i = 1; i < strlen(__cigar); i++)
        {
            if (i == cum_size + indel_offset && idx + 2 < curr_align->exons.size())
            {
                snprintf (cig_buf, (size_t) 255, "%d", count) ;
                for (ii=0; ii < strlen(cig_buf); ii++)
                    cigar[pos + ii] = cig_buf[ii] ;
                pos += strlen(cig_buf) ;
                cigar[pos++] = last ;
                count = 0 ;
                last = ' ' ;
                snprintf (cig_buf, (size_t) 255, "%d", curr_align->exons[idx + 2] - curr_align->exons[idx + 1]) ;
                for (ii=0; ii < strlen(cig_buf); ii++)
                    cigar[pos + ii] = cig_buf[ii] ;
                pos += strlen(cig_buf) ;
                cigar[pos++] = 'N' ;
                idx += 2 ;
                cum_size += (curr_align->exons[idx + 1] - curr_align->exons[idx]) ;
                assert(curr_align->exons[idx+1] - curr_align->exons[idx] >  0) ;
                assert(curr_align->exons[idx+1] - curr_align->exons[idx] <= MAX_EXON_LEN) ;
            }

            if (__cigar[i] != (int)last) 
            {
                if (last != ' ')
                {
                    snprintf (cig_buf, (size_t) 255, "%d", count) ;
                    for(ii=0; ii < strlen(cig_buf); ii++)
                        cigar[pos + ii] = cig_buf[ii] ;
                    pos += strlen(cig_buf) ;
                    cigar[pos++] = last ;
                }
                if (last == 'D')
                    indel_offset -= count;
                if (last == 'I')
                    indel_offset += count;
                count = 1 ;
                last = __cigar[i] ;
            }
            else
                count += 1 ;
        }
        snprintf (cig_buf, (size_t) 255, "%d", count) ;
        for (ii=0; ii < strlen(cig_buf); ii++)
            cigar[pos + ii] = cig_buf[ii] ;
        pos += ii ;
        cigar[pos++] = last ;
        if (last == 'D')
            indel_offset -= count ; 
        else if (last == 'I')
            indel_offset += count ; 

        // handle trimmed reads end
		if ((_config.POLYTRIM_STRATEGY || _config.RTRIM_STRATEGY) && polytrim_cut_end>0)
        {
            snprintf (cig_buf, (size_t) 255, "%d", polytrim_cut_end) ;
            for (ii=0; ii < strlen(cig_buf); ii++)
                cigar[pos + ii] = cig_buf[ii] ;
            pos += strlen(cig_buf) ;
            cigar[pos++] = 'S' ;
        }
        if (cum_size + indel_offset + polytrim_cut_start + polytrim_cut_end != curr_read->length()) 
            fprintf(stdout, "WARNING - block sum does not match readlength: block_sum=%i, readlength=%i, read=%s, read_id=%s \n", cum_size + polytrim_cut_start + polytrim_cut_end + indel_offset, curr_read->length(), curr_read->data(), curr_align->read_id) ;
            //fprintf(stderr, "cum_size %i, trim_start %i, trim_end %i, read_length %i, read %s , indel_offset %i, read anno %s \n", cum_size, curr_align->polytrim_cut_start, curr_align->polytrim_cut_end, curr_read->length(), curr_read->data(), indel_offset, curr_align->read_anno) ;
        //assert(cum_size + indel_offset + curr_align->polytrim_cut_start + curr_align->polytrim_cut_end == curr_read->length()) ;
        cigar[pos] = 0 ;

//		if (curr_align->orientation=='+' || curr_align->exons.size() < 3)
        fprintf(MY_OUT_FP, "\t%s\t*\t0\t0", cigar) ; 
/*		else
		{
			// reverse order of cigar
            char rcigar[500] ;
            uint32_t marker = 0 ;
            for (int k=1; k<strlen(cigar); k++)
            {
                if (cigar[strlen(cigar)-k-1] <= '9')
                    continue ;
                else
                {
                    for (uint32_t kk = 0; kk < k - marker; kk++)
                        rcigar[marker+kk]=cigar[strlen(cigar)-k+kk] ;
                    marker = k;
                }
            }
            for (uint32_t kk = 0; kk < strlen(cigar) - marker; kk++)
                rcigar[marker+kk]=cigar[kk] ;
            
			rcigar[strlen(cigar)]=0 ;
			fprintf(MY_OUT_FP, "\t%s\t*\t0\t0", rcigar) ; 
		}*/

        if (curr_align->orientation=='+')
            fprintf(MY_OUT_FP, "\t%s\t%s", curr_read->data(), curr_read->quality(0)) ;
        else
        {
            // reverse order of quality
            char qual[500] ;
            for (int k=0; k<((int)curr_read->length()); k++)
                qual[k]=(curr_read->quality(0))[((int)(curr_read->length()))-k-1] ;
            qual[((int)(curr_read->length()))]=0 ;
            
            // complementary reverse read 
            char cr_read[500] ;
            for (int i=0; i<((int)curr_read->length()); i++)
                cr_read[i] = get_compl_base(curr_read->data()[((int)(curr_read->length()))-i-1]) ;

            cr_read[((int)(curr_read->length()))]=0 ;

            fprintf(MY_OUT_FP, "\t%s\t%s", cr_read, qual) ;
        }

        fprintf(MY_OUT_FP, "\tNM:i:%i", curr_align->num_mismatches + curr_align->num_gaps) ;
        if (H0 > 0)
            fprintf(MY_OUT_FP, "\tH0:i:%i", H0) ;
        if (H1 > 0)
            fprintf(MY_OUT_FP, "\tH1:i:%i", H1) ;
        if (H2 > 0)
            fprintf(MY_OUT_FP, "\tH2:i:%i", H2) ;

        if (curr_align->spliced)
        {
            fprintf(MY_OUT_FP, "\tXS:A:%c", curr_align->strand) ;
            fprintf(MY_OUT_FP, "\tXe:i:%i", min_exon_len) ;
            fprintf(MY_OUT_FP, "\tXI:i:%i", max_intron_len) ;
            fprintf(MY_OUT_FP, "\tXi:i:%i", min_intron_len) ;
        }
        else if (_config.STRAND > -1) {
            if ((( curr_align->orientation == '+') && _config.STRAND) || ((curr_align->orientation == '-') && ! _config.STRAND))
                fprintf(MY_OUT_FP, "\tXS:A:+") ;
            else
                fprintf(MY_OUT_FP, "\tXS:A:-") ;
        }

        fprintf(MY_OUT_FP, "\tXN:i:%i", (int)curr_align->exons.size()/2) ;
        fprintf(MY_OUT_FP, "\tZS:f:%2.3f", curr_align->qpalma_score) ;
        fprintf(MY_OUT_FP, "\tAS:i:%i", (int)(100*curr_align->qpalma_score)) ;
        fprintf(MY_OUT_FP, "\tHI:i:%i", j) ;
        fprintf(MY_OUT_FP, "\n") ;

	}

	return top_alignments.size() ;
}


