#pragma once

#include <pthread.h>

struct alignment_t ;

class TopAlignments
{
public:
    TopAlignments() ;

	u_int8_t report_unspliced_hit(HIT *hit)  ;
	int construct_aligned_string(HIT *hit, int *num_gaps_p, int *num_mismatches_p, int *num_matches_p);
	alignment_t *gen_alignment_from_hit(HIT *best_hit) ;

	void start_best_alignment_record()  ;
	void check_alignment(struct alignment_t * alignment) ;
	void end_best_alignment_record(int rtrim_cut, int polytrim_start, int polytrim_end) ;
	void add_alignment_record(alignment_t *alignment, int num_alignments) ;

	size_t size()
	{
		return top_alignments.size() ;
	}
	alignment_t * get_alignment(int idx)
	{
		return top_alignments[idx] ;
	}
	
protected:

	int32_t compare_score(alignment_t *a1, alignment_t *a2) ;

	std::vector<alignment_t *> top_alignments;
	int num_spliced_alignments;
	int num_unspliced_alignments;
	const int verbosity ;
	static const int MAX_EXON_LEN ;

	pthread_mutex_t top_mutex;
} ;
