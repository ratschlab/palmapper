#pragma once

#include "palmapper.h"
#include "GenomeMaps.h"
#include <pthread.h>

class Mapper ;
class GenomeMaps ;
class QPalma ;
struct HIT;

typedef struct alignment_t {
  double qpalma_score;
  uint32_t num_matches;
  uint32_t num_gaps;
  char read_anno[4*Config::MAX_READ_LENGTH] ;
  std::vector<int> exons;
  Chromosome const *chromosome;
  char orientation;
  char strand ;
  char read_id[Config::MAX_READ_ID_LENGTH];
  int min_exon_len ;
  int max_intron_len ;
  bool spliced ;
  //bool rescued ;
  //int rescue_start, rescue_end ;
  int rtrim_cut ;
  int polytrim_cut_start ;
  int polytrim_cut_end ;
  char from_gm ;

  HIT* hit ;
  int num ;
  std::vector<std::string> intron_consensus ;

} ALIGNMENT;

class TopAlignments
{
public:
    TopAlignments(GenomeMaps* genomemaps_) ;
	~TopAlignments()
	{
		free(ALIGNSEQ);
		clean_top_alignment_record() ;
	}

	u_int8_t report_unspliced_hit(Read const &read, HIT *hit, int num, QPalma const * qpalma)  ;
	int construct_aligned_string(Read const &read, HIT *hit, int *num_gaps_p, int *num_mismatches_p, int *num_matches_p);
	alignment_t *gen_alignment_from_hit(Read const &read, HIT *best_hit, QPalma const *qpalma) ;

	void clean_top_alignment_record()  ;
	void start_top_alignment_record()  ;
	void check_alignment(struct alignment_t * alignment) ;
	void end_top_alignment_record(Read const &read, FILE *OUT_FP, FILE *SP_OUT_FP, int rtrim_cut, int polytrim_start, int polytrim_end) ;
	void add_alignment_record(Read const &read, alignment_t *alignment, int num_alignments) ;

	int print_top_alignment_records(Read const &read, FILE *OUT_FP, FILE *SP_OUT_FP) ;
	int print_top_alignment_records_bedx(Read const &read, FILE *OUT_FP, FILE *SP_OUT_FP);
	int print_top_alignment_records_shorebed(Read const &read, FILE *OUT_FP) ;
	int print_alignment_shorebed(Read const &read, FILE *OUT_FP, HIT* hit, unsigned int num)  ;
	int print_top_alignment_records_sam(Read const &read, FILE *OUT_FP) ;

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
	bool alignment_is_equal(alignment_t *a1, alignment_t *a2) ;

	std::vector<alignment_t *> top_alignments;
	int num_spliced_alignments;
	int num_unspliced_alignments;
	const int verbosity ;
	int MAX_EXON_LEN ;
	char *ALIGNSEQ;

	pthread_mutex_t top_mutex;

	static pthread_mutex_t num_mutex ;
	static bool num_mutex_initialized ;
	static int num_unspliced_best ;
	static int num_unspliced_suboptimal ;
	static int num_spliced_best ;
	static int num_spliced_suboptimal ;

	GenomeMaps* genomemaps ;
} ;
