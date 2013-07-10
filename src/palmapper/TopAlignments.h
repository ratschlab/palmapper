#pragma once

#include <pthread.h>
#include <sstream>
#include <palmapper/Config.h>
#include <palmapper/GenomeMaps.h>
#include <palmapper/JunctionMap.h>
#include <palmapper/VariantMap.h>
#include <map>
#include <set>

class Mapper ;
class GenomeMaps ;
class QPalma ;
struct HIT;

typedef struct alignment_t {
	double qpalma_score;
	double sort_key;
	
	uint32_t num_matches_ref;
	uint32_t num_mismatches_ref;
	uint32_t qual_mismatches_ref;
	uint32_t num_gaps_ref;

	bool considered_variants ;
	uint32_t num_mismatches_var;
	uint32_t num_gaps_var;

	std::string read_anno ;
	//char read_anno[4*Config::MAX_READ_LENGTH] ;
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
	bool passed_filters;
	bool remapped ;
	std::vector<Variant> found_variants ; //Existing variants used to obtain a better alignment than on the reference
	std::map<int,int> variant_positions;
	std::vector<Variant> align_variants ; //Variants discovered from an alignment
	std::vector<int> aligned_positions ;

	HIT* hit ;
	int num ;
	std::vector<std::string> intron_consensus ;
	std::vector<bool> non_consensus_intron ;
	bool non_consensus_alignment ;
	
} ALIGNMENT;

class TopAlignments {

public:
    TopAlignments(GenomeMaps* genomemaps_) ;
	~TopAlignments()
	{
		free(ALIGNSEQ);
		clean_top_alignment_record() ;
	}

	u_int8_t report_unspliced_hit(Read const &read, HIT *hit, int num, QPalma const * qpalma)  ;
	int construct_aligned_string(Read const &read, HIT *hit, int *num_gaps_p, int *num_mismatches_p, int *qual_mismatches_p, int *num_matches_p);
	alignment_t *gen_alignment_from_hit(Read const &read, HIT *best_hit, QPalma const *qpalma) ;
	int construct_aligned_string(HIT *hit, int *num_gaps_p, int *num_mismatches_p, int *num_matches_p);
	
	inline void free_alignment_record(alignment_t *&alignment) 
	{
		alignment->variant_positions.clear() ;
		alignment->align_variants.clear() ;
		alignment->found_variants.clear() ;
		alignment->aligned_positions.clear() ;
		alignment->non_consensus_intron.clear() ;
		alignment->intron_consensus.clear() ;
		alignment->exons.clear() ;

		delete alignment ;
		alignment=NULL ;
	}

	void init_top_alignment_index()  ;
	void update_top_alignment_index()  ;
	alignment_t * add_alignment_record(alignment_t *alignment, int num_alignments);
	

	void clean_top_alignment_record()  ;
	void start_top_alignment_record()  ;
	void check_alignment(struct alignment_t * alignment) ;
	void end_top_alignment_record(Read const &read, std::ostream *OUT_FP, std::ostream *SP_OUT_FP,std::ostream *VARIANTS_FP, int rtrim_cut, int polytrim_start, int polytrim_end, 
								  JunctionMap &junctionmap, VariantMap & variants) ;

	int print_top_alignment_records(Read const &read, std::ostream *OUT_FP, std::ostream *SP_OUT_FP) ;
	int print_top_alignment_records_bedx(Read const &read, std::ostream *OUT_FP, std::ostream *SP_OUT_FP);
	int print_top_alignment_records_shorebed(Read const &read, std::ostream *OUT_FP, std::ostream *SP_OUT_FP) ;
	int print_alignment_shorebed(Read const &read, std::ostream *OUT_FP, alignment_t* align, unsigned int num)  ;
	int print_top_alignment_records_sam(Read const &read, std::ostream *OUT_FP, std::ostream *SP_OUT_FP) ;
	int print_top_alignment_variants(int rank, int total, std::ostream *OUT_FP, std::vector<Variant> v, const char * chr);
	
	static void print_bam_header(Genome& genome, FILE *OUT_FP);
	static FILE* open_bam_pipe(std::string & out_fname) ;
	static int close_bam_pipe(FILE * FP) ;

	size_t size()
	{
		return top_alignments.size() ;
	}
	alignment_t * get_alignment(int idx)
	{
		return top_alignments[idx] ;
	}

	bool stop_aligning() ;
	void update_max_editops()  ;
	int get_max_editops()  ;
	int get_max_mismatches() ;
	int get_max_gaps() ;
	

protected:

	int32_t compare_score(alignment_t *a1, alignment_t *a2) ;
	bool alignment_is_equal(alignment_t *a1, alignment_t *a2) ;
	int alignment_is_opposite(alignment_t *a1, alignment_t *a2) ;

	int spliced_is_overlapping(alignment_t *a1, alignment_t *a2) ;
	int non_consensus_overlaps_consensus(alignment_t *a1, alignment_t *a2) ;
	void sort_top_alignment_list()  ;  
	void qsort_top_alignments(alignment_t** output, int size);
	
	void determine_transcription_direction(char strand,char orientation, int side, char &transcription, char &read_forward);
    bool overlap(alignment_t* alignment1, alignment_t* alignment2); 
	
	std::vector<alignment_t *> top_alignments;
	int num_spliced_alignments;
	int num_unspliced_alignments;
	const int verbosity ;
	char *ALIGNSEQ;

	int num_filtered;
	int current_ind;
	int temp_ind;
	int max_editops ;
	bool _stop_aligning ;
	bool _drop_alignments ;
	
	//pthread_mutex_t top_mutex;

	GenomeMaps* genomemaps ;
} ;
