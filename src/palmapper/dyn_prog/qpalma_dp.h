// Authors: Bettina Hepp, Uta Schulze, Cheng Soon Ong, Fabio De Bona, Gunnar Raetsch, Geraldine Jean, Soeren Sonnenburg
// Copyright (C) 2005-2010 by Friedrich Miescher Laboratory, Tuebingen, Germany

#ifndef _QPALMA_DP_H_
#define _QPALMA_DP_H_

#include <stdint.h>
#include "penalty_info_dp.h"
#include "fill_matrix.h"
#include "debug_tools.h"
#include "stdint.h"
#include <palmapper/Genome.h>
#include <palmapper/VariantMap.h>

bool fast_result_align(const std::vector<SeedElem*>& seed_matrix_left, const std::vector<SeedElem*>& seed_matrix_right, int z, int est_len, int dna_len, int* result_length_ptr, 
			char* est, char* dna, double* prb, int* s_align, int* e_align, int* mparam, double* alignmentscores, int* max_score_positions, 
					   penalty_struct* qparam, mode currentMode, double score_seed,std::vector<int>& vfound);




//extern void print_align(Pre_score* matrix, int length_est,  int length_dna, Align_pair* vektor, int result_length, int print_matrix);


/* 
 * The function myalign_ calculates a new alignment given the scoring scheme
 * Its input arguments are:
 *
 * num_path(id)               ->    an integer specifying which alignment has to be done
 * dna                        ->    
 * est                        ->
 * h                          ->
 * mmatrix                    -> 
 * donor                      ->
 * acceptor                   ->
 * remove_duplicate_scores    -> a boolean
 * print_matrix               -> a boolean
 *
 * [SpliceAlign, EstAlign, weightMatch, Gesamtscores, dnaest] = myalign_local(...
 *
 * the idea of the qualityScores array is as follows
 *
 * consider a matrix of 24 plifs
 * 
 * -> row major
 *
 */

class Alignment {

   private:
      int* splice_align;
      int* est_align;
      int* mmatrix_param;
      double* alignmentscores;
      struct penalty_struct** qualityFeaturesAllPaths;

      int dna_len;
      int est_len;
      int mlen;
      int nr_paths;

      int result_len;
      int* DNA_ARRAY;
      int* EST_ARRAY;

      int numQualSuppPoints;
      int numPlifs;
      bool use_quality_scores;

      uint32_t splice_align_size ;
      uint32_t est_align_size ;
      uint32_t mmatrix_param_size ;
      uint32_t alignmentscores_size ;
      uint32_t numPathsPlifs ;

	  int verbosity ;
	  
      INT len;
      REAL *limits;
      REAL *penalties;
      INT max_len;
      INT min_len;
      REAL *cache;
      enum ETransformType transform ;
      INT id;
      char * name;
      INT use_svm;
	  
	  std::vector<int> fvariants;
	  

   public:
      Alignment(int numQPlifs ,int numq, bool use_qscores, int verbosity_);
      ~Alignment() {
		  cleanup() ;
      }

      void cleanup() ;
	  double init_seed_position (int hit_read, int hit_dna, int hit_len, int& seed_i, int &seed_j, char* read, int read_len, char* dna, int d_len, struct penalty_struct* qualityScores, double* matchmatrix,int mm_len,double* prb);

	  template<bool use_variants, bool snp_merged>
      void myalign_fast(char strand, Chromosome const &chr, int nr_paths_p, char* dna, int dna_len_p, char* est, int est_len_p, double* prb, struct penalty_struct h, double* matchmatrix, int mm_len,
			 double* donor, int d_len, double* acceptor, int a_len, struct penalty_struct* qualityScores, 
			 bool remove_duplicate_scores, int hit_read, int hit_dna, double best_match, int max_number_introns, 
						int max_gap, int max_mism, int max_edit_op, int min_match,bool remapping, bool variant_mapping, int no_gap_end,int min_exom_len,int min_intron_len, const std::vector<variant_cache_t *> &variant_cache);

      void getDNAEST();
      void getAlignmentResults(int* s_align, int* e_align,
			       int* mmatrix_p, double* alignscores, double* qScores);
      int getResultLength() { return result_len; }
	  std::vector<int> getVariants() { return fvariants; }
      void getAlignmentArrays(int* dna_align, int* est_align);
      double scoreUnsplicedAlignment(const char * align_seq, double * prb, int read_length, struct penalty_struct* qualityScores, double* matchmatrix, char strand) ;
};

#endif  // _QPALMA_DP_H_
