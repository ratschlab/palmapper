// Authors: Bettina Hepp, Uta Schulze, Cheng Soon Ong, Fabio De Bona, Gunnar Raetsch, Geraldine Jean, Soeren Sonnenburg
// Copyright (C) 2005-2010 by Friedrich Miescher Laboratory, Tuebingen, Germany

/*********************************************************************************************/
// fast_fill_matrix 
// Not fills the entire matrix: based on a banded semi-global alignment algorithm adapted to spliced alignment
// Applies twice the algorithm from a "seed position" computed from a long hit of GenomeMapper
// First unspliced alignment and in unspliced alignment, first in the diagonal from seed position 
/*********************************************************************************************/

/*
  matchmatrix: columnwise
    - A C G T N (dna)
  - x x x x x x
  A x x x x x x	
  C x x x x x x
  G x x x x x x
  T x x x x x x
  N x x x x x x
  (est)
*/

/*
  alignment matrix: columnwise
  |DNA . . .
  -+------------
  R|00 ...
  E|0.
  A|. .
  D|.  .
  .|.
  .|
  .|
*/

#include "fill_matrix.h"

static const float variant_penalty = 0.0001 ; // prefer non-variant alignments


template
void fast_fill_matrix<true,true>(int nr_paths_par, int*max_score_positions, int read_len, int dna_len, char* read, char* dna, double* prb, penalty_struct* functions, 
								 double* matchmatrix, penalty_struct* qualityScores, double* donor, double* acceptor, bool remove_duplicate_scores,int seed_i, int seed_j, 
								 std::vector<SeedElem *>& seed_matrix_left, std::vector<SeedElem *>& seed_matrix_right, int max_number_introns, int max_number_deletions,
								 int max_gap, int max_mism, int max_edit_op, int min_match, int verbosity,mode currentMode, bool remapping, 
								 int no_gap_end,int min_exon_len,int min_intron_len, const std::vector<variant_cache_t *> &variant_cache);



template
void fast_fill_matrix<true,false>(int nr_paths_par, int*max_score_positions, int read_len, int dna_len, char* read, char* dna, double* prb, penalty_struct* functions, 
								  double* matchmatrix, penalty_struct* qualityScores, double* donor, double* acceptor, bool remove_duplicate_scores,int seed_i, int seed_j, 
								  std::vector<SeedElem *>& seed_matrix_left, std::vector<SeedElem *>& seed_matrix_right, int max_number_introns, int max_number_deletions,
								  int max_gap, int max_mism, int max_edit_op, int min_match, int verbosity,mode currentMode, bool remapping, 
								  int no_gap_end,int min_exon_len,int min_intron_len, const std::vector<variant_cache_t *> &variant_cache);


template
void fast_fill_matrix<false,false>(int nr_paths_par, int*max_score_positions, int read_len, int dna_len, char* read, char* dna, double* prb, penalty_struct* functions, 
								   double* matchmatrix, penalty_struct* qualityScores, double* donor, double* acceptor, bool remove_duplicate_scores,int seed_i, int seed_j, 
								   std::vector<SeedElem *>& seed_matrix_left, std::vector<SeedElem *>& seed_matrix_right, int max_number_introns, int max_number_deletions, 
								   int max_gap, int max_mism, int max_edit_op, int min_match, int verbosity,mode currentMode, bool remapping,
								   int no_gap_end,int min_exon_len,int min_intron_len, const std::vector<variant_cache_t *> &variant_cache);

int number_fill_matrix;

static const int MAX_SPLICE_MISMATCH_QUAL_SINGLE=40 ;
static const int MAX_SPLICE_MISMATCH_QUAL_TOTAL=80 ;
static const int MAX_SPLICE_MISMATCH_QUAL_EXTEND=3 ;
static const int MAX_SPLICE_MISMATCH_NUM_N=3 ;

#include "fill_matrix.h"
#include "debug_tools.h"

#define D_USE_QUALITY_SCORES 1

#define PERFORM_EXTRA_CHECKS false

inline bool myisfinite(double x)
{
	if (x<=-ALMOST_INFINITY || x>=ALMOST_INFINITY)
		return false ;
	return true ;
}

inline bool isnotminusinf(double x)
{
	if (x<=-ALMOST_INFINITY)
		return false ;
	return true ;
}


/* IUPAC code
 * A -> 1
 * B (C/G/T) -> 6
 * C -> 2
 * D (A/G/T) -> 7
 * G -> 3
 * H (A/C/T) -> 8
 * K (G/T) -> 9
 * M (A/C) -> 10
 * N -> 5
 * R (A/G) -> 11
 * S (G/C) -> 12
 * T -> 4
 * V (A/C/G) -> 13
 * W (A/T) -> 14
 * Y (C/T) -> 15
*/


int char_map[133]={-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 6, 2, 7, -1, -1, 3, 8, -1, -1, 9, -1, 10, 5, -1, -1, -1, 11, 12, 4, -1, 13, 14, -1, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1} ;



	

/* system of Palma. Otherwise scoring_functions should point to an array of
 * plifs scoring the respective pairs of characters together with the EST
 * quality score.
 */

//There is match if they are identical bases or one of the dna/read is a 'N'
inline bool is_a_match(int dna, int read)
{
	if (dna==read)
		return true;
	if (dna==5 || read == 5)
		return true;
	return false;
}




template<bool use_variants>
std::vector<FoundVariant> getDeletionsfromVariants(int position, std::vector<int>& endpositions, const std::vector<variant_cache_t*> & cache, bool right_side){
	std::vector<FoundVariant> dels;

	
	if (!use_variants || position <0 || position >= (int)cache.size() || cache[position]==NULL || cache[position]->end_positions.empty())
		return dels ;
	
	if (PERFORM_EXTRA_CHECKS)
		assert(cache[position]->end_positions.size() == cache[position]->id_dels.size());
	

	//Add new positions to the end of endpositions vector
	//fprintf(stdout,"start copy\n");
	
	for (unsigned int i=0; i<cache[position]->end_positions.size();i++){
		//	fprintf(stdout,"position: %i, end pos: %i id %i\n",position, cache[position]->end_positions[i],cache[position]->id_dels[i]);
		
		endpositions.push_back(cache[position]->end_positions[i]);
		FoundVariant fv;
		if (right_side){
			fv.start_pos=position;
			fv.end_pos=cache[position]->end_positions[i];
		}
		else{
			fv.start_pos=cache[position]->end_positions[i];
			fv.end_pos=position;
		}
		
		fv.id=cache[position]->id_dels[i];
		dels.push_back(fv);
	}
//	fprintf(stdout,"end copy\n");
	// if (endpositions.empty())
	// 	endpositions.assign(cache[position]->end_positions.begin(),cache[position]->end_positions.end());
	// else
	// 	endpositions.insert(endpositions.end(),cache[position]->end_positions.begin(),cache[position]->end_positions.end());

	// dels.assign(cache[position]->id_dels.begin(),cache[position]->id_dels.end());	

	return dels;	
}

template <bool use_variants, bool snp_merged>
char getBestDnaChar(char readChar, char dnaChar, int position, const std::vector<variant_cache_t*> & variant_cache )
{ 

	int readInt = check_char(readChar);	
	int temp[4]={0,0,0,0};

	if (!snp_merged)
	{
		if ( !use_variants)
			return dnaChar;

		if (PERFORM_EXTRA_CHECKS)
			assert(position >=0 && position < (int)variant_cache.size());
		
		if (variant_cache[position]==NULL || readInt == 5)
			return dnaChar;
		
		//Take variants into account and output the best score among the different possibilities (match first I guess)

	
		for(unsigned int j=0; j<variant_cache[position]->snps.size();j++){
			int variantInt= check_char(variant_cache[position]->snps[j]);
			//fprintf(stdout,"a variant for %i with %i (read=%i)\n",dnaInt,variantInt,readInt);
			
			switch ( variantInt )
			{		
			case 1:
			case 2:
			case 3:
			case 4:
				temp[variantInt-1]=1;
				break;
			case 5:
				return readChar;
			case 6:
				temp[1]=1;
				temp[2]=1;
				temp[3]=1;
				break;
			case 7:
				temp[0]=1;
				temp[2]=1;
				temp[3]=1;
				break;
			case 8:
				temp[0]=1;
				temp[1]=1;
				temp[3]=1;
				break;
			case 9:
				temp[2]=1;
				temp[3]=1;
				break;
			case 10:
				temp[0]=1;
				temp[1]=1;
				break;
			case 11:
				temp[0]=1;
				temp[2]=1;
				break;
			case 12:
				temp[1]=1;
				temp[2]=1;
				break;
			case 13:
				temp[0]=1;
				temp[1]=1;
				temp[2]=1;
				break;
			case 14:
				temp[0]=1;
				temp[3]=1;
				break;
			case 15:
				temp[1]=1;
				temp[3]=1;
				break;
			}
			
			if (temp[readInt-1]==1)
				return readChar;
		}
	
		return dnaChar;
	}
	else
	{
		int dnaInt= check_char(dnaChar);
		
		switch ( dnaInt )
		{		
		case 1:
		case 2:
		case 3:
		case 4:
			temp[dnaInt-1]=1;
			break;
		case 5:
			return readChar;
		case 6:
			temp[1]=1;
			temp[2]=1;
			temp[3]=1;
			break;
		case 7:
			temp[0]=1;
			temp[2]=1;
			temp[3]=1;
			break;
		case 8:
			temp[0]=1;
			temp[1]=1;
			temp[3]=1;
			break;
		case 9:
			temp[2]=1;
			temp[3]=1;
			break;
		case 10:
			temp[0]=1;
			temp[1]=1;
			break;
		case 11:
			temp[0]=1;
			temp[2]=1;
			break;
		case 12:
			temp[1]=1;
			temp[2]=1;
			break;
		case 13:
			temp[0]=1;
			temp[1]=1;
			temp[2]=1;
			break;
		case 14:
			temp[0]=1;
			temp[3]=1;
			break;
		case 15:
			temp[1]=1;
			temp[3]=1;
			break;
		}
		
		if (temp[readInt-1]==1)
			return readChar;
		else{
			const char map[8]="-ACGTN*";
			
			for (int i=0;i<4;i++){
				if(temp[i]==1)
					return map[i+1];
			}
		}
		
		return dnaChar;
	}
}

template <bool use_variants, bool snp_merged>
double getBestScoreWithVariants(mode currentMode, double* matchmatrix, penalty_struct* qualityScores,int mlen, char dnaChar, char readChar, double baseScore, int position, int &dnaValue, int &snp_id, const std::vector<variant_cache_t*> & variant_cache )
{
	double score;
	int dnaInt= check_char(dnaChar);
	dnaValue=dnaInt;
	snp_id=-1;
	int readInt = check_char(readChar);
	
	if (!snp_merged)
	{
		if (PERFORM_EXTRA_CHECKS && use_variants)
		{
			if (position <0 || position >= (int)variant_cache.size())
			{
				fprintf(stderr,"ERROR: pos %i cache size %i\n", position,  (int)variant_cache.size()); // BUG-TODO
				return -ALMOST_INFINITY ;
			}
		}
		
		if (currentMode == USE_QUALITY_SCORES){
			score = getScore(qualityScores,mlen,dnaInt,readInt,baseScore);
		}
		else{
			score = (matchmatrix[mlen* dnaInt +readInt]);
		}	
	
		if ( !use_variants)
			return score;
 
		assert(position >=0 && position < (int)variant_cache.size());
	
		if ( variant_cache[position]==NULL)
			return score;
		//Take variants into account and output the best score among the different possibilities (match first I guess)
	
		for(unsigned int j=0; j<variant_cache[position]->snps.size();j++){
			int dnaValuetmp=-1;
			double tmpscore= getScoreIupac(currentMode,matchmatrix,qualityScores,baseScore,mlen, variant_cache[position]->snps[j], readChar, dnaValuetmp) - variant_penalty ;
		
			if (tmpscore > score)
			{
				dnaValue=dnaValuetmp;
				snp_id=variant_cache[position]->id_snps[j];
				score=tmpscore;	
			}
		}
	
		return score;
	}
	else
	{
		//SNP merged with DNA
		score= getScoreIupac(currentMode,matchmatrix,qualityScores,baseScore,mlen, dnaChar, readChar, dnaValue);
		return score;
	}
}

template<bool use_variants, bool snp_merged>
double getBestGapWithVariants(mode currentMode, double* matchmatrix, penalty_struct* qualityScores,int mlen, char dnaChar, int position, int &dnaValue, int &snp_id, const std::vector<variant_cache_t*> & variant_cache )
{
	double score;
	int dnaInt= check_char(dnaChar);
	dnaValue=dnaInt;
	snp_id=-1;


	if (! snp_merged)
	{
		
		if (currentMode == USE_QUALITY_SCORES){
			score = matchmatrix[dnaInt];
		}
		else{
			score = matchmatrix[mlen* dnaInt];
		}	
		if (!use_variants)
			return score;

		if (PERFORM_EXTRA_CHECKS)
			assert(position >=0 && position < (int)variant_cache.size());
		if (variant_cache[position]==NULL)
			return score;
	
		//Take variants into account and output the best score among the different possibilities (match first I guess)
		for(unsigned int j=0; j<variant_cache[position]->snps.size();j++){
			int dnaValuetmp=-1;
			double tmpscore= getGapIupac(currentMode,matchmatrix,qualityScores,mlen, variant_cache[position]->snps[j],  dnaValuetmp) - variant_penalty ;
		
			if (tmpscore > score)
			{
				dnaValue=dnaValuetmp;
				snp_id=variant_cache[position]->id_snps[j];
				score=tmpscore;	
			}
		}
	
		return score;
	}
	else
	{
		//SNP merged with DNA
		score= getGapIupac(currentMode,matchmatrix,qualityScores,mlen, dnaChar,  dnaValue);
		return score;
	}
}



void initScoreCache(struct penalty_struct* qualityScores, int mlen) 
{
	const int max_qual = 100 ;
	for (int estChar=1; estChar<6; estChar++)
		for (int dnaChar=0; dnaChar<6; dnaChar++)
		{
			int currentPos = (estChar-1)*6+dnaChar;
			struct penalty_struct *currentPen = &qualityScores[currentPos];
			if (currentPen->cache)
				continue ;

			double * cache = (double*)malloc(sizeof(double)*max_qual) ;
			for (int qual=0; qual<max_qual; qual++)
				cache[qual]=getScore(qualityScores, mlen, dnaChar, estChar, qual) ;

			currentPen->cache = cache ;
		}
}


void clean_seed_matrix_vector(std::vector<SeedElem*> &matrix, int nr_paths){

	//  fprintf(stdout,"Clean seed matrix...\n");
	for (int n = 0; n < (int)matrix.size(); n++){
		if (matrix[n]!=NULL){
			for (int p=nr_paths-1;p>=0;p--){
				delete[] matrix[n]->matrices[p];
				if (matrix[n]->best_score_pos[p]!=NULL)
					matrix[n]->best_score_pos[p]->next_seed=NULL;
				delete matrix[n]->best_score_pos[p];
			}
			delete[] matrix[n]->matrices;
			delete[] matrix[n]->best_score_pos;
			delete[] matrix[n]->best_scores;
			delete matrix[n];
		} 
	}
	matrix.clear();
	// fprintf(stdout,"Clean seed matrix...END\n");
}


// Sort best_scores and best_score_pos for the input seed_matrix
// TODO: replace "bubble sort" with a more efficient routine
// (not sure this is correct this routine is correct ... )
void sort_best_scores(SeedElem* seed, int nr_paths){

	double last_score=seed->best_scores[nr_paths-1];
  
	PosScore* temp_posscore=NULL;
	int z=nr_paths-2;

	while(z>=0 && last_score>seed->best_scores[z]){
    
		temp_posscore=seed->best_score_pos[z+1];
		seed->best_score_pos[z+1]=seed->best_score_pos[z];
		seed->best_score_pos[z]=temp_posscore;
		seed->best_scores[z+1]=seed->best_scores[z];
		seed->best_scores[z]=last_score;
		z--;
	}
	temp_posscore=NULL;

}


// Sort paths for a given position in the matrix according to scores
// TODO: replace "bubble sort" with a more efficient routine
// (not sure this is correct this routine is correct ... )
void sort_best_paths(Prev_score*matrices[], int nr_paths,int matrix_position){

	double last_score=((Prev_score*)matrices[nr_paths-1]+matrix_position)->value;
  
	Prev_score temp_prevscore;
	int z=nr_paths-2;

	while(z>=0 && last_score>((Prev_score*)matrices[z]+matrix_position)->value){
		temp_prevscore.value=last_score;
		temp_prevscore.prev_i=((Prev_score*)matrices[z+1]+matrix_position)->prev_i;
		temp_prevscore.prev_j=((Prev_score*)matrices[z+1]+matrix_position)->prev_j;
		temp_prevscore.prev_matrix_no=((Prev_score*)matrices[z+1]+matrix_position)->prev_matrix_no;
		temp_prevscore.num_matches=((Prev_score*)matrices[z+1]+matrix_position)->num_matches;
		temp_prevscore.num_mismatches=((Prev_score*)matrices[z+1]+matrix_position)->num_mismatches;
		temp_prevscore.num_gaps=((Prev_score*)matrices[z+1]+matrix_position)->num_gaps;
		temp_prevscore.snp_id=((Prev_score*)matrices[z+1]+matrix_position)->snp_id;
		temp_prevscore.snp_int=((Prev_score*)matrices[z+1]+matrix_position)->snp_int;

		((Prev_score*)matrices[z+1]+matrix_position)->value= ((Prev_score*)matrices[z]+matrix_position)->value;
		((Prev_score*)matrices[z+1]+matrix_position)->prev_i= ((Prev_score*)matrices[z]+matrix_position)->prev_i;
		((Prev_score*)matrices[z+1]+matrix_position)->prev_j= ((Prev_score*)matrices[z]+matrix_position)->prev_j;
		((Prev_score*)matrices[z+1]+matrix_position)->prev_matrix_no=((Prev_score*)matrices[z]+matrix_position)->prev_matrix_no;
		((Prev_score*)matrices[z+1]+matrix_position)->num_matches=((Prev_score*)matrices[z]+matrix_position)->num_matches;
		((Prev_score*)matrices[z+1]+matrix_position)->num_mismatches=((Prev_score*)matrices[z]+matrix_position)->num_mismatches;
		((Prev_score*)matrices[z+1]+matrix_position)->num_gaps=((Prev_score*)matrices[z]+matrix_position)->num_gaps;
		((Prev_score*)matrices[z+1]+matrix_position)->snp_id=((Prev_score*)matrices[z]+matrix_position)->snp_id;
		((Prev_score*)matrices[z+1]+matrix_position)->snp_int=((Prev_score*)matrices[z]+matrix_position)->snp_int;

		((Prev_score*)matrices[z]+matrix_position)->value=temp_prevscore.value;
		((Prev_score*)matrices[z]+matrix_position)->prev_i=temp_prevscore.prev_i;
		((Prev_score*)matrices[z]+matrix_position)->prev_j=temp_prevscore.prev_j;
		((Prev_score*)matrices[z]+matrix_position)->prev_matrix_no=temp_prevscore.prev_matrix_no;
		((Prev_score*)matrices[z]+matrix_position)->num_matches=temp_prevscore.num_matches;
		((Prev_score*)matrices[z]+matrix_position)->num_mismatches=temp_prevscore.num_mismatches;
		((Prev_score*)matrices[z]+matrix_position)->num_gaps=temp_prevscore.num_gaps;
		((Prev_score*)matrices[z]+matrix_position)->snp_id=temp_prevscore.snp_id;
		((Prev_score*)matrices[z]+matrix_position)->snp_int=temp_prevscore.snp_int;

		z--;
	}
}

template <bool use_variants, bool snp_merged>
int check_min_matches(SeedElem* seed, int nr_paths, int matrix_position, int min_matches, int*matrices, int i, int j, int prev_shift, char* read, int read_len, char* dna, int dna_len, double* read_scores, int verbosity, const std::vector<variant_cache_t*> & variant_cache)
{
	double prevValue;
	int num;
	int matrix_number=0;
  
	for (int z=0;z<nr_paths;z++)
    {
		int diff_i = min_matches ;

		prevValue = ((Prev_score*)seed->matrices[z] +matrix_position)->value ; 	    
		num = ((Prev_score*)seed->matrices[z] +matrix_position)->num_matches;

		bool conserved_seq=true ;
		int num_N=0 ;
		int num_mismatch_score=0 ;
		int num_mismatch_extend=0 ;
		int pos_i=i + prev_shift;
		int pos_j=j + prev_shift;
		for (int mn=1; mn<=diff_i; mn++) // BUG-TODO -> leads to these warnings: Warning: _config.NUM_MISMATCHES==0 && alignment_mm=1 -> setting to 0 
			// couldn't figure out how to fix it
			// either one end of the exon had a mismatch too much or another end
		{
			// fprintf(stdout,"[check_min_matches] From position %i-%i with %c-%c\n",pos_i,pos_j,read[pos_i],dna[pos_j]);
			if (pos_i<0 || pos_j<0 || pos_i>=read_len || pos_j>=dna_len)
			{
				fprintf(stdout,"ERROR: [check_min_matches] out of bounds\n");
				conserved_seq=false ;
				break ;
			}
			//Take possible SNP into account
			char dna_char= getBestDnaChar<use_variants,snp_merged>(read[pos_i], dna[pos_j], pos_j, variant_cache);
			if (read[pos_i]!=dna_char)
			{
				if (check_char(read[pos_i])!=5 && check_char(dna_char)!=5 && read_scores[pos_i]>=MAX_SPLICE_MISMATCH_QUAL_SINGLE)
				{
					conserved_seq=false ;
					//fprintf(stdout,"[check_min_matches] mismatch with high quality\n");
					break ;
				}
				else
				{
					diff_i++ ;
					if (check_char(read[pos_i])==5 || check_char(dna_char)==5)
						num_N++ ;
					else
					{
						num_mismatch_score+=read_scores[pos_i] ;
						num_mismatch_extend++ ;
					}
					if (num_N>MAX_SPLICE_MISMATCH_NUM_N || num_mismatch_score>MAX_SPLICE_MISMATCH_QUAL_TOTAL || num_mismatch_extend>MAX_SPLICE_MISMATCH_QUAL_EXTEND) 
					{
						//fprintf(stdout,"[check_min_matches] too many N or too many mismatches or num_mismatch_score too high\n");
						conserved_seq=false;
						break;  
					}
				}
			}
			pos_i+=prev_shift ;
			pos_j+=prev_shift ;
		}

		if (!conserved_seq && num>=min_matches+num_N && num_N<=MAX_SPLICE_MISMATCH_NUM_N){
			//if (verbosity>0)
				fprintf(stdout,"ERROR: [check_min_matches] PROBLEM with number of matches from position %i-%i (%i for %i) with %i N\n",i,j,num,min_matches,num_N);
		}
		if (isnotminusinf(prevValue) && conserved_seq)
		{
			//fprintf(stdout,"Possible matrix: %i with value %f\n",z,prevValue);
			matrices[matrix_number]=z;
			matrix_number++;
		}
    }
	return matrix_number;
}

void print_restricted_matrix(Prev_score* matrices[],int nr_paths, int matrix_len, int row_len){

	for(int z=0; z<nr_paths;z++){
		fprintf(stdout,"Matrix %i:\n",z);
		for (int i=0; i<matrix_len;i++){
			fprintf(stdout,"%f(%i,%i,%i) ", ((Prev_score*)matrices[z]+i)->value,((Prev_score*)matrices[z]+i)->num_gaps,((Prev_score*)matrices[z]+i)->num_mismatches,((Prev_score*)matrices[z]+i)->num_matches);
			if ((i+1)%row_len==0)
				fprintf(stdout,"\n");
      
		}
	}
	fprintf(stdout,"\n\n");
}


template <bool use_variants, bool snp_merged, bool right_side,int nr_paths_par>
void fast_fill_side_unspliced_first(std::vector<SeedElem*> &seed_matrix, int read_len, int dna_len, char* read, char* dna, double* prb, penalty_struct* functions, double* matchmatrix, penalty_struct* qualityScores, double* main_site_scores, double* comp_site_scores, std::vector<int>& comp_sites,	int seed_read, int seed_dna, double* best_match_scores, bool first_seed,int max_number_introns,	int max_number_deletions, int max_gap, int max_mism, int max_edit_op, int min_match, int verbosity, mode currentMode, bool remapping, std::vector<FoundVariant> deletion_id, int no_gap_end,int min_exon_len,int min_intron_len, const std::vector<variant_cache_t*> & variant_cache)
{

//	fprintf(stdout,"START: Fill %s side of the matrix from position %i-%i (%i-%i) (%i,%i,%i,num_intron=%i)...\n",(right_side==true)?"right":"left",seed_read, seed_dna,read_len, dna_len,max_gap,max_mism,max_edit_op, max_number_introns);
  
	/***************************************************/
	/*Initialization */
	/***************************************************/
	const int mlen = 6; // length of matchmatrix
	int readChar, dnaChar; // 0:'-', 1:'a', 2:'c', 3:'g', 4:'t', 5:'n'
  
	double baseScore=0;
	double *read_scores = prb ;

	double prevValue;
	double tempValue;
	double putativeValue;
	double globalValue;
	int prevGaps;
	int prevMism;

	int matrix_position;
	int matrix_prev_position;
  
	int prev_shift;
	if(right_side)
		prev_shift=-1;
	else
		prev_shift=1;

	// Boolean that indicates if complementary splice sites had been scanned before
	bool comp_sites_filled=false;
  
	// Possible paths with enough matches before a given splice site
	int * possible_matrices= new int[nr_paths_par];
	// Possible splice sites for starting a new alignment
	std::vector<splice_pos *> possible_sites;
  
	if (currentMode==USE_QUALITY_SCORES)
		assert((int)functions->limits[functions->len-1]==2);
  
	//width of the band according to the number of gaps
	int row_len=max_gap*2+1;

	//Length of read to align from seed position and of dna from seed_position to beginning or end
	int i_len;
	int j_len;
	if (right_side){
		i_len=read_len-seed_read;
		j_len=dna_len-seed_dna;
	}
	else{
		i_len=seed_read+1;
		j_len=seed_dna+1;
	}
	int matrix_len=row_len*i_len;

	//Init score matrices
	Prev_score** matrices = new Prev_score*[nr_paths_par];
	for (int z=0; z<nr_paths_par; z++) {
		matrices[z] = new Prev_score[matrix_len]; //Allocation only for the band around the seed position
		memset(matrices[z],0,sizeof(Prev_score)*matrix_len);
	}

	// Diagonal activation: first, only activating of diagonal from seed position
	bool* disabled_diagonal=new bool[row_len];
	for (int d=0;d<max_gap;d++){
		disabled_diagonal[d]=true;
		disabled_diagonal[row_len-1-d]=true;
	}
	disabled_diagonal[max_gap]=false;

	//Init seed position score in matrices
	//Best matrix
	int dnaInt;
	int snp_id;

	if (PERFORM_EXTRA_CHECKS && use_variants)
	{
		if (seed_dna <0 || seed_dna >= (int)variant_cache.size())
		{
			fprintf(stderr,"ERROR: pos %i cache size %i %s:%i\n", seed_dna,  (int)variant_cache.size(), __FILE__, __LINE__); // BUG-TODO
			assert(0) ;
			return ;
		}
	}
	((Prev_score*)matrices[0]+ max_gap)->value = getBestScoreWithVariants<use_variants,snp_merged>(currentMode, matchmatrix, qualityScores, mlen, dna[seed_dna], read[seed_read], read_scores[seed_read],
																								   seed_dna,dnaInt,snp_id, variant_cache);
	// if (snp_id != -1){
	// 	fprintf(stdout,"Found better score with snp id %i\n",snp_id);
	// }
	
	((Prev_score*)matrices[0]+ max_gap)->prev_i = seed_read+prev_shift;
	((Prev_score*)matrices[0]+ max_gap)->prev_j = seed_dna+prev_shift;
	((Prev_score*)matrices[0]+ max_gap)->prev_matrix_no = 0;
	((Prev_score*)matrices[0]+ max_gap)->num_matches = (is_a_match(dnaInt,check_char(read[seed_read]))? 1:0);
	((Prev_score*)matrices[0]+ max_gap)->num_mismatches = 1-((Prev_score*)matrices[0]+ max_gap)->num_matches ;
	((Prev_score*)matrices[0]+ max_gap)->num_gaps =0;
  	((Prev_score*)matrices[0]+ max_gap)->snp_id =snp_id;
  	((Prev_score*)matrices[0]+ max_gap)->snp_int =dnaInt;
	//Other matrices: score at -INF
	for(int z=1; z<nr_paths_par;z++){
		((Prev_score*)matrices[z]+ max_gap)->value = -ALMOST_INFINITY; // -inf
		((Prev_score*)matrices[z]+ max_gap)->prev_i = 0; //seed_read+prev_shift;
		((Prev_score*)matrices[z]+ max_gap)->prev_j = 0; //seed_dna+prev_shift;
		((Prev_score*)matrices[z]+ max_gap)->prev_matrix_no = 0;
		((Prev_score*)matrices[z]+ max_gap)->num_matches = 0;
		((Prev_score*)matrices[z]+ max_gap)->num_mismatches = 0;
		((Prev_score*)matrices[z]+ max_gap)->num_gaps = 0;
		((Prev_score*)matrices[z]+ max_gap)->snp_id =-1;
		((Prev_score*)matrices[z]+ max_gap)->snp_int =-1;
	}

	//Alignment from the seed with (seed_read, seed_dna) position
	SeedElem* current_seed= new SeedElem();
	current_seed->read_pos=seed_read;
	current_seed->dna_pos=seed_dna;
	current_seed->max_gaps=max_gap;
	current_seed->max_mm=max_mism;
	current_seed->max_introns=max_number_introns;
	current_seed->deletion_id=deletion_id;
	current_seed->best_scores= new double[nr_paths_par];
	current_seed->best_prev_score= -ALMOST_INFINITY;
	current_seed->best_score_pos=new PosScore*[nr_paths_par];
	for(int z =0; z<nr_paths_par;z++){
		current_seed->best_scores[z]=-ALMOST_INFINITY;
		PosScore * pscore= new PosScore();    
		pscore->read_pos=0;
		pscore->dna_pos=0;
		pscore->num_gaps=0;
		pscore->num_mm=0;
		pscore->num_introns=0;
		pscore->next_seed=NULL;
		pscore->path_number=0; //path number in next matrix
		pscore->path_number_matrices=0;
		pscore->partial_score=0; //score from seed position to (read_pos,dna_pos)
		current_seed->best_score_pos[z]=pscore;
	}
	current_seed->matrices=matrices;
	seed_matrix.push_back(current_seed);
	//if (seed_matrix.size()==100000)
	//	fprintf(stdout, "seed_matrix.size()=100000\n") ;
	if (seed_matrix.size()==1000000)
	{
		fprintf(stdout, "Warning: seed_matrix.size()=1000000, aborting alignment\n") ; // BUG-TODO
		throw std::bad_alloc() ;
	}
	//if (seed_matrix.size()==10000000)
	//	fprintf(stdout, "seed_matrix.size()=10000000\n") ;
	//if (seed_matrix.size()==100000000)
	//	fprintf(stdout, "seed_matrix.size()=100000000\n") ;
	

	//Can align all read nucleotides: all read nucleotides are covered by the band
	if (j_len+max_gap >= i_len){
    
		/***************************************************/
		/*banded SW algorithm                              */
		/***************************************************/
		for(int ni=0; ni<i_len;ni++){
      
			//Position in the read
			int i= seed_read-prev_shift*ni;

			//dna interval for this row
			int left_bound =seed_dna+(i-seed_read)-max_gap;
			int right_bound=seed_dna+(i-seed_read)+max_gap;
      
			for (int nj=0; nj<row_len;nj++){
				int j;
				if (right_side)
					j=left_bound+nj;
				else
					j=right_bound-nj;

				matrix_position= ni*row_len+nj;


				//Diagonal not active or position out of bounds
				if(disabled_diagonal[nj] || (right_side && (j<seed_dna || j>=dna_len))||(!right_side &&(j>seed_dna || j<0))){
					for(int z=0;z<nr_paths_par;z++){
						((Prev_score*)matrices[z] + matrix_position)->value = -ALMOST_INFINITY ;
						((Prev_score*)matrices[z] + matrix_position)->prev_i = 0;
						((Prev_score*)matrices[z] + matrix_position)->prev_j = 0;
						((Prev_score*)matrices[z] + matrix_position)->prev_matrix_no = 0;
						((Prev_score*)matrices[z] + matrix_position)->num_matches = 0;
						((Prev_score*)matrices[z] + matrix_position)->num_mismatches = 0;
						((Prev_score*)matrices[z] + matrix_position)->num_gaps = 0;
						((Prev_score*)matrices[z] + matrix_position)->snp_id =-1;
						((Prev_score*)matrices[z] + matrix_position)->snp_int =-1;
					}
					continue;
				}


				// seed position filled before
				if(i==seed_read && j==seed_dna){
					if(ni==i_len-1){
						current_seed->best_scores[0]=((Prev_score*)current_seed->matrices[0] +matrix_position)->value; 	    
						current_seed->best_score_pos[0]->read_pos=i;
						current_seed->best_score_pos[0]->dna_pos=j;
						current_seed->best_score_pos[0]->num_gaps=0;
						current_seed->best_score_pos[0]->num_mm=((Prev_score*)current_seed->matrices[0] +matrix_position)->num_mismatches;
						current_seed->best_score_pos[0]->num_introns=0;
						current_seed->best_score_pos[0]->next_seed=NULL;
						current_seed->best_score_pos[0]->path_number=0;
						current_seed->best_score_pos[0]->path_number_matrices=0;
						current_seed->best_score_pos[0]->partial_score=((Prev_score*)current_seed->matrices[0] +matrix_position)->value;
					}
					continue;
				}

	
				/********************************************************/
				/* 1. DETECTION OF POSSIBLE SPLICE SITE                 */
				/********************************************************/
				if (main_site_scores[j]>-ALMOST_INFINITY && main_site_scores[j]<ALMOST_INFINITY){	
					if ((max_number_introns>0) && (matrix_position-row_len>=0)){

						//Compute minimum number of matches needed with relaxing the constraint close to the first seed position
						int diff_i;
						if (first_seed){
							diff_i=prev_shift*(seed_read-i);
							if (diff_i>min_match)
								diff_i=min_match;
						}
						else
							diff_i=min_match;

						//Number of paths that has at least diff_i consecutive matches that lead to the previous position in the diagonal
						int number;
						number=check_min_matches<use_variants,snp_merged>(current_seed,nr_paths_par,matrix_position-row_len, diff_i,possible_matrices, i, j, prev_shift, read, read_len, dna, dna_len, 
															   read_scores, verbosity, variant_cache);
						if (number>0){
							splice_pos* sp= new splice_pos();
							sp->site=j; //position of the splice site
							sp->i=i+prev_shift; //position of the previous position in the diagonal
							sp->matrix_pos=matrix_position-row_len; //position in the matrix of the previous position in the diagonal
							sp->number=number;
							sp->matrices=new int[number];
							for(int e=0;e<number;e++)
								sp->matrices[e]=possible_matrices[e];
							possible_sites.push_back(sp);
						}
					}
				}
				/*1. END */
	


				/*************************************************************/
				/* 2. FILLING THE DIAGONAL (MATCH/MISMATCH)                  */
				/* right side: (i-1,j-1)->(i,j); left side: (i+1,j+1)->(i,j) */
				/*************************************************************/
				dnaChar = check_char(dna[j]);
				readChar = check_char(read[i]);
				
				assert(dnaChar!=-1 && readChar!=-1);

				if (currentMode == USE_QUALITY_SCORES)
					baseScore = read_scores[i];
				else
					baseScore =0;
				
				// Best score of what it leaves to align
				if (ni<i_len-1){
					if(right_side)
						putativeValue= best_match_scores[i+1];
					else
						putativeValue= best_match_scores[0]-best_match_scores[i];
				}
				else
					putativeValue=0;

				if (PERFORM_EXTRA_CHECKS && use_variants)
				{
					if (j <0 || j >= (int)variant_cache.size())
					{
						fprintf(stderr, "ERROR: pos %i cache size %i %s:%i\n", j,  (int)variant_cache.size(), __FILE__, __LINE__); // BUG-TODO
						assert(0) ;
						return ;
					}
				}

				//Least good best global score known
				globalValue= current_seed->best_scores[nr_paths_par-1];
	
				//Information about previous position in the diagonal
				matrix_prev_position= matrix_position-row_len;
				int num_unfilled=0;
				for(int z=0;z<nr_paths_par;z++){
					Prev_score* actMatrix = (Prev_score*)matrices[z]; 
					prevValue = ((Prev_score*)actMatrix +matrix_prev_position)->value ; 
					prevGaps=((Prev_score*)actMatrix +matrix_prev_position)->num_gaps;
					prevMism=((Prev_score*)actMatrix +matrix_prev_position)->num_mismatches;
	
					if (isnotminusinf(prevValue)){

						tempValue = prevValue + getBestScoreWithVariants<use_variants,snp_merged>(currentMode, matchmatrix, qualityScores, mlen, dna[j], read[i], baseScore,  
																								  j, dnaInt,snp_id, variant_cache);
						// if (snp_id != -1){
						// 	//fprintf(stdout,"Found better score with snp id %i\n",snp_id);
						// }
						
						bool isMatching = is_a_match(dnaInt,readChar);
						if (isnotminusinf(tempValue) && (isMatching || (prevMism+1<=max_mism && prevMism+prevGaps+1<=max_edit_op)) && tempValue+putativeValue>globalValue){
							((Prev_score*)actMatrix + matrix_position)->value = tempValue;
							((Prev_score*)actMatrix + matrix_position)->prev_i = i+prev_shift; 
							((Prev_score*)actMatrix + matrix_position)->prev_j = j+prev_shift; 
							((Prev_score*)actMatrix + matrix_position)->prev_matrix_no = z;
							if (isMatching) {
								((Prev_score*)actMatrix + matrix_position)->num_matches = ((Prev_score*)actMatrix +matrix_prev_position)->num_matches+1; 
								((Prev_score*)actMatrix + matrix_position)->num_mismatches = prevMism;
							}
							else{
								((Prev_score*)actMatrix + matrix_position)->num_matches = 0 ;
								((Prev_score*)actMatrix + matrix_position)->num_mismatches = prevMism+1;
							}	      
							((Prev_score*)actMatrix + matrix_position)->num_gaps = prevGaps;
							((Prev_score*)actMatrix + matrix_position)->snp_id = snp_id;
							((Prev_score*)actMatrix + matrix_position)->snp_int = dnaInt;

							//Last row to fill: check now if unspliced alignment is better
							if(ni==i_len-1){
								if (tempValue>current_seed->best_scores[nr_paths_par-1]){
									current_seed->best_scores[nr_paths_par-1]=tempValue;
									current_seed->best_score_pos[nr_paths_par-1]->read_pos=i;
									current_seed->best_score_pos[nr_paths_par-1]->dna_pos=j;
									current_seed->best_score_pos[nr_paths_par-1]->num_gaps=((Prev_score*)actMatrix + matrix_position)->num_gaps;									
									current_seed->best_score_pos[nr_paths_par-1]->num_mm=((Prev_score*)actMatrix + matrix_position)->num_mismatches;
									current_seed->best_score_pos[nr_paths_par-1]->num_introns=0;
									current_seed->best_score_pos[nr_paths_par-1]->next_seed=NULL;
									current_seed->best_score_pos[nr_paths_par-1]->path_number=0;
									current_seed->best_score_pos[nr_paths_par-1]->path_number_matrices=z;
									current_seed->best_score_pos[nr_paths_par-1]->partial_score=tempValue;
									// Resort best_score_pos and best_scores
									sort_best_scores(current_seed, nr_paths_par);
								}
							}
						}
						//Too much mismatches/edit operations or better alignment with intron: desactivation of the diagonal
						else{
							((Prev_score*)actMatrix + matrix_position)->value = -ALMOST_INFINITY ;
							((Prev_score*)actMatrix + matrix_position)->prev_i = 0;
							((Prev_score*)actMatrix + matrix_position)->prev_j = 0;
							((Prev_score*)actMatrix + matrix_position)->prev_matrix_no = 0;
							((Prev_score*)actMatrix + matrix_position)->num_matches = 0;
							((Prev_score*)actMatrix + matrix_position)->num_mismatches = 0;
							((Prev_score*)actMatrix + matrix_position)->num_gaps = 0;
							((Prev_score*)actMatrix + matrix_position)->snp_id =-1;
							((Prev_score*)actMatrix + matrix_position)->snp_int =-1;

							num_unfilled++;
						}
					}
					else{
						((Prev_score*)actMatrix + matrix_position)->value = -ALMOST_INFINITY ;
						((Prev_score*)actMatrix + matrix_position)->prev_i = 0;
						((Prev_score*)actMatrix + matrix_position)->prev_j = 0;
						((Prev_score*)actMatrix + matrix_position)->prev_matrix_no = 0;
						((Prev_score*)actMatrix + matrix_position)->num_matches = 0;
						((Prev_score*)actMatrix + matrix_position)->num_mismatches = 0;
						((Prev_score*)actMatrix + matrix_position)->num_gaps = 0;
						((Prev_score*)actMatrix + matrix_position)->snp_id =-1;
						((Prev_score*)actMatrix + matrix_position)->snp_int =-1;
						num_unfilled++;
					}
				}
				// If no solution for all paths, desactivate this diagonal
				if(num_unfilled==nr_paths_par)
					disabled_diagonal[nj]=true;
				/*2. END */
				
				
                /******************************************************************************************/
				/* DELETION VARIANTS                                                                      */
				/* Search for deletions from the next position                                            */
				/* Call recursively like a particular intron                                              */
				/******************************************************************************************/
				
				//Partial alignment to this position still valid
				//Some read still need to be aligned
				if (!disabled_diagonal[nj] && ((right_side && i-prev_shift <read_len) || (!right_side && i-prev_shift >=0))){
					
					//Search for deletions from next DNA position
					std::vector<int> endpositions;
					std::vector<FoundVariant> idsdeletions;
					//fprintf(stdout,"position after match: %i\n",j-prev_shift);

					if (max_number_deletions>=1)
						idsdeletions = getDeletionsfromVariants<use_variants>(j-prev_shift,endpositions, variant_cache,right_side);
					
					std::vector<std::vector<FoundVariant> > tabidsdeletions;
					//Prepare to allow to jump over several consecutive deletions
					for (int d=0;d<(int)idsdeletions.size();d++){
						std::vector<FoundVariant> vtemp;
						vtemp.push_back(idsdeletions[d]);
						tabidsdeletions.push_back(vtemp);
						vtemp.clear();
					}
					idsdeletions.clear();
					while(endpositions.size()>0){
						
                        //Start DNA position for next seed is one after/before the end position of the deletion
						int jj=endpositions[endpositions.size()-1] -prev_shift;
						endpositions.pop_back();
						std::vector<FoundVariant> current_variant_ids;
						current_variant_ids=tabidsdeletions[tabidsdeletions.size()-1];
						tabidsdeletions.pop_back();
						
						if ( (right_side && jj>=dna_len) || (!right_side && jj<0))
							continue;
						
						if (PERFORM_EXTRA_CHECKS)
						{
							assert(jj>=0 && jj<dna_len);
							assert(i-prev_shift>=0 && i-prev_shift<read_len);
						}
						
						if (right_side && PERFORM_EXTRA_CHECKS){
							assert (jj>j-prev_shift);
						}
						
						if (max_number_deletions>=1)
							idsdeletions = getDeletionsfromVariants<use_variants>(jj,endpositions, variant_cache,right_side);
						
						for (int d=0;d<(int)idsdeletions.size();d++){
							std::vector<FoundVariant> vtemp=current_variant_ids;
							vtemp.push_back(idsdeletions[d]);
							tabidsdeletions.push_back(vtemp);
							vtemp.clear();
						}	
						idsdeletions.clear();						
						//Maybe here look if the seed was already filled but not sure...

						
//						fprintf(stdout, "New fill matrix from deletion (%i): %i-%i\n",(int)seed_matrix.size(),j-prev_shift,jj+prev_shift);
						
						prevMism=((Prev_score*)matrices[0] + matrix_position)->num_mismatches;
						prevGaps=((Prev_score*)matrices[0] + matrix_position)->num_gaps;
						int seed_index=seed_matrix.size();
						fast_fill_side_unspliced_first<use_variants,snp_merged,right_side,nr_paths_par>(seed_matrix, read_len, dna_len, read, dna, prb, functions, matchmatrix, qualityScores, 
																						   main_site_scores, comp_site_scores, comp_sites, i-prev_shift,
																						   jj, best_match_scores, false,
																						   max_number_introns,max_number_deletions-1,max_gap-prevGaps,max_mism-prevMism,
																						   max_edit_op-(prevGaps+prevMism),min_match, verbosity,currentMode, 
																						   remapping, current_variant_ids, no_gap_end,min_exon_len,min_intron_len, variant_cache);
						
						
						
						//Keep best scores						
						for (int zz=0; zz<nr_paths_par;zz++){								
							for (int z=0; z<nr_paths_par;z++){								
								double priorScore= ((Prev_score*)current_seed->matrices[z] +matrix_position)->value ; 	    
								if(priorScore+seed_matrix[seed_index]->best_scores[zz] - variant_penalty > current_seed->best_scores[nr_paths_par-1]){											
									//fprintf(stdout, "Found better score with deletion: %i-%i (id=%i) Previous:%i-%i Next:%i-%i\n",j-prev_shift,endpositions[d],idsdeletions[d],i,j,i-prev_shift,jj);
									current_seed->best_scores[nr_paths_par-1]=priorScore+seed_matrix[seed_index]->best_scores[zz] - variant_penalty;
									current_seed->best_score_pos[nr_paths_par-1]->read_pos=i;
									current_seed->best_score_pos[nr_paths_par-1]->dna_pos=j;
									current_seed->best_score_pos[nr_paths_par-1]->num_gaps=((Prev_score*)current_seed->matrices[z] +matrix_position)->num_gaps + seed_matrix[seed_index]->best_score_pos[zz]->num_gaps;
									current_seed->best_score_pos[nr_paths_par-1]->num_mm=((Prev_score*)current_seed->matrices[z] +matrix_position)->num_mismatches+ seed_matrix[seed_index]->best_score_pos[zz]->num_mm;
									current_seed->best_score_pos[nr_paths_par-1]->num_introns=seed_matrix[seed_index]->best_score_pos[zz]->num_introns;					
									current_seed->best_score_pos[nr_paths_par-1]->next_seed=seed_matrix[seed_index];
									current_seed->best_score_pos[nr_paths_par-1]->path_number=zz;
									current_seed->best_score_pos[nr_paths_par-1]->path_number_matrices=z;
									current_seed->best_score_pos[nr_paths_par-1]->partial_score=priorScore - variant_penalty ;
									// Resort best_score_pos and best_scores
									sort_best_scores(current_seed, nr_paths_par);
								}
								else{
									break;
								}
							}
						}
						
					}
					endpositions.clear();
					tabidsdeletions.clear();

					
				}/* END VARIANT DELETION */

				

				/******************************************************************************************/
				/* 3. POSITION FILLED WITH A MISMATCH OR NOT FILLED                                       */
				/* previous position in the diagonal can lead to a gap on read or dna instead of mismatch */
				/* Also look at possible splice sites stored from this diagonal                           */
				/******************************************************************************************/
				if (disabled_diagonal[nj] || !is_a_match(dnaInt,readChar)){
					//if (possible_sites.size()>10)
					//fprintf(stdout, "possible_sites.size()=%i\n", (int)possible_sites.size()) ;
					/* A. Explore splice sites */

					if(!remapping){
						for(int ss=(int)possible_sites.size()-1;ss>=0;ss--){
	    
							//Position of the splice site 'G' of 'GT/C' or 'G' of 'AG' on DNA
							int posj=((splice_pos*)possible_sites[ss])->site;

							//Last position before splice site
							int posi=possible_sites[ss]->i;
							int matrix_pos=((splice_pos*)possible_sites[ss])->matrix_pos;

							int number=possible_sites[ss]->number;
							int* pmatrices=possible_sites[ss]->matrices;
	    
							
							if (matrix_pos%row_len==nj){

								prevGaps=((Prev_score*)matrices[0] + matrix_pos)->num_gaps;
								prevMism=((Prev_score*)matrices[0] + matrix_pos)->num_mismatches;
								
								int left_i;
								if(right_side)
									left_i=read_len-1-posi;
								else
									left_i=posi;
								if (left_i+(max_gap-prevGaps)<min_exon_len)
									continue;
								
								//Figure out possible complementary splice sites if not done before
								if (first_seed && !comp_sites_filled){		
									comp_sites_filled=true;
									int first_ss=((splice_pos*)possible_sites[0])->site;
									if(right_side){
										for(int comp_ss=first_ss+2;comp_ss<dna_len-1;comp_ss++){
											if (comp_site_scores[comp_ss]>-ALMOST_INFINITY && comp_site_scores[comp_ss]<ALMOST_INFINITY)
												comp_sites.push_back(comp_ss);
										}
									}
									else{
										for(int comp_ss=first_ss-2;comp_ss>0;comp_ss--){
											if (comp_site_scores[comp_ss]>-ALMOST_INFINITY && comp_site_scores[comp_ss]<ALMOST_INFINITY)
												comp_sites.push_back(comp_ss);
										}
									}
								}


								//Compute minimum number of matches needed with relaxing the constraint at the ends of the read
								int diff_g_i;
								if(right_side)
								{
									diff_g_i=read_len-1-posi;
									if (posi==read_len-1 || diff_g_i>min_match)
										diff_g_i=min_match;
								}
								else
								{
									diff_g_i=posi;
									if (posi==0 || posi>min_match)
										diff_g_i=min_match;
								}

								//Alignment from each complementary splice site
								for(int comp_ss=0; comp_ss<(int)comp_sites.size();comp_ss++){
									
									
									
									// Next complementary splice site relative to the current splice site
									int jj=comp_sites[comp_ss]; 
									if (-prev_shift*jj<=posj*-prev_shift)
										continue;
									
									//Too short intron
									if (((right_side && (jj-posj+1<min_intron_len)) || (!right_side && (posj-jj+1<min_intron_len))))
										continue;

									//Check at least diff_i consecutive matches after this complementary splice site
									int diff_i=diff_g_i;
									bool conserved_seq=true;	      
									int conserved_seq_mismatches=0 ;
									int ii=posi;
									int num_N=0 ;
									int num_mismatch_score=0 ;
									int num_mismatch_extend=0 ;
									for(int nm=1; nm<=diff_i;nm++)
									{
										ii=ii-prev_shift; 
										jj=jj-prev_shift; // One position after 'G' of 'AG' or before 'G' of 'GT/C'
										if ( (right_side && (ii>=read_len || jj>=dna_len)) || (!right_side && (ii<0 || jj<0)))
										{
											conserved_seq_mismatches++ ;
											conserved_seq=false;
											break;		   
										}
										else{
											//Take possible SNP into account
											char dna_char= getBestDnaChar<use_variants,snp_merged>(read[ii], dna[jj], jj, variant_cache);
											if (read[ii]!=dna_char) 
											{
												if (check_char(read[ii])!=5 && check_char(dna_char)!=5 && read_scores[ii]>=MAX_SPLICE_MISMATCH_QUAL_SINGLE)
												{
													conserved_seq=false;
													break;
												}
												else
												{
													diff_i++ ;
													if (check_char(read[ii])==5 || check_char(dna_char)==5)
														num_N++ ;
													else
													{
														num_mismatch_score+=read_scores[ii] ;
														num_mismatch_extend++ ;
													}
													if (num_N>MAX_SPLICE_MISMATCH_NUM_N || num_mismatch_score>MAX_SPLICE_MISMATCH_QUAL_TOTAL || num_mismatch_extend>MAX_SPLICE_MISMATCH_QUAL_EXTEND)
													{
														conserved_seq=false;
														break;  
													}
												}
											}
										}
									}

									if (conserved_seq){

										int seed_already_filled=-1;
										bool continue_searching=true;
										double tempSplicedScore = main_site_scores[posj]+comp_site_scores[comp_sites[comp_ss]];
										double currentScore= ((Prev_score*)current_seed->matrices[pmatrices[0]] +matrix_pos)->value ; 	    

										for(int num_seed=0; num_seed< (int)seed_matrix.size();num_seed++){
										
											if(seed_matrix[num_seed]->read_pos==posi-prev_shift && 
											   seed_matrix[num_seed]->dna_pos==comp_sites[comp_ss]-prev_shift && seed_matrix[num_seed]->max_introns==max_number_introns-1 && seed_matrix[num_seed]->deletion_id.empty())
											{
								
												if ((seed_matrix[num_seed]->best_score_pos[0]->num_gaps>max_gap-prevGaps && seed_matrix[num_seed]->best_score_pos[0]->num_mm>=max_mism-prevMism)||
													(seed_matrix[num_seed]->best_score_pos[0]->num_gaps>=max_gap-prevGaps && seed_matrix[num_seed]->best_score_pos[0]->num_mm>max_mism-prevMism)){
													continue_searching=false;
													continue;												
												}
											
												if ((seed_matrix[num_seed]->max_gaps>=max_gap-prevGaps && seed_matrix[num_seed]->max_mm>=max_mism-prevMism)&&
													(seed_matrix[num_seed]->best_score_pos[0]->num_gaps<=max_gap-prevGaps && seed_matrix[num_seed]->best_score_pos[0]->num_mm<=max_mism-prevMism)){
												
													if (seed_matrix[num_seed]->best_scores[0] > -ALMOST_INFINITY && tempSplicedScore+currentScore >= seed_matrix[num_seed]->best_prev_score){
														continue_searching=true;
														seed_already_filled=num_seed;
														break;
													}
													else{													
														continue_searching=false;
														break;
													}												
												}

												if (seed_matrix[num_seed]->max_gaps < max_gap-prevGaps && seed_matrix[num_seed]->max_mm < max_mism-prevMism){		
													if (tempSplicedScore+currentScore >= seed_matrix[num_seed]->best_prev_score){
														continue_searching=true;
														continue;
													}
													else{													
														continue_searching=false;
														break;
													}												
												}
											}
										}
									
														  
										// if (!continue_searching)									   
										// 	fprintf(stdout,"**1** Don't fill this seed matrix and don't compare results\n");
										// else{
										// 	if (seed_already_filled!=-1)
										// 		fprintf(stdout,"**2** Use already filled matrix\n");
										// 	//else
										// 		//fprintf(stdout,"**3** New seed matrix\n");
										// }
									
										if (continue_searching){
										
											if(seed_already_filled==-1){
												//Number of this seedElem in the vector seedMatrix (because of recursive calls that add new seedElem)
												seed_already_filled=seed_matrix.size();
												std::vector<FoundVariant> vtemp;
												fast_fill_side_unspliced_first<use_variants,snp_merged,right_side,nr_paths_par>(seed_matrix, read_len, dna_len, read, dna, prb, functions, matchmatrix, qualityScores, 
																												   main_site_scores, comp_site_scores, comp_sites, posi-prev_shift,
																												   comp_sites[comp_ss]-prev_shift, best_match_scores, false,
																												   max_number_introns-1,max_number_deletions,max_gap-prevGaps,
																												   max_mism-prevMism,max_edit_op-(prevGaps+prevMism),min_match, 
																												   verbosity,currentMode, remapping,vtemp, no_gap_end,min_exon_len,min_intron_len,
																												   variant_cache);
											}

											//Keep best scores
											//double tempSplicedScore = main_site_scores[posj]+comp_site_scores[comp_sites[comp_ss]];
											//fprintf(stdout,"Intron score %f %f \n",main_site_scores[posj],comp_site_scores[comp_sites[comp_ss]]);
		  
											for (int z=0; z<nr_paths_par;z++){
												for (int zz=0; zz<number;zz++){
													double priorScore= ((Prev_score*)current_seed->matrices[pmatrices[zz]] +matrix_pos)->value ; 	    
													//fprintf(stdout,"prior score %f new score %f splice score %f\n",priorScore,seed_matrix[seed_already_filled]->best_scores[z],tempSplicedScore);
												
													if(tempSplicedScore +priorScore>seed_matrix[seed_already_filled]->best_prev_score)
														seed_matrix[seed_already_filled]->best_prev_score=tempSplicedScore +priorScore;
												
													if(tempSplicedScore +priorScore+seed_matrix[seed_already_filled]->best_scores[z] > current_seed->best_scores[nr_paths_par-1]){											
														current_seed->best_scores[nr_paths_par-1]=tempSplicedScore+priorScore+seed_matrix[seed_already_filled]->best_scores[z];
														current_seed->best_score_pos[nr_paths_par-1]->read_pos=posi;
														current_seed->best_score_pos[nr_paths_par-1]->dna_pos=posj+prev_shift;
														current_seed->best_score_pos[nr_paths_par-1]->num_gaps=((Prev_score*)current_seed->matrices[pmatrices[zz]] +matrix_pos)->num_gaps + seed_matrix[seed_already_filled]->best_score_pos[z]->num_gaps;
														current_seed->best_score_pos[nr_paths_par-1]->num_mm=((Prev_score*)current_seed->matrices[pmatrices[zz]] +matrix_pos)->num_mismatches+ seed_matrix[seed_already_filled]->best_score_pos[z]->num_mm;
														current_seed->best_score_pos[nr_paths_par-1]->num_introns=seed_matrix[seed_already_filled]->best_score_pos[z]->num_introns +1;					
														current_seed->best_score_pos[nr_paths_par-1]->next_seed=seed_matrix[seed_already_filled];
														current_seed->best_score_pos[nr_paths_par-1]->path_number=z;
														current_seed->best_score_pos[nr_paths_par-1]->path_number_matrices=pmatrices[zz];
														current_seed->best_score_pos[nr_paths_par-1]->partial_score=priorScore+tempSplicedScore;
														// Resort best_score_pos and best_scores
														sort_best_scores(current_seed, nr_paths_par);
													}
													else{
														break;
													}
												}
											}
										}
									}
								
								}//End look through complementary splice sites

								delete[] possible_sites[ss]->matrices;
								possible_sites[ss]->matrices=NULL;
								delete possible_sites[ss];
								if (ss!=(int)possible_sites.size()-1){
									possible_sites[ss]=possible_sites[possible_sites.size()-1];
								}
								possible_sites.pop_back();
							}
						}
						/* A.END */
					}//remapping
				
	  
					/* B. Allow gaps from the previous position in the diagonal */
					matrix_prev_position= matrix_position-row_len;
	  
					
				   
					/* B.1. Gap on DNA sequence */
					int leftover_i;
					if (right_side){
						leftover_i= read_len -i;
					}
					else{
						leftover_i= i+1;
					}
					if (nj!=0 && matrix_prev_position>=0 && leftover_i > no_gap_end){
						dnaChar = check_char(dna[j+prev_shift]);
						readChar = check_char(read[i]);
						assert(dnaChar!=-1 && readChar!=-1);

						if (currentMode == USE_QUALITY_SCORES)
							baseScore = read_scores[i];
	    
						// Best score of what it leaves to align
						if (ni<i_len-1){
							if(right_side)
								putativeValue= best_match_scores[i+1];
							else
								putativeValue= best_match_scores[0]-best_match_scores[i];
						}
						else
							putativeValue=0;
	    
						//Worst best global score known
						globalValue= current_seed->best_scores[nr_paths_par-1];
						//fprintf(stdout,"Putative value: %f Best score: %f\n",putativeValue,globalValue);
	    

						for(int z=0;z<nr_paths_par;z++){
							Prev_score* actMatrix = (Prev_score*)matrices[z]; 

							//Information about the previous position for path z
							prevGaps=((Prev_score*)actMatrix + matrix_prev_position)->num_gaps;
							prevMism=((Prev_score*)actMatrix + matrix_prev_position)->num_mismatches;
							prevValue = ((Prev_score*)actMatrix +matrix_prev_position)->value ;

							//Gap possible according to the number of gaps and mismatches at matrix_prev_position
							if (prevGaps<max_gap && prevGaps+prevMism<max_edit_op){

								if (currentMode == USE_QUALITY_SCORES)
									tempValue = prevValue + getScore(qualityScores,mlen,0,readChar,baseScore);
								else
									tempValue = prevValue +(matchmatrix[readChar]); /* score(READ,gap) */

								if (isnotminusinf(tempValue)&& tempValue> ((Prev_score*)matrices[nr_paths_par-1] +matrix_position-1)->value && tempValue+putativeValue>globalValue){
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-1)->value = tempValue;
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-1)->prev_i = i+prev_shift; /* predecessor */
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-1)->prev_j = j+prev_shift; /* predecessor */
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-1)->prev_matrix_no = z;
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-1)->num_matches = 0;
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-1)->num_mismatches = prevMism;
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-1)->num_gaps = prevGaps+1;
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-1)->snp_id =-1;
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-1)->snp_int =-1;
									disabled_diagonal[nj-1]=false;
									sort_best_paths(matrices,nr_paths_par,matrix_position-1);

									//Last row to fill: check now if unspliced alignment is better
									if(ni==i_len-1){
										if (tempValue>current_seed->best_scores[nr_paths_par-1]){
											current_seed->best_scores[nr_paths_par-1]=tempValue;
											current_seed->best_score_pos[nr_paths_par-1]->read_pos=i;
											current_seed->best_score_pos[nr_paths_par-1]->dna_pos=j+prev_shift;
											current_seed->best_score_pos[nr_paths_par-1]->num_gaps=((Prev_score*)actMatrix + matrix_position-1)->num_gaps;									
											current_seed->best_score_pos[nr_paths_par-1]->num_mm=((Prev_score*)actMatrix + matrix_position-1)->num_mismatches;
											current_seed->best_score_pos[nr_paths_par-1]->num_introns=0;
											current_seed->best_score_pos[nr_paths_par-1]->next_seed=NULL;
											current_seed->best_score_pos[nr_paths_par-1]->path_number=0;
											current_seed->best_score_pos[nr_paths_par-1]->path_number_matrices=z;
											current_seed->best_score_pos[nr_paths_par-1]->partial_score=tempValue;
											// Resort best_score_pos and best_scores
											sort_best_scores(current_seed, nr_paths_par);
											//worst best global score known
											globalValue= current_seed->best_scores[nr_paths_par-1];
										}
									}
									
									//Deletion variants from this position
									//Some read still need to be aligned
									if ((right_side && i-prev_shift <read_len) || (!right_side && i-prev_shift >=0)){
										
										//Search for deletions from next DNA position
										std::vector<int> endpositions;
										std::vector<FoundVariant> idsdeletions;
//										fprintf(stdout,"position after gap dna: %i\n",j);
										if (max_number_deletions>=1)
											idsdeletions = getDeletionsfromVariants<use_variants>(j,endpositions, variant_cache,right_side);
										std::vector<std::vector<FoundVariant> > tabidsdeletions;
										//Prepare to allow to jump over several consecutive deletions
										for (int d=0;d<(int)idsdeletions.size();d++){
											std::vector<FoundVariant> vtemp;
											vtemp.push_back(idsdeletions[d]);
											tabidsdeletions.push_back(vtemp);
											vtemp.clear();
										}
										idsdeletions.clear();
										while(endpositions.size()>0){
						
											//Start DNA position for next seed is one after/before the end position of the deletion
											int jj=endpositions[endpositions.size()-1] -prev_shift;
											endpositions.pop_back();
											std::vector<FoundVariant> current_variant_ids;
											current_variant_ids=tabidsdeletions[tabidsdeletions.size()-1];
											tabidsdeletions.pop_back();

											if ( (right_side && jj>=dna_len) || (!right_side && jj<0))
												continue;
											
											if (PERFORM_EXTRA_CHECKS)
											{
												assert(jj>=0 && jj<dna_len);
												assert(i-prev_shift>=0 && i-prev_shift<read_len);
												if (right_side)
													assert (jj>j);
											}
											
											if (max_number_deletions>=1)
												idsdeletions = getDeletionsfromVariants<use_variants>(jj,endpositions, variant_cache,right_side);
											for (int d=0;d<(int)idsdeletions.size();d++){
												std::vector<FoundVariant> vtemp=current_variant_ids;
												vtemp.push_back(idsdeletions[d]);
												tabidsdeletions.push_back(vtemp);
												vtemp.clear();
											}	
											idsdeletions.clear();						
											//Maybe here look if the seed was already filled but not sure...
//											fprintf(stdout, "New fill matrix from deletion (from dna gap)(%i): %i-%i\n",(int)seed_matrix.size(),j-prev_shift,jj+prev_shift);
											
											prevMism=((Prev_score*)actMatrix + matrix_position-1)->num_mismatches;
											prevGaps=((Prev_score*)actMatrix + matrix_position-1)->num_gaps;
											int seed_index=seed_matrix.size();
											fast_fill_side_unspliced_first<use_variants,snp_merged,right_side,nr_paths_par>(seed_matrix, read_len, dna_len, read, dna, prb, functions, matchmatrix, qualityScores, 
																											   main_site_scores, comp_site_scores, comp_sites, i-prev_shift,
																											   jj, best_match_scores, false,
																											   max_number_introns,max_number_deletions-1,max_gap-prevGaps,
																											   max_mism-prevMism,max_edit_op-(prevGaps+prevMism),min_match, verbosity,currentMode, 
																											   remapping,current_variant_ids, no_gap_end,min_exon_len,min_intron_len, variant_cache);
											
											
											//Keep best scores						
											for (int zz=0; zz<nr_paths_par;zz++){								
												double priorScore= ((Prev_score*)actMatrix +matrix_position-1)->value ; 	    
												if(priorScore+seed_matrix[seed_index]->best_scores[zz] - variant_penalty > current_seed->best_scores[nr_paths_par-1]){											
													//fprintf(stdout, "Found better score with deletion: %i-%i (id=%i) Previous:%i-%i Next:%i-%i\n",j-prev_shift,endpositions[d],idsdeletions[d],i,j,i-prev_shift,jj);
													current_seed->best_scores[nr_paths_par-1]=priorScore+seed_matrix[seed_index]->best_scores[zz] - variant_penalty ;
													current_seed->best_score_pos[nr_paths_par-1]->read_pos=i;
													current_seed->best_score_pos[nr_paths_par-1]->dna_pos=j+prev_shift;
													current_seed->best_score_pos[nr_paths_par-1]->num_gaps=prevGaps + seed_matrix[seed_index]->best_score_pos[zz]->num_gaps;
													current_seed->best_score_pos[nr_paths_par-1]->num_mm=prevMism+ seed_matrix[seed_index]->best_score_pos[zz]->num_mm;
													current_seed->best_score_pos[nr_paths_par-1]->num_introns=seed_matrix[seed_index]->best_score_pos[zz]->num_introns;					
													current_seed->best_score_pos[nr_paths_par-1]->next_seed=seed_matrix[seed_index];
													current_seed->best_score_pos[nr_paths_par-1]->path_number=zz;
													current_seed->best_score_pos[nr_paths_par-1]->path_number_matrices=z;
													current_seed->best_score_pos[nr_paths_par-1]->partial_score=priorScore - variant_penalty ;
													// Resort best_score_pos and best_scores
													sort_best_scores(current_seed, nr_paths_par);
												}
												else{
													break;
												}
											}
										}
										endpositions.clear();
										tabidsdeletions.clear();
										
										
					
									}/* END VARIANT DELETION */
								}
							}
						}
					}/* B.1. END */
	  
					/* B.2. Gap on READ sequence */
					if (right_side){
						leftover_i= read_len -i+1;
					}
					else{
						leftover_i= i;
					}
					if (nj!=row_len-1 && matrix_prev_position>=0 && leftover_i>no_gap_end){
						dnaChar = check_char(dna[j]);
						readChar = check_char(read[i+prev_shift]);
						assert(dnaChar!=-1 && readChar!=-1);

						if (currentMode == USE_QUALITY_SCORES)
							baseScore = read_scores[i+prev_shift];

	    
						// Best score of what it leaves to align
						if (ni<i_len){
							if(right_side)
								putativeValue= best_match_scores[i];
							else
								putativeValue= best_match_scores[0]-best_match_scores[i+prev_shift];
						}
						else
							putativeValue=0;

						//Worst best global score known
						globalValue= current_seed->best_scores[nr_paths_par-1];

						for(int z=0;z<nr_paths_par;z++){
							Prev_score* actMatrix = (Prev_score*)matrices[z]; 

							//Information about the previous position for path z
							prevGaps=((Prev_score*)actMatrix + matrix_prev_position)->num_gaps;
							prevMism=((Prev_score*)actMatrix + matrix_prev_position)->num_mismatches;
							prevValue = ((Prev_score*)actMatrix +matrix_prev_position)->value ;

							//Gap possible according to the number of gaps and mismatches at matrix_prev_position
							if (prevGaps<max_gap && prevGaps+prevMism<max_edit_op){
	    
								tempValue = prevValue + getBestGapWithVariants<use_variants,snp_merged>(currentMode, matchmatrix, qualityScores, mlen, dna[j],j, dnaInt,snp_id, variant_cache);
								
								if (isnotminusinf(tempValue)&& tempValue> ((Prev_score*)matrices[nr_paths_par-1] +matrix_position-row_len+1)->value && tempValue+putativeValue>globalValue){
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-row_len+1)->value = tempValue;
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-row_len+1)->prev_i = i+prev_shift; /* predecessor */
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-row_len+1)->prev_j = j+prev_shift; /* predecessor */
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-row_len+1)->prev_matrix_no = z;
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-row_len+1)->num_matches = 0;
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-row_len+1)->num_mismatches = prevMism;
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-row_len+1)->num_gaps = prevGaps+1;
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-row_len+1)->snp_id =snp_id;
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-row_len+1)->snp_int =dnaInt;
									disabled_diagonal[nj+1]=false;
									sort_best_paths(matrices,nr_paths_par,matrix_position-row_len+1);

									//Deletion variants from this position
									//Some read still need to be aligned
									if ((right_side && j-prev_shift <dna_len) || (!right_side && j-prev_shift >=0)){
										
										//Search for deletions from next DNA position
										std::vector<int> endpositions;
										std::vector<FoundVariant> idsdeletions;
//										fprintf(stdout,"position after gap read: %i\n",j-prev_shift);
										if (max_number_deletions>=1)
											idsdeletions = getDeletionsfromVariants<use_variants>(j-prev_shift,endpositions, variant_cache,right_side);
										std::vector<std::vector<FoundVariant> > tabidsdeletions;
										//Prepare to allow to jump over several consecutive deletions
										for (int d=0;d<(int)idsdeletions.size();d++){
											std::vector<FoundVariant> vtemp;
											vtemp.push_back(idsdeletions[d]);
											tabidsdeletions.push_back(vtemp);
											vtemp.clear();
										}
										idsdeletions.clear();
										while(endpositions.size()>0){
						
											//Start DNA position for next seed is one after/before the end position of the deletion
											int jj=endpositions[endpositions.size()-1] -prev_shift;
											endpositions.pop_back();
											std::vector<FoundVariant> current_variant_ids;
											current_variant_ids=tabidsdeletions[tabidsdeletions.size()-1];
											tabidsdeletions.pop_back();

											if ( (right_side && jj>=dna_len) || (!right_side && jj<0))
												continue;

											if (PERFORM_EXTRA_CHECKS)
											{
												assert(jj>=0 && jj<dna_len);
												assert(i>=0 && i<read_len);
												if (right_side)
													assert (jj>j-prev_shift);
											}

											if (max_number_deletions>=1)
												idsdeletions = getDeletionsfromVariants<use_variants>(jj,endpositions, variant_cache,right_side);
											for (int d=0;d<(int)idsdeletions.size();d++){
												std::vector<FoundVariant> vtemp=current_variant_ids;
												vtemp.push_back(idsdeletions[d]);
												tabidsdeletions.push_back(vtemp);
												vtemp.clear();
											}	
											idsdeletions.clear();						
											//Maybe here look if the seed was already filled but not sure...
//											fprintf(stdout, "New fill matrix from deletion (from dna gap)(%i): %i-%i\n",(int)seed_matrix.size(),j-prev_shift,jj+prev_shift);
											
											prevMism=((Prev_score*)actMatrix + matrix_position-row_len+1)->num_mismatches;
											prevGaps=((Prev_score*)actMatrix + matrix_position-row_len+1)->num_gaps;
											int seed_index=seed_matrix.size();
											fast_fill_side_unspliced_first<use_variants,snp_merged,right_side,nr_paths_par>(seed_matrix, read_len, dna_len, read, dna, prb, functions, matchmatrix, qualityScores, 
																											   main_site_scores, comp_site_scores, comp_sites, i,
																											   jj, best_match_scores, false,
																											   max_number_introns,max_number_deletions-1,max_gap-prevGaps,
																											   max_mism-prevMism,max_edit_op-(prevGaps+prevMism),min_match, verbosity,currentMode, 
																											   remapping,current_variant_ids,no_gap_end,min_exon_len,min_intron_len, variant_cache);
											
											
											//Keep best scores						
											for (int zz=0; zz<nr_paths_par;zz++){								
												double priorScore= ((Prev_score*)actMatrix +matrix_position-row_len+1)->value ; 	    
												if(priorScore+seed_matrix[seed_index]->best_scores[zz] - variant_penalty > current_seed->best_scores[nr_paths_par-1]){											
													//fprintf(stdout, "Found better score with deletion: %i-%i (id=%i) Previous:%i-%i Next:%i-%i\n",j-prev_shift,endpositions[d],idsdeletions[d],i,j,i-prev_shift,jj);
													current_seed->best_scores[nr_paths_par-1]=priorScore+seed_matrix[seed_index]->best_scores[zz] - variant_penalty;
													current_seed->best_score_pos[nr_paths_par-1]->read_pos=i+prev_shift;
													current_seed->best_score_pos[nr_paths_par-1]->dna_pos=j;
													current_seed->best_score_pos[nr_paths_par-1]->num_gaps=prevGaps + seed_matrix[seed_index]->best_score_pos[zz]->num_gaps;
													current_seed->best_score_pos[nr_paths_par-1]->num_mm=prevMism+ seed_matrix[seed_index]->best_score_pos[zz]->num_mm;
													current_seed->best_score_pos[nr_paths_par-1]->num_introns=seed_matrix[seed_index]->best_score_pos[zz]->num_introns;					
													current_seed->best_score_pos[nr_paths_par-1]->next_seed=seed_matrix[seed_index];
													current_seed->best_score_pos[nr_paths_par-1]->path_number=zz;
													current_seed->best_score_pos[nr_paths_par-1]->path_number_matrices=z;
													current_seed->best_score_pos[nr_paths_par-1]->partial_score=priorScore - variant_penalty;
													// Resort best_score_pos and best_scores
													sort_best_scores(current_seed, nr_paths_par);
												}
												else{
													break;
												}
											}
										}
										
										endpositions.clear();
										tabidsdeletions.clear();
									}/* END VARIANT DELETION */
									

								}
							}
						}
					}/* B.2. END */
				}/* 3. END */

			}//end j
	
			bool is_bad_row=true;
			for(int t=0;t<row_len;t++)
				is_bad_row=is_bad_row && (disabled_diagonal[t]==1);
			if (is_bad_row){
				break;
			}
		}//end i
	}
  
	////Print best scores
	//for (int z=0; z<nr_paths_par; z++)
	//   fprintf(stdout,"%i best score: %f\n",z+1,current_seed->best_scores[z]);
	//if(first_seed)
	//   print_restricted_matrix(matrices,nr_paths_par,matrix_len,row_len);
  

	// Free memory
	current_seed=NULL;
	for(int s=0;s<(int)possible_sites.size();s++){
		delete[] possible_sites[s]->matrices;
		possible_sites[s]->matrices=NULL;
		delete possible_sites[s];
	}
	possible_sites.clear();
	delete[] possible_matrices;  
	delete[] disabled_diagonal;

	//fprintf(stdout,"END: Fill %s side of the matrix from position %i-%i (%i,%i,%i,num_intron=%i)...\n",(right_side==true)?"right":"left",seed_read, seed_dna,max_gap,max_mism,max_edit_op, max_number_introns);

	//  fprintf(stdout,"Fill a %s side of the matrix from position %i-%i...END\n",right_side?"right":"left",seed_read, seed_dna);
}


template<bool use_variants, bool snp_merged>
void fast_fill_matrix(int nr_paths_par, int*max_score_positions, int read_len, int dna_len, char* read, char* dna, double* prb, penalty_struct* functions, 
					  double* matchmatrix, penalty_struct* qualityScores, double* donor, double* acceptor, bool remove_duplicate_scores,int seed_i, int seed_j, 
					  std::vector<SeedElem *>& seed_matrix_left, std::vector<SeedElem *>& seed_matrix_right, int max_number_introns, int max_number_deletions,
					  int max_gap, int max_mism, int max_edit_op, int min_match, int verbosity,mode currentMode, bool remapping,
					  int no_gap_end,int min_exon_len,int min_intron_len, const std::vector<variant_cache_t *> &variant_cache)
{
  
	const int MMATRIX_LEN = 6; // length of matchmatrix
  
	//fprintf(stdout,"Max number of exons %i\n",max_number_introns);

	//printf("Entering fill_matrix...\n");
  
	/*********************************************************************************************/
	/*Best score for a matching alignment from a position i in read sequence   */
	/*********************************************************************************************/
	// fprintf(stdout,"Best match scores...\n");
	double* best_match_scores= new double[read_len];
	double temp_best=0;
	double *read_scores = prb;
  
	for(int i=read_len-1; i >=0;i--){
		
		double score;
		
		if (currentMode == USE_QUALITY_SCORES)
			score=getScore(qualityScores,MMATRIX_LEN,check_char(read[i]),check_char(read[i]),read_scores[i]);
		else
			score=matchmatrix[MMATRIX_LEN*check_char(read[i])+check_char(read[i])];
		//fprintf(stdout,"match score at position %i: %f\n",i,score);
		temp_best+=score;
		best_match_scores[i]=temp_best;
    
	}
  
	//fprintf(stdout,"Best match scores...END\n");
	
	
	/*********************************************************************************************/
	/*Left and right alignments */
	/*********************************************************************************************/
	//fprintf(stdout,"Left and right alignments from %i-%i...\n",seed_i,seed_j);

	assert(nr_paths_par==1) ;

	std::vector<int> comp_sites;
	std::vector<FoundVariant> vtemp;
	fast_fill_side_unspliced_first<use_variants,snp_merged,true,1>(seed_matrix_right,read_len,dna_len, read, dna, prb,functions, matchmatrix,qualityScores, donor,acceptor,comp_sites,seed_i, 
																 seed_j,best_match_scores,true,max_number_introns,max_number_deletions, max_gap,max_mism,max_edit_op,min_match, 
																 verbosity,currentMode,remapping, vtemp,
																 no_gap_end,min_exon_len,min_intron_len, variant_cache);
	//fprintf(stdout,"%ld right sides of the matrix filled...\n",seed_matrix_right.size());
	// for(int n=0;n<seed_matrix_right.size();n++){
	// 	if (((SeedElem*)seed_matrix_right[n])!=NULL)
	// 		fprintf(stdout,"seed position %i %i %f\n",((SeedElem*)seed_matrix_right[n])->read_pos,((SeedElem*)seed_matrix_right[n])->dna_pos,((SeedElem*)seed_matrix_right[n])->best_scores[0]);
	// }
 
	comp_sites.clear();
	fast_fill_side_unspliced_first<use_variants,snp_merged,false,1>(seed_matrix_left,read_len,dna_len, read, dna, prb,functions, matchmatrix,qualityScores, acceptor,donor,comp_sites,seed_i, 
																  seed_j,best_match_scores,true,max_number_introns,max_number_deletions, 
																  max_gap,max_mism,max_edit_op,min_match, verbosity,currentMode,remapping, vtemp,
																  no_gap_end,min_exon_len,min_intron_len, variant_cache);
	comp_sites.clear();
	//fprintf(stdout,"%ld left sides of the matrix filled...\n",seed_matrix_left.size());
  // for(int n=0;n<seed_matrix_left.size();n++){
  //   if (((SeedElem*)seed_matrix_left[n])!=NULL)
  //     fprintf(stdout,"seed position %i %i %f\n",((SeedElem*)seed_matrix_left[n])->read_pos,((SeedElem*)seed_matrix_left[n])->dna_pos,((SeedElem*)seed_matrix_left[n])->best_scores[0]);
  // }

	//fprintf(stdout,"Left and right alignments...END\n");

 
	/*********************************************************************************************/
	/* Find out the nr_paths_par best combinations */
	/*********************************************************************************************/

	memset(max_score_positions,0,sizeof(int)*2*nr_paths_par);
	if(seed_matrix_left.size()>=1 && seed_matrix_left[0]!=NULL && seed_matrix_right.size()>=1 && seed_matrix_right[0]!=NULL)
	{ 
		int best_left_left=1;
		int best_left_right=0;
		int best_right_left=0;
		int best_right_right=1;

		max_score_positions[0]=0;
		max_score_positions[1]=0;

		for(int z=1;z<nr_paths_par;z++){
			double temp_left= seed_matrix_left[0]->best_scores[best_left_left] + seed_matrix_right[0]->best_scores[best_left_right];
			double temp_right= seed_matrix_left[0]->best_scores[best_right_left] + seed_matrix_right[0]->best_scores[best_right_right];

			if (temp_left > temp_right){
				max_score_positions[2*z]=best_left_left;
				max_score_positions[2*z+1]=best_left_right;
				best_left_left++;
			}
			else{
				max_score_positions[2*z]=best_right_left;
				max_score_positions[2*z+1]=best_right_right;
				best_right_right++;
			}
		}
		//fprintf(stdout, "best_left_left=%i, best_left_right=%i\n", best_left_left, best_left_right) ;
	}
	/*else
		if (seed_matrix_left[0]!=NULL)
			fprintf(stdout, "seed_matrix_left[0]!=NULL\n") ;
		else
		fprintf(stdout, "seed_matrix_right[0]!=NULL\n") ;*/
	
	/*********************************************************************************************/
	/* Display results */
	/*********************************************************************************************/

	// for (int z=0;z<nr_paths_par;z++){

	//   int z_path_left=max_score_positions[2*z] ; //path number for the left seed matrix
	//   int z_path_right=max_score_positions[2*z+1] ; //path number for the right seed matrix

	//   fprintf(stdout,"Align result for %i left matrix and %i right matrix \n",z_path_left,z_path_right);
  
	//   if(seed_matrix_left[0]!=NULL && seed_matrix_right[0]!=NULL && 
	//      seed_matrix_left[0]->best_scores[z_path_left]>-ALMOST_INFINITY && seed_matrix_right[0]->best_scores[z_path_right]>-ALMOST_INFINITY){

	//     fprintf(stdout,"score: %f\n",seed_matrix_left[0]->best_scores[z_path_left] + seed_matrix_right[0]->best_scores[z_path_right]);
	//     fprintf(stdout, "seed score %f\n",getScore(qualityScores,MMATRIX_LEN,check_char(read[seed_i]),check_char(dna[seed_j]),read_scores[seed_i]));

	//     SeedElem *next_seed=seed_matrix_left[0];
	//     int zz=z_path_left;

	//     int rstart;
	//     int dstart;
    
	//     std::vector<char> read_align;
	//     std::vector<char> dna_align;
	//     int max_gaps2;

	//     while (next_seed!=NULL){
	// 		std::vector<char> read_align_temp;
	// 		std::vector<char> dna_align_temp;      
	// 		max_gaps2=next_seed->max_gaps;
	// 		Prev_score** matrix=next_seed->matrices;
	// 		rstart=next_seed->best_score_pos[zz]->read_pos;
	// 		dstart=next_seed->best_score_pos[zz]->dna_pos;
	// 		int rseed=next_seed->read_pos;
	// 		int dseed=next_seed->dna_pos;
	// 		int prev_z= next_seed->best_score_pos[zz]->path_number_matrices;
			
	// 		while(!(rstart==rseed+1 && dstart==dseed+1)){
				
	// 			//fprintf(stdout,"%i-%i:%i-%i (%i)\n",rstart,dstart,rseed,dseed,prev_z);
	// 			int matrix_position= (rseed-rstart)*(max_gaps2*2+1)+(dseed-(rseed-rstart))-dstart+max_gaps2; 
	// 			int prev_i=((Prev_score*)matrix[prev_z]+matrix_position)->prev_i;
	// 			int prev_j=((Prev_score*)matrix[prev_z]+matrix_position)->prev_j;
	// 			prev_z=((Prev_score*)matrix[prev_z]+matrix_position)->prev_matrix_no;
	// 			//fprintf(stdout,"prev %i-%i:%i-%i (%i)\n",rstart,dstart,prev_i,prev_j, prev_z);
	// 			assert(rstart<=prev_i && dstart<=prev_j);
				
	// 			if (prev_i==rstart && prev_j==dstart+1){//read gap
	// 				dna_align_temp.push_back(dna[dstart]);
	// 				read_align_temp.push_back('-');	    
	// 			}
				
				
	// 			else if(prev_i==rstart+1 && prev_j==dstart){//dna gap
	// 				dna_align_temp.push_back('-');
	// 				read_align_temp.push_back(read[rstart]);
	// 			}
				
				
	// 			else if(prev_i==rstart+1 && prev_j==dstart+1){//match/mismatch
	// 				dna_align_temp.push_back(dna[dstart]);
	// 				read_align_temp.push_back(read[rstart]);
	// 			}
				
	// 			rstart=prev_i;
	// 			dstart=prev_j;
				
	// 		}
			
     
	// 		rstart=next_seed->best_score_pos[zz]->read_pos;
	// 		dstart=next_seed->best_score_pos[zz]->dna_pos;
	
	// 		int tmp=next_seed->best_score_pos[zz]->path_number;
	//     	next_seed=next_seed->best_score_pos[zz]->next_seed;
	// 		zz=tmp;
	
	// 		std::vector<char>::reverse_iterator rit;
	// 		for ( rit= dna_align_temp.rbegin() ; rit < dna_align_temp.rend(); ++rit )
	// 			dna_align.push_back(*rit);
	// 		for ( rit= read_align_temp.rbegin() ; rit < read_align_temp.rend(); ++rit )
	// 			read_align.push_back(*rit);
			
	// 		dna_align_temp.clear();
	// 		read_align_temp.clear();
	
	// 		if (next_seed!=NULL){
	// 			if (next_seed->deletion_id.empty()){
	// 				for(int n=dstart-1;n>=next_seed->dna_pos+1;n--){
	// 					dna_align.push_back(dna[n]);
	// 					read_align.push_back('*');
	// 				}
	// 			}
				
	// 			else{
	// 				for(int n=dstart-1;n>=next_seed->dna_pos+1;n--){
	// 					dna_align.push_back(dna[n]);
	// 					read_align.push_back('d');
	// 				}
	// 			}
				
	// 		} 	   
	// 	}
	  
   
	//     std::vector<char>::reverse_iterator rit;
	//     fprintf(stdout,"DNA  ALIGN: ");
	//     for ( rit= dna_align.rbegin() ; rit < dna_align.rend(); ++rit )
	// 		fprintf(stdout,"%c",*rit);
	//     fprintf(stdout,"\nREAD ALIGN: ");
	//     for ( rit= read_align.rbegin() ; rit < read_align.rend(); ++rit )
	// 		fprintf(stdout,"%c",*rit);
	//     fprintf(stdout,"\n");
      
	//     dna_align.clear();
	//     read_align.clear();
      
	//     next_seed=seed_matrix_right[0];
	//     zz=z_path_right;
		
	//     while (next_seed!=NULL){
	// 		std::vector<char> read_align_temp;
	// 		std::vector<char> dna_align_temp;      
	// 		max_gaps2=next_seed->max_gaps;      
	// 		Prev_score** matrix=next_seed->matrices;
	// 		rstart=next_seed->best_score_pos[zz]->read_pos;
	// 		dstart=next_seed->best_score_pos[zz]->dna_pos;
	// 		int rseed=next_seed->read_pos;
	// 		int dseed=next_seed->dna_pos;
	// 		int prev_z= next_seed->best_score_pos[zz]->path_number_matrices;
			
	// 		while(!(rstart==rseed-1 && dstart==dseed-1)){
				
	// 			int matrix_position= (rstart-rseed)*(max_gaps2*2+1)+dstart-(dseed+rstart-rseed)+max_gaps2;
	// 			//	fprintf(stdout,"(%i-%i),(%i,%i) %i\n",rstart,dstart,rseed,dseed,matrix_position);
	// 			int prev_i=((Prev_score*)matrix[prev_z]+matrix_position)->prev_i;
	// 			int prev_j=((Prev_score*)matrix[prev_z]+matrix_position)->prev_j;
	// 			prev_z=((Prev_score*)matrix[prev_z]+matrix_position)->prev_matrix_no;
	// 			//fprintf(stdout,"prev %i-%i:%i-%i (%i)\n",rstart,dstart,prev_i,prev_j, prev_z);
				
	// 			if (prev_i==rstart && prev_j==dstart-1){//read gap
	// 				dna_align_temp.push_back(dna[dstart]);
	// 				read_align_temp.push_back('-');
	// 			}
	// 			else if(prev_i==rstart-1 && prev_j==dstart){//dna gap
	// 				dna_align_temp.push_back('-');
	// 				read_align_temp.push_back(read[rstart]);
	// 			}
	// 			else if(prev_i==rstart-1 && prev_j==dstart-1){//match/mismatch
	// 				dna_align_temp.push_back(dna[dstart]);
	// 				read_align_temp.push_back(read[rstart]);
	// 			}
				
				
	// 			rstart=prev_i;
	// 			dstart=prev_j;
	// 		}
			
			
	// 		rstart=next_seed->best_score_pos[zz]->read_pos;
	// 		dstart=next_seed->best_score_pos[zz]->dna_pos;
			
	// 		int tmp=next_seed->best_score_pos[zz]->path_number;
	// 		next_seed=next_seed->best_score_pos[zz]->next_seed;
	// 		zz=tmp;
			
			
	// 		for ( rit= dna_align_temp.rbegin() ; rit < dna_align_temp.rend(); ++rit )
	// 			dna_align.push_back(*rit);
	// 		for ( rit= read_align_temp.rbegin() ; rit < read_align_temp.rend(); ++rit )
	// 			read_align.push_back(*rit);
			
	// 		dna_align_temp.clear();
	// 		read_align_temp.clear();
			
	// 		if (next_seed!=NULL){
	// 			if (next_seed->deletion_id.empty()){
	// 				for(int n=dstart+1;n<=next_seed->dna_pos-1;n++){
	// 					dna_align.push_back(dna[n]);
	// 					read_align.push_back('*');
	// 				} 
	// 			}
	// 			else{
	// 				for(int n=dstart+1;n<=next_seed->dna_pos-1;n++){
	// 					dna_align.push_back(dna[n]);
	// 					read_align.push_back('d');
	// 				} 
	// 			}
				
				
	// 		}
	//     }


	//     std::vector<char>::iterator it;
	//     fprintf(stdout,"DNA  ALIGN: ");
	//     for ( it= dna_align.begin() ; it < dna_align.end(); it++ )
	// 		fprintf(stdout,"%c",*it);
	//     fprintf(stdout,"\nREAD ALIGN: ");
	//     for ( it= read_align.begin() ; it < read_align.end(); it++ )
	// 		fprintf(stdout,"%c",*it);
	//     fprintf(stdout,"\n");
	//     dna_align.clear();
	//     read_align.clear();
	//   }
	// }
  
	// SeedElem* next_seed= seed_matrix_left[0];
	// while(next_seed!=NULL){
	// 	fprintf(stdout,"seed %i-%i\n",next_seed->read_pos,next_seed->dna_pos);
	// 	next_seed=next_seed->best_score_pos[0]->next_seed;
	// }
	

	/*********************************************************************************************/
	/* Clean structures */
	/*********************************************************************************************/
	//  fprintf(stdout,"Clean structures...\n");
 
	delete[] best_match_scores;
	best_match_scores=NULL;

	
	//  fprintf(stdout,"Clean structures...END\n");
}



