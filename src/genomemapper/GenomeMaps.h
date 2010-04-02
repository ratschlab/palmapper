#pragma once

// ##############################################################
// ####### GenomeMaps ###########################################
// ##############################################################

//#define CHR_MAP_DNAARRAY
//#define CHR_MAP_DNAARRAY_CLASS CDNAArray4
// define this if CDNAArray4 is chosen
//#define CHR_MAP_DNAARRAY_2BIT


#ifndef CHR_MAP_DNAARRAY
#define MASK_MAPPED_READ_BEST      1
#define MASK_MAPPED_READ           2
#define MASK_SPLICED_READ_BEST      4
#define MASK_SPLICED_READ           8
#else
#ifdef CHR_MAP_DNAARRAY_2BIT
#define MASK_MAPPED_READ_BEST      1
#define MASK_MAPPED_READ           4
#define MASK_SPLICED_READ_BEST      2
#define MASK_SPLICED_READ_BEST_ORIG      4
#define MASK_SPLICED_READ           8
#else
#define MASK_MAPPED_READ_BEST      1
#define MASK_MAPPED_READ           2
#define MASK_SPLICED_READ_BEST      4
#define MASK_SPLICED_READ           8
#endif
#endif
#define MASK_MAPPED_REGION         16
#define MASK_REPETITIVE_SEED       32
#define MASK_REPETITIVE_SEED_MANY1  64
#define MASK_REPETITIVE_SEED_MANY2  128


class GenomeMaps
{
public:
	GenomeMaps() ;
	~GenomeMaps() ;

	inline unsigned char CHR_MAP(Chromosome const &chr, size_t index)
	{
#ifdef CHR_MAP_DNAARRAY
		return CHR_MAP_a[chr.nr()]->get_elem(index) ;
#else
		assert(CHR_MAP_c!=NULL) ;
		
		/*if (CHR_MAP_c[chr][index] != CHR_MAP_a[chr]->get_elem(index))
		  {
		  fprintf(stdout, "get: chr=%i, index=%i: %i != %i\n", (int)chr, (int)index,  (int)CHR_MAP_c[chr][index], (int)CHR_MAP_a[chr]->get_elem(index)) ;
		  exit(-1) ;
		  }*/
		
		return CHR_MAP_c[chr.nr()][index] ;
#endif
	}
	
	inline void CHR_MAP_set(Chromosome const &chr, size_t index, unsigned char c)
	{
#ifdef CHR_MAP_DNAARRAY
#ifdef CHR_MAP_DNAARRAY_2BIT
		assert(c<4) ;
#else // CHR_MAP_DNAARRAY_2BIT
		assert(c<16) ;
#endif // CHR_MAP_DNAARRAY_2BIT
		CHR_MAP_a[chr.nr()]->set_elem(index, c) ;
#else // CHR_MAP_DNAARRAY
		CHR_MAP_c[chr.nr()][index]=c ;
		//CHR_MAP_a[chr]->set_elem(index, c) ;
		
		/*if (CHR_MAP_c[chr][index] != CHR_MAP_a[chr]->get_elem(index))
		  {
		  fprintf(stdout, "set: chr=%i, index=%i: %i != %i\n", (int)chr, (int)index,  (int)CHR_MAP_c[chr][index], (int)CHR_MAP_a[chr]->get_elem(index)) ;
		  exit(-1) ;
		  }*/
#endif // CHR_MAP_DNAARRAY
	}
	

	int init_reporting() ;
	int report_repetitive_seed(Chromosome const &chr, int chr_start, int count) ;
	int report_mapped_region(Chromosome const &chr, int chr_start, int chr_end, int num_matches)  ;
	int report_mapped_read(Chromosome const &chr, int start, int end, int num_matches, int nbest_hit) ;
	int report_spliced_read(Chromosome const &chr, std::vector<int> & exons, int num_matches, int nbest_hit) ;
    //extern int report_splice_site(int chr, int pos, char strand, char type) ;
	int do_reporting(int force=0) ;
	int read_reporting() ;
	int write_reporting() ;
	int clean_reporting() ;
	
protected:
	unsigned char **CHR_MAP_c ;
#ifdef CHR_MAP_DNAARRAY
	std::vector<CHR_MAP_DNAARRAY_CLASS*> CHR_MAP_a  ;
#endif
	void to_dnaarray(int chr=-1) ;
	void from_dnaarray(int chr = -1) ;

	int reported_repetitive_seeds  ;
	int reported_mapped_regions  ;
	int reported_mapped_reads  ;
	int reported_spliced_reads  ;
    //int reported_splice_sites  ;
	
	int covered_mapped_read_positions  ;
	int covered_mapped_read_positions_best  ;
	int covered_spliced_read_positions  ;
	int covered_spliced_read_positions_best  ;
	int covered_repetitive_seed_positions  ;
	int covered_repetitive_seed_positions_many1  ;
	int covered_repetitive_seed_positions_many2  ;
	int covered_mapped_region_positions  ;
    //int covered_splice_site_positions  ;

	static clock_t last_report ;

public:
	int REPORT_REPETITIVE_SEED_DEPTH_EXTRA ;

} ;

