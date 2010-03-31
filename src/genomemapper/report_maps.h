// ##############################################################
// ####### REPORT_MAPS ##########################################
// ##############################################################

//#define CHR_MAP_DNAARRAY
//#define CHR_MAP_DNAARRAY_CLASS CDNAArray4
// define this if CDNAArray4 is chosen
//#define CHR_MAP_DNAARRAY_2BIT

extern unsigned char **CHR_MAP_c ;
#ifdef CHR_MAP_DNAARRAY
extern std::vector<CHR_MAP_DNAARRAY_CLASS*> CHR_MAP_a  ;
#endif

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

extern int init_reporting() ;
extern int report_repetitive_seed(Chromosome const &chr, int chr_start, int count) ;
extern int report_mapped_region(Chromosome const &chr, int chr_start, int chr_end, int num_matches)  ;
extern int report_mapped_read(Chromosome const &chr, int start, int end, int num_matches, int nbest_hit) ;
extern int report_spliced_read(Chromosome const &chr, std::vector<int> & exons, int num_matches, int nbest_hit) ;
//extern int report_splice_site(int chr, int pos, char strand, char type) ;
extern int do_reporting(int force=0) ;
extern int read_reporting() ;
extern int write_reporting() ;
extern int clean_reporting() ;
