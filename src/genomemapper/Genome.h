#pragma once

#include <stdio.h>
#include <string>
#include <vector>

#include <genomemapper/Chromosome.h>


#define CHR_DESC_LENGTH 50

struct STORAGE_ENTRY{
	unsigned char id[4];
};

struct INDEX_ENTRY {
	unsigned int num;
	//STORAGE_ENTRY *last_entry;
	unsigned int offset;
};

typedef struct meta_idx_file_entry {
	int slot;
	unsigned int num;
} META_INDEX_ENTRY;

typedef struct position_structure {
	unsigned int pos;
	unsigned int chr;
} POS;

class Hits ;

class Genome {
public:
	Genome();
	Genome(bool dontCare) {hits=NULL ;}

	void set_hits(Hits* hits_)
	{
		hits=hits_ ;
	}
	
	Chromosome &chromosome(int index) {
		return _chromosomes[index];
	}
	unsigned int nrChromosomes() const {
		return NUM_CHROMOSOMES;
	}

	int BINARY_CODE[4];
	unsigned int LONGEST_CHROMOSOME;

#ifndef BinaryStream_MAP
	void index_pre_buffer(STORAGE_ENTRY* index_mmap, STORAGE_ENTRY* buffer, long index, long size) ;
#endif

private:
	int init_constants();
	int alloc_index_memory() ;
	int build_index();
	int load_genome();
	int read_meta_index_header(FILE *META_INDEX_FP);
	int read_meta_index(FILE *META_INDEX_FP);
	int read_chr_index(FILE *CHR_INDEX_FP);
	void mmap_indices();
	int gm_mmap(size_t length, int prot, int flags, int fd, off_t offset, void *map, const char* path);
	int mmap_full_file(const char *path, void **map, size_t * size_p);

	unsigned int NUM_CHROMOSOMES;
	Chromosome* _chromosomes;

	static bool initClass();
	static bool classInitialized;

	static int valid_char[256];
	static char compl_char[256];
	static char upper_char[256];

	friend char unique_base(char c);
	friend char mytoupper(char c); // considerably faster
	int is_valid_char(char cc) {
		return Genome::valid_char[(int)cc] ;
	}
	friend char get_compl_base(char c);

public:
	size_t INDEX_SIZE ;


	POS *BLOCK_TABLE;
	unsigned int BLOCK_TABLE_SIZE;
	
#if 1 // dd
	INDEX_ENTRY *INDEX;
	INDEX_ENTRY *INDEX_REV;
	
#ifndef BinaryStream_MAP
	STORAGE_ENTRY *INDEX_REV_MMAP;
	STORAGE_ENTRY *INDEX_FWD_MMAP;
#else
	CBinaryStream<STORAGE_ENTRY>* INDEX_REV_MMAP;
	CBinaryStream<STORAGE_ENTRY>* INDEX_FWD_MMAP;
#endif
#endif // dd

	unsigned long int MAX_POSITIONS;

	Hits* hits ;
};

inline char unique_base(char c)
{
	return (c == 'A' || c == 'C' || c == 'G' || c == 'T');
}

inline char mytoupper(char c) // considerably faster
{
	return Genome::upper_char[(int)c] ;
}

inline char get_compl_base(char c)
{
	return Genome::compl_char[(int)c] ;
}
