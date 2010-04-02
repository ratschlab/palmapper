#pragma once

#include <stdio.h>
#include <string>
#include <vector>

#include <genomemapper/Chromosome.h>

#define CHR_DESC_LENGTH 50

class Genome {
public:
	Genome();
	Genome(bool dontCare) {}
	Chromosome &chromosome(int index) {
		return _chromosomes[index];
	}
	unsigned int nrChromosomes() const {
		return NUM_CHROMOSOMES;
	}


private:
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
