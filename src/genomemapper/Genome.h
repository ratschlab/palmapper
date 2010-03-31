#pragma once

#include <stdio.h>
#include <string>
#include <vector>

#define CHR_DESC_LENGTH 50

#define USE_CHR_BIN
#define USE_CHR_BIN_CLASS CDNAArray4

class Chromosome {
	friend class Genome;
public:
	Chromosome();

	char operator[](unsigned int index) const {
//		assert(index < _length);
#ifdef USE_CHR_BIN
		return CHR_SEQ_dd->get_char(index);
#else
		return _data[index];
#endif
	}

	unsigned int length() const {
		return _length;
	}

	unsigned int nr() const {
		return _nr;
	}

	char const *desc() const {
		return _desc.c_str();
	}
private:
#ifdef USE_CHR_BIN
	USE_CHR_BIN_CLASS* CHR_SEQ_dd;
#endif
	std::string _desc;
	unsigned int _nr;
	unsigned int _length;
	char *_data;
};

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

//#ifdef USE_CHR_BIN
//extern std::vector<USE_CHR_BIN_CLASS*> CHR_SEQ_a ;
//#endif
//extern char** CHR_SEQ_c ;

//inline char CHR_SEQ(size_t chr, size_t index)
//{
//#ifdef USE_CHR_BIN
//	return CHR_SEQ_a[(chr)]->get_char(index) ;
//#else
//	return CHR_SEQ_c[(chr)][(index)] ;
//#endif
//}

//extern unsigned int* CHR_LENGTH;
//extern char** CHR_DESC;

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
