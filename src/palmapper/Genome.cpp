// Authors: Korbinian Schneeberger and Joerg Hagmann 
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>


#include <palmapper/Genome.h>
#include <palmapper/Util.h>

#include <palmapper/palmapper.h>

#ifdef PMINDEX
#undef VERSION
#include <pmindex/pmindex.h>

#endif

#ifndef PMINDEX
#include <bwa/bwtaln.h>
#include <bwa/bwtmyaln.h>
#endif


inline char get_compl_base_(char c)
{
	switch (c)
	{
	case 'A':
		return 'T';
	case 'C':
		return 'G';
	case 'G':
		return 'C';
	case 'T':
		return 'A';
	default:
		return c;
	}
	return 0;
}

bool Genome::classInitialized = initClass();
bool Genome::initClass() {
	for (int cc=0; cc<255; cc++)
	{
		upper_char[cc] = (cc >= 'a' && cc <= 'z') ? cc - 32 : cc;
		valid_char[cc] = ::strchr("ACGTNRYMKWSBDHV", upper_char[cc]) != NULL ? upper_char[cc] : 0;
		compl_char[cc] = get_compl_base_(cc);
	}
	return true;
}

int Genome::valid_char[256];
char Genome::compl_char[256];
char Genome::upper_char[256];

Genome::Genome(int numChromosomes) {
	MAX_POSITIONS = 0;
	
	NUM_CHROMOSOMES = numChromosomes ;

#ifndef PMINDEX 
	INDEX_SIZE = Config::INDEX_SIZE_16 ;
#endif

	INDEX=NULL;
#ifndef PMINDEX 
	BLOCK_TABLE=NULL;
#endif

#ifndef PMINDEX 
	assert(numChromosomes==0) ;
 	if (_config.VERBOSE) { printf("Reading in indices\n"); }
	build_index();

	if (_config.VERBOSE) printf("Reading in genome\n");
#endif
#ifdef PMINDEX
	if (NUM_CHROMOSOMES>0)
		_chromosomes = new Chromosome[NUM_CHROMOSOMES];
	else
		_chromosomes = NULL ;
	
#endif

#ifndef PMINDEX 
	load_genome();
	init_constants(); // updated
#endif
}

Genome::~Genome() {
	delete[] _chromosomes;
	free(INDEX);
	free(BLOCK_TABLE);

#ifndef PMINDEX 
	if (_config.BWA_INDEX == 1)
		bwa_seed2genome_destroy() ;
#endif
}

int Genome::alloc_index_memory()
{

        //if ( (MEM_MGR.first_entry = (STORAGE_ENTRY *) malloc (NUM_POS * (1+_config.MAP_REVERSE) * sizeof(STORAGE_ENTRY)) ) == NULL) {
        //        fprintf(stderr, "ERROR : not enough memory for mem_master (1)\n");
        //        exit(1);
        //}

	//MEM_MGR.num_bins = 0;
	//MEM_MGR.next_unused_entry = MEM_MGR.first_entry;

#ifndef PMINDEX
	INDEX_SIZE=_config.INDEX_SIZE_16 ;
	if (_config.INDEX_DEPTH==15)
		INDEX_SIZE = _config.INDEX_SIZE_15 ;
	if (_config.INDEX_DEPTH==14)
		INDEX_SIZE = _config.INDEX_SIZE_14 ;
	if (_config.INDEX_DEPTH==13)
		INDEX_SIZE = _config.INDEX_SIZE_13 ;
	if (_config.INDEX_DEPTH==12)
		INDEX_SIZE = _config.INDEX_SIZE_12 ;
	if (_config.VERBOSE>0)
		fprintf(stdout, "_config.INDEX_DEPTH=%i, INDEX_SIZE=%ld, sizeof(INDEX_ENTRY)=%ld, index size=%ld\n",  (int)_config.INDEX_DEPTH, (long int)INDEX_SIZE, (long int)sizeof(INDEX_ENTRY), (long int)sizeof(INDEX_ENTRY)*INDEX_SIZE) ;
#endif

	if ( (INDEX = (INDEX_ENTRY *) calloc (INDEX_SIZE, sizeof(INDEX_ENTRY)) ) == NULL) {
		fprintf(stderr, "ERROR : not enough memory for mem_master (2)\n");
		exit(1);
	}

	/*if ( _config.MAP_REVERSE && ((INDEX_REV = (INDEX_ENTRY *) calloc (INDEX_SIZE, sizeof(INDEX_ENTRY))) == NULL) ) {
		fprintf(stderr, "ERROR : not enough memory for mem_master (3)\n");
		exit(1);
	}*/

	return(0);
}

#ifdef PMINDEX
int Genome::load_genome()
{	
	FILE *genome_fp = fopen(GENOME_FILE_NAME, "r");
	if (genome_fp == NULL) {
		fprintf(stderr, "ERROR : Couldn't open genome file %s\n",
				GENOME_FILE_NAME);
		exit(1);
	}

	char line[513];
	unsigned int fp = ftell(genome_fp);
	unsigned int linelen;

	//std::vector<Chromosome> chromosomes ;
	_chromosomes = new Chromosome[10000] ;//NUM_CHROMOSOMES];

#ifdef USE_CHR_BIN
	fprintf(stdout, "using compact genome representation (class %s)\n", USE_CHR_BIN_CLASS::get_class_name()) ;
#endif	
	unsigned int i ;
	for (i=0; ; ++i) { 
		fprintf(stdout, ".") ;
		//chromosomes.push_back(Chromosome()) ; 
		Chromosome & chr = _chromosomes[i];
		chr._nr = i;
		std::string data ;
		data.reserve(200000000) ;
		char desc[100] ;//, rest[1000] ;
		
		//unsigned int pos = 0;
		
		line[0] = '\0';
		fseek(genome_fp, fp, SEEK_SET);
		
		while (line[0] != '>')
			if (fgets(line, 512, genome_fp) == 0) 
			{
				break;
			}
		if (line[0]=='>')
		{
			int ret = sscanf(line, ">%100s", desc) ;
			assert(ret>=1) ;
			//fprintf(stdout, "found chromosome %s\n", desc) ;
		}
		fflush(stdout) ; 
	
		if (fgets(line, 512, genome_fp) == NULL || line[0] == '>') 
		{
			break ;
			//fprintf(stderr, "ERROR: cannot find sequence \"%s\"!\n",chr.desc());
			//exit(1);
		}
		while (line[0] != '>') 
		{
			linelen = strcspn(line, " \n\t");
			if (linelen > 0 && (line[linelen] == '\t' || line[linelen] == ' ')) {
				fprintf(stderr, "ERROR: white space character unequal to newline found in genome input file '%s' in chromosome '%s'!\n", GENOME_FILE_NAME, chr.desc());
				exit(0);
			} 
			char buf[linelen+1] ;
			for (unsigned int j=0; j!=linelen; j++)
			{
				if (is_valid_char((int)line[j]))
				{
					buf[j]= is_valid_char((int)line[j]) ;
				}
				else 
				{
					char c=mytoupper(line[j]) ;
					fprintf(stderr,"ERROR: Character '%c' encountered in chromosome '%s'! Only IUPAC-code is accepted!\n", c, chr.desc());
					exit(0);
				} 
			}
			buf[linelen]=0 ;
			//if (data.size()==0)
			//fprintf(stdout, "%s\n", buf) ; 
			data+=buf ;
			
			fp = ftell(genome_fp);
			if (fgets(line, 512, genome_fp) == NULL) break;
		}
		
		//data[chr.length()] = '\0';
		chr._length= data.size() ;
		
		if (chr.length() != data.size()) {
			fprintf(stderr, "ERROR: Idx file seems to be corrupted. Chromosome %d has %d characters at the end! (%d %d)\n",i+1, (int)data.size()-chr.length(), (int)data.size(), chr.length());
			exit(1);
		}
		char *mbuf=new char[data.size()+1] ;
		strcpy(mbuf, data.c_str()) ;
		
		chr.data(mbuf, chr.length(), desc);
	}
	NUM_CHROMOSOMES=i ;
	fprintf(stdout, "read %i chromosomes\n", NUM_CHROMOSOMES) ;

	//delete[] _chromosomes ;
	//_chromosomes = new Chromosome[NUM_CHROMOSOMES];
	//for (unsigned int ii=0; ii<NUM_CHROMOSOMES; ii++)
	//	_chromosomes[ii]=chromosomes[ii] ;
	
	fclose(genome_fp);
	
	return 0;	
}
#endif

#ifndef PMINDEX
int Genome::load_genome()
{	

	if (_config.NO_LOAD_GENOME)
	{
		fprintf(stdout, "Not loading genome sequence...\n") ;
		return 0 ;
	}
	
	FILE *genome_fp = fopen(_config.GENOME_FILE_NAME.c_str(), "r");
	if (genome_fp == NULL) {
		fprintf(stderr, "ERROR : Couldn't open genome file %s\n",
				_config.GENOME_FILE_NAME.c_str());
		exit(1);
	}

//	if ((CHR_SEQ_c = (char**) malloc (NUM_CHROMOSOMES * sizeof(char**))) == NULL) {
//		fprintf(stderr, "ERROR : not enough memory for genome\n");
//		exit(1);
//	}

	char line[513];
	unsigned int fp = ftell(genome_fp);
	unsigned int linelen;

#ifdef USE_CHR_BIN
		fprintf(stdout, "using compact genome representation (class %s)\n", USE_CHR_BIN_CLASS::get_class_name()) ;
#endif	

	for (unsigned int i=0; i!=NUM_CHROMOSOMES; ++i) 
	{
		Chromosome &chr = _chromosomes[i];
		chr._nr = i;
		char *data = new char[chr.length() + 1];
//		if ((chr._data = (char*) malloc ((chr._length + 1) * sizeof(char))) == NULL) {
//			fprintf(stderr, "ERROR : not enough memory for genome\n");
//			exit(1);
//		}
		
		unsigned int pos = 0;
		
		line[0] = '\0';
		fseek(genome_fp, fp, SEEK_SET);
		
		while (line[0] != '>')
			if (fgets(line, 512, genome_fp) == 0) {}
		
		if (fgets(line, 512, genome_fp) == NULL || line[0] == '>') {
			fprintf(stderr, "ERROR: cannot find sequence \"%s\"!\n",chr.desc());
			exit(1);
		}
		while (line[0] != '>') 
		{
			linelen = strcspn(line, " \n\t");
			if (linelen > 0 && (line[linelen] == '\t' || line[linelen] == ' ')) {
				fprintf(stderr, "ERROR: white space character unequal to newline found in genome input file '%s' in chromosome '%s'!\n", _config.GENOME_FILE_NAME.c_str(), chr.desc());
				exit(0);
			}
			for (unsigned int j=0; j!=linelen; j++)
			{
				/*if (c=='A' || c=='C' ||
				    c=='G' || c=='T' ||
				    c=='N' || c=='R' ||
				    c=='Y' || c=='M' ||
				    c=='K' || c=='W' ||
				    c=='S' || c=='B' ||
				    c=='D' || c=='H' ||
				    c=='V' ) */
				if (is_valid_char((int)line[j]))
				{
					data[pos++] = is_valid_char((int)line[j]) ;
				}
				else 
				{
					char c=mytoupper(line[j]) ;
					fprintf(stderr,"ERROR: Character '%c' encountered in chromosome '%s'! Only IUPAC-code is accepted!\n", c, chr.desc());
					exit(0);
				} 
			}
			
			fp = ftell(genome_fp);
			if (fgets(line, 512, genome_fp) == NULL) break;
		}
		
		data[chr.length()] = '\0';

		if (chr.length() != strlen(data)) {
			fprintf(stderr, "ERROR: Idx file seems to be corrupted. Chromosome %d has %d characters at the end! (%d %d)\n",i+1, (int)strlen(data)-chr.length(), (int)strlen(data), chr.length());
			exit(1);
		}

		chr.data(data, chr.length());
	}
	
	fclose(genome_fp);
	
	return 0;	
}
#endif

#ifndef PMINDEX
int Genome::build_index()
{

	if (_config.BWA_INDEX==0)
	{	
		FILE *META_INDEX_FP = Util::openFile(_config.META_INDEX_FILE_NAME, "r");
		FILE *CHR_INDEX_FP = Util::openFile(_config.CHR_INDEX_FILE_NAME, "r");
		// handle meta information
		read_meta_index_header(META_INDEX_FP);  // updated
		alloc_index_memory(); // updated		
		if (!_config.NO_LOAD_GENOME)
			read_meta_index(META_INDEX_FP); // updated
		
		// mmap map files into memory
		if (!_config.NO_LOAD_GENOME)
			mmap_indices(); // updated
		
		// handle index information
		read_chr_index(CHR_INDEX_FP);

	}
	else
	{

		if (_config.VERBOSE) { printf("\tIndex depth is %d\n", _config.INDEX_DEPTH); }

		if (_config.HITLEN_LIMIT == 0) _config.HITLEN_LIMIT = _config.INDEX_DEPTH;
		else if (_config.HITLEN_LIMIT < _config.INDEX_DEPTH) {
			fprintf(stderr, "\n!!! WARNING: Hitlength limit is smaller than seedlength, it will be set to seedlength!\n\n");
			_config.HITLEN_LIMIT = _config.INDEX_DEPTH;
		}

		fprintf(stdout, "bwa_seed2genome_init ... ") ;
		bwa_seed2genome_init(_config.GENOME_FILE_NAME.c_str(), NULL) ;
		read_chr_bwa();
		fprintf(stdout, "done.\n") ;
	}
 
	return(0);
}

int Genome::read_meta_index_header(FILE *META_INDEX_FP)
{

	if (_config.VERBOSE) { printf("Reading in meta index\n"); }

	////////////////////////////////////////////////////////////////////////////////////
	// Get rev index flag (only for downcompatibility to previous GM versions)
	char dummy;
	if (fread(&dummy, sizeof(char), 1, META_INDEX_FP) == 0) {
                fprintf(stderr, "ERROR: cant read meta index file\n");
                exit(1);
        }

	// Index depth
	if (fread(&_config.INDEX_DEPTH, sizeof(int), 1, META_INDEX_FP) == 0) {
                fprintf(stderr, "ERROR: cant read meta index file\n");
                exit(0);
        }
	assert(_config.INDEX_DEPTH<=MAX_INDEX_DEPTH) ;

	if (_config.VERBOSE) { printf("\tIndex depth is %d\n", _config.INDEX_DEPTH); }

	if (_config.HITLEN_LIMIT == 0) _config.HITLEN_LIMIT = _config.INDEX_DEPTH;
	else if (_config.HITLEN_LIMIT < _config.INDEX_DEPTH) {
		fprintf(stderr, "\n!!! WARNING: Hitlength limit is smaller than seedlength, it will be set to seedlength!\n\n");
		_config.HITLEN_LIMIT = _config.INDEX_DEPTH;
	}

	////////////////////////////////////////////////////////////////////////////////////
	// Get number of chromosomes
	if (fread(&NUM_CHROMOSOMES, sizeof(int), 1, META_INDEX_FP) == 0) {
                fprintf(stderr, "ERROR: cant read meta index file\n");
                exit(0);
        }
	if (_config.VERBOSE) { printf("\tNb of chromosomes is %d\n", NUM_CHROMOSOMES); }

	_chromosomes = new Chromosome[NUM_CHROMOSOMES];

//	// alloc space for chomosome size
//	if ((CHR_LENGTH = (unsigned int *) malloc (NUM_CHROMOSOMES * sizeof(unsigned int*))) == NULL)
//	{
//		fprintf(stderr, "[read_meta_index_header] ERROR Could not allocate memory for genome memory\n");
//		exit(1);
//	}
//	// and descriptions
//	if ((CHR_DESC = (char**) malloc (NUM_CHROMOSOMES * sizeof(char**))) == NULL) {
//		fprintf(stderr, "ERROR : not enough memory for genome description\n");
//		exit(1);
//	}

	// initialize hit region mapping


	////////////////////////////////////////////////////////////////////////////////////
	// Get number of positions in index?
	int NUM_POS ;
	if (fread(&NUM_POS, sizeof(int), 1, META_INDEX_FP) == 0) {
                fprintf(stderr, "ERROR: cant read meta index file\n");
                exit(0);
        }

	////////////////////////////////////////////////////////////////////////////////////
	// Size of longest chromosome
	if (fread(&LONGEST_CHROMOSOME, sizeof(int), 1, META_INDEX_FP) == 0) {
                fprintf(stderr, "ERROR: cant read meta index file\n");
                exit(0);
        }

	//////////////////////////////////////////////////////////////
	// read block table:
	//////////////////////////////////////////////////////////////
	unsigned int blocks;
	if (fread(&blocks, sizeof(unsigned int), 1, META_INDEX_FP) == 0) {
                fprintf(stderr, "ERROR: cant read meta index file\n");
                exit(0);
        }

	//TODO: dd introduce an overloaded operator new?
	BLOCK_TABLE = (POS *) malloc (blocks * sizeof(POS));
	if (BLOCK_TABLE == NULL)
	{
		fprintf(stderr, "[read_meta_index_header] Could not allocate memory\n");
		exit(1);
	}

	 if (fread(BLOCK_TABLE, sizeof(POS), blocks, META_INDEX_FP) == 0)
	 {
                fprintf(stderr, "ERROR: cant read meta index file (blocktable)\n");
                exit(0);
        }

	// print block table (debugging)
/*	printf("-------------------------------------------\n");
	int i;
	for (i=0; i!=blocks; ++i) {
		printf("| block %5d | pos %9d | chr %5d |\n", i, BLOCK_TABLE[i].pos, BLOCK_TABLE[i].chr+1);
	}
	printf("-------------------------------------------\n");*/

	return(0);
}

int Genome::read_meta_index(FILE *META_INDEX_FP)
{
 	META_INDEX_ENTRY file_entry;
	int used_slots = 0;
	unsigned int index_offset = 0;//, index_rev_offset = 0;

	while (fread(&file_entry, sizeof(file_entry), 1, META_INDEX_FP) == 1)
	{
		if (file_entry.slot < 0) continue;	//downcompatibility to previous GM versions

		used_slots++;

		index_offset += file_entry.num;

		INDEX[file_entry.slot].num = file_entry.num;
		INDEX[file_entry.slot].offset = index_offset ;

		if (file_entry.num > MAX_POSITIONS)
			MAX_POSITIONS = file_entry.num;
	}

	fclose(META_INDEX_FP);

  	return 0;
}


int Genome::read_chr_bwa()
{
	unsigned int chr_num = 0;
	int chrlen;
	int numchr;
	char *chr_desc=NULL;
	unsigned int max_len=0;
	
	if (_config.VERBOSE || true) { printf("Reading in index\n"); }

	bwa_seed2genome_numchr(&numchr);
	NUM_CHROMOSOMES=numchr;
	
	fprintf(stdout, "\tnum chromosome is %i, ", NUM_CHROMOSOMES); 
	_chromosomes = new Chromosome[NUM_CHROMOSOMES];

	while (chr_num != NUM_CHROMOSOMES) {

		//chromosome
		if (_config.VERBOSE) { printf("\tchromosome ID is %d, ", chr_num+1); }

		//chromosome length
		bwa_seed2genome_lenchr(chr_num,&chrlen);
		_chromosomes[chr_num]._length = chrlen;
		if ((int)chrlen> (int)max_len)
			max_len=chrlen;
		
		//chromosome description
		bwa_seed2genome_descchr(chr_num,chr_desc);
		_chromosomes[chr_num].desc(chr_desc);
		if (_config.VERBOSE) { printf("description is %s\n", _chromosomes[chr_num].desc()); }

		chr_num++;
	} // for every chromosome

	if (_config.VERBOSE) { printf("Finished parsing index\n"); }
	
	LONGEST_CHROMOSOME=max_len;
	
  	return(0);
}

int Genome::read_chr_index(FILE *CHR_INDEX_FP)
{
	unsigned int chr_num = 0;
	unsigned int chr;
	unsigned int chrlen;
	char chr_desc[CHR_DESC_LENGTH];

	if (_config.VERBOSE) { printf("Reading in index\n"); }

	while (chr_num != NUM_CHROMOSOMES) {

		chr_num++;

		//HEADER OF CHROMOSOME ENTRY

		//chromosome
		if (fread(&chr, sizeof(unsigned int), 1, CHR_INDEX_FP) != 1) {
			printf("Early stop in index file (1).\nCorrupted file?\n");
			exit(1);
		}
		if ( _config.VERBOSE) { printf("\tchromosome ID is %d, ", chr+1); }

		//chromosome length
		if (fread(&chrlen, sizeof(unsigned int), 1, CHR_INDEX_FP) != 1) {
			printf("Early stop in index file (2).\nCorrupted file?\n");
			exit(1);
		}
		_chromosomes[chr]._length = chrlen;

		//chromosome description
		if (fread(&chr_desc, sizeof(char), CHR_DESC_LENGTH, CHR_INDEX_FP) != CHR_DESC_LENGTH) {
			printf("Early stop in index file (3).\nCorrupted file?\n");
			exit(1);
		}
		_chromosomes[chr].desc(chr_desc);
		if (_config.VERBOSE) { printf("description is %s\n", _chromosomes[chr].desc()); }

	} // for every chromosome

	fclose(CHR_INDEX_FP);

	if (_config.VERBOSE) { printf("Finished parsing index\n"); }

  	return(0);
}

////////////////////////////////////////////////////////
// Code originally written by Andre Noll.
// Code adapted from Paraslash.
////////////////////////////////////////////////////////

#ifndef BinaryStream_MAP

void Genome::index_pre_buffer(STORAGE_ENTRY* index_mmap, STORAGE_ENTRY* buffer, long index, long size) const
{
	assert(index>=0);
	assert(size>=0);

	for (long i=0; i<size; i++)
		buffer[i]=index_mmap[index+i] ;
}

void Genome::mmap_indices()
{
	void *fwd;
	size_t fwd_size=0 ;

	// Map fwd index
	if (mmap_full_file(_config.INDEX_FWD_FILE_NAME.c_str(), &fwd, &fwd_size) != 0) {
		printf("ERROR: Could not get file status\n");
                exit(1);
	}

	INDEX_FWD_MMAP = (STORAGE_ENTRY *)fwd;

    if (_config.INDEX_PRECACHE)
    {
        STORAGE_ENTRY buffer[1024] ;

		clock_t last_time = clock() ;
        for (size_t i=0; i<(fwd_size/sizeof(STORAGE_ENTRY))-1024; i+=1024)
		{
            index_pre_buffer(INDEX_FWD_MMAP, buffer, i, 1024) ;
			
			if (i%1000==0 && ((1.0*clock()-last_time)/CLOCKS_PER_SEC>0.01)) // the timing is tricky. Clock only measures user time, but this statement mostly uses system or io time
			{
				fprintf(stdout, "Linearly reading index files to fill caches: %3.1f%%\r", (100.0*i)/(fwd_size/sizeof(STORAGE_ENTRY))) ;
				last_time = clock() ;
			}
		}
		fprintf(stdout, "Linearly reading index files to fill caches: %3.1f%%\r", 100.0) ;
        fprintf(stdout, "\n") ;
		}
	return;
}

int Genome::mmap_full_file(const char *path, void **map, size_t * size_p)
{

	int fd, ret, mmap_prot, mmap_flags, open_mode;
	struct stat file_status;
	size_t size;

	// Set modus
	open_mode = O_RDONLY;
	mmap_prot = PROT_READ;
	mmap_flags = MAP_PRIVATE ;
#ifndef __APPLE__
	if (_config.INDEX_PRECACHE)
		mmap_flags |= MAP_POPULATE ;
#endif

	// Open file to get file size
	ret = open(path, open_mode, 0);
	if (ret < 0) {
		printf("ERROR: Could not open file\n");
                exit(1);
	}
	fd = ret;
        if (fstat(fd, &file_status) < 0) {
		printf("ERROR: Could not get file status\n");
		exit(1);
	}

	// Map file
	size = file_status.st_size;
	*size_p=size ;
	ret = gm_mmap(size, mmap_prot, mmap_flags, fd, 0, map, path);

	*size_p = size ;

	return ret;

}

int Genome::gm_mmap(size_t length, int prot, int flags, int fd, off_t offset, void *map, const char* path) {

	void **m = (void**)map;

	*m = mmap(NULL, length, prot, flags, fd, offset);

	if (*m == MAP_FAILED) {
	  printf("ERROR: Could not map file: %s\n", path);
	  perror("mmap");
                exit(1);
	}

	return 0;
}

#else

void Genome::mmap_indices()
{

	// Map fwd index
	INDEX_FWD_MMAP = new CBinaryStream<STORAGE_ENTRY>(_config.INDEX_FWD_FILE_NAME);

	if (_config.INDEX_PRECACHE)
	{
	//	fprintf(stdout, "Linearly reading index files to fill caches: 1/2 ") ;

		INDEX_FWD_MMAP->read_and_forget() ;

	//	fprintf(stdout, "Done\n") ;
	}

//		printf("ERROR: Could not get file status\n");
//                exit(1);
	/*for (int i=0; i<1000; i++)
	{
		STORAGE_ENTRY se=(*INDEX_FWD_MMAP)[i*10];
		fprintf(stderr, "se[%d]=%0x%0x%0x%0x\n",i, se.id[0], se.id[1], se.id[2], se.id[3]);
	}
	exit(1);*/
}

#endif // BinaryStream_MAP



int Genome::init_constants() 
{

	if (_config.BWA_INDEX==0){	
		if (_config.INDEX_DEPTH == 5) {
			BINARY_CODE[0] = 0; //binary number: 0000 0000 0000 0000 0000 0000
			BINARY_CODE[1] = 256; //binary number: 0000 0000 0000 0001 0000 0000
			BINARY_CODE[2] = 512; //binary number: 0000 0000 0000 0010 0000 0000
			BINARY_CODE[3] = 768; //binary number: 0000 0000 0000 0011 0000 0000
			BINARY_CODE[4] = 1023;     //binary number: 0000 0000 0000 0000 0000 0011 1111 1111
		}
		if (_config.INDEX_DEPTH == 6) {
			BINARY_CODE[0] = 0; //binary number: 0000 0000 0000 0000 0000 0000
			BINARY_CODE[1] = 1024; //binary number: 0000 0000 0000 0100 0000 0000
			BINARY_CODE[2] = 2048; //binary number: 0000 0000 0000 1000 0000 0000
			BINARY_CODE[3] = 3072; //binary number: 0000 0000 0000 1100 0000 0000
			BINARY_CODE[4] = 4095;    //binary number: 0000 0000 0000 0000 0000 1111 1111 1111
		}
		if (_config.INDEX_DEPTH == 7) {
			BINARY_CODE[0] = 0; //binary number: 0000 0000 0000 0000 0000 0000
			BINARY_CODE[1] = 4096; //binary number: 0000 0000 0001 0000 0000 0000
			BINARY_CODE[2] = 8192; //binary number: 0000 0000 0010 0000 0000 0000
			BINARY_CODE[3] = 12288; //binary number: 0000 0000 0011 0000 0000 0000
			BINARY_CODE[4] = 16383;    //binary number: 0000 0000 0000 0000 0011 1111 1111 1111
		}
		if (_config.INDEX_DEPTH == 8) {
			BINARY_CODE[0] = 0; //binary number: 0000 0000 0000 0000 0000 0000
			BINARY_CODE[1] = 16384; //binary number: 0000 0000 0100 0000 0000 0000
			BINARY_CODE[2] = 32768; //binary number: 0000 0000 1000 0000 0000 0000
			BINARY_CODE[3] = 49152; //binary number: 0000 0000 1100 0000 0000 0000
			BINARY_CODE[4] = 65535;    //binary number: 0000 0000 0000 0000 1111 1111 1111 1111
		}
		if (_config.INDEX_DEPTH == 9) {
			BINARY_CODE[0] = 0; //binary number: 0000 0000 0000 0000 0000 0000
			BINARY_CODE[1] = 65536; //binary number: 0000 0001 0000 0000 0000 0000
			BINARY_CODE[2] = 131072; //binary number: 0000 0010 0000 0000 0000 0000
			BINARY_CODE[3] = 196608; //binary number: 0000 0011 0000 0000 0000 0000
			BINARY_CODE[4] = 262143;    //binary number: 0000 0000 0000 0011 1111 1111 1111 1111
		}
		if (_config.INDEX_DEPTH == 10) {
			BINARY_CODE[0] = 0; //binary number: 0000 0000 0000 0000 0000 0000
			BINARY_CODE[1] = 262144; //binary number: 0000 0100 0000 0000 0000 0000
			BINARY_CODE[2] = 524288; //binary number: 0000 1000 0000 0000 0000 0000
			BINARY_CODE[3] = 786432; //binary number: 0000 1100 0000 0000 0000 0000
			BINARY_CODE[4] = 1048575;    //binary number: 0000 0000 0000 1111 1111 1111 1111 1111
		}
		if (_config.INDEX_DEPTH == 11) {
			BINARY_CODE[0] = 0; //binary number: 0000 0000 0000 0000 0000 0000
			BINARY_CODE[1] = 1048576; //binary number: 0001 0000 0000 0000 0000 0000
			BINARY_CODE[2] = 2097152; //binary number: 0010 0000 0000 0000 0000 0000
			BINARY_CODE[3] = 3145728; //binary number: 0011 0000 0000 0000 0000 0000
			BINARY_CODE[4] = 4194303;    //binary number: 0000 0000 0011 1111 1111 1111 1111 1111
		}
		if (_config.INDEX_DEPTH == 12) {
			BINARY_CODE[0] = 0; //binary number: 0000 0000 0000 0000 0000 0000
			BINARY_CODE[1] = 4194304; //binary number: 0100 0000 0000 0000 0000 0000
			BINARY_CODE[2] = 8388608; //binary number: 1000 0000 0000 0000 0000 0000
			BINARY_CODE[3] = 12582912; //binary number: 1100 0000 0000 0000 0000 0000
			BINARY_CODE[4] = 16777215;    //binary number: 0000 0000 1111 1111 1111 1111 1111 1111
		}
		if (_config.INDEX_DEPTH == 13) {
			BINARY_CODE[0] = 0; //binary number: 0000 0000 0000 0000 0000 0000 0000
			BINARY_CODE[1] = 16777216; //binary number: 0001 0000 0000 0000 0000 0000 0000
			BINARY_CODE[2] = 33554432; //binary number: 0010 0000 0000 0000 0000 0000 0000
			BINARY_CODE[3] = 50331648; //binary number: 0011 0000 0000 0000 0000 0000 0000
			BINARY_CODE[4] = 67108863;    //binary number: 0000 0011 1111 1111 1111 1111 1111 1111
		}
		if (_config.INDEX_DEPTH == 14) {
			BINARY_CODE[0] = 0; //binary number: 0000 0000 0000 0000 0000 0000 0000
			BINARY_CODE[1] = 67108864; //binary number: 0001 0000 0000 0000 0000 0000 0000
			BINARY_CODE[2] = 134217728; //binary number: 0010 0000 0000 0000 0000 0000 0000
			BINARY_CODE[3] = 201326592; //binary number: 0011 0000 0000 0000 0000 0000 0000
			BINARY_CODE[4] = 268435455;    //binary number: 0000 1111 1111 1111 1111 1111 1111 1111
		}
		if (_config.INDEX_DEPTH == 15) {
		  BINARY_CODE[0] = 0;             //binary number: 0000 0000 0000 0000 0000 0000 0000 0000
		  BINARY_CODE[1] = 268435456;     //binary number: 0001 0000 0000 0000 0000 0000 0000 0000
		  BINARY_CODE[2] = 536870912;     //binary number: 0010 0000 0000 0000 0000 0000 0000 0000
		  BINARY_CODE[3] = 805306368;     //binary number: 0011 0000 0000 0000 0000 0000 0000 0000
		  BINARY_CODE[4] = 1073741823; //binary number: 0011 1111 1111 1111 1111 1111 1111 1111
		}
		if (_config.INDEX_DEPTH == 16) {
		  BINARY_CODE[0] = 0;             //binary number: 0000 0000 0000 0000 0000 0000 0000 0000
		  BINARY_CODE[1] = 1073741824;     //binary number: 0100 0000 0000 0000 0000 0000 0000 0000
		  BINARY_CODE[2] = 2147483648;     //binary number: 1000 0000 0000 0000 0000 0000 0000 0000
		  BINARY_CODE[3] = 3221225472;     //binary number: 1100 0000 0000 0000 0000 0000 0000 0000
		  BINARY_CODE[4] = 4294967295; //binary number: 1111 1111 1111 1111 1111 1111 1111 1111
		}
		if (_config.INDEX_DEPTH>16 || _config.INDEX_DEPTH<5)
		{
			fprintf(stderr, "ERROR: _config.INDEX_DEPTH out of range\n") ;
			//exit(1) ;
		}
	}
	else{
		// Initialize binary codes to 0. Not used for bwt-based index
		BINARY_CODE[0] = 0;             
		BINARY_CODE[1] = 0;     
		BINARY_CODE[2] = 0;     
		BINARY_CODE[3] = 0;     
		BINARY_CODE[4] = 0;
		
	}
	
	return (0);
}
#endif
