#pragma once

#include <genomemapper/Config.h>
#include <genomemapper/dyn_prog/qpalma_dp.h>

struct region_t {
	int32_t start;
	int32_t end;
	bool erased ;
	bool from_map ;
	char orientation ;
	bool* read_map ;
	//int32_t chromosome ;
	//char strand ;
};

class QPalma ;

struct perform_alignment_t
{
	std::string read_string ;
	std::string read_quality ;
	std::string dna ; 
	std::vector<region_t *> current_regions ;
	std::vector<int> positions ;
	Chromosome const *contig_idx ;
	char strand ;
	int ori ;
	int num_reported ;
	int ret ;
	pthread_t thread ;
	int hit_read;
	int hit_dna;
	int hit_length;
	QPalma * qpalma ;
	bool joined ;
} ;

struct alignment_parameter_struct {
	struct penalty_struct h, a, d, *qualityPlifs;
	int num_qualityPlifs;
	double *matchmatrix;
	int matchmatrix_dim[2];
	int quality_offset ;
} ;


const int num_filter_reasons=4  ;

class QPalma
{
    // initialization
	////////////////////

public:

	QPalma(Genome* genome_, Hits * hits, TopAlignments* topalignments_, GenomeMaps* genomemaps_, int verbosity_=2) ;
	~QPalma() ;

	int map_splice_sites(std::string file_template, char type, float &splice_site_threshold, bool estimate_thresh, bool do_report) ;
protected:
	int init_spliced_align(const char *fname, struct penalty_struct &h,
						   struct penalty_struct &a, struct penalty_struct &d,
						   struct penalty_struct *&qualityPlifs, int &num_qualityPlifs,
						   double*&matchmatrix, int dims[2], int &quality_offset) ;
	int check_splice_files(std::string file_template) ;
	void skip_comment_lines(FILE* fd) ;
	int read_matrix(FILE* fd, double *& matrix, char*& name, int dims[2]) ;
	int read_plif(FILE *fd, struct penalty_struct &plif) ;
	int clean_alignment_parameters() ;
	int init_alignment_parameters(std::string qpalma_file) ;
	int compare_double(const void *a, const void *b) ;
	
    // qpalma filtering
	////////////////////
	
public:
	int qpalma_filter(struct alignment_t *ali, int num_N) ;
	void qpalma_filter_stat(bool spliced) ;
	void qpalma_filter_stat_report() ;

protected:
	int get_num_splicesites(std::string file_template, const char* type, Chromosome const &chr, char strand, int start, int end, float thresh) ;
	
	
    // qpalma alignment
	////////////////////

public:
	int capture_hits();
	int perform_alignment(std::string &read_string, std::string &read_quality, std::string &dna, std::vector<region_t *> &regions, std::vector<int> &positions,
						  Chromosome const &contig_id, char strand, int ori, int & num_reported,int hit_read, int hit_dna, int hit_length) ;
	float score_unspliced(const char * read_anno) ;
	void capture_hits_timing(int read_count=-1, float this_read=-1.0) ;
	
protected:
	
	int get_splicesite_positions(std::string file_template, const char *type, Chromosome const &chr, char strand, int start, int end, float thresh, bool store_pos,
								 std::vector<int> &positions) ;
	
	
	int get_string_from_region(Chromosome const &chrN, region_t *region, std::string &str) ;
	void add_buffer_to_region(int ori, Chromosome const &chrN, int32_t nregion) ;
	void qsort(region_t** output, int size) ;
	void recover_long_regions(std::vector<region_t*> &long_regions_output, std::vector<region_t*> long_regions, std::vector<region_t*> current_regions) ;
	int convert_dna_position(int real_position, size_t* cum_length, const std::vector<region_t *> &current_regions) ;
	int get_first_read_map(bool* read_map) ;
	void print_hit(HIT *hit) ;
	void print_region(region_t *region, const char * bla)  ;
	void print_map(bool* read_map, const char *name) ;

	int perform_alignment_starter(std::string read_string, std::string read_quality, std::string dna, std::vector<region_t *> current_regions, std::vector<int> positions, Chromosome const &contig_idx, char strand, int ori, int hit_read_position, int hit_dna_position, int hit_length) ;
	int perform_alignment_wait(int & num_reported) ;

	void delete_regions() ;
	void delete_long_regions(std::vector<std::vector<region_t *> > *long_regions) ;
	int rescue_alignment(std::string & read_anno, int ori, int &num_A, int &num_T, int &num) ;


	// inline helpers
	////////////////////

	inline int ori_map(char c) 
	{
		if (c == '+')
			return 0;
		if (c == '-')
			return 1;
		fprintf(stdout, "ori: %c\n", c);
		
		assert(0);
	};

	inline std::string reverse(std::string str) 
	{
		for (int i = 0; i < (int)str.length() / 2; i++) 
		{
			char c = str[i];
			str[i] = str[str.length() - i - 1];
			str[str.length() - i - 1] = c;
		}
		
		return str;
	}
	
	inline std::vector<int> reverse(std::vector<int> vec) 
	{
		for (int i = 0; i < (int)vec.size() / 2; i++) 
		{
			int c = vec[i];
			vec[i] = vec[vec.size() - i - 1];
			vec[vec.size() - i - 1] = c;
		}
	
		return vec;
	}
	
	inline void reverse(double *vec, int len) 
	{
		for (int i = 0; i < (int)len / 2; i++) {
			double c = vec[i];
			vec[i] = vec[len - i - 1];
			vec[len - i - 1] = c;
		}
	}
	
	inline std::string complement(std::string str) {
		for (int i = 0; i < (int)str.length(); i++) {
			char c = str[i];
			switch (c) {
			case 'a':
				str[i] = 't';
				break;
			case 'c':
				str[i] = 'g';
				break;
			case 'g':
				str[i] = 'c';
				break;
			case 't':
				str[i] = 'a';
				break;
			case 'A':
				str[i] = 'T';
				break;
			case 'C':
				str[i] = 'G';
				break;
			case 'G':
				str[i] = 'C';
				break;
			case 'T':
				str[i] = 'A';
				break;
			case '[':
				str[i] = ']';
				break;
			case ']':
				str[i] = '[';
				break;
			case '-':
				str[i] = '-';
				break;
			default:
				if (c >= 'a' && c <= 'z')
					str[i] = 'n';
				else if (c >= 'A' && c <= 'Z')
					str[i] = 'N';
				else
					assert(0);
			}
		}
		
		return str;
	}

public:
	std::vector<std::vector<region_t *> > regions[2];

protected:
	
	struct alignment_parameter_struct *alignment_parameters;

	std::vector<perform_alignment_t*> thread_data ;

	clock_t region_align_time ;
	clock_t region1_time ;
	clock_t align_time ;
	int read_count ;
	
	long int total_dna_length ;
	long int total_alignments ;
	
	static clock_t last_timing_report ;
	
	const int verbosity ;
	const int MIN_NUM_MATCHES ;
	
	int total_num_threads ;
	int total_num_thread_tasks  ;


	static clock_t last_filter_report ;

	int qpalma_filter_reason  ;
	int qpalma_filter_stat_spliced[num_filter_reasons] ;
	int qpalma_filter_stat_unspliced[num_filter_reasons] ;

	Genome * genome ;
	Hits * hits ;
	TopAlignments* topalignments ;
	GenomeMaps* genomemaps ;
	
} ;
