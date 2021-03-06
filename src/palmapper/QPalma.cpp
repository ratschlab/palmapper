// Authors: Geraldine Jean, Gunnar Raetsch, Lisa Thalheim, Dominik Diesch
// Copyright (C) 2009-2011 by Friedrich Miescher Laboratory, Tuebingen, Germany

#include <assert.h>
#include <string>
#include <sys/time.h>
#include <time.h>

#include "palmapper.h"
#include "IntervalQuery.h"
#include "dyn_prog/qpalma_dp.h"

#include <palmapper/Genome.h>
#include <palmapper/QPalma.h>
#include <palmapper/JunctionMap.h>

const float QPalma::NON_CONSENSUS_SCORE = -123456;

static const bool perform_extra_checks = false ;
static const std::string verbose_read_id = "C0HTRACXX120406:3:1301:16294:65782"  ;
static const int verbose_read_level = 0 ;

#define TIME_CODE(x)

inline void get_vector_IUPAC(char c, std::vector<int> &l)
{
	c=toupper(c);

	switch (c)
	{		
	case 'A':
		l[0]=1;
		break;
	case 'B':
		l[1]=1;
		l[2]=1;
		l[3]=1;
		break;
	case 'C':
		l[1]=1;
		break;
	case 'D':
		l[0]=1;
		l[2]=1;
		l[3]=1;
		break;
	case 'G':
		l[2]=1;
		break;
	case 'H':
		l[0]=1;
		l[1]=1;
		l[3]=1;
		break;
	case 'K':
		l[2]=1;
		l[3]=1;
		break;
	case 'M':
		l[0]=1;
		l[1]=1;
		break;
	case 'N':
		l[0]=1;
		l[1]=1;
		l[2]=1;
		l[3]=1;
		break;
	case 'R':
		l[0]=1;
		l[2]=1;
		break;
	case 'S':
		l[1]=1;
		l[2]=1;
		break;
	case 'T':
		l[3]=1;
		break;
	case 'V':
		l[0]=1;
		l[1]=1;
		l[2]=1;
		break;
	case 'W':
		l[0]=1;
		l[3]=1;
		break;
	case 'Y':
		l[1]=1;
		l[3]=1;
		break;
	}
}

inline char get_IUPAC_code(char c1, char c2)
{
	
	std::vector<int> temp(4,0);
	get_vector_IUPAC(c1,temp);
	get_vector_IUPAC(c2,temp);
	
	if (temp[0]==1 && temp[1]==1 && temp[2]==1 && temp[3]==1)
		return 'N';
	if  (temp[0]==1 && temp[1]==1 && temp[2]==1)
		return 'V';
	if  (temp[0]==1 && temp[1]==1 && temp[3]==1)
		return 'H';
	if  (temp[0]==1 && temp[2]==1 && temp[3]==1)
		return 'D';
	if  (temp[1]==1 && temp[2]==1 && temp[3]==1)
		return 'B';
	if  (temp[0]==1 && temp[1]==1)
		return 'M';
	if  (temp[0]==1 && temp[2]==1)
		return 'R';
	if  (temp[0]==1 && temp[3]==1)
		return 'W';
	if  (temp[1]==1 && temp[2]==1)
		return 'S';
	if  (temp[1]==1 && temp[3]==1)
		return 'Y';
	if  (temp[2]==1 && temp[3]==1)
		return 'K';
	if  (temp[0]==1)
		return 'A';
	if  (temp[1]==1)
		return 'C';
	if  (temp[2]==1)
		return 'G';
	if  (temp[3]==1)
		return 'T';

	assert(false);
	
	return c1;
	
}

bool IUPACMatch(char c1, char c2)
{
	if ((c1=='a' || c1=='A' || c1=='c' || c1=='C' || c1=='g' || c1=='G' || c1=='t' || c1=='T') && (c2=='a' || c2=='A' || c2=='c' || c2=='C' || c2=='g' || c2=='G' || c2=='t' || c2=='T'))
		return c1==c2 ;
	
	std::vector<int> l1,l2 ;
	l1.resize(4,0) ;
	l2.resize(4,0) ;	
	get_vector_IUPAC(c1, l1) ;
	get_vector_IUPAC(c1, l2) ;
	for (unsigned int i=0; i<l1.size() && i<l2.size(); i++)
		if (l1[i]==1 && l2[i]==1)
			return true ;
	return false ;
}


int map_back(char c)
{
	switch (c){
	case 'A':
	case 'a':
		return 1;
	case 'c':
	case 'C':
		return 2;
	case 'g':
	case 'G':
		return 3;
	case 't':
	case 'T':
		return 4;
	case 'n':
	case 'N':
		return 5;
	case '-':
		return 0;
	case '*':
		return 6;
	}
	return 7;
}

bool alignment_pass_filters (int min_intron_len, int max_intron_len,  int mm, int gaps, int num_exons, int min_exon_len, bool remapping)
{
	//fprintf(stdout,"%i<intron<%i, %i<exon\n",_config.SPLICED_SHORTEST_INTRON_LENGTH,_config.SPLICED_LONGEST_INTRON_LENGTH,_config.SPLICED_MIN_SEGMENT_LENGTH);

	return (max_intron_len<_config.SPLICED_LONGEST_INTRON_LENGTH || remapping) && (min_intron_len>=_config.SPLICED_SHORTEST_INTRON_LENGTH || remapping) && (mm <= _config.NUM_MISMATCHES && gaps <= _config.NUM_GAPS && mm+gaps <= _config.NUM_EDIT_OPS) &&(num_exons ==1 ||(num_exons >= 2 && (num_exons <= _config.SPLICED_MAX_INTRONS+1 || remapping) && (min_exon_len >= _config.SPLICED_MIN_SEGMENT_LENGTH || remapping))) ;
}

void get_annotated_splice_positions( std::vector<int> &pos, JunctionMap &annotatedjunctions,const char * type, int start, int end, int chr, char strand)
{
	// find a lower bound on the index with binary search
	annotatedjunctions.lock() ;
	//std::deque<Junction>::iterator it = annotatedjunctions.junctionlist[chrN].begin();
	std::deque<Junction>::iterator it = my_lower_bound_by_start(annotatedjunctions.junctionlist_by_start[chr].begin(), annotatedjunctions.junctionlist_by_start[chr].end(), start-100) ;

	while (it != annotatedjunctions.junctionlist_by_start[chr].end() && (*it).start <= end)
	{

		if ((*it).strand!=strand){
			it++;
			continue;
		}
		//fprintf(stdout,"Splice junction: %i-%i %c (%i-%i %c)\n", (*it).start, (*it).end, (*it).strand,start, end,strand );		
		if ((*it).start>=start && (*it).start< end){
			if (strcmp(type, "don")==0 && strand=='+') 
				pos.push_back((*it).start);
			if (strcmp(type, "acc")==0 && strand=='-') 
				pos.push_back((*it).start);
		}
		if ((*it).end>=start && (*it).end< end){
			if (strcmp(type, "acc")==0 && strand=='+') 
				pos.push_back((*it).end);
			if (strcmp(type, "don")==0 && strand=='-') 
				pos.push_back((*it).end);
		}
		it++;
		
	}
	
	annotatedjunctions.unlock();
}

void get_annotated_splice_sites(std::vector<int> &acc_pos, std::vector<int> &don_pos, JunctionMap &annotatedjunctions, int start, int end, int chr, char strand)
{
	//fprintf(stdout,"Get splice sites\n");
	// find a lower bound on the index with binary search
	annotatedjunctions.lock() ;
	//std::deque<Junction>::iterator it = annotatedjunctions.junctionlist[chrN].begin();
	std::deque<Junction>::iterator it = my_lower_bound_by_start(annotatedjunctions.junctionlist_by_start[chr].begin(), annotatedjunctions.junctionlist_by_start[chr].end(), start-100) ;

	while (it != annotatedjunctions.junctionlist_by_start[chr].end() && (*it).start <= end)
	{


		if ((*it).strand!=strand){
			it++;
			continue;
		}

		//fprintf(stdout,"Splice junction: %i-%i %c (%i-%i %c) cons=%i\n", (*it).start, (*it).end, (*it).strand,start, end,strand, (*it).consensus );		
		
		if ((*it).start>=start && (*it).start< end){
			if (strand=='+'){
				don_pos.push_back((*it).start);
			}
			else
				acc_pos.push_back((*it).start);
		}
		if ((*it).end>=start && (*it).end< end){
			if (strand=='+'){
				acc_pos.push_back((*it).end);
			}
			
			else
				don_pos.push_back((*it).end);
		}
		it++;
		

			
	}
	
	annotatedjunctions.unlock();
}


int QPalma::get_transcription_direction(int side,int orientation) const
{

	//Return 1 if the transcription is on the forward direction
	//Return -1 if the transcription is on the negative direction
	//Return 0 if we don't know

	//No strand specific information: need to try all combinations to look for splice sites
	if (side<=-1 || _config.PROTOCOL == PROTOCOL_UNSTRANDED)
		return 0;
	
	
	//Left reads: 
	//same orientation than the transcription for first strand protocol
	//opposite orientation than the transcription for second strand protocol
	if (side==1){
		if (_config.PROTOCOL == PROTOCOL_FIRST)
			return (orientation == 0)?1:-1;
		else
			return (orientation == 0)?-1:1;
	}
	
	//Right reads:
	//opposite orientation than the transcription for first strand protocol
	//same orientation than the transcription for second strand protocol
	if (_config.PROTOCOL == PROTOCOL_FIRST)
		return (orientation == 0)?-1:1;
	else
		return (orientation == 0)?1:-1;
		   
}

QPalma::QPalma(Genome* genome_, GenomeMaps* genomemaps_, int verbosity_)
//  : 	verbosity(4), MIN_NUM_MATCHES(_config.QPALMA_MIN_NUM_MATCHES)
	: 	verbosity(verbosity_), MIN_NUM_MATCHES(_config.QPALMA_MIN_NUM_MATCHES)  
{
	genome=genome_ ;
	genomemaps = genomemaps_ ;
	alignment_parameters = NULL;
	
	if (_config.SPLICED_HITS)
	{
		if (_config.QPALMA_FILE.length()>0 || _config.NO_QPALMA)
		{
			int ret=init_alignment_parameters(_config.QPALMA_FILE) ;
			if (ret!=0)
			{
				fprintf(stderr, "init_alignment_parameters failed\n") ;
				exit(ret) ;
			}
		}
		else
		{
				fprintf(stderr, "alignment parameter file not provided (use -qpalma option)\n") ;
				exit(-1) ;
		}
		
		if (_config.ACC_FILES.length()>0)
		{
			int ret=check_splice_files(_config.ACC_FILES) ;
			if (ret!=0)
			{
				fprintf(stderr, "check_splice_file failed\n") ;
				exit(ret) ;
			}
		}
		else
		{
			if (!_config.NO_SPLICE_PREDICTIONS && !_config.ANNOTATED_SPLICE_SITES_FILE.length()>0)
			{
				fprintf(stderr, "acceptor splice site predictions file not provided (use -acc option)\n") ;
				exit(-1) ;
			}
		}
		
		if (_config.DON_FILES.length()>0)
		{
			int ret=check_splice_files(_config.DON_FILES) ;
			if (ret!=0)
			{
				fprintf(stderr, "check_splice_file failed\n") ;
				exit(ret) ;
			}
		}
		else
		{
			if (!_config.NO_SPLICE_PREDICTIONS && !_config.ANNOTATED_SPLICE_SITES_FILE.length()>0)
			{
				fprintf(stderr, "donor splice site predictions file not provided (use -acc option)\n") ;
				exit(-1) ;
			}
		}
	}
}

QPalma::Result::Result(Read const &read, QPalma const &qpalma)
	:	_qpalma(qpalma), _read(read)
{
	if (&qpalma == NULL)
		return;
	qpalma_filter_reason = -1 ;
	for (int ori=0; ori<2; ori++)
		for (uint32_t i = 0; i < qpalma.genome->nrChromosomes(); i++)
		{
			std::vector<region_t *> r; 
			regions[ori].push_back(r);
		}
}


QPalma::~QPalma()
{
	clean_alignment_parameters() ;
}


int QPalma::check_splice_files(std::string file_template)
{
	char basename[1000] ;
	
	for (int chr=1; chr<(int)genome->nrChromosomes(); chr++)
	{
		char strand='+' ;
		
		while (1)
		{
			sprintf(basename, file_template.c_str(), chr, strand) ;
		    char posname[1000] ;
			sprintf(posname, "%s.pos", basename) ;
		    char confcumname[1000] ;
			sprintf(confcumname, "%s.Conf_cum", basename) ;

			FILE *fd=fopen(posname, "r") ;
			if (fd==NULL)
			{
				fprintf(stderr, "%s does not exist\n", posname) ;
				return -1 ;
			}
			int ret = fseek(fd, -4, SEEK_END) ;
			if (ret==0) // file may be empty
			{
				int buf;
				int retu = fread(&buf, sizeof(int), 1, fd) ;
				if (retu!=1)
				{
					fprintf(stderr, "fread error on %s\n", posname) ;
					return -1 ;
				}
				if (buf>(int)genome->chromosome(chr-1).length())
				{
					fprintf(stderr, "splice site position (%i,%s) greater than chromosome length (%i, %s)\n", buf, posname, genome->chromosome(chr-1).length(),genome->chromosome(chr-1).desc()) ;
					return -1 ;
				}
				//else
				//	fprintf(stdout, "last_pos=%i, size=%i\n", buf, genome->chromosome(chr-1).length()) ;
			}
			//else
			//fprintf(stderr, "fseek %s failed\n", posname) ;
			fclose(fd) ;

			fd=fopen(confcumname, "r") ;
			if (fd==NULL)
			{
				fprintf(stderr, "%s does not exist\n", confcumname) ;
				return -1 ;
			}
			fclose(fd) ;
			
			struct stat pos_stat ;
			int statret = stat(posname, & pos_stat) ;
			if (statret!=0)
			{
			    fprintf(stderr, "stat on %s failed\n", posname) ;
			    return -1 ;
			}

			struct stat confcum_stat ;
			statret = stat(confcumname, & confcum_stat) ;
			if (statret!=0)
			{
			    fprintf(stderr, "stat on %s failed\n", confcumname) ;
			    return -1 ;
			}
			
			if (pos_stat.st_size!=confcum_stat.st_size)
			    fprintf(stdout, "**The two files %s and %s have different sizes (%ld, %ld).**\n**Continuing anyway.**\n\n", posname, confcumname, (long int)pos_stat.st_size, (long int)confcum_stat.st_size);

			if (strand=='+')
				strand='-' ;
			else
				break ;
		}
	}
	return 0 ;
}

int QPalma::compare_double(const void *a, const void *b) 
{
	double *ia = (double *) a;
	double *ib = (double *) b;

	return (int) (*ia - *ib);
}

double QPalma::estimate_splice_score_threshold(std::string file_template, double splice_site_threshold) 
{
	char basename[1000] ;

	size_t max_len=0 ;
	int max_chr = 0 ;
	char max_strand='+' ;
	
	for (int chr=0; chr<(int)genome->nrChromosomes(); chr++)
	{
		char strand='+' ;
		
		while (1)
		{
			sprintf(basename, file_template.c_str(), chr+1, strand) ;
		    char posname[1000] ;
			sprintf(posname, "%s.pos", basename) ;
		    char confcumname[1000] ;
			sprintf(confcumname, "%s.Conf_cum", basename) ;

			FILE *fd=fopen(posname, "r") ;
			if (fd==NULL)
			{
				fprintf(stderr, "%s does not exist\n", posname) ;
				return -1 ;
			}
			fseek(fd, 0, SEEK_END) ;
			size_t s=ftell(fd) ;
			if (s>max_len)
			{
				max_len=s ;
				max_chr=chr;
				max_strand=strand;
			}
			
			fclose(fd) ;

			if (strand=='+')
				strand='-' ;
			else
				break ;
		}
	}
	//fprintf(stdout, "Using chr=%i strand=%c for threshold estimates\n", max_chr, max_strand) ;

	IntervalQuery iq;
	const int num_scores = 1;
	const char *score_names[num_scores] = { "Conf_cum" };
	const int num_intervals = 1 ;
	int interval_matrix[num_intervals * 2];
	int cum_length[num_intervals + 1];
	cum_length[0] = 0;
	interval_matrix[0] = 0 ;
	interval_matrix[1] = genome->chromosome(max_chr).length();
	
	// acc
	sprintf(basename, file_template.c_str(), max_chr+1, max_strand) ;
	
	int acc_size = iq.query(basename, (char**) score_names, num_scores,
							interval_matrix, num_intervals);
	assert(acc_size>0) ;
	
	int* acc_pos = NULL ;
	int* acc_index = NULL ;
	double* acc_score = NULL ;
	try
	{
		acc_pos = new int[acc_size] ;
		acc_index = new int[acc_size] ;
		acc_score = new double[acc_size] ;
	}
	catch (std::bad_alloc&)
	{
		fprintf(stderr, "[map_splice_sites] Could not allocate memory\n");
		delete[] acc_pos ;
		delete[] acc_index ;
		delete[] acc_score ;
		return -1;
	}
	
	iq.getResults(acc_pos, acc_index, acc_score);
	iq.cleanup();
	
	try
	{
		std::vector<double> acc_score2(acc_size, 0.0) ;

		for (int i=0; i<acc_size; i++)
			acc_score2[i]=acc_score[i] ;

		std::sort(acc_score2.begin(), acc_score2.end()) ; 
		
		if (_config.VERBOSE)
			fprintf(stdout, "estimate_threshold: min=%1.3f  max=%1.3f %1.2f%% percentile=%1.3f\n", acc_score2[0], acc_score2[acc_size-1], splice_site_threshold*100, acc_score2[acc_size-1-(int)(acc_size*splice_site_threshold)]) ;
		
		return acc_score2[acc_size-1-(int)(acc_size*splice_site_threshold)] ;
	}
	catch (std::bad_alloc&)
	{
		fprintf(stderr, "[map_splice_sites] Could not allocate memory\n");
		delete[] acc_pos ;
		delete[] acc_index ;
		delete[] acc_score ;
		return -1;
	}
	

}

int QPalma::map_splice_sites(std::string file_template, char type, float &splice_site_threshold, bool estimate_thresh, bool do_report)
{
	
	if (estimate_thresh)
	{
		splice_site_threshold=estimate_splice_score_threshold(file_template, splice_site_threshold) ;
		if (splice_site_threshold==-1)
			return -1 ;
	}

	/*
	  This part is disabled. The splice site map is not used anymore. (Is this efficient????)

	char basename[1000] ;

	if (!do_report)
		return 0 ;
	
	for (int chr=0; chr<(int)genome->nrChromosomes(); chr++)
	{
		char strand='+' ;
		
		while (1)
		{
			sprintf(basename, file_template.c_str(), chr+1, strand) ;
			char posname[1000] ;
			sprintf(posname, "%s.pos", basename) ;
			char confcumname[1000] ;
			sprintf(confcumname, "%s.Conf_cum", basename) ;
			
			FILE *fd=fopen(posname, "r") ;
			if (fd==NULL)
			{
				fprintf(stderr, "%s does not exist\n", posname) ;
				return -1 ;
			}
			fclose(fd) ;
			
			fd=fopen(confcumname, "r") ;
			if (fd==NULL)
			{
				fprintf(stderr, "%s does not exist\n", confcumname) ;
				return -1 ;
			}
			fclose(fd) ;
			
			{
				IntervalQuery iq;
				const int num_scores = 1;
				const char *score_names[num_scores] = { "Conf_cum" };
				const int num_intervals = 1 ;
				int interval_matrix[num_intervals * 2];
				int cum_length[num_intervals + 1];
				cum_length[0] = 0;
				interval_matrix[0] = 0 ;
				interval_matrix[1] = genome->chromosome(chr).length();
				
				// acc
				sprintf(basename, file_template.c_str(), chr+1, strand) ;
				
				int acc_size = iq.query(basename, (char**) score_names, num_scores,
										interval_matrix, num_intervals);
				if (acc_size < 0) 
					return -1;
				
				int* acc_pos = NULL ;
				int* acc_index = NULL ;
				double* acc_score = NULL ;
				try
				{
					acc_pos = new int[acc_size] ;
					acc_index = new int[acc_size] ;
					acc_score = new double[acc_size] ;
				}
				catch (std::bad_alloc&)
				{
					fprintf(stderr, "[map_splice_sites] Could not allocate memory\n");
					delete[] acc_pos ;
					delete[] acc_index ;
					delete[] acc_score ;
					return -1;
				}
				
				iq.getResults(acc_pos, acc_index, acc_score);
				iq.cleanup();
				
				int num_reported = 0 ;
				if (do_report)
					for (int i=0; i<acc_size; i++)
					{
						if (acc_score[i]>splice_site_threshold)
						{
							num_reported++ ;
							report_splice_site(chr, acc_pos[i], strand, type) ;
						}
					}
				
				delete[] acc_pos ;
				delete[] acc_index ;
				delete[] acc_score ;
				
				//fprintf(stderr, "reported %i splice sites for %i %c\n", num_reported, chr, strand) ;
			}
			
			if (strand=='+')
				strand='-' ;
			else
				break ;
		}
	}
	*/

	return 0 ;
	
}

int QPalma::init_alignment_parameters(std::string qpalma_file)
{
	// initialize parameters for alignment
	if (alignment_parameters == NULL)
	{
		alignment_parameters = (struct alignment_parameter_struct*) malloc(
			sizeof(struct alignment_parameter_struct));
		
		if (alignment_parameters == NULL) 
		{
			fprintf(stderr, "[init_alignment_parameters] Could not allocate memory\n");
			return -1;
		}
		
		int ret = init_spliced_align(qpalma_file.c_str(),
									 alignment_parameters->h, alignment_parameters->a,
									 alignment_parameters->d, alignment_parameters->qualityPlifs,
									 alignment_parameters->num_qualityPlifs,
									 alignment_parameters->matchmatrix,
									 alignment_parameters->matchmatrix_dim,
									 alignment_parameters->quality_offset);
		if (ret < 0)
			return ret;
		return 0 ;
	} 
	else
		assert(0) ;
	return -1 ;
}

int QPalma::clean_alignment_parameters()
{
	if (alignment_parameters==NULL)
		return -1 ;
	
	free(alignment_parameters->h.name) ;
	free(alignment_parameters->h.limits) ;
	free(alignment_parameters->h.penalties) ;
	free(alignment_parameters->a.name) ;
	free(alignment_parameters->a.limits) ;
	free(alignment_parameters->a.penalties) ;
	free(alignment_parameters->d.name) ;
	free(alignment_parameters->d.limits) ;
	free(alignment_parameters->d.penalties) ;
	for (int i=0; i<alignment_parameters->num_qualityPlifs; i++)
	{
		free(alignment_parameters->qualityPlifs[i].name) ;
		free(alignment_parameters->qualityPlifs[i].limits) ;
		free(alignment_parameters->qualityPlifs[i].penalties) ;
	}
	free(alignment_parameters->qualityPlifs) ;
	free(alignment_parameters->matchmatrix) ;
	free(alignment_parameters) ;
	alignment_parameters = NULL ;

	return 0 ;
}


	
int QPalma::read_plif(FILE *fd, struct penalty_struct &plif) {
	char buf[1000] ;

	int narg = fscanf(fd, "%1000s:\t", buf);
	assert(narg==1);
	if (strlen(buf) > 1 && buf[strlen(buf) - 1] == ':')
		buf[strlen(buf) - 1] = 0;
	plif.name = strdup(buf);
	assert(plif.name!=NULL);

	narg = fscanf(fd, "%i\t%i\t\t", &plif.min_len, &plif.max_len);
	//fprintf(stdout, "narg=%i\n", narg) ;
	assert(narg==2);
	//assert(strlen(buf)==0) ;
	plif.transform = T_LINEAR;

	double values[100];
	int values_cnt = 0;
	while (1) {
		narg = fscanf(fd, "%lf,", &values[values_cnt]);
		if (narg != 1)
			break;
		//fprintf(stdout, "%i:%f  ", values_cnt, values[values_cnt]) ;
		values_cnt++;
		assert(values_cnt<100);
	}

	plif.len = values_cnt / 2;
	plif.limits = (double*) malloc(sizeof(double) * plif.len);
	plif.penalties = (double*) malloc(sizeof(double) * plif.len);

	if (plif.limits == NULL || plif.penalties == NULL) {
		fprintf(stderr,	"[read_plif] Could not allocate memory\n");
		free(plif.limits);
		free(plif.penalties);
		return -1;
	}

	for (int i = 0; i < plif.len; i++) {
		plif.limits[i] = values[i];
		plif.penalties[i] = values[plif.len + i];
	}
	plif.cache = NULL;
	plif.use_svm = 0;
	plif.next_pen = NULL;
	return 0;
}

int QPalma::read_matrix(FILE* fd, double *& matrix, char*& name, int dims[2]) {
	char buf[1000];

	int narg = fscanf(fd, "%1000s:\t", buf);
	assert(narg==1);
	if (buf[strlen(buf) - 1] == ':')
		buf[strlen(buf) - 1] = '\0';
	name = strdup(buf);
	assert(name!=NULL);

	narg = fscanf(fd, "%i\t%i\t", &dims[0], &dims[1]);
	assert(narg==2);

	matrix = (double*) malloc(sizeof(double) * dims[0] * dims[1]);

	if (matrix == NULL) {
		fprintf(stderr,
				"[read_matrix] Could not allocate memory (name = %s)\n", name);
		return -1;
	}

	int cnt = 0;
	while (1) {
		narg = fscanf(fd, "%lf,", &matrix[cnt]);
		if (narg != 1)
			break;
		assert(cnt<dims[0]*dims[1]);
		cnt++;
	}
	return 0;
}

void QPalma::add_const_to_plif(struct penalty_struct & p, double INDEL_PENALTY)
{
	for (int i=0; i<p.len; i++)
		p.penalties[i]+=INDEL_PENALTY ;
}

int check_qpalma_scores(struct penalty_struct *&qualityPlifs, int &num_qualityPlifs,
						double*&matchmatrix, int verbosity)
{
		int num_qpalma_score_fixes=0 ;
		for (int i=0; i<qualityPlifs[0].len; i++)
			for (int rna=1; rna<6; rna++)
			{
				double & m_score=qualityPlifs[(rna-1)*6+rna].penalties[i] ;
				double max_mm_score=-1000 ;
				double min_mm_score=1000 ;
				double max_ins_score=matchmatrix[1] ;
				for (int dna=1; dna<5; dna++)
				{
					double mm_score=qualityPlifs[(rna-1)*6+dna].penalties[i] ;
					if (dna!=rna && max_mm_score<mm_score)
						max_mm_score=mm_score ;
					if (dna!=rna && min_mm_score>mm_score)
						min_mm_score=mm_score ;

					double & ins_score=matchmatrix[dna] ;
					if (ins_score>max_ins_score)
						max_ins_score=ins_score ;
				}
				double & del_score=qualityPlifs[(rna-1)*6].penalties[i] ;
				const double min_gap = 1e-4 ;
				
				if (m_score<max_mm_score+min_gap || m_score<del_score+min_gap || m_score<max_ins_score+min_gap)
				{
					char acgtn[7]="-ACGTN" ;
					if (acgtn[rna]!='N' && verbosity>0)
						fprintf(stdout, "Warning: unfavorable qpalma match parameters for %i:%c: m_score=%2.3f max_mm_score=%2.3f del_score=%2.3f max_ins_score=%2.3f\n", i, acgtn[rna], m_score, max_mm_score, del_score, max_ins_score) ;
					m_score=max(max_mm_score,max(del_score, max_ins_score))+min_gap ;
					num_qpalma_score_fixes++ ;
				}
				if (min_mm_score<del_score+min_gap)
				{
					char acgtn[7]="-ACGTN" ;
					//if (verbosity >= (rna==5))
					if (acgtn[rna]!='N' && verbosity>0)
						fprintf(stdout, "Warning: unfavorable qpalma del parameters for %i:%c: min_mm_score=%2.3f del_score=%2.3f max_ins_score=%2.3f\n", i, acgtn[rna], min_mm_score, del_score, max_ins_score) ;
					del_score=min_mm_score-min_gap ;
					num_qpalma_score_fixes++ ;
				}
				if (min_mm_score<max_ins_score+min_gap)
				{
					char acgtn[7]="-ACGTN" ;
					//if (verbosity >= (rna==5))
					if (acgtn[rna]!='N' && verbosity>0)
						fprintf(stdout, "Warning: unfavorable qpalma ins parameters for %i:%c: min_mm_score=%2.3f del_score=%2.3f max_ins_score=%2.3f\n", i, acgtn[rna], min_mm_score, del_score, max_ins_score) ;
					for (int dna=1; dna<5; dna++)
						if (min_mm_score<matchmatrix[dna]+min_gap)
						{
							matchmatrix[dna]=min_mm_score-min_gap ;
							num_qpalma_score_fixes++ ;
						}
				}
			}
		return num_qpalma_score_fixes ;
}


int QPalma::init_spliced_align(const char *fname, struct penalty_struct &h,
							   struct penalty_struct &a, struct penalty_struct &d,
							   struct penalty_struct *&qualityPlifs, int &num_qualityPlifs,
							   double*&matchmatrix, int dims[2], int &quality_offset) 
{

	if (_config.NO_QPALMA){
		
		//No scoring function for intron length (h), acceptor sites (a), donor sites (d), edit operations
		init_penalty_struct(h);
		init_penalty_struct(a);
		init_penalty_struct(d);
		qualityPlifs=NULL;
		num_qualityPlifs=0;

		// No quality-based score
		quality_offset=33; //Sanger Fastq format
		
		//Fill matchmatrix from mismatch, gap and match penalties
		//   - A C G T N
		// - x x x x x x
		// A x x x x x x
		// C x x x x x x
		// G x x x x x x
		// T x x x x x x
		// N x x x x x x
		dims[0]=6;
		dims[1]=6;
		
		matchmatrix = (double*) malloc(sizeof(double) * dims[0]*dims[1]);

		if (matchmatrix == NULL) {
			fprintf(stderr,
					"[init_spliced_align] Could not allocate memory\n");
			return -1;
		}
		
		
		for (int i=0;i<dims[0];i++){ //DNA
			for (int j=0;j<dims[1];j++){ //Read

				if (i==0){ // Gap on DNA
					matchmatrix[i*dims[0]+j] = -_config.GAP_SCORE;
					continue;
				}
				if (j == 0){ //Gap on Read
					matchmatrix[i*dims[0]+j] = -_config.GAP_SCORE;
					continue;
				}
				
				if (i==j){ //Match
					matchmatrix[i*dims[0]+j] = -_config.M_SCORE;
					continue;
				}
				
				//Mismatch
				matchmatrix[i*dims[0]+j] = -_config.MM_SCORE;
				
			}
		}
	}
	else
	{
		num_qualityPlifs = 30;
		qualityPlifs = (struct penalty_struct*) malloc(
			sizeof(struct penalty_struct) * num_qualityPlifs);
	
		if (qualityPlifs == NULL) {
			fprintf(stderr, "[init_spliced_align] Could not allocate memory\n");
			return -1;
		}
		char *matrix_name = NULL;
	
		FILE * fd = Util::openFile(fname, "r");
		Util::skip_comment_lines(fd) ;
		int ret = read_plif(fd, h);
		if (ret < 0)
			return ret;

		assert(strcmp(h.name, "h")==0);
		h.max_len = _config.SPLICED_LONGEST_INTRON_LENGTH ;
		for (int i = 0; i < h.len; i++)
			if (h.limits[i] > 2)
				h.limits[i] = 2;

		Util::skip_comment_lines(fd) ;
		ret = read_plif(fd, d);
		if (ret < 0)
			return ret;
		assert(strcmp(d.name, "d")==0);
		d.use_svm = 1;

		Util::skip_comment_lines(fd) ;
		ret = read_plif(fd, a);
		assert(strcmp(a.name, "a")==0);
		if (ret < 0)
			return ret;
		a.use_svm = 1;
		for (int i = 0; i < num_qualityPlifs; i++)
		{
			Util::skip_comment_lines(fd) ;
			read_plif(fd, qualityPlifs[i]);
		}
		Util::skip_comment_lines(fd) ;
		ret = read_matrix(fd, matchmatrix, matrix_name, dims);
		if (ret < 0)
		{
			fprintf(stderr, "init_spliced_align: reading mmatrix failed\n") ;
			return ret;
		}
		assert(strcmp(matrix_name, "mmatrix")==0);
		free(matrix_name);
		matrix_name=NULL ;

		if (_config.QPALMA_INDEL_PENALTY!=0.0)
		{
			fprintf(stdout, "adding %lf to indel penalties\n", _config.QPALMA_INDEL_PENALTY) ;
			
			add_const_to_plif(qualityPlifs[0], -_config.QPALMA_INDEL_PENALTY) ;
			add_const_to_plif(qualityPlifs[6], -_config.QPALMA_INDEL_PENALTY) ;
			add_const_to_plif(qualityPlifs[12], -_config.QPALMA_INDEL_PENALTY) ;
			add_const_to_plif(qualityPlifs[18], -_config.QPALMA_INDEL_PENALTY) ;
			add_const_to_plif(qualityPlifs[24], -_config.QPALMA_INDEL_PENALTY) ;
			assert(dims[0]*dims[1]==6) ;
			for (int i=0; i<6; i++)
				matchmatrix[i] -= _config.QPALMA_INDEL_PENALTY ;
		}

		/*fprintf(stderr, "Warning: hack") ;
		add_const_to_plif(qualityPlifs[0+1], 10) ;
		add_const_to_plif(qualityPlifs[6+2], 10) ;
		add_const_to_plif(qualityPlifs[12+3], 10) ;
		add_const_to_plif(qualityPlifs[18+4], 10) ;
		add_const_to_plif(qualityPlifs[24+5], 10) ;*/
		
		int num_qpalma_score_fixes=check_qpalma_scores(qualityPlifs, num_qualityPlifs, matchmatrix, verbosity) ;
		if (num_qpalma_score_fixes>0 && verbosity>=0)
			fprintf(stdout, "Warning: fixed %i qpalma score elements\n", num_qpalma_score_fixes) ;
		assert(check_qpalma_scores(qualityPlifs, num_qualityPlifs, matchmatrix, 0)==0) ;

		if (!_config.NO_QPALMA_SCORE_FIX)
		{
			double s=0 ;

			for (int i=1; i<6; i++)
				s+=matchmatrix[i] ;
			s+=qualityPlifs[0].penalties[qualityPlifs[0].len-1] ;
			s+=qualityPlifs[6].penalties[qualityPlifs[6].len-1] ;
			s+=qualityPlifs[12].penalties[qualityPlifs[12].len-1] ;
			s+=qualityPlifs[18].penalties[qualityPlifs[18].len-1] ;
			s+=qualityPlifs[24].penalties[qualityPlifs[24].len-1] ;
			_config.GAP_SCORE = s/10 ;

			s=0 ;
			s+=qualityPlifs[0+1].penalties[qualityPlifs[0+1].len-1] ;
			s+=qualityPlifs[6+2].penalties[qualityPlifs[6+2].len-1] ;
			s+=qualityPlifs[12+3].penalties[qualityPlifs[12+3].len-1] ;
			s+=qualityPlifs[18+4].penalties[qualityPlifs[18+4].len-1] ;
			s+=qualityPlifs[24+5].penalties[qualityPlifs[24+5].len-1] ;
			_config.M_SCORE = s/5 ;

			s=0 ;
			s+=qualityPlifs[0+2].penalties[qualityPlifs[0+2].len-1] ;
			s+=qualityPlifs[0+3].penalties[qualityPlifs[0+3].len-1] ;
			s+=qualityPlifs[0+4].penalties[qualityPlifs[0+4].len-1] ;
			s+=qualityPlifs[0+5].penalties[qualityPlifs[0+5].len-1] ;

			s+=qualityPlifs[6+1].penalties[qualityPlifs[6+1].len-1] ;
			s+=qualityPlifs[6+3].penalties[qualityPlifs[6+3].len-1] ;
			s+=qualityPlifs[6+4].penalties[qualityPlifs[6+4].len-1] ;
			s+=qualityPlifs[6+5].penalties[qualityPlifs[6+5].len-1] ;

			s+=qualityPlifs[12+1].penalties[qualityPlifs[12+1].len-1] ;
			s+=qualityPlifs[12+2].penalties[qualityPlifs[12+2].len-1] ;
			s+=qualityPlifs[12+4].penalties[qualityPlifs[12+4].len-1] ;
			s+=qualityPlifs[12+5].penalties[qualityPlifs[12+5].len-1] ;

			s+=qualityPlifs[18+1].penalties[qualityPlifs[18+1].len-1] ;
			s+=qualityPlifs[18+2].penalties[qualityPlifs[18+2].len-1] ;
			s+=qualityPlifs[18+3].penalties[qualityPlifs[18+3].len-1] ;
			s+=qualityPlifs[18+5].penalties[qualityPlifs[18+5].len-1] ;

			s+=qualityPlifs[24+1].penalties[qualityPlifs[24+1].len-1] ;
			s+=qualityPlifs[24+2].penalties[qualityPlifs[24+2].len-1] ;
			s+=qualityPlifs[24+3].penalties[qualityPlifs[24+3].len-1] ;
			s+=qualityPlifs[24+4].penalties[qualityPlifs[24+4].len-1] ;
			_config.MM_SCORE = s/20 ;

			_config.MM_SCORE = _config.M_SCORE - _config.MM_SCORE ;
			_config.GAP_SCORE = _config.M_SCORE - _config.GAP_SCORE ;
			_config.M_SCORE = 0 ;

			assert(_config.MM_SCORE>=0) ;
			assert(_config.GAP_SCORE>=0) ;
			if (_config.GAP_SCORE<_config.MM_SCORE)
				fprintf(stdout, "Warning: QPALMA parameters specify unfavourable combination of scores for mismatches and gaps\n") ;
			
			//fprintf(stdout, "GAP=%f, M=%f, MM=%f\n", _config.GAP_SCORE, _config.M_SCORE, _config.MM_SCORE) ;
		}

		double *quality_offset_matrix ;
		int quality_offset_dims[2] ;
		Util::skip_comment_lines(fd) ;
		ret = read_matrix(fd, quality_offset_matrix, matrix_name, quality_offset_dims);
		if (ret < 0)
		{
			fprintf(stderr, "init_spliced_align: reading quality_offset failed\n") ;
			return ret;
		}
		assert(strcmp(matrix_name, "prb_offset")==0);
		free(matrix_name);
		assert(quality_offset_dims[0]==1) ;
		assert(quality_offset_dims[1]==1) ;
		quality_offset=(int)quality_offset_matrix[0] ;
		free(quality_offset_matrix) ;

		fclose(fd) ;
	}

	return 0;
}



int QPalma::get_splicesite_positions(std::string file_template, const char * type, Chromosome const &chr, char strand, int start, int end, float thresh, bool store_pos,
									 std::vector<int> &positions) const
{
	char basename[1000] ;

	// some bounds checking/corrections
	if (start<0)
		start=0 ;
	if (start>(int)chr.length())
		start=chr.length()-1 ;
	
	if (end<0)
		end=0 ;
	if (end>(int)chr.length())
		end=chr.length()-1 ;

	int num=0 ;
	if (_config.NO_SPLICE_PREDICTIONS)
	{
		if (strcmp(type, "acc")==0)
		{
			if (strand=='+')
		    {
				for (int i=start; i<end; i++)
					for (unsigned int j=0; j<_config.ACC_CONSENSUS.size(); j++)
						if (i>start && IUPACMatch(chr[i-1], _config.ACC_CONSENSUS[j][0] && chr[i]==_config.ACC_CONSENSUS[j][1]))
						{
							num++ ;
							if (store_pos)
								positions.push_back(i) ;
						}
		    } 
			else
		    {
				for (int i=start; i<end; i++)
					for (unsigned int j=0; j<_config.ACC_CONSENSUS_REV.size(); j++)
						if (i>start && IUPACMatch(chr[i-1], _config.ACC_CONSENSUS_REV[j][0]) && IUPACMatch(chr[i], _config.ACC_CONSENSUS_REV[j][1]))
						{
							num++ ;
							if (store_pos)
								positions.push_back(i) ;
						}
		    } 
		}
		if (strcmp(type, "don")==0)
		{
			if (strand=='+')
		    {
				for (int i=start; i<end; i++)
					for (unsigned int j=0; j<_config.DON_CONSENSUS.size(); j++)
						if (i<end-1 && IUPACMatch(chr[i],_config.DON_CONSENSUS[j][0]) && IUPACMatch(chr[i+1],_config.DON_CONSENSUS[j][1]))
						{
							num++ ;
							if (store_pos)
								positions.push_back(i) ;
						}
		    } 
			else
		    {
				for (int i=start; i<end; i++)
					for (unsigned int j=0; j<_config.DON_CONSENSUS_REV.size(); j++)
						if (i<end-1 && IUPACMatch(chr[i], _config.DON_CONSENSUS_REV[j][0]) && IUPACMatch(chr[i+1],_config.DON_CONSENSUS_REV[j][1]))
						{
							num++ ;
							if (store_pos)
								positions.push_back(i) ;
						}
		    } 
		}
	}
	else
	{
		IntervalQuery iq;
		const int num_scores = 1;
		const char *score_names[num_scores] = { "Conf_cum" };
		const int num_intervals = 1 ;
		int interval_matrix[num_intervals * 2];
		interval_matrix[0] = start ;
		interval_matrix[1] = end ;
		
		sprintf(basename, file_template.c_str(), chr.nr() + 1, strand) ;
		
		int ss_size = iq.query(basename, (char**) score_names, num_scores,
							   interval_matrix, num_intervals);
		if (ss_size < 0) 
			return -1;
		
		int* ss_pos = NULL ;
		int* ss_index = NULL ;
		double* ss_score = NULL ;
		try
		{
			//int* ss_pos = (int*) malloc(sizeof(int) * ss_size);
			ss_pos = new int[ss_size] ;
			//int* ss_index = (int*) malloc(sizeof(int) * ss_size);
			ss_index = new int[ss_size] ;
			
			//double* ss_score = (double*) malloc(sizeof(double) * ss_size);
			ss_score = new double[ss_size] ;
		}
		catch (std::bad_alloc&)
		{
			if (ss_pos == NULL || ss_index == NULL || ss_score == NULL) 
			{
				fprintf(stderr, "[map_splice_sites] Could not allocate memory\n");
				delete[] ss_pos;
				delete[] ss_index;
				delete[] ss_score;
				return -1;
			}
		}
		
		
		iq.getResults(ss_pos, ss_index, ss_score);
		iq.cleanup();
		
		for (int i=0; i<ss_size; i++)
		{
			if (ss_score[i]>=thresh)
			{
				num++ ;
				if (store_pos)
				{
					positions.push_back(ss_pos[i]) ;
					//scores.push_back(ss_score[i]) ;
				}
			}
		}
		
		delete[] ss_pos ;
		delete[] ss_index ;
		delete[] ss_score ;
	}
	
	return num ;
}

int QPalma::get_string_from_region(Chromosome const &chrN, region_t *region, std::string &str) const
{

	assert(region->end > region->start);

	int32_t length = region->end - region->start;

	str.assign("") ;
	for (int i=0; i<length; i++)
		str.push_back(chrN[region->start+i]) ;
	//str.assign(CHR_SEQ[chrN] + region->start, length);

	// add (up to) three N letters on each side
	if (str.size()>0)
	{
		str[0]='N' ;
		str[str.size()-1]='N' ;
	}
	if (str.size()>1)
	{
		str[1]='N' ;
		str[str.size()-2]='N' ;
	}
	if (str.size()>2)
	{
		str[2]='N' ;
		str[str.size()-3]='N' ;
	}
	return 0 ;
}

void QPalma::Result::add_buffer_to_region(int ori, Chromosome const &chrN, int32_t nregion, int read_len) 
{
	region_t *region = regions[ori][chrN.nr()][nregion];
	int left=_config.SPLICED_CLUSTER_TOLERANCE ;
	int right=_config.SPLICED_CLUSTER_TOLERANCE ;

	int read_pos = region->read_pos ;
	int hit_len = region->hit_len ;
	
	if (read_pos>=0 && hit_len>0)
	{
		//if (ori==0)  // works slower!
		{
			if (region->read_pos+_config.SPLICED_CLUSTER_TOLERANCE > left)
				left = region->read_pos+_config.SPLICED_CLUSTER_TOLERANCE ;
			if (read_len - region->read_pos - region->hit_len + _config.SPLICED_CLUSTER_TOLERANCE > right)
				right = read_len - region->read_pos - region->hit_len + _config.SPLICED_CLUSTER_TOLERANCE ;
		}
		//if (ori==1) // works slower!
		{
			if (region->read_pos+_config.SPLICED_CLUSTER_TOLERANCE > right)
				right = region->read_pos+_config.SPLICED_CLUSTER_TOLERANCE ;
			if (read_len - region->read_pos - region->hit_len + _config.SPLICED_CLUSTER_TOLERANCE > left)
				left = read_len - region->read_pos - region->hit_len + _config.SPLICED_CLUSTER_TOLERANCE ;
		}
	}
	else
	{
		left+=read_len ;
		right+=read_len ;
	}
	
	assert(left>=_config.SPLICED_CLUSTER_TOLERANCE) ;
	assert(right>=_config.SPLICED_CLUSTER_TOLERANCE) ;
	
	if (region->start <= left)
		region->start = 0 ;
	else
		region->start -= left;

	if (region->end < (int) chrN.length())
	{
		if (region->end + right >= (int)chrN.length())
			region->end = chrN.length() - 1; // chr length because end bound is excluded
		else
			region->end += right;
	}
}

/** performs a quicksort on an array output of length size
 * it is sorted from in ascending (for type T) */
void QPalma::qsort(region_t** output, int size) const 
{
	if (size == 2) 
	{
		if (output[0]->start > output[1]->start) {
			region_t *c = output[0];
			output[0] = output[1];
			output[1] = c;
		}
		return;
	}
	//T split=output[random(0,size-1)];
	region_t *split = output[size / 2];

	int32_t left = 0;
	int32_t right = size - 1;

	while (left <= right) {
		while (output[left]->start < split->start)
			left++;
		while (output[right]->start > split->start)
			right--;

		if (left <= right) {
			region_t *c = output[left];
			output[left] = output[right];
			output[right] = c;
			left++;
			right--;
		}
	}

	if (right + 1 > 1)
		qsort(output, right + 1);

	if (size - left > 1)
		qsort(&output[left], size - left);
}

void QPalma::qsort_length(region_t** output, int size) const 
{
	if (size == 2) 
	{
		if (output[0]->end-output[0]->start < output[1]->end-output[1]->start) {
			region_t *c = output[0];
			output[0] = output[1];
			output[1] = c;
		}
		return;
	}
	region_t *split = output[size / 2];

	int32_t left = 0;
	int32_t right = size - 1;

	while (left <= right) {
		while (output[left]->end-output[left]->start > split->end-split->start)
			left++;
		while (output[right]->end-output[right]->start < split->end-split->start)
			right--;

		if (left <= right) {
			region_t *c = output[left];
			output[left] = output[right];
			output[right] = c;
			left++;
			right--;
		}
	}

	if (right + 1 > 1)
		qsort_length(output, right + 1);

	if (size - left > 1)
		qsort_length(&output[left], size - left);
}


void QPalma::print_map(Read const &read, bool* read_map, const char *name)
{
	fprintf(stdout, "# read_map (%s): ", name);
	for (size_t i = 0; i < read.length(); i++)
		if (read_map[i])
			fprintf(stdout, "1");
		else
			fprintf(stdout, "0");
	fprintf(stdout, "\n");
}


////////////////////////////////////////////////////////////////////////////////
// filter


int QPalma::get_num_splicesites(std::string file_template, const char* type, Chromosome const &chr, char strand, int start, int end, float thresh) const
{
	std::vector<int> positions ;
	return get_splicesite_positions(file_template, type, chr, strand, start, end, thresh, false, positions) ;
}


int QPalma::qpalma_filter(Result &result, struct alignment_t *ali, int num_N,JunctionMap &annotatedjunctions) const
{
	assert(ali->exons.size()<=2) ;
	static bool use_ss_for_filtering = true ;
	
	Chromosome const &chr = *ali -> chromosome ;
	int start = ali -> exons[0] ;
	int end = ali -> exons[1] ;
	
	unsigned int num_gaps = ali -> num_gaps_ref ;
	unsigned int num_matches = ali -> num_matches_ref ;

	if (num_gaps > _config.FILTER_BY_MAX_GAPS || result._read.length()-num_matches > _config.FILTER_BY_MAX_MISMATCHES+num_N)
	{
		if (verbosity>=1)
			fprintf(stdout, "filter decides YES: num_gaps=%i, num_mismatches=%i, num_N=%i\n", num_gaps, result._read.length()-num_matches, num_N) ;
		
		result.qpalma_filter_reason=2 ; // cheap positive filter
			
		return 1 ;
	}

	float thresh_acc = _config.FILTER_BY_SPLICE_SITES_THRESH_ACC ;
	float thresh_don = _config.FILTER_BY_SPLICE_SITES_THRESH_DON ;
	int region = _config.FILTER_BY_SPLICE_SITES_REGION ;

	if (region==-1 || !_config.FILTER_BY_SPLICE_SITES)
	{
		result.qpalma_filter_reason=0 ; // cheap negative filter
		return 0 ;
	}

	if ( result._read.length()-num_matches < _config.FILTER_BY_SPLICE_SITES_EDIT_MIN+num_N )
	{
		result.qpalma_filter_reason=0 ; // cheap negative filter
		return 0 ;
	}
	
	if (use_ss_for_filtering)
		try
		{
			int num_ss = 0;
			
			if (_config.ACC_FILES.length()>0 && _config.DON_FILES.length()>0){
				

				num_ss+= get_num_splicesites(_config.ACC_FILES, "acc", chr, '+', start-region, end+region, thresh_acc) + 
				get_num_splicesites(_config.ACC_FILES, "acc", chr, '-', start-region, end+region, thresh_acc) +
				get_num_splicesites(_config.DON_FILES, "don", chr, '+', start-region, end+region, thresh_don) +
				get_num_splicesites(_config.DON_FILES, "don", chr, '-', start-region, end+region, thresh_don) ;
			}
			if (_config.SCORE_ANNOTATED_SPLICE_SITES){
				std::vector<int> positions;
				int chrN = chr.nr();
				
				get_annotated_splice_positions(positions,annotatedjunctions,"acc",start-region, end+region,chrN,'+');
				get_annotated_splice_positions(positions,annotatedjunctions,"don",start-region, end+region,chrN,'+');
				get_annotated_splice_positions(positions,annotatedjunctions,"acc",start-region, end+region,chrN,'-');
				get_annotated_splice_positions(positions,annotatedjunctions,"don",start-region, end+region,chrN,'-');
				num_ss += positions.size();
				//fprintf(stdout,"Positions %i\n",positions.size());
				
				positions.clear();
				

			}
			
			
			if ( num_ss>0 )
			{
				if (verbosity>=1)
					fprintf(stdout, "filter decides YES: num_ss=%i, gaps=%i, num_mismatches=%i, num_N=%i\n", num_ss, num_gaps, result._read.length()-num_matches, num_N) ;
				
				result.qpalma_filter_reason=3 ; // positive splice site filter
				
				return 1 ;
			}
		}
		catch (IntervalQueryException & e)
		{
			fprintf(stdout, "Warning: do not use splice site files for triggering qpalma alignments\n") ;
			use_ss_for_filtering=false ;
		}
	
	if (verbosity>=1)
		fprintf(stdout, "filter decides NO: num_gaps=%i, num_mismatches=%i, num_N=%i\n", num_gaps, result._read.length()-num_matches, num_N) ;
	
	result.qpalma_filter_reason=1 ; // negative splice site filter

	return 0 ;
}


bool QPalma::sort_regions(std::vector<region_t*> & regions, enum sortorder_enum sortorder) const
{
	// Sort regions by start position
	region_t ** arr = NULL ;
	size_t nbr_regions=regions.size();
	try 
    {
		arr = new region_t*[nbr_regions] ;			
    }
	catch (std::bad_alloc&)
    {
		fprintf(stderr, "[sort_regions] ERROR Could not allocate memory\n");
		return false ;
    }
	for (int i = 0; i < (int)nbr_regions; i++)
		arr[i] = regions[i];  
	if (sortorder==sortorder_startpos)
		qsort(arr, nbr_regions);
	else
	{
		assert(sortorder==sortorder_length) ;
		qsort_length(arr, nbr_regions);
	}
	
	for (int i=0; i<(int)nbr_regions; i++)
		regions[i]=arr[i] ;
	delete[] arr ;

	if (perform_extra_checks && sortorder==sortorder_startpos)
		for (int i = 0; i + 1 < (int)nbr_regions; i++)
			assert(regions[i]->start<=regions[i+1]->start);
	return true ;
}


////////////////////////////////////////////////////////////////////////////////
// alignment

/** find long regions included in the set of current regions */
void QPalma::recover_long_regions(Read const &read, std::vector<region_t*> &long_regions_output, std::vector<region_t*> long_regions, std::vector<region_t*> current_regions) const
{
	unsigned int reg=0;
	unsigned int lr=0;
  
	//fprintf(stdout,"**start** %i-%i\n",current_regions.size(),nbr_long_regions);
	while (reg < current_regions.size() && lr < long_regions.size())
	{
		//fprintf(stdout,"current region: %i-%i and current long region: %i-%i\n", current_regions[reg]->start,current_regions[reg]->end, long_regions[lr]->start, long_regions[lr]->end);
    
		if (long_regions[lr]->start >= current_regions[reg]->start && long_regions[lr]->end <= current_regions[reg]->end)
		{
			long_regions_output.push_back(long_regions[lr]);
			lr++;
		}
		else if (long_regions[lr]->start < current_regions[reg]->start)
			lr++;
		else // current long region coordinates after current region
			reg++;
	}
}

/** Gives the relative position on dna sequence to align  */
int QPalma::convert_dna_position(int real_position, const std::vector<size_t> & cum_length, const std::vector<region_t *> &current_regions) const
{
	for (size_t j = 0; j < current_regions.size(); j++)
		if (real_position >= current_regions[j]->start && real_position< current_regions[j]->end)
			return real_position - current_regions[j]->start + cum_length[j];
  
	return -1;
}

void QPalma::reverse_variants(std::vector<Variant> &variants, int dna_len) const
{
	for (unsigned int i = 0; i < variants.size(); i++)
	{
			int start_tmp=variants[i].position;
			variants[i].position=dna_len-(variants[i].end_position);
			variants[i].end_position=dna_len-(start_tmp);
			variants[i].ref_str=reverse(complement(variants[i].ref_str));
			variants[i].variant_str=reverse(complement(variants[i].variant_str));
	}
}

int QPalma::get_first_read_map(Read const &read, bool* read_map) const
{
	if (read_map)
		for(unsigned int i=0; i < read.length(); i++)
			if (read_map[i])
				return i;
	
	return -1;
}

// GenomeMapper's handling of spliced reads is very simplistic. It just finds
// exact 12-mers and returns these as spliced hits. It does _not_ do an
// alignment of any kind. The following routine takes a pair of such
// spliced hits and aligns them against the region of the genome where the
// hits were reported.


void QPalma::print_hit(HIT *hit) {
	fprintf(stdout, "# hit= {start=%i, end=%i, readpos=%i}\n", hit->start,
			hit->end, hit->readpos);
}

void QPalma::print_region(region_t *region, const char * bla)  const
{
	fprintf(stdout, "# region %s = {start=%i, end=%i, ori=%c, erased=%i, from_map=%i}\n", bla,
			region->start, region->end, region->orientation ? '1' : '0', region->erased, region->from_map);
}


void QPalma::Result::delete_regions()
{
	for (int ori = 0; ori < 2; ori++)
		for (int32_t chrN = 0; chrN < (int)regions[ori].size(); chrN++)
			for (int32_t nregion = 0; nregion < (int)regions[ori][chrN].size(); nregion++) 
			{
				if (regions[ori][chrN][nregion]==NULL)
					continue ;
				delete[] regions[ori][chrN][nregion]->read_map;
				regions[ori][chrN][nregion]->read_map = NULL;
				delete regions[ori][chrN][nregion];
				regions[ori][chrN][nregion]=NULL ;
			}
}

void QPalma::delete_long_regions(std::vector<std::vector<region_t *> > *long_regions) const
{ 
	//long_regions need to be deleted: deep copy of initial regions into long_regions
	for (int ori = 0; ori < 2; ori++)
		for (int32_t chrN = 0; chrN < (int)long_regions[ori].size(); chrN++)
			for (int32_t nregion = 0; nregion < (int)long_regions[ori][chrN].size(); nregion++) 
			{
				if (long_regions[ori][chrN][nregion]==NULL)
					continue ;
				delete[] long_regions[ori][chrN][nregion]->read_map;
				long_regions[ori][chrN][nregion]->read_map = NULL;
				delete long_regions[ori][chrN][nregion];
				long_regions[ori][chrN][nregion]=NULL ;
			}
}


int QPalma::recapture_hits(Hits &hits, Result &result, bool const non_consensus_search,JunctionMap &annotatedjunctions, VariantMap & variants) const
{
	TIME_CODE(clock_t start_time = clock();)
	
	int myverbosity=verbosity ;
	
	if (myverbosity!=verbose_read_level && std::string(result._read.id())==verbose_read_id)
		myverbosity=verbose_read_level ;

	Read const &read(hits.getRead());

	result.cleanup();

	int32_t num_hits = 0; // TODO debugging only
	int32_t num_hits_dropped = 0; // TODO debugging only
  
	// Examine all hits and construct a list of region where these hits map.
	// regions is a list of region clusters per chromosome, sorted by start position
  
	result.long_regions[0] = result.regions[0];
	result.long_regions[1] = result.regions[1];
	std::vector<std::vector<region_t *> > * const regions = result.regions;
	std::vector<std::vector<region_t *> > * const long_regions = result.long_regions;

	for (unsigned int nump=0; nump<read.prealigned_info.size(); nump++)
	{
		struct t_prealigned & prealigned_info = *read.prealigned_info[nump] ;
		// clean up data generated for the previous read
		
		int best_start=-1, best_end=-1, max_len=0 ;
		for (unsigned int i=0; i<prealigned_info.exons.size(); i+=2)
		{
			if (prealigned_info.exons[i+1]-prealigned_info.exons[i]>max_len)
			{
				max_len = prealigned_info.exons[i+1]-prealigned_info.exons[i] ;
				best_start=prealigned_info.exons[i] ;
				best_end = prealigned_info.exons[i+1] ;
			}
		}
		assert(best_start>=0 && best_end>=0) ;
		
		region_t *new_lregion=NULL;
		try {
			new_lregion = new region_t();
			new_lregion->read_map = new bool[read.length()];
		} catch (std::bad_alloc&) 
		{
			fprintf(stderr, "[recapture_hits] allocating memory for read_map failed\n");
			result.delete_regions();
			delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
			return -1;
		}  
		int delta=-1 ;
		new_lregion->start = best_start + delta ;
		new_lregion->end = best_end+1 + delta ;
		new_lregion->from_map = false ;
		new_lregion->read_pos = -1 ; 
		for (size_t ii=0; ii < read.length(); ii++)
			if (prealigned_info.aligned_positions[ii]==best_start)
			{
				new_lregion->read_pos = ii ;
				break ;
			}
		assert(new_lregion->read_pos>=0) ;
		
		new_lregion->hit_len = new_lregion->end-new_lregion->start ;
		
		for (size_t ii = 0; ii < read.length(); ii++)
			if (prealigned_info.aligned_positions[ii]>=best_start && prealigned_info.aligned_positions[ii]<=best_end)
				new_lregion->read_map[ii] = true ;
			else
				new_lregion->read_map[ii] = false ;
		
		int contig_no=genome->find_desc(prealigned_info.contig.c_str()) ;
		if (contig_no<0)
		{
			fprintf(stderr, "ERROR: Contig name %s not found\n", prealigned_info.contig.c_str()) ;
			return -1 ;
		}
		long_regions[ori_map(prealigned_info.orientation)][contig_no].push_back(new_lregion);
		
		for (unsigned int i=0; i<prealigned_info.exons.size(); i+=2)
		{
			// Create first region for this hit
			region_t *new_region = NULL;
			try {
				new_region = new region_t();
				new_region->read_map = new bool[read.length()];
			} catch (std::bad_alloc&) 
			{
				fprintf(stderr, "[recapture_hits] allocating memory for read_map failed\n");
				result.delete_regions();
				delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
				return -1;
			}  
			
			new_region->start = prealigned_info.exons[i] ;
			new_region->end = prealigned_info.exons[i+1]+1;
			
			new_region->from_map = false ;
			for (size_t ii = 0; ii < read.length(); ii++)
				new_region->read_map[ii] = true ;
			
			for (size_t ii = 0; ii < read.length(); ii++)
				if (prealigned_info.aligned_positions[ii]>=new_region->start && prealigned_info.aligned_positions[ii]<=new_region->end)
					new_region->read_map[ii] = true ;
				else
					new_region->read_map[ii] = false ;
			
			new_region->read_pos = -1 ; 
			for (size_t ii=0; ii < read.length(); ii++)
				if (prealigned_info.aligned_positions[ii]==new_region->start)
				{
					new_region->read_pos = ii ;
					break ;
				}
			//assert(new_region->read_pos>=0) ;
			
			new_region->hit_len = new_region->end-new_region->start ;
			
			regions[ori_map(prealigned_info.orientation)][contig_no].push_back(new_region);
		}
	}
	
	if (myverbosity >= 1)
		fprintf(stdout,	"# [capture_hits] Captured %d hits, dropped %i hits for read %s\n",
				num_hits - num_hits_dropped, num_hits_dropped, read.id());
  
	// sort regions by starting position
	for (int ori = 0; ori < 2; ori++)
		for (int32_t chrN = 0; chrN < (int)regions[ori].size(); chrN++) 
		{
			if (regions[ori][chrN].size() == 0)
				continue;
			bool ret=sort_regions(regions[ori][chrN], sortorder_startpos) ;
			ret = ret && sort_regions(long_regions[ori][chrN], sortorder_startpos) ;
			if (!ret)
			{
				fprintf(stderr, "[capture_hits] allocating memory for sort_regions failed\n");
				result.delete_regions();
				delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
				return -1;
			}
		}


	// TODO: deduplicate code from capture_hits

	// Surround the regions with a buffer of size _config.SPLICED_CLUSTER_TOLERANCE.
	// If any of them overlap before/after extension, merge them to
	// avoid duplicating subsequences.
	//
    for (int ori = 0; ori < 2; ori++){
		for (int32_t chrN = 0; chrN < (int)regions[ori].size(); chrN++)
		{
			Chromosome const &chromosome = genome->chromosome(chrN);
			size_t nbr_regions=regions[ori][chrN].size();
			if (nbr_regions == 0)
				continue;
			if (nbr_regions == 1) 
			{
				// Nothing to merge, just extend the one existing region to include a buffer.
				result.add_buffer_to_region(ori, chromosome, 0, read.length());
				continue;
			}
	
			for (int i = 0; i < (int)nbr_regions; i++) 
				result.add_buffer_to_region(ori, chromosome, i, read.length());
	
			for (int nregion = 0; nregion < (int)nbr_regions - 1; nregion++) 
			{
				if (regions[ori][chrN][nregion]->erased)
					continue;
				size_t next = nregion + 1;
				while (next < nbr_regions && regions[ori][chrN][next]->erased)
					next++;
				if (next >= nbr_regions || regions[ori][chrN][next]->erased)
					continue;
	    
				if ((int32_t) ((regions[ori][chrN][nregion])->end) >= (int32_t) ((regions[ori][chrN][next])->start))
				{
					if (myverbosity>=2)
					{
						print_region(regions[ori][chrN][nregion], "region 1 for merge")  ;
						print_region(regions[ori][chrN][next], "region 2 for merge")  ;
					} ;

					// regions[ori] are overlapping or adjacent. Merge them into one.
					if ((regions[ori][chrN][nregion])->end < (regions[ori][chrN][next])->end)
						(regions[ori][chrN][nregion])->end = (regions[ori][chrN][next])->end;
		
					// merge the two read_maps
					for (size_t i = 0; i < read.length(); i++)
						regions[ori][chrN][nregion]->read_map[i] = regions[ori][chrN][nregion]->read_map[i]
							|| regions[ori][chrN][next]->read_map[i];
		
					// Get rid of the extra region.
					regions[ori][chrN][next]->erased = true;
					delete[] regions[ori][chrN][next]->read_map; 
					regions[ori][chrN][next]->read_map = NULL;
		
					// consider this as a map-region, only if any both were generated from a map (important for downstream filtering)
					regions[ori][chrN][nregion]->from_map = regions[ori][chrN][nregion]->from_map && regions[ori][chrN][next]->from_map ;

					if (myverbosity>=2)
						print_region(regions[ori][chrN][nregion], "region merged")  ;
		
					// make sure this item is looked at again to merge it with the next item if necessary
					nregion-- ;
				}
			}
		}
	}

	TIME_CODE(_stats.qpalma_region_align_time += clock() - start_time;)
	return 0;
}

region_t *gen_new_region(size_t start_pos, size_t end_pos, size_t read_pos, size_t read_len)
{
	// Create first region for this hit
	region_t *new_region = NULL;
	try {
		new_region = new region_t();
		new_region->read_map = new bool[read_len];
	} catch (std::bad_alloc&) 
	{
		fprintf(stderr, "[capture_hits] allocating memory for read_map failed\n");
		return NULL;
	}  
	new_region->erased = false ;
	new_region->start = start_pos ;
	new_region->end = start_pos ;
	new_region->from_map = false ;
	for (size_t ii = 0; ii < read_len; ii++)
		new_region->read_map[ii] = false ;
	
	new_region->read_pos = read_pos ;
	new_region->hit_len = 0 ;
	
	for (size_t ii = 0; ii < end_pos - start_pos && read_pos + ii < read_len; ii++, new_region->hit_len++, new_region->end++)
			 new_region->read_map[read_pos + ii] = true;
    return new_region ;
}
	
int QPalma::capture_hits(Hits &hits, Result &result, bool const non_consensus_search,JunctionMap &annotatedjunctions, VariantMap & variants) const
{
	Read const &read(hits.getRead());

	TIME_CODE(clock_t start_time = clock();)

	int myverbosity=verbosity ;
	
	if (myverbosity!=verbose_read_level && std::string(result._read.id())==verbose_read_id)
		myverbosity=verbose_read_level ;

	// clean up data generated for the previous read

	result.cleanup();

	HIT const *hit;
	int32_t num_hits = 0, num_long_hits=0; // TODO debugging only
	int32_t num_hits_dropped = 0; // TODO debugging only
  
	// Examine all hits and construct a list of region where these hits map.
	// regions is a list of region clusters per chromosome, sorted by start position
  
	result.long_regions[0] = result.regions[0];
	result.long_regions[1] = result.regions[1];
	std::vector<std::vector<region_t *> > * const regions = result.regions;
	std::vector<std::vector<region_t *> > * const long_regions = result.long_regions;
  

	//TODO: Real length of a hit is i-1
	for (int32_t i = read.length(); i >= _config.SPLICED_HIT_MIN_LENGTH_SHORT; i--) 
	{
		hit = hits.HIT_LISTS_OPERATOR[i];
		
		while (hit != NULL) 
		{
		  //fprintf(stdout, "long hit: %i  ", i) ; 
			//TODO: Real length of a hit is i-1
			bool consider_as_long_hit = (i >= _config.SPLICED_HIT_MIN_LENGTH_LONG) || hit->aligned ;
			if (consider_as_long_hit && num_long_hits >= _config.SPLICED_MAX_NUM_ALIGNMENTS) 
			{
				consider_as_long_hit = false;
				if (myverbosity >= 2)
					fprintf(stdout, "# ignoring long hits from now on\n");
				break ;
			}
			if (consider_as_long_hit)
			{
				num_long_hits++;

				region_t *new_lregion=NULL;
				try {
					new_lregion = new region_t();
					new_lregion->read_map = new bool[read.length()];
				} catch (std::bad_alloc&) 
				{
					fprintf(stderr, "[capture_hits] allocating memory for read_map failed\n");
					result.delete_regions();
					delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
					return -1;
				}  
				
				new_lregion->start = hit->start;
				new_lregion->end = hit->start;
				new_lregion->from_map = false ;
				new_lregion->read_pos = hit->readpos ;
				new_lregion->hit_len = 0 ;
				new_lregion->erased = false ;
				
				for (size_t ii = 0; ii < read.length(); ii++)
					new_lregion->read_map[ii] = false ;
				for (size_t ii = 0; ii < hit->end - hit->start && hit->readpos + ii< read.length(); ii++,  new_lregion->hit_len++, new_lregion->end++)
					new_lregion->read_map[hit->readpos + ii] = true;

                new_lregion->end = (new_lregion->end > hit->chromosome->length())?hit->chromosome->length():new_lregion->end;
				long_regions[ori_map(hit->orientation)][hit->chromosome->nr()].push_back(new_lregion);


				if (false) 
				{
					const size_t region_len = 100000 ;
					char buf[region_len+1] ;
					const Chromosome & chr=*hit->chromosome ;
					for (size_t ii=0; hit->end+ii<hit->chromosome->length() && ii<region_len; ii++)
					{
						buf[ii]=chr[hit->end+ii] ;
						buf[ii+1]=0 ;
					}
					std::string sbuf=buf ;
					char read_str[read.length()+1] ;
					strcpy(read_str, read.data()) ;
					const int k=6 ;
					for (size_t read_pos=hit->readpos+(hit->end-hit->start); read_pos<read.length()-k; read_pos+=k/2)
					{
						size_t found = sbuf.find(&read_str[read_pos], 0, k) ;
						if (found!=std::string::npos)
						{
							//fprintf(stdout, "found kmer match: %ld %ld\n", read_pos, found) ;
							region_t *new_region = gen_new_region(hit->end+found, hit->end+found+k, read_pos, read.length()) ;
							regions[ori_map(hit->orientation)][hit->chromosome->nr()].push_back(new_region);
						}
					}
				}
			}
			
			hit = hit->next;
		}
	}

	//TODO: Real length of a hit is i-1
	for (int32_t i = read.length(); i >= _config.SPLICED_HIT_MIN_LENGTH_SHORT; i--)
	{
		hit = hits.HIT_LISTS_OPERATOR[i];
    
		while (hit != NULL) 
		{
			num_hits++;

			// first check whether it is close enough to a long enough hit
			bool found = false;
			for (int32_t nregion = 0; nregion < (int)long_regions[ori_map(hit->orientation)][hit->chromosome->nr()].size(); nregion++) 
			{
				//if (perform_extra_checks)
				//	assert(long_regions[ori_map(hit->orientation)][hit->chromosome->nr()][nregion]->islong) ;
				
				int32_t rs = long_regions[ori_map(hit->orientation)][hit->chromosome->nr()][nregion]->start - _config.SPLICED_LONGEST_INTRON_LENGTH*_config.SPLICED_MAX_INTRONS - read.length() ;
				if (rs < 0)
					rs = 0;
				int32_t re = long_regions[ori_map(hit->orientation)][hit->chromosome->nr()][nregion]->end + _config.SPLICED_LONGEST_INTRON_LENGTH*_config.SPLICED_MAX_INTRONS + read.length();
				if (((int)hit->start >= rs) && ((int)hit->end <= re)) 
				{
					found = true;
					break;
				}
			}
			if (!found) {
				num_hits_dropped++;
				//	    fprintf(stdout, "dropped %i-%i (%i)\n", hit->start, hit->end,hit->readpos) ;
				hit = hit->next;
				continue;
			}
			if (myverbosity >= 2)
				fprintf(stdout, "# captured %i-%i  (%i,%i) %c real length: %i mism-gaps(%i-%i)\n", hit->start,
						hit->end, hit->readpos, hit->end - hit->start, 
						hit->orientation,i,hit->mismatches,hit->gaps);

			// Create first region for this hit
			region_t *new_region = NULL;
			try {
				new_region = new region_t();
				new_region->read_map = new bool[read.length()];
			} catch (std::bad_alloc&) 
			{
				fprintf(stderr, "[capture_hits] allocating memory for read_map failed\n");
				result.delete_regions();
				delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
				return -1;
			}  
			new_region->erased = false ;
			new_region->start = hit->start;
			new_region->end = hit->start;
			new_region->from_map = false ;
			for (size_t ii = 0; ii < read.length(); ii++)
				new_region->read_map[ii] = hit->aligned ; // fake read_map for reads that could be aligned by genomemapper
			assert(hit->end >= hit->start) ;
			
			new_region->read_pos = hit->readpos ;
			new_region->hit_len = 0 ;
			
			for (size_t ii = 0; ii < hit->end - hit->start && hit->readpos + ii < read.length(); ii++, new_region->hit_len++, new_region->end++)
				new_region->read_map[hit->readpos + ii] = true;
			
            new_region->end = (new_region->end > hit->chromosome->length())?hit->chromosome->length():new_region->end;

			regions[ori_map(hit->orientation)][hit->chromosome->nr()].push_back(new_region);

			hit = hit->next;
		}
	}

	if (myverbosity >= 1)
		fprintf(stdout,	"# [capture_hits] Captured %d hits, dropped %i hits for read %s\n",
				num_hits - num_hits_dropped, num_hits_dropped, read.id());
  
	//find short hits in close vicinity
	//find_short_hits();
  
	if (_config.QPALMA_USE_MAP || _config.QPALMA_USE_SPLICE_SITES)
    {
		// add regions from reporting map
		//const int take_report_map = MASK_MAPPED_READ_BEST | MASK_MAPPED_READ | MASK_SPLICED_READ_BEST | MASK_SPLICED_READ | MASK_MAPPED_REGION ;
		//const int take_report_map = 0 ;
      
		//const int take_report_map = MASK_MAPPED_READ_BEST | MASK_SPLICED_READ_BEST ;
		const int take_report_map = MASK_MAPPED_READ_BEST | MASK_SPLICED_READ_BEST ; 
      
		int added_map_regions = 0 ; 
		int added_map_total_len = 0 ;
      
		for (int ori = 0; ori < 2; ori++)
		{
			for (int32_t chrN = 0; chrN < (int)long_regions[ori].size(); chrN++)
			{
				Chromosome const &chromosome = genome->chromosome(chrN);
				if (long_regions[ori][chrN].size() == 0)
					continue;
				for (int i=0; i<(int)long_regions[ori][chrN].size(); i++)
				{
					int start = long_regions[ori][chrN][i]->start - _config.SPLICED_LONGEST_INTRON_LENGTH*_config.SPLICED_MAX_INTRONS - read.length();
					if (start<0)
						start=0 ;
					int end = long_regions[ori][chrN][i]->end + _config.SPLICED_LONGEST_INTRON_LENGTH*_config.SPLICED_MAX_INTRONS + read.length();
					if (end>(int)chromosome.length())
						end=chromosome.length();
					int midpoint = (end+start)/2 ;
		  
					
					if (_config.QPALMA_USE_SPLICE_SITES && !_config.NO_SPLICE_PREDICTIONS)

					{
						std::vector<int> positions;
						int num_acc=0, num_don=0 ;
						
						// find all acceptor splice sites downstream and donor splice site upstream
						if (ori==0)
						{
							if (_config.ACC_FILES.length()>0 && _config.DON_FILES.length()>0){
								num_acc = get_splicesite_positions(_config.ACC_FILES, "acc", chromosome, '+', midpoint, end, _config.QPALMA_USE_SPLICE_SITES_THRESH_ACC, true, positions) ;
								num_don = get_splicesite_positions(_config.DON_FILES, "don", chromosome, '+', start, midpoint, _config.QPALMA_USE_SPLICE_SITES_THRESH_DON, true, positions) ;
							}
							
							if (_config.SCORE_ANNOTATED_SPLICE_SITES){
								//fprintf(stdout,"Annotated splice site for pseudo chromosome seq\n");
								
								get_annotated_splice_positions(positions,annotatedjunctions,"acc",midpoint,end,chrN,'+');
								get_annotated_splice_positions(positions,annotatedjunctions,"don",start,midpoint,chrN,'+');
							}
						}
						
						else
						{
							if (_config.ACC_FILES.length()>0 && _config.DON_FILES.length()>0){
								num_don = get_splicesite_positions(_config.DON_FILES, "don", chromosome, '-', midpoint, end, _config.QPALMA_USE_SPLICE_SITES_THRESH_DON, true, positions) ;
								num_acc = get_splicesite_positions(_config.ACC_FILES, "acc", chromosome, '-', start, midpoint, _config.QPALMA_USE_SPLICE_SITES_THRESH_ACC, true, positions) ;
								}
							if (_config.SCORE_ANNOTATED_SPLICE_SITES){
								get_annotated_splice_positions(positions,annotatedjunctions,"acc",midpoint,end,chrN,'-');
								get_annotated_splice_positions(positions,annotatedjunctions,"don",start,midpoint,chrN,'-');
							}	
							
							
						}
						
						if (myverbosity>=3)
							fprintf(stdout, "splice site map: chr=%i, ori=%i: num_acc=%i\t num_don=%i\n", chrN, ori, num_acc, num_don) ;

						for (size_t jj=0; jj<positions.size(); jj++)
						{
							region_t *new_region = NULL;
							try {
								new_region = new region_t();
								new_region->read_map = new bool[read.length()];
							} catch (std::bad_alloc&) {
								fprintf(stderr, "[capture_hits] allocating memory for read_map failed\n");
								result.delete_regions();
								delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
								return -1;
							}
							new_region->erased=false ;
							new_region->start = positions[jj] - _config.INDEX_DEPTH ;
							new_region->end = positions[jj] + _config.INDEX_DEPTH ;
                            new_region->end = (new_region->end > chromosome.length())?chromosome.length():new_region->end;

							for (size_t ii = 0; ii < read.length(); ii++)
								new_region->read_map[ii] = false;
							
							//new_region->strand = (ori == 0) ? '+' : '-' ;
							//new_region->chromosome = chrN ;
							new_region->from_map = true ;
							
							new_region->read_pos = -111 ; // TODO
							new_region->hit_len = -111 ;  // TODO					

							regions[ori][chrN].push_back(new_region);
							assert(new_region->end > new_region->start) ;
			  
							//fprintf(stdout, "ss region %i-%i (seed %i-%i)\n", new_region->start, new_region->end, start, end) ;
			  
							added_map_regions++ ;
							added_map_total_len += new_region->end - new_region->start; 
						}
					}
		  
					if (_config.QPALMA_USE_MAP)
					{
		      
						int region_start = -1 ;
						int region_end = -1 ;
						for (int p=start; p<end; p+=Config::QPALMA_USE_MAP_WINDOW)//p+=1)
						{
							int val = genomemaps->CHR_MAP(chromosome,p) ;
			  
							if ((val & take_report_map)!=0 && end-p>Config::QPALMA_USE_MAP_WINDOW)//p<end-1)
							{
								if (region_start==-1){
									region_start=p;
								}
								region_end=p;
							}
							else
							{
								if (region_start!=-1)
								{ // region finished
				  
									region_t *new_region = NULL;
									try {
										new_region = new region_t();
										new_region->read_map = new bool[read.length()];
									} catch (std::bad_alloc&) {
										fprintf(stderr, "[capture_hits] allocating memory for read_map failed\n");
										result.delete_regions();
										delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
										return -1;
									}
									new_region->erased = false ;
									new_region->start = region_start - _config.INDEX_DEPTH ;
									if (end-p<=Config::QPALMA_USE_MAP_WINDOW)
										new_region->end = end-1 + _config.INDEX_DEPTH ;
									else
										new_region->end = region_end + _config.INDEX_DEPTH ;
									for (size_t ii = 0; ii < read.length(); ii++)
										new_region->read_map[ii] = false;

                                    new_region->end = (new_region->end > chromosome.length())?chromosome.length():new_region->end;

									//new_region->strand = (ori == 0) ? '+' : '-' ;
									//new_region->chromosome = chrN ;
									new_region->from_map = true ;
				  
									new_region->read_pos = -111 ; // TODO
									new_region->hit_len = -111 ;  // TODO					

									regions[ori][chrN].push_back(new_region);
									assert(new_region->end > new_region->start) ;

									added_map_regions++ ;
									added_map_total_len += new_region->end - new_region->start; 
								}
								region_start = -1 ;
								region_end = -1 ;
							}
						}						
					}
				}
			}
		}
		if (myverbosity>=1)
			fprintf(stdout, "# added %i regions from reporting map (total len=%i, average len=%i)\n", added_map_regions, added_map_total_len, added_map_total_len/(added_map_regions+1)) ;
    }
  
	// sort regions by starting position
	for (int ori = 0; ori < 2; ori++)
		for (int32_t chrN = 0; chrN < (int)regions[ori].size(); chrN++) 
		{
			if (regions[ori][chrN].size() == 0)
				continue;

			Chromosome const &chromosome = genome->chromosome(chrN);
			size_t nbr_regions=regions[ori][chrN].size();
			for (int i = 0; i < (int)nbr_regions; i++) 
				result.add_buffer_to_region(ori, chromosome, i, read.length());

			bool ret=sort_regions(regions[ori][chrN], sortorder_startpos) ;
			ret = ret && sort_regions(long_regions[ori][chrN], sortorder_startpos) ;
			if (!ret)
			{
				fprintf(stderr, "[capture_hits] allocating memory for sort_regions failed\n");
				result.delete_regions();
				delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
				return -1;
			}
		}

	// Surround the regions with a buffer of size _config.SPLICED_CLUSTER_TOLERANCE.
	// If any of them overlap before/after extension, merge them to
	// avoid duplicating subsequences.
	//
    for (int ori = 0; ori < 2; ori++){
		for (int32_t chrN = 0; chrN < (int)regions[ori].size(); chrN++)
		{
			size_t nbr_regions=regions[ori][chrN].size();
			if (nbr_regions == 0)
				continue;
	
			for (int nregion = 0; nregion < (int)nbr_regions - 1; nregion++) 
			{
				if (regions[ori][chrN][nregion]->erased)
					continue;
				size_t next = nregion + 1;
				while (next < nbr_regions && regions[ori][chrN][next]->erased)
					next++;
				if (next >= nbr_regions || regions[ori][chrN][next]->erased)
					continue;
	    
				if ((int32_t) ((regions[ori][chrN][nregion])->end) >= (int32_t) ((regions[ori][chrN][next])->start))
				{
					if (myverbosity>=2)
					{
						print_region(regions[ori][chrN][nregion], "region 1 for merge")  ;
						print_region(regions[ori][chrN][next], "region 2 for merge")  ;
					} ;

					// regions[ori] are overlapping or adjacent. Merge them into one.
					if ((regions[ori][chrN][nregion])->end < (regions[ori][chrN][next])->end)
						(regions[ori][chrN][nregion])->end = (regions[ori][chrN][next])->end;
					
					if (perform_extra_checks)
					{
						if (!(regions[ori][chrN][nregion]->start <= regions[ori][chrN][next]->start))
							fprintf(stdout, "warning: Failed: regions[ori][chrN][nregion]->start = %i <= %i = regions[ori][chrN][next]->start\n", regions[ori][chrN][nregion]->start, regions[ori][chrN][next]->start) ;
						//assert(regions[ori][chrN][nregion]->start <= regions[ori][chrN][next]->start) ;
					}
					
					// merge the two read_maps
					for (size_t i = 0; i < read.length(); i++)
						regions[ori][chrN][nregion]->read_map[i] = regions[ori][chrN][nregion]->read_map[i]
							|| regions[ori][chrN][next]->read_map[i];
		
					// Get rid of the extra region.
					regions[ori][chrN][next]->erased = true;
					delete[] regions[ori][chrN][next]->read_map; 
					regions[ori][chrN][next]->read_map = NULL;
		
					// consider this as a map-region, only if any both were generated from a map (important for downstream filtering)
					regions[ori][chrN][nregion]->from_map = regions[ori][chrN][nregion]->from_map && regions[ori][chrN][next]->from_map ;

					if (myverbosity>=2)
						print_region(regions[ori][chrN][nregion], "region merged")  ;
		
					// make sure this item is looked at again to merge it with the next item if necessary
					nregion-- ;
				}
			}
		}
	}
  
	//TODO: mt region1_time += clock() - start_time;
  
	
	/*printf("Clusters after consolidation: \n");
	  for (int ori=0; ori<1; ori++)
	  for (int32_t chrN = 0; chrN < regions[ori].size(); chrN++)
	  for (int32_t tmp = 0; tmp < regions[ori][chrN].size(); tmp++)
	  {
	  if (regions[ori][chrN][tmp]->erased)
	  continue ;
	  printf("\t%d - %d\n", regions[ori][chrN][tmp]->start, regions[ori][chrN][tmp]->end);
	  assert(regions[ori][chrN][tmp]->start<regions[ori][chrN][tmp]->end) ;
	  if (tmp>0)
	  assert(regions[ori][chrN][tmp-1]->end<regions[ori][chrN][tmp]->start) ;
	  }
	  fprintf(stdout, "_config.SPLICED_CLUSTER_TOLERANCE=%i\n", _config.SPLICED_CLUSTER_TOLERANCE) ;*/
  
	int num_merged = 0;
	for (int ori = 0; ori < 2; ori++)
		for (size_t chrN = 0; chrN < regions[ori].size(); chrN++)
			for (size_t tmp = 0; tmp < regions[ori][chrN].size(); tmp++)
			{
				if (perform_extra_checks)
					assert(regions[ori][chrN][tmp]->end>regions[ori][chrN][tmp]->start) ;
				if (!regions[ori][chrN][tmp]->erased)
					num_merged++;
			}
	if (myverbosity >= 2)
		fprintf(stdout, "# Merged hit list has %i items\n", num_merged);
	TIME_CODE(_stats.qpalma_region_align_time += clock() - start_time;)
	return 0;
}


int QPalma::find_regions_for_long_regions(const region_t *long_region, const std::vector<region_t*> & regions, std::vector<region_t*> & current_regions, const Read & read, int myverbosity) const
{
	int num_read_map=0 ;

	int midpoint = (long_region->start+long_region->end)/2 ;

	// look at the regions downstream. only keep those which lead to an alignment string shorter than _config.QPALMA_USE_MAP_MAX_SIZE
	int right_covered_region=0 ;
	int right_removed = 0 ;
	int left_covered_region=0 ;
	int left_removed=0 ;
	int added_map_regions=0 ;
	int added_map_total_len=0 ;

	// initialize read_map
	bool read_map[read.length()] ;
	for (size_t i = 0; i < read.length(); i++)
		read_map[i] = false ;
	
	for (size_t nregion=0; nregion<regions.size(); nregion++)
	{
		if (regions[nregion]->erased)
			continue;

		if (regions[nregion]->start > midpoint)
		{
			if (regions[nregion]->start-midpoint > _config.SPLICED_LONGEST_INTRON_LENGTH*_config.SPLICED_MAX_INTRONS + (int)read.length())
				continue ;

			if (right_covered_region>_config.QPALMA_USE_MAP_MAX_SIZE)
			{
				if (!regions[nregion]->from_map)
				{
					if (perform_extra_checks)
						assert(regions[nregion]->start-midpoint <= _config.SPLICED_LONGEST_INTRON_LENGTH*_config.SPLICED_MAX_INTRONS + (int)read.length());
					
					added_map_regions++ ;
					added_map_total_len += regions[nregion]->end - regions[nregion]->start;
					current_regions.push_back(regions[nregion]) ;
					if (myverbosity>=2)
					{
						char c=' ' ;
						if (midpoint>=regions[nregion]->start && midpoint<=regions[nregion]->end)
							c='*' ;
						
						fprintf(stdout, "[find_regions_for_long_regions] region %i: %i - %i (from_map=%i)\t%c\n", (int)nregion, regions[nregion]->start, regions[nregion]->end, 
								regions[nregion]->from_map, c) ;
					}

					for (size_t i = 0; i < read.length(); i++)
						if (regions[nregion]->read_map[i])
							read_map[i] = true ;
				}
				else
					right_removed++ ;
			}
			else
			{
				added_map_regions++ ;
				added_map_total_len += regions[nregion]->end - regions[nregion]->start;
				current_regions.push_back(regions[nregion]) ;
				if (myverbosity>=2)
				{
					char c=' ' ;
					if (midpoint>=regions[nregion]->start && midpoint<=regions[nregion]->end)
						c='*' ;
					
					fprintf(stdout, "[find_regions_for_long_regions] region %i: %i - %i (from_map=%i)\t%c\n", (int)nregion, regions[nregion]->start, regions[nregion]->end, 
							regions[nregion]->from_map, c) ;
				}
				for (size_t i = 0; i < read.length(); i++)
					if (regions[nregion]->read_map[i])
						read_map[i] = true ;
			}
			if (regions[nregion]->from_map)
				right_covered_region+= regions[nregion]->end - regions[nregion]->start + 2*_config.SPLICED_CLUSTER_TOLERANCE ;
		}
		if (regions[nregion]->start <= midpoint)
		{
			if (midpoint - regions[nregion]->end > _config.SPLICED_LONGEST_INTRON_LENGTH*_config.SPLICED_MAX_INTRONS + (int)read.length()) 
				continue ;
			
			if (left_covered_region > _config.QPALMA_USE_MAP_MAX_SIZE)
			{
				if (!regions[nregion]->from_map)
				{
					if (perform_extra_checks)
						assert(midpoint - regions[nregion]->end <= _config.SPLICED_LONGEST_INTRON_LENGTH*_config.SPLICED_MAX_INTRONS + (int)read.length()) ;
					added_map_regions++ ;
					added_map_total_len += regions[nregion]->end - regions[nregion]->start;
					current_regions.push_back(regions[nregion]) ;
					if (myverbosity>=2)
					{
						char c=' ' ;
						if (midpoint>=regions[nregion]->start && midpoint<=regions[nregion]->end)
							c='*' ;
						
						fprintf(stdout, "[find_regions_for_long_regions] region %i: %i - %i (from_map=%i)\t%c\n", (int)nregion, regions[nregion]->start, regions[nregion]->end, 
								regions[nregion]->from_map, c) ;
					}
					for (size_t i = 0; i < read.length(); i++)
						if (regions[nregion]->read_map[i])
							read_map[i] = true ;
				}
				else
					left_removed++ ;
			}
			else
			{
				added_map_regions++ ;
				added_map_total_len += regions[nregion]->end - regions[nregion]->start;
				current_regions.push_back(regions[nregion]) ;
				if (myverbosity>=2)
				{
					char c=' ' ;
					if (midpoint>=regions[nregion]->start && midpoint<=regions[nregion]->end)
						c='*' ;
					
					fprintf(stdout, "[find_regions_for_long_regions] region %i: %i - %i (from_map=%i)\t%c\n", (int)nregion, regions[nregion]->start, regions[nregion]->end, 
							regions[nregion]->from_map, c) ;
				}
				for (size_t i = 0; i < read.length(); i++)
					if (regions[nregion]->read_map[i])
						read_map[i] = true ;
			}
			if (regions[nregion]->from_map)
				left_covered_region += (regions[nregion]->end - regions[nregion]->start) + 2*_config.SPLICED_CLUSTER_TOLERANCE ;
		}

	}

	if (myverbosity>=2)
		fprintf(stdout, "cover: %i, %i\tregions removed: %i, %i\n", left_covered_region, right_covered_region, left_removed, right_removed) ;

	if (myverbosity>=1)
		fprintf(stdout, "# took  %i regions (total len=%i, average len=%i)\n", added_map_regions, added_map_total_len, added_map_total_len/(added_map_regions+1)) ;

	for (size_t i = 0; i < read.length(); i++)
		if (read_map[i])
			num_read_map++;

	return num_read_map ;
}


int QPalma::capture_hits_2(Hits &hits, Result &result, bool non_consensus_search, JunctionMap &annotatedjunctions, VariantMap &variants) const 
{
	// Dominik marker
	// std::vector<std::vector<region_t *> > long_regions_other[2] ;
	// if (paired_read_processing)
	//      	long_regions_other = get_long_regions_from_other_end(long_regions, semaphore_whatever_is_necessary) ;

	int myverbosity=verbosity ;
	
	if (myverbosity!=verbose_read_level && std::string(result._read.id())==verbose_read_id)
		myverbosity=verbose_read_level ;

	Read const &read(hits.getRead());
	TIME_CODE(clock_t start_time = clock();)

	std::vector<std::vector<region_t *> > *regions = result.regions;
	std::vector<std::vector<region_t *> > *long_regions = result.long_regions;

	// Construct a pseudo genomic sequence from the regions by stringing
	// them together. A new sequence is started when two adjacent regions
	// are farther apart than the maximum assumed intron length.
	// Also, we build up a map of positions on the generated strings to
	// their corresponding positions in the genomic sequence. This is
	// so we can map the results of the alignment back onto the actual
	// genomic sequence later.

	int num_alignments_reported = 0 ;
  
	std::string read_seq[2];
	read_seq[0] = std::string(read.data(), read.length());
	read_seq[1] = reverse(complement(read_seq[0]));
	if (myverbosity >= 3) {
		fprintf(stdout, "# read[0]: %s\n", read_seq[0].c_str());
		fprintf(stdout, "# read[1]: %s\n", read_seq[1].c_str());
	}
  
	std::string read_quality[2];
	read_quality[0] = std::string(read.quality(0), read.length());
	read_quality[1] = reverse(read_quality[0]);
	if (myverbosity >= 3)
		fprintf(stdout, "# readqual[0]: %s\n", read_quality[0].c_str());
  
	for (int ori = 0; ori < 2; ori++)
	{

		hits.topAlignments().init_top_alignment_index();
	
		for (size_t chrN = 0; chrN < regions[ori].size(); chrN++) 
		{
			Chromosome const &chr = genome->chromosome(chrN);
			std::vector<region_t *> current_regions;
			current_regions.clear();
			if (regions[ori][chrN].size() == 0)
				continue;

			if (myverbosity>=1)
				fprintf(stdout, "start performing alignment for %lu long regions (%s)\n", long_regions[ori][chrN].size(), read.id()) ;

			for (int i=0; i<(int)long_regions[ori][chrN].size(); i++)
			{
				current_regions.clear();
				int num_read_map=find_regions_for_long_regions(long_regions[ori][chrN][i], regions[ori][chrN], current_regions, read, myverbosity) ; 
				
				if (myverbosity >= 1)
					fprintf(stdout, "# region list covering %i bases\n", num_read_map);
				if (num_read_map >= _config.SPLICED_HIT_MIN_LENGTH_COMB)
				{
					int ret=capture_hits_3(chr, hits, result, read, read_seq, read_quality, ori, current_regions, long_regions[ori][chrN][i], 
										   non_consensus_search, annotatedjunctions, variants, myverbosity) ;
					if (ret<0)
					{
						if (myverbosity>=1)
							fprintf(stdout, "Warning: alignment of read %s failed\n", read.id()) ;
						result.delete_regions();
						delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
						return ret ;
					}

					if (hits.topAlignments().stop_aligning())
					  {
						if (myverbosity>=1)
							fprintf(stdout, "Warning: stopped alignments of read %s\n", read.id()) ;
						result.delete_regions();
						delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements

						return 0 ;
					  }

					num_alignments_reported+=ret ;
					hits.topAlignments().update_top_alignment_index();
				}
				else
					if (myverbosity >= 1)
						fprintf(stdout,	"# dropped region list covering only %i bases\n", num_read_map);
			}
		}
	}
	
	if (!_config.MAP_JUNCTIONS){
		result.delete_regions();
		delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
	}

	TIME_CODE(if (myverbosity >= 1 || ((clock()-_stats.qpalma_last_timing_report)/CLOCKS_PER_SEC>=10))
    {
		_stats.qpalma_last_timing_report = clock() ;
		_stats.qpalma_timing(((float) clock() - start_time) / CLOCKS_PER_SEC);
    })

	return num_alignments_reported ;

}

int QPalma::capture_hits_3(Chromosome const &chr, Hits &hits, Result &result, Read const &read, 
						   std::string read_seq[2], std::string read_quality[2], int ori, 
						   std::vector<region_t*> & regions, region_t* long_region,
						   bool non_consensus_search, JunctionMap &annotatedjunctions, VariantMap &variants,
						   int myverbosity) const 
{
	int num_alignments_reported=0 ;

	if (myverbosity!=verbose_read_level && std::string(result._read.id())==verbose_read_id)
		myverbosity=verbose_read_level ;

	//Iterate all long regions and start alignment
	int transcription_direction = get_transcription_direction(_config.STRAND, ori) ;

	int hit_read_position = get_first_read_map(read, long_region->read_map);
	if (perform_extra_checks)
		assert( hit_read_position == long_region->read_pos ) ;
	
	int hit_len= long_region->end-long_region->start;
	if (perform_extra_checks)
	{
		int this_num_read_map = 0;
		for (size_t i = 0; i < read.length(); i++)
			if (long_region->read_map[i])
				this_num_read_map++;
		assert(this_num_read_map==hit_len) ;
	}
	
	if(ori==1)
		hit_read_position = read.length()-hit_len-hit_read_position+1;
	if (perform_extra_checks)
		assert (hit_read_position>=0 && hit_len >0);
	
	std::vector<region_t*> & current_regions = regions ;
	std::string current_seq ;
	std::vector<int> current_positions;
	
	for (int nregion=0; nregion < (int)current_regions.size(); nregion++)
	{
		std::string str ;
		int ret = get_string_from_region(chr, regions[nregion], str);
		if (ret < 0)
			return ret ;
		current_seq+=str ;
		
		for (int p = 0; p < regions[nregion]->end - regions[nregion]->start; p++) 
		{
			//Fake positions for N blocks
			if (p<3 || p>= regions[nregion]->end - regions[nregion]->start -3)
				current_positions.push_back(-2);
			else
				current_positions.push_back(regions[nregion]->start + p);
		}
	}
	
	if (myverbosity>=1)
	{
		fprintf(stdout,"read id %s curr len %lu\n", read.id(), current_positions.size());
		fprintf(stdout,	"# Starting point for alignments: read %i, dna %i, len %i\n",hit_read_position, long_region->start, hit_len);					  
		fprintf(stdout,	"# Number of current regions %lu\n", current_regions.size());
	}
	
	if (transcription_direction >=0)
	{
		//fprintf(stdout,	"1)hit read position %i\n",hit_read_position);					  
		int ret = perform_alignment_starter_variant(result, hits, read_seq[ori], read_quality[ori], current_seq, current_regions,
													current_positions, chr, '+', ori, hit_read_position,
													long_region->start, hit_len, non_consensus_search, num_alignments_reported, false,
													annotatedjunctions, variants, myverbosity);
		if (ret < 0)
			return ret;
	}
	
	if (transcription_direction <=0)
	{
		//fprintf(stdout,	"# Starting point for alignments: read %i, dna %i, len %i\n",read.length()-(hit_read_position+hit_len),
		//      corres_long_regions[lr]->end, hit_len);					  
		//fprintf(stdout,	"# Number of current regions %i\n",(int)current_regions.size());					  
		//fprintf(stdout,	"2)hit read position %i\n", read.length()-(hit_read_position+hit_len));					  
		
		int ret = perform_alignment_starter_variant(result, hits, read_seq[1 - ori], read_quality[1 - ori], current_seq, current_regions, 
													current_positions, chr, '-', 1-ori, read.length()-(hit_read_position+hit_len),
													long_region->end-1, hit_len, non_consensus_search, num_alignments_reported, false,
													annotatedjunctions, variants, myverbosity); // end nucleotide in dna not included
		if (ret < 0)
			return ret;
	}

	return num_alignments_reported ;
}


template<int myverbosity, bool discover_variants, bool remapping>
void *perform_alignment_wrapper(QPalma::perform_alignment_t *data)
{
	try
	{
		assert(data->qpalma!=NULL) ;
		const bool do_timing=false ;
		
		data->ret = data->qpalma->perform_alignment<myverbosity, discover_variants, remapping, do_timing>
			(*data->result, *data->readMappings, data->read_string, data->read_quality, data->dna,
			 data->current_regions, data->positions, *data->contig_idx,
			 data->strand, data->ori, data->hit_read,
			 data->hit_dna,data->hit_length, data->non_consensus_search, data->aln, 
			 *data->annotatedjunctions, *data->variants, *data->variant_list, *data->variant_positions) ;
	}
	catch (std::bad_alloc&)
	{
		fprintf(stderr, "Warning: Alignment of read %s failed (out of memory exception caught)\n", data->result->_read.id()) ;
		data->ret = -1 ; // failed
	}

	return data ;
}

template<int myverbosity, bool discover_variants>
void QPalma::perform_alignment_wrapper2(QPalma::perform_alignment_t *data) const
{
	if (data->remapping)
		perform_alignment_wrapper<myverbosity,discover_variants,true>(data);
	else
		perform_alignment_wrapper<myverbosity,discover_variants,false>(data);
}


void QPalma::perform_alignment_wrapper1(QPalma::perform_alignment_t *data, int myverbosity) const
{
	bool discover_variants = _config.DISCOVER_VARIANTS ;
	if (data->remapping)
		discover_variants = false ;

	if (myverbosity!=verbose_read_level && std::string(data->result->_read.id())==verbose_read_id)
		myverbosity=verbose_read_level ;
	
	switch(myverbosity)
	{
	case 0:
		if (discover_variants)
			perform_alignment_wrapper2<0,true>(data);
		else
			perform_alignment_wrapper2<0,false>(data);
		break ;
	case 1:
		if (discover_variants)
			perform_alignment_wrapper2<1,true>(data);
		else
			perform_alignment_wrapper2<1,false>(data);
		break ;
	case 2:
		if (discover_variants)
			perform_alignment_wrapper2<2,true>(data);
		else
			perform_alignment_wrapper2<2,false>(data);
		break ;
	default:
		if (discover_variants)
			perform_alignment_wrapper2<3,true>(data);
		else
			perform_alignment_wrapper2<3,false>(data);
		break ;
	}
}

int QPalma::junctions_remapping(Hits &hits, Result &result, JunctionMap &junctionmap, int nb_spliced_alignments, JunctionMap &annotatedjunctions, VariantMap & variants) const 
{

	std::vector<std::vector<region_t *> > *long_regions = result.long_regions;
	// if (nb_spliced_alignments > 0){
	// 	result.delete_regions();
	// 	delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
	// 	return 0;
	// }
	int junction_tol = 10 ; // the hit may overlap by this length with the junction
	
	int myverbosity=verbosity ;
	if (myverbosity!=verbose_read_level && std::string(result._read.id())==verbose_read_id)
		myverbosity=verbose_read_level ;
	if (myverbosity>=2)
		fprintf(stdout, "remapping read %s\n", result._read.id()) ;

	Read const &read(hits.getRead());

	// Current regions correspond to regions around the junction
	std::vector<region_t *> current_regions;
	// Positions are real positions corresponding to current regions
	std::vector<int> current_positions;
	
	int num_alignments_reported = 0 ;
	
 	std::string read_seq[2];
	read_seq[0] = std::string(read.data(), read.length());
	read_seq[1] = reverse(complement(read_seq[0]));
	if (myverbosity >= 3) {
		fprintf(stdout, "# read[0]: %s\n", read_seq[0].c_str());
		fprintf(stdout, "# read[1]: %s\n", read_seq[1].c_str());
	}
  
	std::string read_quality[2];
	read_quality[0] = std::string(read.quality(0), read.length());
	read_quality[1] = reverse(read_quality[0]);
	if (myverbosity >= 3)
		fprintf(stdout, "# readqual[0]: %s\n", read_quality[0].c_str());
  

	for (int ori = 0; ori < 2; ori++){

		int transcription_direction = get_transcription_direction(_config.STRAND,ori);
						
		hits.topAlignments().init_top_alignment_index();
	
		for (size_t chrN = 0; chrN < long_regions[ori].size(); chrN++) 
		{
			std::vector<region_t *> & this_long_regions = long_regions[ori][chrN] ;

			if (this_long_regions.size() == 0)
				continue;

			current_regions.clear();
			current_positions.clear();
			
			Chromosome const &chr = genome->chromosome(chrN);			

			for (size_t nregion = 0; nregion < this_long_regions.size(); nregion++)
			{
				if (this_long_regions[nregion]==NULL)
					continue ;
				int rstart_ = this_long_regions[nregion]->start;
				int rend_ = this_long_regions[nregion]->end; // half-open
				
                int available_left = read.length();
                int available_right = read.length();

                // reset junction tolerance to at most half the long region length
                junction_tol = 10 < ((rend_ - rstart_) / 2) ? 10 : (rend_ - rstart_) / 2;

				if (myverbosity>=3) // was4
					fprintf(stdout, "long region: %i-%i (%ld,%i)\n", rstart_, rend_, chrN, ori);

                std::vector<std::vector<Junction> > combinations;
                std::vector<unsigned int> cum_spans_right;
                std::vector<unsigned int> cum_spans_left;

				// find a lower bound on the indexes with binary search (first hit that fullfils < x)
				junctionmap.lock() ;
				std::deque<Junction>::iterator it_s = my_lower_bound_by_start(junctionmap.junctionlist_by_start[chrN].begin(), junctionmap.junctionlist_by_start[chrN].end(), rend_) ;
				std::deque<Junction>::iterator it_e = my_lower_bound_by_end(junctionmap.junctionlist_by_end[chrN].begin(), junctionmap.junctionlist_by_end[chrN].end(), rstart_ - 1);

                // iterate over left junctions
                // junctions have closed intervals
                while (true) {
                    
                    if (it_e == junctionmap.junctionlist_by_end[chrN].end())
                        break;

                    //fprintf(stdout, "[left juncs]: %i-%i\n", it_e->start, it_e->end);
                    if (it_e->end >= rstart_ + junction_tol) {
                        if (it_e == junctionmap.junctionlist_by_end[chrN].begin())
                            break;
                        it_e--;
                        continue;
                    }

                    //Continue only if the strand with the splice junction is consistent with the transcription direction
                    if ((transcription_direction==1 && it_e->strand=='-') || (transcription_direction==-1 && it_e->strand=='+') )		
                    {
                        if (myverbosity>=3)
                            fprintf(stdout, "junction_remapping: skipping junction at %ld:%i-%i on strand %c with transcription direction %i\n", chrN, it_e->start, it_e->end, it_e->strand, transcription_direction) ;
                        if (it_e == junctionmap.junctionlist_by_end[chrN].begin())
                            break;
                        it_e--;
                        continue ;
                    }
                    
                    // augment existing combinations
                    size_t num_combs = combinations.size();
                    for (size_t i = 0; i < num_combs; i++) {
                        if (combinations.size() <= _config.JUNCTION_MAX_NUM_ALIGNMENTS && (it_e->end < combinations.at(i).back().start - 1) && 
                           (cum_spans_left.at(i) + (combinations.at(i).back().start - it_e->end - 1) < available_left) &&
                           combinations.at(i).size() < _config.JUNCTION_REMAP_MAX_INTRONS) { // less than available as we need at least one nt before the junction
                            std::vector<Junction> tmp = combinations.at(i);
                            tmp.push_back(*it_e);
                            combinations.push_back(tmp);
                            cum_spans_right.push_back(0);
                            cum_spans_left.push_back(cum_spans_left.at(i) + combinations.at(i).back().start - it_e->end - 1);
                        }
                    }
                    // try junction alone
                    if (combinations.size() <= _config.JUNCTION_MAX_NUM_ALIGNMENTS && (it_e->end > (rstart_ - available_left)) && (it_e->end < (rstart_ + junction_tol))) {
                        std::vector<Junction> tmp;
                        tmp.push_back(*it_e);
                        combinations.push_back(tmp);
                        cum_spans_right.push_back(0);
                        if (it_e->end < rstart_)
                            cum_spans_left.push_back(rstart_ - it_e->end - 1);
                        else
                            cum_spans_left.push_back(0);
                    } else {
                        break;
                    }
                    if (it_e > junctionmap.junctionlist_by_end[chrN].begin())
                        it_e--;
                    else
                        break;
                } 

                // try additional right junctions
                if (combinations.size() <= 2*_config.JUNCTION_MAX_NUM_ALIGNMENTS)
                {
                    // iterate over right junctions
                    while (it_s != junctionmap.junctionlist_by_start[chrN].end()) {
                        // fprintf(stdout, "[right juncs]: %i-%i\n", it_s->start, it_s->end);
                        
                        if (it_s->start < (rend_ - junction_tol)) {
                            it_s++;
                            continue;
                        }

                        //Continue only if the strand with the splice junction is consistent with the transcription direction
                        if ((transcription_direction==1 && it_s->strand=='-') || (transcription_direction==-1 && it_s->strand=='+') )		
                        {
                            if (myverbosity>=3)
                                fprintf(stdout, "junction_remapping: skipping junction at %ld:%i-%i on strand %c with transcription direction %i\n", chrN, it_s->start, it_s->end, it_s->strand, transcription_direction) ;
                            it_s++;
                            continue ;
                        }

                        // augment existing combinations
                        size_t num_combs = combinations.size();
                        for (size_t i = 0; i < num_combs; i++) {
                            if ((combinations.size() <= 2*_config.JUNCTION_MAX_NUM_ALIGNMENTS) && (it_s->start > combinations.at(i).back().end + 1) &&  
                                (cum_spans_right.at(i) + (it_s->start - combinations.at(i).back().end - 1) < available_right) &&
                                combinations.at(i).size() < _config.JUNCTION_REMAP_MAX_INTRONS) {
                                std::vector<Junction> tmp = combinations.at(i);
                                tmp.push_back(*it_s);
                                combinations.push_back(tmp);
                                cum_spans_left.push_back(0);
                                cum_spans_right.push_back(cum_spans_right.at(i) + it_s->start - combinations.at(i).back().end - 1);
                            }
                        }
                        // try junction alone
                        if (combinations.size() <= 2*_config.JUNCTION_MAX_NUM_ALIGNMENTS && (it_s->start < (available_right + rend_)) && (it_s->start >= (rend_ - junction_tol))) {
                            std::vector<Junction> tmp;
                            tmp.push_back(*it_s);
                            combinations.push_back(tmp);
                            cum_spans_left.push_back(0);
                            if (it_s->start > rend_)
                                cum_spans_right.push_back(it_s->start - rend_);
                            else
                                cum_spans_right.push_back(0);
                        } else {
                            break;
                        }
                        it_s++;
                    }
                }
				junctionmap.unlock() ;

                if (myverbosity >=3 ) {
                    fprintf(stdout, "found %i possible junction combinations\n", combinations.size());
                }

                // iterate over found combinations
                for (size_t c_idx = 0; c_idx < combinations.size(); c_idx++) 
                {
                        if (myverbosity >= 3) {
                            fprintf(stdout, " %i junctions in current combination\n", combinations.at(c_idx).size());
                            for (size_t j =  0; j < combinations.at(c_idx).size(); j++)
                                fprintf(stdout, "\t(%i - %i) [ %c %c ... %c %c ]\n", combinations.at(c_idx).at(j).start, combinations.at(c_idx).at(j).end, 
                                                                                   chr[combinations.at(c_idx).at(j).start], chr[combinations.at(c_idx).at(j).start + 1],
                                                                                   chr[combinations.at(c_idx).at(j).end - 1], chr[combinations.at(c_idx).at(j).end]);
                        }
					 	/*if (myverbosity>4)
						{
							fprintf(stdout, "try to align\n") ;
							fprintf(stdout, "    junction: (%i-%i)-(%i-%i) (%ld,%c)\n",int1_start,int1_end,int2_start,int2_end,chrN,strand);
							fprintf(stdout, "    junction %c%c-%c%c\n",chr[(*it).start],chr[(*it).start+1],chr[(*it).end-1],chr[(*it).end]);
							fprintf(stdout, "long region: %i-%i (%ld,%i)\n", rstart, rend, chrN, ori);
						}*/

                        int rstart = rstart_;
                        int rend = rend_;

                        char strand = combinations.front().front().strand ;

                        std::string current_seq;
                        current_seq.assign("") ;
                        // go through junction list to find leftmost junction (left junction with highest index)
                        size_t leftmost = combinations.at(c_idx).size();
                        for (size_t j = 0; j < combinations.at(c_idx).size(); j++)
                        {
                            // fprintf(stdout, "curr comb (%i) end = %i; rstart=%i; junction_tol=%i\n", j, combinations.at(c_idx).at(j).end, rstart, junction_tol);
                            if (combinations.at(c_idx).at(j).end < (rstart + junction_tol)) 
                                leftmost = j;
                            else
                                break;
                        }
                        if (myverbosity >= 3) { 
                            fprintf(stdout, "try to align\n") ;
                            fprintf(stdout, "\tcomb regions left: ") ;
                        }
                        // handle left junctions
                        if (leftmost < combinations.at(c_idx).size())
                        {
                            for (int j = leftmost; j >= 0; j--)
                            {
                                //Create regions from junction
                                region_t* new_region1 = new region_t; 
                                bool read_map1[read.length()] ;
                                new_region1->read_map = read_map1 ;
                                if (j == leftmost)
                                    new_region1->start = combinations.at(c_idx).at(j).start - (available_left - cum_spans_left.at(c_idx));
                                else
                                    new_region1->start = combinations.at(c_idx).at(j+1).end + 1;
                                new_region1->end = combinations.at(c_idx).at(j).start; // interval open
                                new_region1->from_map = true ;
                                for (size_t ii = 0; ii < read.length(); ii++)
                                    new_region1->read_map[ii] = false;
                                new_region1->read_pos = -111 ; // TODO
                                new_region1->hit_len = -111 ;  // TODO					
                                // check chromsome boundaries 
                                new_region1->start = new_region1->start < 0 ? 0 : new_region1->start;
                                new_region1->end = new_region1->end > chr.length() ? chr.length() : new_region1->end;

                                current_regions.push_back(new_region1);
            
                                for (int p = new_region1->start; p < new_region1->end; p++) 
                                {
                                    current_positions.push_back(p);
                                    current_seq.push_back(chr[p]) ;							
                                }
                                if (myverbosity >= 3)
                                    fprintf(stdout, "(%i-%i)%c%c-%c%c", new_region1->start, new_region1->end, chr[combinations.at(c_idx).at(j).start], chr[combinations.at(c_idx).at(j).start + 1],
                                            chr[combinations.at(c_idx).at(j).end - 1], chr[combinations.at(c_idx).at(j).end]);
                            }
                            // adapt rstart according to junction_tol
                            if (rstart <= combinations.at(c_idx).front().end)
                            {
                                rstart = combinations.at(c_idx).front().end + 2;
                            } else {
                                region_t* new_region1 = new region_t ; 
                                bool read_map1[read.length()] ;
                                new_region1->read_map = read_map1 ;
                                new_region1->start = combinations.at(c_idx).front().end + 1;
                                new_region1->end = rstart; // interval open
                                for (size_t ii = 0; ii < read.length(); ii++)
                                    new_region1->read_map[ii] = false;
                                new_region1->read_pos = -111 ; // TODO
                                new_region1->hit_len = -111 ;  // TODO					
                                // check chromsome boundaries 
                                new_region1->start = new_region1->start < 0 ? 0 : new_region1->start;
                                new_region1->end = new_region1->end > chr.length() ? chr.length() : new_region1->end;

                                current_regions.push_back(new_region1);
            
                                for (int p = new_region1->start; p < new_region1->end; p++) 
                                {
                                    current_positions.push_back(p);
                                    current_seq.push_back(chr[p]) ;							
                                }
                                if (myverbosity >= 3) 
                                    fprintf(stdout, "(%i-%i) (%ld,%c)\n", new_region1->start, new_region1->end, chrN, strand);
                            }
                        } else { 
                            // no junction on the left, add whole left part as region
                            region_t* new_region1 = new region_t ; 
                            bool read_map1[read.length()] ;
                            new_region1->read_map = read_map1 ;
                            new_region1->start = rstart - available_left;
                            new_region1->end = rstart;
                            new_region1->from_map = true ;
                            for (size_t ii = 0; ii < read.length(); ii++)
                                new_region1->read_map[ii] = false;
                            new_region1->read_pos = -111 ; // TODO
                            new_region1->hit_len = -111 ;  // TODO					
                            // check chromsome boundaries 
                            new_region1->start = new_region1->start < 0 ? 0 : new_region1->start;
                            new_region1->end = new_region1->end > chr.length() ? chr.length() : new_region1->end;
        
                            current_regions.push_back(new_region1);

                            for (int p = new_region1->start; p < new_region1->end; p++) 
                            {
                                current_positions.push_back(p);
                                current_seq.push_back(chr[p]) ;							
                            }
                            if (myverbosity >= 3)
                                fprintf(stdout, "(%i-%i) (%ld,%c)\n", new_region1->start, new_region1->end, chrN, strand);
                        }

                        // adapt rend, if necessary
                        if ((leftmost < combinations.at(c_idx).size() - 1) || (leftmost == combinations.at(c_idx).size()))
                        {
                            size_t start_idx = leftmost + 1;
                            if (leftmost == combinations.at(c_idx).size())
                                start_idx = 0;
                            if (combinations.at(c_idx).at(start_idx).start < rend)
                                rend = combinations.at(c_idx).at(start_idx).start;
                        }

                        // add long region as region
                        region_t* new_region1 = new region_t;
                        bool read_map_hit[read.length()] ;
                        new_region1->read_map = read_map_hit ;
                        new_region1->start = rstart;
                        new_region1->end = rend; // interval open
                        for (size_t ii = 0; ii < read.length(); ii++)
                            new_region1->read_map[ii] = false;
                        new_region1->read_pos = -111 ; // TODO
                        new_region1->hit_len = -111 ;  // TODO					
                        // check chromsome boundaries 
                        new_region1->start = new_region1->start < 0 ? 0 : new_region1->start;
                        new_region1->end = new_region1->end > chr.length() ? chr.length() : new_region1->end;

                        current_regions.push_back(new_region1);
    
                        for (int p = new_region1->start; p < new_region1->end; p++) 
                        {
                            current_positions.push_back(p);
                            current_seq.push_back(chr[p]) ;							
                        }
                        if (myverbosity >= 3) {
                            fprintf(stdout, "\tlong region: %i-%i (chr:%ld,ori:%i)\n", rstart, rend, chrN, ori);
                            fprintf(stdout, "\tcomb regions right: ") ;
                        }
                        
                        // there is a at least one junction on the right
                        if ((leftmost < combinations.at(c_idx).size() - 1) || (leftmost == combinations.at(c_idx).size()))
                        {
                            size_t start_idx = leftmost + 1;
                            if (leftmost == combinations.at(c_idx).size())
                                start_idx = 0;

                            // iterate over remaining right junctions of this combination
                            for (size_t j = start_idx; j < combinations.at(c_idx).size(); j++) 
                            {
                                if (combinations.at(c_idx).at(j).start == rend)
                                    continue;
                                //Create regions from junction
                                region_t* new_region1 = new region_t; 
                                bool read_map1[read.length()] ;
                                new_region1->read_map = read_map1 ;
                                if (j == start_idx)
                                    new_region1->start = rend ; // interval open
                                else
                                    new_region1->start = combinations.at(c_idx).at(j-1).end + 1;
                                new_region1->end = combinations.at(c_idx).at(j).start; // interval open
                                new_region1->from_map = true ;
                                for (size_t ii = 0; ii < read.length(); ii++)
                                    new_region1->read_map[ii] = false;
                                new_region1->read_pos = -111 ; // TODO
                                new_region1->hit_len = -111 ;  // TODO					
                                // check chromsome boundaries 
                                new_region1->start = new_region1->start < 0 ? 0 : new_region1->start;
                                new_region1->end = new_region1->end > chr.length() ? chr.length() : new_region1->end;

                                current_regions.push_back(new_region1);
            
                                for (int p = new_region1->start; p < new_region1->end; p++) 
                                {
                                    current_positions.push_back(p);
                                    current_seq.push_back(chr[p]) ;							
                                }
                                if (myverbosity >= 3)
                                    fprintf(stdout, "(%i-%i)%c%c-%c%c", new_region1->start, new_region1->end, chr[combinations.at(c_idx).at(j).start], chr[combinations.at(c_idx).at(j).start + 1],
                                            chr[combinations.at(c_idx).at(j).end - 1], chr[combinations.at(c_idx).at(j).end]);
                            }
                            // add last region for remaining read
                            region_t* new_region1 = new region_t ; 
                            bool read_map1[read.length()] ;
                            new_region1->read_map = read_map1 ;
                            new_region1->start = combinations.at(c_idx).back().end + 1;
                            new_region1->end = combinations.at(c_idx).back().end + available_left - cum_spans_left.at(c_idx) + 1; // interval open
                            new_region1->from_map = true ;
                            for (size_t ii = 0; ii < read.length(); ii++)
                                new_region1->read_map[ii] = false;
                            new_region1->read_pos = -111 ; // TODO
                            new_region1->hit_len = -111 ;  // TODO					
                            // check chromsome boundaries 
                            new_region1->start = new_region1->start < 0 ? 0 : new_region1->start;
                            new_region1->end = new_region1->end > chr.length() ? chr.length() : new_region1->end;
        
                            current_regions.push_back(new_region1);

                            for (int p = new_region1->start; p < new_region1->end; p++) 
                            {
                                current_positions.push_back(p);
                                current_seq.push_back(chr[p]) ;							
                            }
                            if (myverbosity >= 3)
                                fprintf(stdout, "(%i-%i) (%ld,%c)\n", new_region1->start, new_region1->end, chrN, strand);
                        } else {
                            // add whole right side without junction
                            region_t* new_region1 = new region_t ; 
                            bool read_map1[read.length()] ;
                            new_region1->read_map = read_map1 ;
                            new_region1->start = rend;
                            new_region1->end = rend + available_right;
                            new_region1->from_map = true ;
                            for (size_t ii = 0; ii < read.length(); ii++)
                                new_region1->read_map[ii] = false;
                            new_region1->read_pos = -111 ; // TODO
                            new_region1->hit_len = -111 ;  // TODO					
                            // check chromsome boundaries 
                            new_region1->start = new_region1->start < 0 ? 0 : new_region1->start;
                            new_region1->end = new_region1->end > chr.length() ? chr.length() : new_region1->end;
        
                            current_regions.push_back(new_region1);

                            for (int p = new_region1->start; p < new_region1->end; p++) 
                            {
                                current_positions.push_back(p);
                                current_seq.push_back(chr[p]) ;							
                            }
                            if (myverbosity >= 3)
                                fprintf(stdout, "(%i-%i) (%ld,%c)\n", new_region1->start, new_region1->end, chrN, strand);
                        }
						
                        //Take the first long region  to start alignment
                        int hit_read_position = get_first_read_map(read, this_long_regions[nregion]->read_map);
                        if (ori==0)
                            hit_read_position += (rstart - rstart_);
                        else
                            hit_read_position += (rend_ - rend);
                        
                        int hit_len = rend - rstart ; // rend is open interval
                        
                        if(ori==1){
                            hit_read_position = read.length() - hit_len - hit_read_position + 1;
                        }
                        
                        if (perform_extra_checks)
                        {
                            assert (hit_read_position>=0) ;
                            if (!(hit_len >0))
                            {
                                fprintf(stderr, "ERROR: hitlen=%i\n", hit_len) ; // BUG-TODO
                                hit_len=0 ;
                                result.delete_regions();
                                delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
                                return -1 ;
                            }
                        }
                        if (myverbosity >= 3)
                        {
                            fprintf(stdout,"read id %s curr len %i\n",read.id(), (int)current_positions.size());
                            fprintf(stdout,	"# Starting point for alignments: read %i, dna %i, len %i\n",hit_read_position, rstart /*this_long_regions[nregion]->start*/, hit_len);					  
                            fprintf(stdout,	"# Number of current regions %i\n",(int)current_regions.size());
                            
                            fprintf(stdout,"DNA:%s\n",current_seq.c_str());
                            fprintf(stdout,"READ:%s\n",read_seq[ori].c_str());
                        }
                        
                        
						int ret = 0 ;

                        if (strand == '+')
                        {
                            ret = perform_alignment_starter_variant(result, hits, read_seq[ori], read_quality[ori], 
                                                                    current_seq, current_regions, current_positions, 
                                                                    chr, '+', ori, hit_read_position, 
                                                                    rstart /*this_long_regions[nregion]->start*/, 
                                                                    hit_len, false, num_alignments_reported, true,
                                                                    annotatedjunctions, variants, myverbosity);
                        }
                        else
                        {
                            if (perform_extra_checks)
                                assert (read.length()-(hit_read_position+hit_len)>=0);
                            
                            ret = perform_alignment_starter_variant(result, hits, read_seq[1 - ori], read_quality[1 - ori], 
                                                                    current_seq, current_regions, current_positions, 
                                                                    chr, '-', 1-ori, read.length()-(hit_read_position+hit_len),
                                                                    rend /*this_long_regions[nregion]->end*/ -1, 
                                                                    hit_len, false, num_alignments_reported, true, 
                                                                    annotatedjunctions, variants, myverbosity);
                        }
						
                        for (size_t i = 0; i < current_regions.size(); i++) {
                            delete current_regions.at(i);
                        }
						current_regions.clear();
						current_positions.clear();
				  						
						if (ret < 0)
						{
							delete_long_regions(long_regions); 
							return ret;
						}
						if (hits.topAlignments().stop_aligning())
						  {
							delete_long_regions(long_regions); 
							return 0 ;
						  }

						num_alignments_reported += ret;						
						hits.topAlignments().update_top_alignment_index();

                }// iterate over combinations
			}//loop through long regions for a chromosome and ori 
		}// loop chr
	}// loop ori
	
			
	result.delete_regions();
	delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
	
	
	//fprintf(stdout, "num_alignments_reported=%i\n", num_alignments_reported) ;
	return num_alignments_reported ;
	
}

template <bool do_timing>	
std::vector<Variant> QPalma::identify_variants(std::string dna, std::vector<int> positions, 
											   Chromosome const &contig_idx, VariantMap & variants, std::map<int,int> & variantpos) const
{

	int k=0;
	while (positions[k]==-2)
		k++;
	int start_pos = positions[k] ;
	k=positions.size()-1;
	while (positions[k]==-2)
		k--;
	int end_pos = positions[k] ;
	int chr=contig_idx.nr() ;
	std::vector<int> map(end_pos-start_pos+1, -1) ; 
	for (unsigned int i=0; i<positions.size(); i++)
	{
		if (positions[i]!=-2)
		{//(dna[i]!='N')
			if (perform_extra_checks)
			{		
				assert(positions[i]-start_pos>=0) ;
				assert(positions[i]-start_pos<=end_pos-start_pos) ;
			}
			map[positions[i]-start_pos] = i ;
		}
	}

	variants.lock() ;

	
	TIME_CODE(clock_t start_time ;
	if (do_timing)
		start_time = clock();)

	std::vector<Variant>::iterator it = my_lower_bound(variants.variantlist[chr].begin(), variants.variantlist[chr].end(), start_pos-100) ;
	TIME_CODE(if (do_timing)
		_stats.variant_lower_bound_time += clock() - start_time;
	if (do_timing)
		start_time = clock();)

	int num_found_variants = 0 ;
	int num_found_SNP_variants = 0 ;
	int num_checked_variants = 0 ;
	std::vector<Variant> variant_list ;
	
	//position of the first variant pointed by it
	int index_position= it - variants.variantlist[chr].begin() ;
	
	while (it != variants.variantlist[chr].end() && (*it).position <= end_pos)
	{
		bool found = false ;
		std::vector<int> match_pos ;
		
		Variant v ;
		v.position = -1 ;
		
		if ((*it).type == pt_SNP)
		{
			assert((*it).ref_len==1 && (*it).variant_len==1) ;
			if ((*it).position-start_pos>=0 && (*it).position-start_pos<end_pos-start_pos && map[(*it).position-start_pos]>=0)
			{
				v=*it ;
				v.position = map[(*it).position-start_pos] ; 
				v.end_position = v.position+v.ref_len ;
				if (perform_extra_checks)
					assert(v.position>=0) ;
				
				if (dna[v.position]!='N')
				{
					//if (v.ref_str[0]!=dna[v.position])
					//	fprintf(stdout, "%c \t %c (%i) \t %s \n", contig_idx[(*it).position], dna[v.position], v.position, v.ref_str.c_str()) ;
					if (perform_extra_checks)
					{
						if (contig_idx[(*it).position]!=dna[v.position])
						{
							fprintf(stderr, "ERROR: sequence mismatch at %i position: %c|%c|%c != %c  (1)\n", v.position, contig_idx[(*it).position-1], contig_idx[(*it).position],contig_idx[(*it).position+1], dna[v.position]) ; 
							assert(contig_idx[(*it).position]==dna[v.position]) ;
						}
						if (v.ref_str[0]=='A' || v.ref_str[0]=='C' || v.ref_str[0]=='G' || v.ref_str[0]=='T')
						{
							if (contig_idx[(*it).position]!=v.ref_str[0])
							{
								fprintf(stderr, "ERROR: sequence mismatch: %c != %c  (2)\n", contig_idx[(*it).position], v.ref_str[0]) ; // BUG-TODO
								assert(contig_idx[(*it).position]==v.ref_str[0]) ;
							}
						}
					}
					found = true ;
				}
			}
		}
		if ((*it).type == pt_insertion && _config.DP_MAX_DELETIONS>0)
		{
			assert((*it).variant_len>0 && (*it).ref_len==0) ;
			if ((*it).position-start_pos>0 && (*it).position-start_pos<=end_pos-start_pos && map[(*it).position-start_pos]>=0)
			{	
				v=*it ;
				v.position = map[(*it).position-start_pos] ;
				//	fprintf(stdout,"found insertion:v.position=%i=map[%i] with start_pos=%i\n",v.position,(*it).position,start_pos);
				if (perform_extra_checks)
					assert(v.position>=0) ;
				v.end_position = v.position + v.ref_len ;
				found = true ;
			}
		}
		if (((*it).type == pt_deletion || (*it).type == pt_substitution) && _config.DP_MAX_DELETIONS>0)
		{
			assert((*it).ref_len>0) ;
			v=*it ;
			v.position = -1 ;
			v.ref_str = "" ;
			unsigned int first_stretch = 0, last_not_ok = 0 ;
			std::vector<int> deleted_positions  ;
			for (int i=0; i<(*it).ref_len; i++)
			{
				if ((*it).position+i-start_pos<0)
					continue ;
				if ((signed)((*it).position+i-start_pos)>=(signed)(end_pos-start_pos))
				{
					last_not_ok = (*it).ref_len-1 ;
					break ;
				}
				if (map[(*it).position+i-start_pos]>=0)
				{
					deleted_positions.push_back(map[(*it).position+i-start_pos]) ;
					if (i==(signed)first_stretch) first_stretch++ ;
				} else
					last_not_ok = i ;
			}
			if (perform_extra_checks)
				assert((*it).ref_len - 1 - last_not_ok>=0) ;
			int last_stretch = (*it).ref_len - 1 - last_not_ok ;
			
			if (deleted_positions.size()>0 && (*it).type == pt_deletion)
			{
				assert((*it).ref_len>0 && (*it).variant_len==0) ;
				v.position = deleted_positions[0] ;
				if (perform_extra_checks)
					assert(v.position>=0) ;
				v.ref_len = 0 ;
				v.ref_str = "" ;
				for (int i=deleted_positions[0]; i<=deleted_positions[deleted_positions.size()-1]; i++)
				{
					v.ref_len++ ;
					v.ref_str += dna[i];
				}
				v.end_position = v.position+v.ref_len ;
				found = true ;
				//fprintf(stdout, "%s\t%s\n", (*it).ref_str.c_str(), v.ref_str.c_str()) ;
			}
			if (deleted_positions.size()>0 && (*it).type == pt_substitution)
			{
				assert((*it).ref_len>0 && (*it).variant_len>0) ;
				if (deleted_positions.size()==(unsigned)(*it).ref_len)
				{
					// use the full substitution
					
					v.ref_str = (*it).ref_str ;
					v.position = deleted_positions[0] ;
					if (perform_extra_checks)
						assert(v.position>=0) ;
					v.end_position = v.position+v.ref_len ;
					found = true ;
					//fprintf(stdout, "full subst: %s\t%s\t%s\n", (*it).ref_str.c_str(), v.ref_str.c_str(), v.variant_str.c_str()) ;
				} 
				else 
					if (false) // buggy
					{
						if ((signed)first_stretch>=(signed)last_stretch && first_stretch>0)
						{
							// only use the 5' end of the substitution
							v.position = deleted_positions[0] ;
							if (perform_extra_checks)
								assert(v.position>=0) ;
							v.ref_len = 0 ;
							v.ref_str = "" ;
							v.variant_str = "" ;
							for (unsigned int i=0; i<first_stretch; i++)
							{
								v.ref_len++ ;
								v.ref_str += dna[deleted_positions[i]] ;
								if (i<(*it).variant_str.size())
									v.variant_str += (*it).variant_str[i] ;
							}
							v.end_position = v.position+v.ref_len ;
							found = true ;
							//fprintf(stdout, "partial subst 5': %s -> %s\t%s -> %s\n", (*it).ref_str.c_str(), (*it).variant_str.c_str(), v.ref_str.c_str(), v.variant_str.c_str()) ;
						} 
						if ((signed)last_stretch>(signed)first_stretch && last_stretch>0)
						{
							// only use the 3' end of the substitution
							v.position = deleted_positions[deleted_positions.size()-1]-last_stretch + 1 ; // +/- 1 ? 
							if (perform_extra_checks)
								assert(v.position>=0) ;
							v.ref_len = 0 ;
							v.ref_str = "" ;
							v.variant_str = "" ;
							for (int i=0; i<last_stretch && deleted_positions.size() - i - 1>0; i++)
							{
								v.ref_len++ ;
								v.ref_str = dna[deleted_positions[deleted_positions.size() - i - 1]] + v.ref_str ;
								if (i<(signed)(*it).variant_str.size())
									v.variant_str = (*it).variant_str[(*it).variant_str.size() - i - 1] + v.variant_str;
							}
							v.end_position = v.position+v.ref_len ;
							found = true ;
							//fprintf(stdout, "partial subst 3': %s -> %s\t%s -> %s\n", (*it).ref_str.c_str(), (*it).variant_str.c_str(), v.ref_str.c_str(), v.variant_str.c_str()) ;
						}
					}
			}
		}
		
		num_checked_variants++ ;
		if (found)
		{
			num_found_variants++ ;
			v.non_used_count=0;
			v.used_count=0;
			v.conf_count=0;
			v.non_conf_count=0;
			variant_list.push_back(v) ;
			if ((*it).type == pt_SNP)
				num_found_SNP_variants++ ;
			
			//Store the position of this particular variant from VariantMap object
			variantpos[v.id]=index_position;
			//fprintf(stdout,"Position %i, id %i, position %i\n", index_position,v.id,v.position);
//			fprintf(stdout,"%i\t%i\t%i\t%i\t%i\t%s\t%s\n",
//					(*it).id,(*it).type, chr,(*it).position,(*it).variant_len-(*it).ref_len,(*it).ref_str.c_str(), (*it).variant_str.c_str());
		}
		it++ ;
		index_position++;
		
	}
	//fprintf(stdout, "%i\t%i\t%i found\n", num_found_SNP_variants, num_found_variants, num_checked_variants) ;
	variants.unlock() ;

	TIME_CODE(if (do_timing)
				  _stats.variant_identification_time += clock() - start_time;)
	
	return variant_list ;
}

int insert_variants(std::vector<Variant> & variant_list, std::string & dna, std::vector<region_t *> & current_regions, std::vector<int> &positions, Chromosome const &contig_idx)
{
	//fprintf(stdout, "not implemented yet: %lu variants\n", variant_list.size()) ;
	int n_SNP=0 ;
	for (unsigned int i=0; i<variant_list.size(); i++)
	{
		if (variant_list[i].type == pt_SNP)
		{
			if (perform_extra_checks)
				assert(dna[variant_list[i].position] == variant_list[i].ref_str[0]) ;
			char c = variant_list[i].variant_str[0] ;
			if (c!='A' && c!='C' && c!='G' && c!='T')
				c='N' ;
			dna[variant_list[i].position] = c ;
			n_SNP ++ ;
		}
	}
	fprintf(stdout, "filled in %i SNPs; %lu\n", n_SNP, variant_list.size()-n_SNP) ;

	return 0 ;
}

int QPalma::perform_alignment_starter_variant(Result &result, Hits &readMappings, 
											  std::string read_string, std::string read_quality, 
											  std::string dna, std::vector<region_t *> current_regions, std::vector<int> positions, 
											  Chromosome const &contig_idx, char strand, int ori,
											  int hit_read_position, int hit_dna_position, int hit_length, 
											  bool non_consensus_search, int& num_alignments_reported, bool remapping, 
											  JunctionMap &annotatedjunctions, VariantMap & variants, int myverbosity) const
{

    std::map<int, int> variant_positions;
    std::vector<Variant> variant_list ;
	if (_config.USE_VARIANTS)
    {
		const bool do_timing = false ;
		variant_list = identify_variants<do_timing>(dna, positions, contig_idx, variants, variant_positions) ;
    }
		
    return perform_alignment_starter_single(result, readMappings, 
                                            read_string, read_quality, 
                                            dna, current_regions, positions, 
                                            contig_idx, strand, ori,
                                            hit_read_position, hit_dna_position, hit_length, 
                                            non_consensus_search, num_alignments_reported, remapping, 
                                            annotatedjunctions, variants, variant_list, variant_positions, myverbosity) ;
} 

// TODO: dd remove relicts from multithreading
int QPalma::perform_alignment_starter_single(Result &result, Hits &readMappings, 
											 std::string read_string, std::string read_quality, 
											 std::string dna, std::vector<region_t *> current_regions, std::vector<int> positions, 
											 Chromosome const &contig_idx, char strand, int ori,
											 int hit_read_position, int hit_dna_position, int hit_length, 
											 bool non_consensus_search, int& num_alignments_reported, bool remapping, 
											 JunctionMap &annotatedjunctions, const VariantMap & variants, std::vector<Variant> & variant_list, 
											 std::map<int, int> & variant_pos, int myverbosity) const
{
	struct perform_alignment_t data ;
	ALIGNMENT* non_consensus_alignment = NULL ;
	ALIGNMENT* consensus_alignment = NULL ;

	try
	{
		{
			//data = new struct perform_alignment_t ;
			
			data.result = &result;
			data.readMappings = &readMappings;
			data.read_string=read_string ;
			data.read_quality=read_quality ;
			data.dna = dna ;
			data.current_regions=current_regions ;
			data.positions=positions ;
			data.contig_idx=&contig_idx;
			data.strand =strand ;
			data.ori = ori ;
			data.ret = -1000 ;
			data.hit_read = hit_read_position;
			data.hit_dna = hit_dna_position;
			data.hit_length = hit_length;
			data.qpalma = this ;
			data.joined=false ;
			data.aln=NULL ;
			data.non_consensus_search=false ;
			data.remapping=remapping;
			data.annotatedjunctions=&annotatedjunctions;
			data.variant_list = &variant_list ;
			data.variants = &variants ;
			data.variant_positions = &variant_pos;
			
			perform_alignment_wrapper1(&data, myverbosity);

			consensus_alignment = data.aln ;
			data.aln = NULL ;
			data.joined = true;
			//delete data ;
		}
		
		if (non_consensus_search) 
		{
			//data = new struct perform_alignment_t ;

			data.result = &result;
			data.readMappings = &readMappings;
			data.read_string=read_string ;
			data.read_quality=read_quality ;
			data.dna = dna ;
			data.current_regions=current_regions ;
			data.positions=positions ;
			data.contig_idx=&contig_idx;
			data.strand =strand ;
			data.ori = ori ;
			data.ret = -1000 ;
			data.hit_read = hit_read_position;
			data.hit_dna = hit_dna_position;
			data.hit_length = hit_length;
			data.qpalma = this ;
			data.joined=false ;
			data.aln=NULL ;
			data.non_consensus_search=true ;
			data.remapping=remapping;
			data.annotatedjunctions=&annotatedjunctions;
			data.variant_list = &variant_list ;
			data.variants = &variants ;
			data.variant_positions = &variant_pos;

			perform_alignment_wrapper1(&data, myverbosity);
				
			non_consensus_alignment=data.aln ;
			data.aln=NULL;

			//delete data ;
		}
		
		if (consensus_alignment && (consensus_alignment->passed_filters || non_consensus_search))
		{
			if (myverbosity>=1)
				fprintf(stdout, "added one consensus alignment\n") ;
			readMappings.topAlignments().add_alignment_record(consensus_alignment, 1) ;
			num_alignments_reported++ ;
			if (non_consensus_alignment)
				delete non_consensus_alignment ;
		}
		else
		{
			if (non_consensus_alignment && non_consensus_alignment->passed_filters)
			{
				if (myverbosity>=1)
					fprintf(stdout, "added one nonconsensus alignment\n") ;
				readMappings.topAlignments().add_alignment_record(non_consensus_alignment, 1) ;
				num_alignments_reported++ ;
			}
			if (consensus_alignment)
				delete consensus_alignment ;
			
		}
		
		
		return 0 ; // in the meaning of perform_alignment this corresponds to a spliced alignment
	}
	catch (std::bad_alloc&)
	{
	    fprintf(stderr, "WARNING: perform_alignment failed with out of memory\n") ;
		delete consensus_alignment ;
		delete non_consensus_alignment ;
	    
	    return -1 ;
	}
}



struct pos_table_str
{
	int pos ;
	char nuc ;
	bool to_be_deleted ;
	double acc ;
	double don ;
	std::vector<pos_table_str *> del_refs ;
	std::vector<int> del_ids ;
	std::vector<Variant*> snps ;
	std::vector<int> snp_ids ;
	std::vector<pos_table_str *> del_origs ;
} ;


inline std::vector<struct pos_table_str *>::iterator  pos_table_lower_bound ( std::vector<struct pos_table_str *>::iterator first, std::vector<struct pos_table_str *>::iterator last, const int& value )
{
	std::vector<struct pos_table_str *>::iterator it;
	long int count, step;
	count = distance(first,last);

	while (count>0)
	{
		it = first;
		step = count/2;
		advance(it,step);

		// find closest non-negative position
		std::vector<struct pos_table_str *>::iterator itf=it, itb=it;
		while ((itf<last && (*itf)!=NULL && (*itf)->pos<0) || (itb>=first && (*itb)!=NULL && itb<last && (*itb)->pos<0))
		{
			if (itf>=last && itb<first)
				break ;
			if (itf<last)
				advance(itf, 1) ;
			if (itb>=first)
				advance(itb, -1) ;
		}
		if (itf<last && (*itf)->pos>=0)
			it=itf ;
		else
			if (itb>=first && (*itf)!=NULL && itb<last && (*itb)->pos>=0)
				it=itb ;
			else
			{
				while (first<last && (*first)->pos<value)
					advance(first, 1) ;
				
				return first;
			}
		
		if ( (*it)->pos < value) 
		{
			first=it; //++it;
			count-=step+1;
		}
		else count=step ;
	}

	while (first<last && (*first)->pos<value)
		advance(first, 1) ;
	
	return first;
}


int find_pos(std::vector< struct pos_table_str *> &pos_table, int position)
{
	// todo: speedup by binary search 
	
	std::vector<struct pos_table_str *>::iterator first = pos_table_lower_bound(pos_table.begin(), pos_table.end(), position) ;
	int p = first - pos_table.begin() ;
	if (perform_extra_checks)
	{
		//fprintf(stderr, "ERROR: Position %i not found \n", position) ; 
		//assert(pos_table[p]->pos==position) ;
	}
	if (p<(int)pos_table.size() && p>=0)
		for (unsigned int i=p; i<pos_table.size() && i-p<10; i++)
			if (pos_table[i]->pos == position)
			return i ;

	for (unsigned int i=0; i<pos_table.size(); i++)
		if (pos_table[i]->pos == position)
		{
			//assert(p==(int)i) ;
			return i ;
		} 

	fprintf(stderr, "ERROR: Position %i not found \n", position) ; 
	assert(0) ;
	
	return -1 ;
}

void change_pos_table_deletion_ends(struct pos_table_str * pos_table_previous_end_p, struct pos_table_str * pos_table_new_end_p)
{

	for (unsigned int i=0; i<pos_table_previous_end_p->del_origs.size();i++){
		struct pos_table_str * origin_p=pos_table_previous_end_p->del_origs[i];
		
		for (unsigned int j=0; j< origin_p->del_refs.size();j++){
			if (origin_p->del_refs[j] == pos_table_previous_end_p){
				origin_p->del_refs[j] = pos_table_new_end_p;	
			}
		}
	}

	pos_table_previous_end_p->del_origs.clear();
}

template <int myverbosity, bool do_timing>	
std::vector<variant_cache_t *> QPalma::create_super_sequence_from_variants(std::vector<Variant> & variants, std::string & dna, double *&acceptor, int &a_len, double *&donor, int &d_len, 
                                                int &hit_dna_pos, std::vector<bool> &ref_map) const 
{
	int seed_ref=hit_dna_pos;

	TIME_CODE(clock_t my_start_time ;
	if (do_timing)
		my_start_time=clock() ;)
	
	struct pos_table_str * pos_table_= new struct pos_table_str[dna.size()] ;
	
	std::vector< struct pos_table_str *> pos_table(dna.size(), NULL) ;
	for (unsigned int i=0; i<dna.size(); i++)
	{
		pos_table[i] = &pos_table_[i] ; //new struct pos_table_str ;
		pos_table[i]->to_be_deleted = false ;
		pos_table[i]->pos = i ;
		pos_table[i]->nuc = dna[i] ;
		pos_table[i]->acc = acceptor[i] ;
		pos_table[i]->don = donor[i] ;
	}
	int nbv_dels=0, nbv_snp=0, nbv_ins=0, nbv_subst=0;
	
	for (unsigned int i=0; i<variants.size(); i++)
	{
		if (variants[i].type == pt_deletion)
		{
			if (variants[i].position-1>=0 && variants[i].position + variants[i].ref_len<(int)pos_table.size()) // cannot encode deletion at the beginning ... how bad is this?
			{
				pos_table[variants[i].position-1]->del_refs.push_back(pos_table[variants[i].position + variants[i].ref_len]) ;
				//pos_table[variants[i].position-1]->del_ids.push_back(variants[i].id) ;
				pos_table[variants[i].position-1]->del_ids.push_back(i) ;
				pos_table[variants[i].position + variants[i].ref_len]->del_origs.push_back(pos_table[variants[i].position -1]) ;
				nbv_dels++ ;
			} 
			else
				if (myverbosity>=2)
					fprintf(stdout, "dropped deletion of length %i at beginning or end of sequence\n", variants[i].ref_len) ;
		}
		
		if (variants[i].type == pt_SNP)
		{
			pos_table[variants[i].position]->snps.push_back(&variants[i]) ;
			pos_table[variants[i].position]->snp_ids.push_back(i) ;
			nbv_snp++ ;
		}
	}	
	TIME_CODE(if (do_timing)
				  _stats.variant_create_super_sequence_from_variants_init_time += clock() - my_start_time;
	if (do_timing)
		my_start_time=clock()) ;

	for (unsigned int i=0; i<variants.size(); i++)
	{
		if (variants[i].type == pt_insertion)
		{
			//Position cannot be 0 because insertion is considered before this position
			int idx = find_pos(pos_table, variants[i].position) ;
			if (idx>0)
			{
				nbv_ins++ ;				
				for (int j=0; j<variants[i].variant_len; j++)
				{
					std::vector<struct pos_table_str*>::iterator it = pos_table.begin() + idx + j ;
					struct pos_table_str *p = new struct pos_table_str ;
					p->nuc = variants[i].variant_str[j] ;
					p->acc = -ALMOST_INFINITY ;
					p->don = -ALMOST_INFINITY ;
					p->pos = -1 ;
					p->to_be_deleted=true ;
					
					pos_table.insert(it, p) ;
				}
				change_pos_table_deletion_ends(pos_table[idx+variants[i].variant_len],pos_table[idx]);
				pos_table[idx-1]->del_refs.push_back(pos_table[idx+variants[i].variant_len]) ;
				pos_table[idx+variants[i].variant_len]->del_origs.push_back(pos_table[idx-1]) ;
				pos_table[idx-1]->del_ids.push_back(i) ; 
				//pos_table[idx-1]->del_ids.push_back(variants[i].id) ; 
			}
			else
				fprintf(stdout, "dropped insertion of length %i at beginning or end of sequence\n", variants[i].variant_len) ;

		}
		
		if (variants[i].type == pt_substitution && false)
		{
			//fprintf(stdout, "%i, %i\n", variants[i].position, variants[i].end_position) ;
			
			int idx_start=find_pos(pos_table, variants[i].position) ;
			int idx_end=find_pos(pos_table, variants[i].end_position) ;

			if (idx_start>0 && idx_end>0)
			{
				nbv_subst++ ;
				const int num_N = 10 ;
				for (int j=0; j<num_N; j++)
				{
					std::vector<struct pos_table_str*>::iterator it = pos_table.begin() + idx_end + 1 + j ;
					struct pos_table_str *p = new struct pos_table_str ;
					p->nuc = 'N' ;
					p->acc = -ALMOST_INFINITY ;
					p->don = -ALMOST_INFINITY ;
					p->pos = -1 ;
					p->to_be_deleted=true ;
					
					pos_table.insert(it, p) ;
				}
				
				for (int j=0; j<variants[i].variant_len; j++)
				{
					std::vector<struct pos_table_str*>::iterator it = pos_table.begin() + idx_end + 1 + j + num_N ;
					struct pos_table_str *p = new struct pos_table_str ;
					p->nuc = variants[i].variant_str[j] ;
					p->acc = -ALMOST_INFINITY ;
					p->don = -ALMOST_INFINITY ;
					p->pos = -1 ;
					p->to_be_deleted=true ;
					
					pos_table.insert(it, p) ;
				}
				// skipping the reference version
				pos_table[idx_start-1]->del_refs.push_back(pos_table[idx_end+1+num_N]) ;
				//pos_table[idx_start-1]->del_ids.push_back(variants[i].id) ;
				pos_table[idx_start-1]->del_ids.push_back(i) ;

				// skipping the variant version
				change_pos_table_deletion_ends(pos_table[idx_end+num_N+variants[i].variant_len],pos_table[idx_end+num_N]); //because of insertion: change ends of former deletion back to before the insertion
				pos_table[idx_end]->del_refs.push_back(pos_table[idx_end+num_N+variants[i].variant_len]) ;
				pos_table[idx_end]->del_ids.push_back(i) ; //-1 because the id is already included in the first deletion
				pos_table[idx_end+num_N+variants[i].variant_len]->del_origs.push_back(pos_table[idx_end]) ;
			}
		}
	}
	TIME_CODE(if (do_timing)
				  _stats.variant_create_super_sequence_from_variants_insertion_init_time += clock() - my_start_time;)

	if (myverbosity>=1)
		fprintf(stdout, "Created supersequence of length %ld from reference sequence of length %ld using %ld variants\n", dna.length(), pos_table.size(), variants.size()) ;
	
	TIME_CODE(if (do_timing)
				  my_start_time=clock() ;)
	
	dna = std::string(pos_table.size(), ' ') ;
	delete[] acceptor ;
	delete[] donor ;
	acceptor=new double[pos_table.size()] ; 
	a_len = pos_table.size() ;
	donor=new double[pos_table.size()] ; 
	d_len = pos_table.size() ;
	//std::vector<bool> ref_map(pos_table.size(), false) ;
	ref_map.assign(pos_table.size(), false) ;

	for (unsigned int i=0; i<pos_table.size(); i++)
	{		
		if (pos_table[i]->pos>=0)
			ref_map[i]=true ;
		//Always compared to the reference
		if (seed_ref == pos_table[i]->pos)
			hit_dna_pos = i ;
		pos_table[i]->pos = i ;	
		acceptor[i] = pos_table[i]->acc ;
		donor[i] = pos_table[i]->don ;
		dna[i] = pos_table[i]->nuc ;
	}
	TIME_CODE(if (do_timing)
				  _stats.variant_create_super_sequence_from_variants_table_init_time += clock() - my_start_time;)

	
	std::vector<variant_cache_t *> variant_cache(d_len,NULL);
	
	int nb_snps = 0 ;
	int nb_dels = 0 ;
	
	TIME_CODE(if (do_timing)
				  my_start_time=clock() ;)
	for (unsigned int i=0; i<pos_table.size(); i++)
	{
		for (unsigned j = 0; j<pos_table[i]->del_refs.size(); j++)
		{

			int start_pos=i +1;
			int end_pos=pos_table[i]->del_refs[j]->pos -1;
			int id=pos_table[i]->del_ids[j] ;
			
			if (perform_extra_checks){
				assert(start_pos < d_len);
				assert(start_pos >=0);
				assert(end_pos < d_len);
				assert(end_pos >=0);
			}
						
			//left side: the end of the deletion is encountered first
			if ((int)end_pos<=hit_dna_pos){
				nb_dels++;
				
				if (variant_cache[end_pos]==NULL){
					variant_cache[end_pos]= new variant_cache_t;
					variant_cache[end_pos]->insertion = -1;
				}
				variant_cache[end_pos]->end_positions.push_back(start_pos);
				variant_cache[end_pos]->id_dels.push_back(id);

				if (variants[id].type == pt_insertion){
					if (variant_cache[start_pos]==NULL){
						variant_cache[start_pos]= new variant_cache_t;
						variant_cache[start_pos]->insertion = -1;
					}	 
					//insertion unique
					if (perform_extra_checks)
						assert(variant_cache[start_pos]->insertion <0);
					variant_cache[start_pos]->insertion = id;
				}
				
			}
			//right side: the beginning of the deletion is encountered first
			else if (start_pos >= hit_dna_pos){
				nb_dels++;

				if (variant_cache[start_pos]==NULL){
					variant_cache[start_pos]= new variant_cache_t;
					variant_cache[start_pos]->insertion = -1;
				}
				variant_cache[start_pos]->end_positions.push_back(end_pos);
				variant_cache[start_pos]->id_dels.push_back(id);	
				if (variants[id].type == pt_insertion){
					//insertion unique
					if (perform_extra_checks)
						assert(variant_cache[start_pos]->insertion <0);
					
					variant_cache[start_pos]->insertion = id;
				}
			}
		}
	}
	TIME_CODE(if (do_timing)
				  _stats.variant_create_super_sequence_from_variants_table_fill1_time += clock() - my_start_time;
	if (do_timing)
		my_start_time=clock() ;)

	if (!_config.IUPAC_SNPS)
	{
		for (unsigned int i=0; i<pos_table.size(); i++)
		{
			for (unsigned j = 0; j<pos_table[i]->snps.size(); j++)
			{
				nb_snps++ ;

				if (variant_cache[i]==NULL)
				{
					variant_cache[i]= new variant_cache_t;
					variant_cache[i]->insertion = -1;
				}
				variant_cache[i]->snps.push_back(pos_table[i]->snps[j]->variant_str[0]);
				variant_cache[i]->id_snps.push_back(pos_table[i]->snp_ids[j]);
			}
		}
	}
	else{
		for (unsigned int i=0; i<pos_table.size(); i++)
		{
			//Merge SNPs and DNA base in dna[i]
			if ((int)pos_table[i]->snps.size() > 0)
			{
				if (variant_cache[i]==NULL){
					variant_cache[i]= new variant_cache_t;
					variant_cache[i]->insertion = -1;
					
				}
				//Put the original DNA value
				variant_cache[i]->snps.push_back(dna[i]);
				variant_cache[i]->id_snps.push_back(pos_table[i]->snp_ids[0]);
				assert(variant_cache[i]->snps.size()==1) ;
				
				char merged_base;
				int start=0;
				//fprintf(stdout,"original base %c\n",dna[i]);				
				//If DNA is 'N', don't include it in merging stuff
				if (dna[i]=='N')
				{
					merged_base = pos_table[i]->snps[0]->variant_str[0];
					start++;
				}
				else
				{
					merged_base = dna[i];
				}
				
				for (unsigned j = start; j<pos_table[i]->snps.size(); j++)
				{
					merged_base = get_IUPAC_code(merged_base, pos_table[i]->snps[j]->variant_str[0]);
					//fprintf(stdout,"SNP base %c -> %c\n",pos_table[i]->snps[j]->variant_str[0],merged_base);
				}
				dna[i]=merged_base;				
			}
		}
	}
	TIME_CODE(if (do_timing)
				  _stats.variant_create_super_sequence_from_variants_table_fill2_time += clock() - my_start_time;)

	for (unsigned int i=0; i<pos_table.size(); i++)
	{
		if (pos_table[i]->to_be_deleted)
			delete pos_table[i] ;
		pos_table[i]=NULL ;
	}
	delete[] pos_table_ ;

	if (myverbosity>=1)
		fprintf(stdout, "found %i variants (%i snps, %i dels from nbv_dels=%i, nbv_snp=%i, nbv_ins=%i, nbv_subst=%i)\n", 
				nb_snps+nb_dels, nb_snps, nb_dels, nbv_dels, nbv_snp, nbv_ins, nbv_subst) ;
	
	return variant_cache ;
}


int report_variant_at_read_pos(Variant & variant, int read_pos)
{

	
	variant.read_pos=read_pos;
	variant.used_count+=1;
	variant.non_used_count=0;
	variant.conf_count=0;
	variant.non_conf_count=0;

	return 1;

}
		
void QPalma::recover_variants_on_ref(Variant &variant,std::vector<int> positions,char strand, int read_len,Chromosome const &contig_idx) const
{
	
	
	if (strand =='+')
	{
		//Get directly original positions
		variant.position=positions[variant.position];
		if (variant.type == pt_insertion)
			variant.end_position=positions[variant.end_position];
		else
			variant.end_position=positions[variant.end_position-1]+1;
	}
	else{	
		//Give sequences according to positive strand
		variant.ref_str=reverse(complement(variant.ref_str));
		variant.variant_str=reverse(complement(variant.variant_str));
		//Give position on read oriented according to positive strand (like in SAM output)
		if (variant.read_pos>=0)
			variant.read_pos=read_len-1-variant.read_pos;
		
		//Get original positions: switch start and end
		//Particular case for insertion where start=end and insertion is before the given position
		//insertion before 3 on +: 1 2 | 3 <=> 3 | 2 1 insertion before 2 on -

		if (variant.type == pt_insertion)
		{
			int start_v=positions[variant.position-1];
			int end_v=positions[variant.end_position-1];
			variant.position=end_v;
			variant.end_position=start_v;
		}
		else
		{
			int start_v=positions[variant.position]+1;
			int end_v=positions[variant.end_position-1];
			variant.position=end_v;
			variant.end_position=start_v;
		}
	}

	if (perform_extra_checks){
		assert (variant.position >=0);
		assert (variant.end_position >=0);
	}
	

	//Case for long deletion or substitution which spans several regions
	if (variant.type == pt_deletion || variant.type == pt_substitution){
		
		if (variant.end_position-variant.position != variant.ref_len){

			int ref_len=0;
			std::string ref_str="";
			
			for (int i=variant.position; i<variant.end_position;i++){
				ref_str+=contig_idx[i];
				ref_len++;
			}
			
			variant.ref_str=ref_str;
			variant.ref_len=ref_len;
			ref_str.clear();
		}
	   		
	}
}



int get_end_position (Variant &variant, bool is_ref, int pos)
{
	if (variant.type == pt_SNP){
		return pos;
	}
	else if (variant.type == pt_deletion){
		return pos + variant.ref_len -1;
	}
	else if (variant.type == pt_substitution){
		if (is_ref)
			return pos + variant.ref_len-1 + 10;
		else
			return pos + variant.variant_len-1 + 10;
	}
	else if (variant.type == pt_insertion){
		return pos + variant.variant_len-1;
	}

	return -1;
	
}

template <int myverbosity>
int QPalma::reconstruct_reference_alignment(std::vector<Variant> & variants, const std::vector<FoundVariant> & found_variants, std::string & dna, const std::vector<bool> & ref_map, 
											int * &s_align, int & s_len, int *&e_align, int & e_len,int *&dna_align,int *&read_align,int &result_length,bool remapping, 
											bool& alignment_passed_filters,	int &alignment_gaps, int &alignment_mm, 
											const std::vector<variant_cache_t *> &variant_cache, bool report_variants, const char * read_id) const
{

	if (myverbosity>=3)
	{
		fprintf(stdout, "DNA: %s\n",(char*)dna.c_str());
		fprintf(stdout,"ref_map: ");	
		for (int j=0; j<(int)ref_map.size();j++){
			fprintf(stdout,"%i",ref_map[j]);		
		}
		fprintf(stdout,"\n");
		fprintf(stdout,"s_align: ");	
		for (int j=0; j<s_len;j++){
			fprintf(stdout,"%i",s_align[j]);		
		}
		fprintf(stdout,"\n");	
		fprintf(stdout,"dna_align:  ");
		for (int j=0; j<result_length;j++){
			fprintf(stdout,"%i",dna_align[j]);		
		}
		fprintf(stdout,"\n");

		fprintf(stdout,"read_align: ")	;
		for (int j=0; j<result_length;j++){
			fprintf(stdout,"%i",read_align[j]);		
		}
		fprintf(stdout,"\n length=%i\n",result_length);	


		for (unsigned int pos=0; pos< variants.size();pos++){
			
			fprintf(stdout,"V[%i] type=%i id=%i pos=%i ref(%i)=%s var(%i)=%s\n",pos,variants[pos].type,variants[pos].id, variants[pos].position,variants[pos].ref_len, (char*)variants[pos].ref_str.c_str(),variants[pos].variant_len,(char *)variants[pos].variant_str.c_str());
		}
		
		fprintf(stdout,"------\n");
		
		for (unsigned int i=0; i< found_variants.size();i++){
			int pos=found_variants[i].id;
			
			fprintf(stdout,"V[%i] type=%i id=%i ref(%i)=%s var(%i)=%s\n",pos,variants[pos].type,variants[pos].id, variants[pos].ref_len, (char*)variants[pos].ref_str.c_str(),variants[pos].variant_len,(char *)variants[pos].variant_str.c_str());
		}
	}

	//Vectors for storing alignment back to the reference
	std::vector<char> dna_back;
	std::vector<int> s_align_back;
	std::vector<int> dna_align_back;
	std::vector<int> read_align_back;

	int align_ind=-1; //Index in alignment string (only aligned nucleotides)
	
	//Current alignment length
	int result_tmp=result_length;

	//Variables to compute intron, exon lengths and edit operation counts for the alignment with variants
	//If no variant used, should be equivalent to the numbers on the reference sequence
	int max_intron_len=0;
	int min_intron_len=10000000;
	int min_exon_len=e_len;
	alignment_mm=0;
	alignment_gaps=0;
	int len_current_exon=-1;
	int len_current_intron=-1;
	int num_exons=0;

	//Variable when spanning a used variant
	int current_variant_end=-1;
	FoundVariant *next_v;
	next_v=NULL;

	int insertion_pos=-1;
	
	int it=-1;
	if (found_variants.size()>0){
		it=0;
		next_v= new FoundVariant;
		*next_v =found_variants[it];
	}
	
	
		
	//Read position to report for a used variant
	int read_pos=-1;

	//Number of used variants to report
	int used_variants=0;
	
	int original_dna=-1;

	//Gap on DNA from alignment: keep them
	while(align_ind+1<result_tmp && dna_align[align_ind+1]==0)
	{
		alignment_gaps++;
		dna_align_back.push_back(dna_align[align_ind+1]);
		read_align_back.push_back(read_align[align_ind+1]);
		align_ind++;
		read_pos++;
	} 

	unsigned int i=0;
	while(i<dna.length())
	{

		//Build reference sequence
		if (ref_map[i])
		{
			//Take the original base if it exists SNP for this position
			if (_config.IUPAC_SNPS)
			{
				if (variant_cache[i] != NULL && !variant_cache[i]->snps.empty())
				{
					assert(variant_cache[i]->snps.size()==1) ; // in this case just one special SNP is added (the reference base)
					dna_back.push_back(variant_cache[i]->snps[0]);
					original_dna=map_back(variant_cache[i]->snps[0]);
				}
				else
				{
					dna_back.push_back(dna[i]);
					original_dna=map_back(dna[i]);
				}
				
			}
			else
			{
				dna_back.push_back(dna[i]);
				original_dna=map_back(dna[i]);
			}
		}	
	
		//Deletion on super sequence was taken
		if (s_align[i]==5)
		{
			align_ind++;
			//Deletion on ref

			//fprintf(stdout,"5: align_ind=%i, ref_map[%i]=%i\n",align_ind,i,ref_map[i]?1:0);

			if(ref_map[i]) 
			{
				s_align_back.push_back(0);

				if (perform_extra_checks)
					assert(align_ind>=0 && align_ind<result_tmp);

				if (_config.IUPAC_SNPS)
					dna_align_back.push_back(original_dna);
				else
				{
					assert(original_dna==dna_align[align_ind]);
					dna_align_back.push_back(dna_align[align_ind]);
				}
				read_align_back.push_back(0);

				if (report_variants)
				{
					//New deletion/substitution start: get end position
					if (current_variant_end == -1){
						if (next_v != NULL)
							current_variant_end=next_v->end_pos;
					}
					//End deletion on ref
					if (current_variant_end == (int)i){

						if (perform_extra_checks){
							assert(next_v !=NULL);
							assert(next_v->id>=0 && next_v->id<(int)variants.size()) ;
						}
						
						used_variants+=report_variant_at_read_pos(variants[next_v->id], read_pos+1);
						
						//Get new variant position
						current_variant_end =-1;
						delete next_v;
						next_v= NULL;
						it++;
						if (perform_extra_checks)
							assert(it>=0) ;
						if (it < (int)found_variants.size()){
							next_v= new FoundVariant;
							*next_v=found_variants[it];
						}
					}
				}
			}
			//Insertion on ref not taken or N part of a imbalanced substitution
			else{
				if (report_variants){
					//End deletion on ref for imbalanced substitution
					if (current_variant_end == (int)i){
						if (perform_extra_checks)
							assert(next_v !=NULL && next_v->id>=0 && next_v->id<(int)variants.size()) ;
						used_variants+=report_variant_at_read_pos(variants[next_v->id],read_pos+1);
						current_variant_end =-1;

						//Get new variant position
						delete next_v;
						next_v=NULL;
						it++;
						if (perform_extra_checks)
							assert(it>=0) ;
						if (it < (int)found_variants.size()){
							next_v=new FoundVariant;
							*next_v=found_variants[it];
						}
					}
				}
				result_length--;
			}
			
			//Gap on DNA from alignment: keep them
			while(align_ind+1<result_tmp && dna_align[align_ind+1]==0){
				alignment_gaps++;
				dna_align_back.push_back(dna_align[align_ind+1]);
				read_align_back.push_back(read_align[align_ind+1]);
				align_ind++;
				read_pos++;
			} 

			i++;
			continue;
		}
		
		
		//Nucleotide on the super sequence DNA is used to align
		if (s_align[i]==0){
			align_ind++;
			if (read_align[align_ind]!=0){
				read_pos++;
			}
			
			//End of an intron
			if (len_current_intron != -1){
				if (len_current_intron>max_intron_len)
					max_intron_len=len_current_intron;
				if (len_current_intron<min_intron_len)
					min_intron_len=len_current_intron;
				len_current_intron=-1;
			}
			
			//Start new exon
			if (len_current_exon ==-1)
				len_current_exon++;
			
			if (perform_extra_checks)
				assert(align_ind>=0 && align_ind<result_tmp);
			
			//Mismatch
			if (dna_align[align_ind]!=read_align[align_ind] 
				&& dna_align[align_ind] != 0 && read_align[align_ind] != 0 &&
				read_align[align_ind]<=4 && dna_align[align_ind]<=4) // BUG-TODO: added this condition to avoid counting SNPs
			{
				alignment_mm++;
			}

			//Gap on read
			if (read_align[align_ind]==0){
				alignment_gaps++;
			}
			
			len_current_exon++;
			
			//This position is from the reference
			if(ref_map[i]){
				s_align_back.push_back(0);

				//SNP: get the original base
				if(_config.IUPAC_SNPS)
				{
					dna_align_back.push_back(original_dna);
				}
				else{
					if (map_back(dna[i])!= dna_align[align_ind]){
						dna_align_back.push_back(map_back(dna[i]));
						
						if (report_variants){
							if (perform_extra_checks)
								assert(next_v != NULL && next_v->id>=0 && next_v->id<(int)variants.size()) ;
							used_variants+=report_variant_at_read_pos(variants[next_v->id], read_pos);
						
							delete next_v;
							next_v=NULL;
							it++;
							if (perform_extra_checks)
								assert(it>=0) ;
							if (it < (int)found_variants.size()){
								next_v=new FoundVariant;
								*next_v=found_variants[it];
							}
						}
					}
					else{
						if (perform_extra_checks)
							assert(align_ind>=0 && align_ind<result_tmp) ;
						dna_align_back.push_back(dna_align[align_ind]);
					}
				}
				
				
				read_align_back.push_back(read_align[align_ind]);
			}
			
			//Insertion on ref used
			else{
				if (report_variants)
				{
					//Start of a new insertion
					if (current_variant_end==-1 && variant_cache[i] != NULL && variant_cache[i]->insertion != -1){
						insertion_pos=variant_cache[i]->insertion;
						current_variant_end=get_end_position(variants[insertion_pos],false,i);
					}
					
					//End of insertion: report it
					if (current_variant_end==(int)i)
					{
						if (perform_extra_checks && insertion_pos==-1)
						{
							fprintf(stderr, "Error: insertion_pos==-1\n") ;
							alignment_passed_filters=false ;
							return 0 ;
						}
						if (perform_extra_checks)
							assert(insertion_pos>=0 && insertion_pos<(int)variants.size()) ; // BUG-TODO
						used_variants+=report_variant_at_read_pos(variants[insertion_pos], read_pos);
						current_variant_end=-1;
						insertion_pos=-1;
					}
				}
				

				//If not a gap on read: keep it as a gap on dna
				if(read_align[align_ind]!=0){
					dna_align_back.push_back(0);
					read_align_back.push_back(read_align[align_ind]);
				}
				//If gap on read: length is actually one less
				else{
					result_length--;
				}
				
			}
			//Gap on DNA from alignment: keep them
			while(align_ind+1<result_tmp && dna_align[align_ind+1]==0){
				alignment_gaps++;
				dna_align_back.push_back(dna_align[align_ind+1]);
				read_align_back.push_back(read_align[align_ind+1]);
				align_ind++;
				read_pos++;
			} 

			i++;
			continue;
		}

		//1,2,3 for introns and 4 for unmapped
		if(ref_map[i])
			s_align_back.push_back(s_align[i]);
		else{
			
			if (report_variants){

				//Start of a new insertion
				if (current_variant_end==-1 && variant_cache[i] != NULL && variant_cache[i]->insertion != -1){
					insertion_pos=variant_cache[i]->insertion;
					current_variant_end=get_end_position(variants[insertion_pos],false,i);
				}
				
				//End of insertion can fall here
				if (current_variant_end==(int)i){
					used_variants+=report_variant_at_read_pos(variants[insertion_pos], read_pos);
					current_variant_end=-1;
					insertion_pos=-1;
				}
			}
			
		}
		
		
		//End of an exon
		if (len_current_exon!=-1){
			num_exons++;
			if(len_current_exon<min_exon_len)
				min_exon_len=len_current_exon;
			len_current_exon=-1;
		}

		if (s_align[i]!=4){
			align_ind++;

			//New intron
			if (ref_map[i]){
				if (_config.IUPAC_SNPS)
					dna_align_back.push_back(original_dna);
				else{
					assert(original_dna==dna_align[align_ind]);
					dna_align_back.push_back(dna_align[align_ind]);
				}
				read_align_back.push_back(read_align[align_ind]);
				if (len_current_intron ==-1)
					len_current_intron++;
				
				len_current_intron++;

			}
			else{
				result_length--;
			}
			
			
		}
	
		//Gap on DNA from alignment: keep them
		while(align_ind+1<result_tmp && dna_align[align_ind+1]==0){
			alignment_gaps++;
			dna_align_back.push_back(dna_align[align_ind+1]);
			read_align_back.push_back(read_align[align_ind+1]);
			align_ind++;
			read_pos++;
		} 	
		i++;

	}
	
	if (len_current_exon!=-1)
	{
		num_exons++;
		if (len_current_exon<min_exon_len)
			min_exon_len=len_current_exon;
		len_current_exon=-1;
	}


	if (next_v!=NULL){
		delete next_v;
		next_v=NULL;
	}
	
		

	if (perform_extra_checks && (result_length!=(int)dna_align_back.size() || result_length!=(int)read_align_back.size()))

	{
		fprintf(stderr, "ERROR: len mismatch %i!=%ld || %i!=%ld\n", result_length, dna_align_back.size(), result_length, read_align_back.size()) ; 
		alignment_passed_filters = false ;
		assert(result_length==(int)dna_align_back.size() && result_length==(int)read_align_back.size()) ;
	}
	else	
	{
		delete[] s_align;
		delete[] dna_align;
		delete[] read_align;
		s_align = new int[s_align_back.size()];
		s_len=s_align_back.size();
		dna_align = new int[dna_align_back.size()];
		read_align = new int[read_align_back.size()];
		dna = std::string(dna_back.size(), ' ') ;
		
		for (unsigned int j=0; j<dna_back.size();j++){
			dna[j]=dna_back[j];
		}
		
		for (unsigned int j=0; j<s_align_back.size();j++){
			s_align[j]=s_align_back[j];
		}
		
		
		for (unsigned int j=0; j<dna_align_back.size();j++){
			dna_align[j]=dna_align_back[j];
		}
		
		for (unsigned int j=0; j<read_align_back.size();j++){
			read_align[j]=read_align_back[j];
		}
		
		s_align_back.clear();
		dna_align_back.clear();
		read_align_back.clear();
		//fprintf(stdout, "alignment_mm=%i, alignment_gaps=%i, max_intron_len=%i, min_intron_len=%i, num_exons=%i\n", alignment_mm, alignment_gaps, max_intron_len, min_intron_len, num_exons) ;

		if (_config.NUM_MISMATCHES==0 && alignment_mm>0)
		{ // in this case the DP found a non-mismatch version, but somehow the subsequent processing determines that there are mismatches
			if (myverbosity>=1 || perform_extra_checks)
				fprintf(stdout, "Warning: _config.NUM_MISMATCHES==0 && alignment_mm=%i -> setting to 0 (%s)\n", alignment_mm, read_id); // BUG-TODO
			alignment_mm=0 ;
		}
		if (_config.NUM_MISMATCHES==1 && alignment_mm>1)
		{ // in this case the DP found a non-mismatch version, but somehow the subsequent processing determines that there are mismatches
			if (myverbosity>=1 || perform_extra_checks)
				fprintf(stdout, "Warning: _config.NUM_MISMATCHES==1 && alignment_mm=%i -> setting to 1 (%s)\n", alignment_mm, read_id); // BUG-TODO
			alignment_mm=1 ;
		}
		if (_config.NUM_MISMATCHES==2 && alignment_mm>2)
		{ // in this case the DP found a non-mismatch version, but somehow the subsequent processing determines that there are mismatches
			if (myverbosity>=1 || perform_extra_checks)
				fprintf(stdout, "Warning: _config.NUM_MISMATCHES==2 && alignment_mm=%i -> setting to 2 (%s)\n", alignment_mm, read_id); // BUG-TODO
			alignment_mm=2 ;
		}
		if (_config.NUM_GAPS==0 && alignment_gaps>0)
		{ // same here for gaps
			if (myverbosity>=1 || perform_extra_checks)
				fprintf(stdout, "Warning: _config.NUM_GAPS==0 && alignment_gaps=%i -> setting to 0 (%s)\n", alignment_gaps, read_id); // BUG-TODO
			alignment_gaps=0 ;
		}

		alignment_passed_filters = alignment_pass_filters(min_intron_len,max_intron_len,alignment_mm,alignment_gaps,num_exons,min_exon_len,remapping);
	}
	
	if (myverbosity>=3)
	{
		fprintf(stdout, "DNA: %s\n",(char*)dna.c_str());
		
		fprintf(stdout,"new s_align: ");	
		for (int j=0; j<s_len;j++){
			fprintf(stdout,"%i",s_align[j]);		
		}
		fprintf(stdout,"\n");	
		fprintf(stdout,"new dna_align:  ");
		for (int j=0; j<result_length;j++){
			fprintf(stdout,"%i",dna_align[j]);		
		}
		fprintf(stdout,"\n");
		fprintf(stdout,"new read_align: ")	;
		
		for (int j=0; j<result_length;j++){
			fprintf(stdout,"%i",read_align[j]);		
		}
		fprintf(stdout,"\n length=%i\n",result_length);	
	}
	
	if (myverbosity>=1)
	{
		
		fprintf(stdout,"Alignment with variants (valid=%i)has: %i mm %i gaps; %i<=intron<=%i; %i<=exon num_exon==%i\n",
				alignment_passed_filters,alignment_mm,alignment_gaps,min_intron_len,max_intron_len,min_exon_len,num_exons);
		
	}
	
	return used_variants ;
}


template <int myverbosity>
bool QPalma::determine_exons(std::vector<int> & exons, const std::string & dna, const std::vector<int> &positions, bool remapping, char strand, const int *s_align, const int *e_align, 
					 int & min_exon_len, int & max_intron_len, int & min_intron_len) const
{
	bool intron_location=false;
		
	bool alignment_valid=true ;
	int exon_start = -1;

	if (perform_extra_checks && !(dna.length()==positions.size()))
	{
		fprintf(stderr, "ERROR: len mismatch dna.length()=%ld != %ld=positions.size()\n", dna.length(), positions.size()) ; 
		assert(dna.length()==positions.size()) ;
	}
	
	for (size_t i = 0; i < dna.length(); i++) 
	{
		//Alignment goes over N sequence
		if (s_align[i] == 0 && positions[i]<0)
		{
			if (myverbosity>=2)
				fprintf(stdout, "Warning: Alignment goes over region boundary. Marking alignment as invalid\n") ;
			return false;
			
			//assert(0) ;
			
			//return false;
		}
		
		if (exon_start == -1 && s_align[i] == 0 && i<dna.length()-1) 
		{
			exon_start = i;
			continue;
		}
		if (exon_start!=-1 && i>0 && s_align[i]==0)
		{
			alignment_valid = remapping || (alignment_valid && ((positions[i-1]+1 == positions[i]) || (positions[i-1] == positions[i]+1))) ;
			if (!alignment_valid)
			{
				if (myverbosity>=3)
					fprintf(stdout, "Warning: alignment positions not consecutive\n") ;
				//return false;
			}
			
			intron_location=remapping && s_align[i]==0 && 
				((strand=='+' && positions[i-1]+1 != positions[i]) || 
				 (strand=='-' && positions[i-1] != positions[i]+1)) ;
			if (myverbosity>=4)
				fprintf(stdout, "pos[%i]=%i   pos[%i]=%i  valid=%i\n", (int)i-1, (int)positions[i-1], (int)i, (int)positions[i], alignment_valid) ;
		}
		
		if (exon_start!=-1 && (s_align[i]!=0 || i==dna.length()-1 || intron_location))
		{
			if (exons.size()>0)
			{
				int intron_len = positions[exon_start]-exons[exons.size()-1] -1 ;
				if (strand=='-')
					intron_len*=-1 ;
				if (intron_len>max_intron_len)
					max_intron_len=intron_len ;
				if (intron_len<min_intron_len)
					min_intron_len=intron_len ;
			}
			if (perform_extra_checks && positions[exon_start]<0)
			{
				fprintf(stderr, "ERROR: positions[exon_start=%i]=%i\n", exon_start, positions[exon_start]) ; // BUG-TODO
				alignment_valid=false ;
				//assert(positions[exon_start]>=0) ;
				return false ;
			}

			exons.push_back(positions[exon_start]) ;
			
			int exon_len;
			if (s_align[i]==0 && !intron_location)
			{
				if (perform_extra_checks && positions[i]<0)
				{
					fprintf(stderr, "ERROR: positions[%i]=%i\n", (int)i, (int)positions[i]) ; 
					alignment_valid=false ;
					assert(positions[i]>=0) ;
				}
				exons.push_back(positions[i]) ;
				//fprintf(stdout,"Exon: %i-%i\n",exon_start,i);
				exon_len=positions[i] - positions[exon_start] ;		
			}				
			else
			{					
				if (perform_extra_checks && positions[i-1]<0)
				{
					fprintf(stderr, "ERROR: positions[%i-1]=%i\n", (int)i, (int)positions[i-1]) ; 
					alignment_valid=false ;
					assert(positions[i-1]>=0) ;
				}
				exons.push_back(positions[i-1]) ;
				//fprintf(stdout,"Exon: %i-%i\n",exon_start,i-1);
				exon_len=positions[i-1] - positions[exon_start] ;		
			}
			
			
			if (strand=='-')
				exon_len*=-1 ;
			exon_len++ ;
			if (perform_extra_checks)
				assert(exon_len>=1) ;
			if (exon_len<min_exon_len)
				min_exon_len = exon_len ;
			if (intron_location)
				exon_start = i ;
			else
				exon_start = -1 ;
			
			intron_location=false;
			
			continue ;
		}
	}
	
	if (perform_extra_checks)
		assert(exon_start==-1) ;

	if (exons.size()==0)
		alignment_valid=false ;
	
	if (myverbosity>=2)
		fprintf(stdout, "alignment valid=%i (0 means that exons went over block boundaries)\n", alignment_valid) ;

	return alignment_valid ;
}

template<int myverbosity, bool discover_variants> 
int QPalma::determine_read_variants(Chromosome const &contig_idx, const int * s_align, const int* e_align, const int *dna_align, const int *est_align, const std::vector<int> & positions, 
									 const VariantMap & variants, std::vector<Variant> & align_variants, std::vector<int> & aligned_positions,
									 Read const & read, const std::string & read_string, const std::string & read_quality, std::string &read_anno, int est_len_p, int result_length, char strand, int ori,
									 int &alignment_matches, int &alignment_gaps, int &alignment_mismatches, int &alignment_qual_mismatches) const
{
	/*fprintf(stdout, "DNA: %d ", s_align[0]);
	  for (int i = 0; i < dna.length(); i++)
	  fprintf(stdout, "%i ", s_align[i]);
	  fprintf(stdout, "\n");
	  
	  fprintf(stdout, "EST: ");
	  for (int i = 0; i < est_len_p; i++)
	  fprintf(stdout, "%i ", e_align[i]);
	  fprintf(stdout, "\n");*/
	
	int start_offset=0 ;
	while (s_align[start_offset]==4)
		start_offset++ ;
	
	//fprintf(stdout, "start_offset=%i\n", start_offset) ;
	
	//int dna_align[result_length];
	//int est_align[result_length];
	char dna_align_str[result_length + 1];
	char est_align_str[result_length + 1];
	
	//alignment.getAlignmentArrays(dna_align, est_align);
	
	
	//	int alignment_matches = 0;
	//int alignment_gaps = 0;
	//int alignment_mismatches = 0 ;
	
	//std::string read_anno = std::string("");
	char map[9] = "-ACGTN*+";

	// check for possibly wrong indels
	bool possibly_wrong_indel=false ;
	if (myverbosity>=1)
		for (int i = 0; i < result_length; i++)
		{
			const int window=3 ;
			if (dna_align[i]==0 && i>=window)
			{
				for (int j=i-window; j<i+window && j<result_length; j++)
					if (est_align[j]==0)
					{
						possibly_wrong_indel=true ;
						fprintf(stdout, "Warning: possibly wrong indel dnapos=%i, rnapos=%i (ori=%i, strand=%c)\n", i, j, ori, strand) ;
					}
			}
		}
	
	if (myverbosity>=1 || possibly_wrong_indel)
	{
		fprintf(stdout, "DNA: ");
		for (int i = 0; i < result_length; i++)
			fprintf(stdout, "%c", map[dna_align[i]]); 
		fprintf(stdout, "\n");
		
		fprintf(stdout, "EST: ");
		for (int i = 0; i < result_length; i++)
			fprintf(stdout, "%c", map[est_align[i]]);
		fprintf(stdout, "\n");
	}
	
	int read_pos=0, dna_pos=0 ;
	int est_gap_start = -1 ;
	int est_gap_end = -1 ;
	int dna_gap_start = -1 ;
	int dna_gap_start_i = -1 ;
	std::string dna_gap = "" ;
	std::string est_gap = "" ;
	
	
	for (int i = 0; i < result_length; i++) 
	{
		if (perform_extra_checks)
			assert(dna_align[i]>=0 && dna_align[i]<=6);
		dna_align_str[i] = map[dna_align[i]];
		
		if (perform_extra_checks)
		{
			if (!(est_align[i]>=0 && est_align[i]<=7))
			{
				fprintf(stderr, "ERROR: est_align[%i]=%i=%c (should be between 0 and 7)\n", i, est_align[i], map[est_align[i]]) ; 
				assert(est_align[i]>=0 && est_align[i]<=7) ;
				return -1 ;
			}
		}
		
		est_align_str[i] = map[est_align[i]];
		
		if (myverbosity>=3)
		{
			if (strand=='+')
				fprintf(stdout, "%i-%i\t%i (%i, %c)\t%i\t%s\n", dna_align[i], est_align[i], dna_pos, 
						positions[start_offset+dna_pos]+1, positions[start_offset+dna_pos]+1>=0?contig_idx[positions[start_offset+dna_pos]+1]:' ', read_pos, read_anno.c_str()) ;
			else
				fprintf(stdout, "%i-%i\t%i (%i, %c)\t%i\t%s\n", dna_align[i], est_align[i], dna_pos, 
						positions[start_offset+dna_pos]+1, positions[start_offset+dna_pos]+1>=0?contig_idx[positions[start_offset+dna_pos]+1]:' ', read_pos, read_anno.c_str()) ;
		}
		
		if (discover_variants)
		{
			if ( est_align[i] !=0 && est_gap_start!=-1)
			{
				if (read_pos>=_config.REPORT_INDEL_TERMINAL_DIST && (int)read_string.length()-read_pos>=_config.REPORT_INDEL_TERMINAL_DIST)
				{
					if (strand=='-')
					{
						est_gap=reverse(complement(est_gap)) ;
					}
					if (myverbosity>=2)
						fprintf(stdout, "deletion: %i:%i (%s) %c/%i\n%s %i\n", positions[start_offset+est_gap_start], positions[start_offset+est_gap_end], est_gap.c_str(), strand, ori, read_anno.c_str(), read_pos) ;
					//assert(abs(positions[start_offset+est_gap_end]-positions[start_offset+est_gap_start])+1==(int)est_gap.size()) ;
					if (positions[start_offset+est_gap_start] < positions[start_offset+est_gap_end])
						variants.report_del_variant(align_variants, contig_idx, positions[start_offset+est_gap_start], est_gap.size(), est_gap, read.id(),
													ori==0?read_pos:est_len_p-1-read_pos, est_len_p) ;
					else
						variants.report_del_variant(align_variants, contig_idx, positions[start_offset+est_gap_end], est_gap.size(), est_gap, read.id(),
													ori==0?read_pos:est_len_p-1-read_pos, est_len_p) ;
				}
				est_gap_start=-1 ;
				est_gap_end=-1 ;
				est_gap="" ;
			}
			if ( dna_align[i] !=0 && dna_gap_start!=-1)
			{
				if (read_pos>=_config.REPORT_INDEL_TERMINAL_DIST && (int)read_string.length()-read_pos>=_config.REPORT_INDEL_TERMINAL_DIST &&
					read_pos-(int)dna_gap.size()>=_config.REPORT_INDEL_TERMINAL_DIST && (int)read_string.length()-(read_pos-(int)dna_gap.size())>=_config.REPORT_INDEL_TERMINAL_DIST)
				{
					std::string flank("NN") ;
					flank[0]=map[dna_align[dna_gap_start_i-1]] ;
					flank[1]=map[dna_align[i]] ;
					
					if (strand=='-')
					{
						dna_gap_start-=2 ; // why is this not symmetric????
						dna_gap=reverse(complement(dna_gap)) ;
						flank=reverse(complement(flank)) ;
					}
					else
						dna_gap_start+=1 ;
					
					
					if (myverbosity>=2)
						fprintf(stdout, "insertion: %i (%s) %c/%i %i (%i) flank=%s\n%s %i\n", positions[start_offset+dna_gap_start], dna_gap.c_str(), strand, ori, 
								i, dna_pos, flank.c_str(), read_anno.c_str(), read_pos) ;
					variants.report_ins_variant(align_variants, contig_idx, positions[start_offset+dna_gap_start]-1, dna_gap.size(), dna_gap, read.id(),
												ori==0?read_pos:est_len_p-1-read_pos, est_len_p, flank.c_str()) ;
				}
				dna_gap_start=-1 ;
				dna_gap="" ;
			}
		}
		
		
		if (est_align[i]!=0 && est_align[i]!=6 && est_align[i]!=7 && dna_align[i]!=0)
		{
			if (est_align[i]==dna_align[i])
			{
				read_anno.push_back(map[est_align[i]]) ;
				if (perform_extra_checks && map[est_align[i]]!='N')
					assert(map[est_align[i]]==read_string[read_pos]) ;

				
				if (start_offset+dna_pos<(int)positions.size())
					aligned_positions.push_back(positions[start_offset+dna_pos]) ;
				else
					fprintf(stderr, "ERROR: start_offset (%i) + dna_pos (%i) >= positions.size() (%ld)\n", start_offset, dna_pos, positions.size()) ;
				
				read_pos++ ;
			}
			else
			{
				read_anno.push_back('[') ;
				if (strand=='+')
				{
					read_anno.push_back(map[dna_align[i]]) ;
					read_anno.push_back(map[est_align[i]]) ;
				} else
				{
					read_anno.push_back(map[est_align[i]]) ;
					read_anno.push_back(map[dna_align[i]]) ;
				}
				read_anno.push_back(']');
				if (map[dna_align[i]]!='N' && est_align[i]!=7 && map[est_align[i]]!='N')
				{ 
					alignment_mismatches++ ;
					alignment_qual_mismatches += read_quality[read_pos]-read.get_quality_offset() ; 
					if (perform_extra_checks)
						assert(map[est_align[i]]==read_string[read_pos]) ;
					
					if (discover_variants && start_offset+dna_pos<(int)positions.size())
						if (read_pos>=_config.REPORT_SNP_TERMINAL_DIST && (int)read_string.length()-read_pos>=_config.REPORT_SNP_TERMINAL_DIST)
						{
							if (myverbosity>=2)
								fprintf(stdout, "SNP: %i (%c->%c) %c/%i\n%s %i\n", positions[start_offset+dna_pos], map[dna_align[i]], map[est_align[i]], strand, ori, read_anno.c_str(), read_pos) ; 
							if (strand=='+')
								variants.report_SNP_variant(align_variants, contig_idx, positions[start_offset+dna_pos], map[dna_align[i]], map[est_align[i]], read.id(),
															ori==0?read_pos:est_len_p-1-read_pos, est_len_p) ;
							else
							{
								std::string dna_letter ;
								dna_letter+=map[dna_align[i]] ;
								dna_letter=complement(dna_letter) ;
								
								std::string est_letter ;
								est_letter+=map[est_align[i]] ;
								est_letter=complement(est_letter) ;
								
								variants.report_SNP_variant(align_variants, contig_idx, positions[start_offset+dna_pos], dna_letter[0], est_letter[0], read.id(),
															ori==0?read_pos:est_len_p-1-read_pos, est_len_p) ;
							}
						}
				}
				
				read_pos++ ;
			}
			alignment_matches += (est_align[i]==dna_align[i]) ;
			dna_pos++ ;
		}
		else if ( est_align[i]==0 )
		{
			if (discover_variants)
			{
				if (est_gap_start==-1)
					est_gap_start = dna_pos ;
				est_gap_end = dna_pos ;
				est_gap+=map[dna_align[i]] ;
			}
			
			if (dna_align[i]!=0)
			{
				read_anno.push_back('[') ;
				if (strand=='+')
				{
					read_anno.push_back(map[dna_align[i]]) ;
					read_anno.push_back('-') ;
				}
				else
				{
					read_anno.push_back('-') ;
					read_anno.push_back(map[dna_align[i]]) ;
				}
				read_anno.push_back(']');
			}
			alignment_gaps++ ;
			dna_pos++ ;
		}
		else if ( dna_align[i]==0)
		{
			if (discover_variants)
			{
				if (dna_gap_start==-1)
				{
					dna_gap_start = dna_pos ;
					dna_gap_start_i = i ;
				}
				dna_gap+=map[est_align[i]] ;
			}
			
			read_anno.push_back('[') ;
			if (strand=='+')
			{
				read_anno.push_back('-') ;
				read_anno.push_back(map[est_align[i]]) ;
			}
			else
			{
				read_anno.push_back(map[est_align[i]]) ;
				read_anno.push_back('-') ;
			}
			read_anno.push_back(']');
			alignment_gaps++;
			if (map[est_align[i]]!='N')
				alignment_qual_mismatches += read_quality[read_pos]-read.get_quality_offset() ; 
			if (perform_extra_checks && map[est_align[i]]!='N')
				assert(map[est_align[i]]==read_string[read_pos]) ;
			read_pos++ ;
		} 
		else if ((est_align[i]==6||est_align[i]==7) && dna_align[i]!=0 && dna_align[i]!=6)
			dna_pos++ ;
		else
		{
			fprintf(stdout, "lost case\n") ;
			assert(0) ;
		}
	}
	if (perform_extra_checks)
	{
		if (!(read_pos==(int)read.length()))
		{
			fprintf(stderr, "ERROR: len mismatch: %i != %i\n", read_pos, (int)read.length()) ;
			assert(read_pos==(int)read.length()) ;
			return -1 ;
		}
		//if (read.length()*4<read_anno.length())
		//fprintf(stderr, "ERROR: read anno longer than 4*read_length: %i*4<%i\n", (int)read.length(), (int)read_anno.length()) ;
	}
	
	dna_align_str[result_length] = 0;
	est_align_str[result_length] = 0;
	
	if (myverbosity >= 3) {
		fprintf(stdout, "# DNA : %s\n", dna_align_str);
		fprintf(stdout, "# Read: %s\n", est_align_str);
	}

	if (strand=='-')
		read_anno=reverse(complement(read_anno)) ;

	return  0 ;
}

template <int myverbosity>
int QPalma::get_splice_predictions(std::vector<region_t *> &current_regions, Chromosome const &contig_idx, Read const &read, JunctionMap &annotatedjunctions, double* prb, char* est,
								   std::vector<size_t> & cum_length, 
								   int &a_len, double *& acceptor, int &d_len, double*& donor, std::string & dna, bool non_consensus_search, char strand, int ori, bool remapping)  const
{
	/* initialize acceptor and donor tables */
	d_len = dna.length() ;
	donor = NULL ;
	try
	{
		donor = new double[d_len] ;
	}
	catch (std::bad_alloc&)
	{
		fprintf(stderr, "[perform_alignment] Could not allocate memory (_read.id() = %s; donor=%i)\n", read.id(), d_len);
		return -1;
	}

	for (int i = 0; i < d_len; i++)
		donor[i] = -ALMOST_INFINITY;

	a_len = dna.length();
	acceptor = NULL ;
	try
	{
		acceptor = new double[a_len] ;
	}
	catch (std::bad_alloc&)
	{
		fprintf(stderr,	"[perform_alignment] Could not allocate memory (_read.id() = %s; acceptor)\n", read.id());
		delete[] donor;
		return -1;
	}
	for (int i = 0; i < a_len; i++)
		acceptor[i] = -ALMOST_INFINITY;

	const int num_intervals = current_regions.size();
	int interval_matrix[num_intervals * 2];
	cum_length.resize(num_intervals + 1, 0);

	cum_length[0] = 0;
	for (int i = 0; i < num_intervals; i++) {
		interval_matrix[i * 2 + 0] = current_regions[i]->start;
		interval_matrix[i * 2 + 1] = current_regions[i]->end;
		cum_length[i + 1] = cum_length[i] + current_regions[i]->end - current_regions[i]->start;
	}

	if (perform_extra_checks && !(cum_length[num_intervals]==dna.length()))
	{
		fprintf(stderr, "ERROR: cum_length[num_intervals]=%ld != %ld=dna.length()\n", cum_length[num_intervals], dna.length()) ; // BUG-TODO
		delete[] acceptor ;
		delete[] donor ;
		return -1 ;
	}
	assert(cum_length[num_intervals]==dna.length()) ;

	if (!_config.NO_SPLICE_PREDICTIONS)
	{
		if (!remapping)
		{
			if (_config.ACC_FILES.length()>0 && _config.DON_FILES.length()>0 &&  !_config.NO_SPLICE_DISCOVERY )
			{
				/* query acceptor and donor scores */
				IntervalQuery iq;
				char basename[1000];
				const int num_scores = 1;
				const char *score_names[num_scores] = { "Conf_cum" };
	  

				// acc
				sprintf(basename, _config.ACC_FILES.c_str(), contig_idx.nr()+1, strand) ;
	    
				int acc_size = iq.query(basename, (char**) score_names, num_scores,
										interval_matrix, num_intervals);
				if (acc_size < 0) 
				{
					delete[] acceptor;
					delete[] donor ;
					return -1;
				}
	    
				int* acc_pos = NULL ;
				int* acc_index = NULL ;
				double* acc_score = NULL ;
				try
				{
					acc_pos = new int[acc_size] ;
					acc_index = new int[acc_size] ;
					acc_score = new double[acc_size] ;
				}
				catch (std::bad_alloc&)
				{
					fprintf(stderr, "[perform_alignment] Could not allocate memory (_read.id() = %s)\n", read.id());
					delete[] donor;
					delete[] acc_pos;
					delete[] acc_index;
					delete[] acc_score;
					return -1;
				}
	    
				iq.getResults(acc_pos, acc_index, acc_score);
				iq.cleanup();

				// don
				sprintf(basename, _config.DON_FILES.c_str(), contig_idx.nr()+1, strand) ;
	    
				int don_size = iq.query(basename, (char**) score_names, num_scores,
										interval_matrix, num_intervals);
				if (don_size < 0) {
					/* cleanup */
					delete[] acceptor;
					delete[] donor;
					delete[] acc_pos;
					delete[] acc_index;
					delete[] acc_score;
					return -1;
				}
	  
				int* don_pos = NULL ;
				int* don_index = NULL ;
				double* don_score = NULL ;
				try
				{
					don_pos = new int[don_size] ;
					don_index = new int[don_size] ;
					don_score = new double[don_size] ;
				}
				catch (std::bad_alloc&)
				{
					fprintf(stderr, "[perform_alignment] Could not allocate memory (_read.id() = %s)\n", read.id());
					/* cleanup */
					delete[] donor ;
					delete[] acc_pos ;
					delete[] acc_index ;
					delete[] acc_score ;
					delete[] don_pos ;
					delete[] don_index ;
					delete[] don_score ;
					return -1;
				}
				iq.getResults(don_pos, don_index, don_score);
				iq.cleanup();


				size_t j;
				size_t nbr_regions=current_regions.size();
				j=0;
	  
				//Splice sites and current regions are sorted by start position
				if (strand=='+'){
					for (int i=0; i<acc_size; i++){
						acc_pos[i]-=2 ; // 1-based -> 0-based  and shift by 1
						//fprintf(stdout,"(%i %i %i %i) ",acc_pos[i],current_regions[j]->start,current_regions[j]->end,j);
						while (j<nbr_regions-1 && acc_pos[i]>=current_regions[j]->end){
							j++;
						}
	      
						//fprintf(stdout,"c(%i %i %i %i) ",acc_pos[i],current_regions[j]->start,current_regions[j]->end,j);
						if (acc_pos[i]>=current_regions[j]->start && acc_pos[i]<current_regions[j]->end)
							acceptor[ acc_pos[i] - current_regions[j]->start + cum_length[j] ] = acc_score[i] ;  
					}
					j=0;
					for (int i = 0; i < don_size; i++) {
						don_pos[i] -= 1; // 1-based -> 0-based
						//fprintf(stdout,"b(%i %i %i %i) ",don_pos[i],current_regions[j]->start,current_regions[j]->end,j);
						while (j<nbr_regions-1 && don_pos[i] >= current_regions[j]->end)
							j++;
						//fprintf(stdout,"(%i %i %i %i) ",don_pos[i],current_regions[j]->start,current_regions[j]->end,j);
						if (don_pos[i]>=current_regions[j]->start && don_pos[i] < current_regions[j]->end){
							donor[don_pos[i] - current_regions[j]->start + cum_length[j]]= don_score[i];
						}
					}
				}
				
				//Negative strand: directly modify the relative coordinate of the splice site in reversed and complemented dna
				else{
					for (int i=0; i<acc_size; i++){
						while (j<nbr_regions-1 && acc_pos[i]>=current_regions[j]->end){
							j++;
						}

						//fprintf(stdout,"(%i %i %i %i) ",acc_pos[i],current_regions[j]->start,current_regions[j]->end,j);
						if (acc_pos[i]>=current_regions[j]->start && acc_pos[i]<current_regions[j]->end){
							acceptor[a_len-1-(acc_pos[i] - current_regions[j]->start + cum_length[j]) ] = acc_score[i] ;
						}
					}
					//	  fprintf(stdout,"\n");
					j=0;
					for (int i = 0; i < don_size; i++) {
						don_pos[i] -= 1; // 1-based -> 0-based
						//fprintf(stdout,"b(%i %i %i %i) ",don_pos[i],current_regions[j]->start,current_regions[j]->end,j);
						while (j<nbr_regions-1 && don_pos[i] >= current_regions[j]->end)
							j++;
						//fprintf(stdout,"(%i %i %i %i) ",don_pos[i],current_regions[j]->start,current_regions[j]->end,j);
						if (don_pos[i]>=current_regions[j]->start && don_pos[i] < current_regions[j]->end){
							donor[d_len-1 -(don_pos[i] - current_regions[j]->start + cum_length[j])]= don_score[i];
						}
					}
				}

				if (myverbosity>=3)
					fprintf(stdout, "got %i acceptors and %i donors\n", acc_size, don_size) ;

				/* cleanup */
				delete[] acc_pos ;
				delete[] acc_index ;
				delete[] acc_score ;
				delete[] don_pos ;
				delete[] don_index ;
				delete[] don_score ;
			}
			
			//Use splice sites from the annotation
			if (_config.SCORE_ANNOTATED_SPLICE_SITES) // this is done, even if _config.NO_SPLICE_DISCOVERY==true
			{
				//fprintf(stdout,"Use annotated splice sites\n");
				std::vector<int> acc_pos_vec;
				std::vector<int> don_pos_vec;

				//Get splice sites for different intervals				
				for (int i = 0; i < num_intervals; i++) {
					get_annotated_splice_sites(acc_pos_vec, don_pos_vec, annotatedjunctions, current_regions[i]->start, current_regions[i]->end, contig_idx.nr(), strand);

					//Transform coordinates
					if (strand=='+'){
						for (int k=0; k<(int)acc_pos_vec.size(); k++){
							acceptor[ acc_pos_vec[k] - current_regions[i]->start + cum_length[i]] = 1.0 ;  
						}
						for (int k=0; k<(int)don_pos_vec.size(); k++) {
							donor[don_pos_vec[k] - current_regions[i]->start + cum_length[i]]= 1.0;
						}
					}

					else{
						for (int k=0; k<(int)acc_pos_vec.size(); k++){
							acceptor[a_len-1-(acc_pos_vec[k] - current_regions[i]->start + cum_length[i])] = 1.0 ;
						}
						for (int k=0; k<(int)don_pos_vec.size(); k++) {
							donor[d_len-1 -(don_pos_vec[k] - current_regions[i]->start + cum_length[i])]= 1.0;
						}
					}
					
					acc_pos_vec.clear();
					don_pos_vec.clear();
					
				}
			}
		}
	}
	return 0 ;
}





int QPalma::postprocess_splice_predictions(Chromosome const &contig_idx,
										   const int a_len, double *acceptor, const int d_len, double* donor, const std::string & dna, bool non_consensus_search, char strand, int ori,
										   bool check_ss, int region_start, int region_end) const
{
	if (_config.NO_SPLICE_DISCOVERY)
		return 0 ;
		
	int match = 0, num = 0;
	for (int i = 2; i < d_len - 2; i++, num++)
	{
		bool is_ss = false ;
		for (unsigned int j=0; j<_config.DON_CONSENSUS.size(); j++)
			if (IUPACMatch(dna[i], _config.DON_CONSENSUS[j][0]) && IUPACMatch(dna[i+1], _config.DON_CONSENSUS[j][1]))
			{
				is_ss = true ;
				break ;
			}
		if (non_consensus_search && !is_ss) 
		{
			if (dna[i]!='N' && dna[i+1]!='N')
				donor[i] = NON_CONSENSUS_SCORE ;
			match += 1 ;//(donor[i] > -ALMOST_INFINITY);
		} 
		else
			if (is_ss) 
			{
				// if no splice predictions, put score 0 for allowed donor splice consensus (e.g. 'GT/C')
				if (_config.NO_SPLICE_PREDICTIONS)
					donor[i] = 0.5 ;
				if ((_config.SCORE_ANNOTATED_SPLICE_SITES||_config.USE_VARIANTS) && donor[i]==-ALMOST_INFINITY)
					donor[i] = NON_CONSENSUS_SCORE ;
				match += (donor[i] > -ALMOST_INFINITY || _config.SCORE_ANNOTATED_SPLICE_SITES || _config.USE_VARIANTS);
			} 
			else 
			{
				match += (_config.SCORE_ANNOTATED_SPLICE_SITES || _config.USE_VARIANTS) || (donor[i] <= -ALMOST_INFINITY)  || dna[i]=='N' || dna[i+1]=='N' ;
				donor[i] = -ALMOST_INFINITY;
			}
	}
	
	if (check_ss && match < num * 0.9)
		fprintf(stdout,
				"Warning: donor predictions do not match genome positions (match=%i  num=%i  strand=%c  ori=%c  chr=%s  start=%i  end=%i)\n",
				match, num, strand, ori == 0 ? '+' : '-', contig_idx.desc(),
				region_start, region_end);
	
	
	match = 0; num = 0;
	for (int i = 2; i < a_len - 2; i++, num++)
	{
		bool is_ss = false ;
		for (unsigned int j=0; j < _config.ACC_CONSENSUS.size(); j++)
			if (i>0 && IUPACMatch(dna[i-1], _config.ACC_CONSENSUS[j][0]) && IUPACMatch(dna[i], _config.ACC_CONSENSUS[j][1]))
			{
				is_ss = true ;
				break ;
			}
		
		if (non_consensus_search && !is_ss)
		{
			if (dna[i]!='N' && dna[i-1]!='N')
				acceptor[i] = NON_CONSENSUS_SCORE ;
			match += 1;
		} 
		else
			if (is_ss)
			{
				if (_config.NO_SPLICE_PREDICTIONS)
					acceptor[i] = 0.5 ;
				if ((_config.SCORE_ANNOTATED_SPLICE_SITES || _config.USE_VARIANTS) && acceptor[i]==-ALMOST_INFINITY)
					acceptor[i] = NON_CONSENSUS_SCORE ;
				match += (acceptor[i] > -ALMOST_INFINITY || _config.SCORE_ANNOTATED_SPLICE_SITES || _config.USE_VARIANTS);
			} 
			else 
			{
				//if (acceptor[i]>-ALMOST_INFINITY && dna[i]!='N' && dna[i-1]!='N')
				//	fprintf(stdout, "acc over %i %f  %c%c\n", i, acceptor[i], dna[i-1], dna[i]) ;
				match += _config.SCORE_ANNOTATED_SPLICE_SITES || _config.USE_VARIANTS || (acceptor[i] <= -ALMOST_INFINITY) || dna[i]=='N' || dna[i-1]=='N' ;
				acceptor[i] = -ALMOST_INFINITY;
			}
	}
	
	if (check_ss && match<num*0.9)
		fprintf(stdout, "Warning: acceptor predictions do not match genome positions (match=%i  num=%i  strand=%c  ori=%c  chr=%s  start=%i  end=%i, %i)\n", 
				match, num, strand, ori==0 ? '+' : '-' , contig_idx.desc(), region_start, region_end, non_consensus_search) ;

	return 0 ;
}

int QPalma::transform_splice_predictions(const int a_len, double *acceptor, const int d_len, double* donor) const
{
	/* apply donor and acceptor plifs (i.e., transform the splice site prediction score to a score that can used in QPalma dynamic program) */
	for (int i = 0; i < d_len; i++)
		if (donor[i] > -ALMOST_INFINITY)
		{
			if (donor[i]==NON_CONSENSUS_SCORE)
			{
				donor[i]=0.1 ; // low donor probability  -> correspondingly low score
				if (not _config.NO_QPALMA)
					donor[i] = lookup_penalty(&alignment_parameters->d, 0, &donor[i]);
			}
			else
			{
				if (not _config.NO_QPALMA)
					donor[i] = lookup_penalty(&alignment_parameters->d, 0, &donor[i]);
			}			
		}
	for (int i = 0; i < a_len; i++)
		if (acceptor[i] > -ALMOST_INFINITY)
		{
			if (acceptor[i]==NON_CONSENSUS_SCORE)
			{
				acceptor[i]=0.1 ; // low acceptor probability  -> correspondingly low score
				if (not _config.NO_QPALMA)
					acceptor[i] = lookup_penalty(&alignment_parameters->a, 0, &acceptor[i]);
			}
			else
			{
				if (not _config.NO_QPALMA)
					acceptor[i] = lookup_penalty(&alignment_parameters->a, 0, &acceptor[i]);
			}
		}
	return 0 ;
}

int QPalma::construct_intron_strings(ALIGNMENT * aln, Chromosome const &contig_idx, const std::vector<int> & exons, 
									 const char strand, const int ori, bool & non_consensus_alignment, const bool remapping, const bool non_consensus_search) const
{
	for (unsigned int ne=0; ne<(exons.size()/2)-1; ne++)
	{
		bool non_consensus_intron = false ;
		
		int istart=exons[ne*2+1] ;
		int istop =exons[(ne+1)*2] ;
		
		if (perform_extra_checks && (istop<4 || istart<2))
		{
			fprintf(stderr, "ERROR: istart=%i, istop=%i\n", istart, istop) ; 
			assert(istop>=4 && istart>=2) ;
			return -1 ;
		}
		
		char buf[1000] ;
		if ((ori==0 && strand=='+') || (ori==1 && strand=='+'))
		{
			sprintf(buf, "%s:%i:%i:%c%c|%c%c%c%c..%c%c%c%c|%c%c", contig_idx.desc(), istart, istop-1, 
					contig_idx[istart-2], contig_idx[istart-1], contig_idx[istart], contig_idx[istart+1],  contig_idx[istart+2], contig_idx[istart+3], 
					contig_idx[istop-4], contig_idx[istop-3], contig_idx[istop-2], contig_idx[istop-1], contig_idx[istop], contig_idx[istop+1]) ;
			//fprintf(stdout, "intron+: %s\n", buf) ;
			
			bool is_consensus_don = false ;
			for (unsigned int j=0; j < _config.DON_CONSENSUS.size(); j++)
				if (IUPACMatch(contig_idx[istart], _config.DON_CONSENSUS[j][0]) && IUPACMatch(contig_idx[istart+1], _config.DON_CONSENSUS[j][1]))
				{
					is_consensus_don = true ;
					break ;
				}
			bool is_consensus_acc = false ;
			for (unsigned int j=0; j < _config.ACC_CONSENSUS.size(); j++)
				if (IUPACMatch(contig_idx[istop-2], _config.ACC_CONSENSUS[j][0]) && IUPACMatch(contig_idx[istop-1], _config.ACC_CONSENSUS[j][1]))
				{
					is_consensus_acc = true ;
					break ;
				}
			if (!is_consensus_don || !is_consensus_acc)
			{
				non_consensus_intron=true ;
				non_consensus_alignment=true ;
			}
			if (!non_consensus_search && !remapping && !_config.USE_VARIANTS)//  // && final_variants.size()==0) // TODO
				if (non_consensus_intron && perform_extra_checks)
					fprintf(stderr, "Warning: acc non_consensus_intron=true\n") ;
			//assert(!non_consensus_intron) ;
		}
		else
		{
			sprintf(buf, "%s:%i:%i:%c%c|%c%c%c%c..%c%c%c%c|%c%c", contig_idx.desc(), istop-1, istart, 
					complement(contig_idx[istop+1]), complement(contig_idx[istop]), complement(contig_idx[istop-1]), complement(contig_idx[istop-2]), complement(contig_idx[istop-3]), complement(contig_idx[istop-4]), 
					complement(contig_idx[istart+3]), complement(contig_idx[istart+2]), complement(contig_idx[istart+1]), complement(contig_idx[istart-0]), complement(contig_idx[istart-1]), complement(contig_idx[istart-2])) ;
			//fprintf(stdout, "intron-: %s\n", buf) ;
			
			bool is_consensus_don = false ;
			for (unsigned int j=0; j < _config.DON_CONSENSUS.size(); j++)
				if (IUPACMatch(complement(contig_idx[istop-1]), _config.DON_CONSENSUS[j][0]) && IUPACMatch(complement(contig_idx[istop-2]), _config.DON_CONSENSUS[j][1]))
				{
					is_consensus_don = true ;
					break ;
				}
			bool is_consensus_acc = false ;
			for (unsigned int j=0; j < _config.ACC_CONSENSUS.size(); j++)
				if (IUPACMatch(complement(contig_idx[istart+1]), _config.ACC_CONSENSUS[j][0]) && IUPACMatch(complement(contig_idx[istart-0]), _config.ACC_CONSENSUS[j][1]))
				{
					is_consensus_acc = true ;
					break ;
				}
			if (!is_consensus_don || !is_consensus_acc)
			{
				non_consensus_intron=true ;
				non_consensus_alignment=true ;
			}
			if (!non_consensus_search && !remapping && !_config.USE_VARIANTS) // && final_variants.size()==0) // TODO
				if (non_consensus_intron && perform_extra_checks) 
					fprintf(stderr, "Warning: acc non_consensus_intron=true\n") ;
			//assert(!non_consensus_intron) ;
			
		}
		aln->non_consensus_intron.push_back(non_consensus_intron) ;
		aln->intron_consensus.push_back(buf) ; 
	}
	return 0 ;
}

void clear_variant_cache(std::vector<variant_cache_t *> & variant_cache)
{
	for (int i=0;i<(int)variant_cache.size();i++){
		if (variant_cache[i] == NULL)
			continue;
		variant_cache[i]->end_positions.clear();
		variant_cache[i]->id_dels.clear();
		variant_cache[i]->id_snps.clear();
		variant_cache[i]->snps.clear();
		delete variant_cache[i];
		variant_cache[i]=NULL;
	}
	
	variant_cache.clear();
}



template<int myverbosity, bool discover_variants, bool remapping, bool do_timing> 
int QPalma::perform_alignment(Result &result, Hits &readMappings, std::string &read_string, std::string &read_quality, std::string &dna, 
							  std::vector<region_t *> &current_regions, std::vector<int> &positions, 
							  Chromosome const &contig_idx, char strand, int ori, int hit_read, int hit_dna, int hit_length, bool non_consensus_search, ALIGNMENT *& aln,  
							  JunctionMap &annotatedjunctions, const VariantMap & variants, std::vector<Variant> & variant_list, std::map<int, int> & variant_pos) const
// ori = read orientation
// strand = dna strand/orientation
//hit_read, hit_dna, hit_length: starting (real) positions and length of the hit for starting the alignment

{
	Read const &read(result._read);
	bool isunspliced = false; // is only true, if a valid alignment was reported and it is unspliced

	if (myverbosity>=2)
	{
		for (size_t i=0; i<current_regions.size(); i++)
		{
			char c=' ' ;
			if (hit_dna>=current_regions[i]->start && hit_dna<=current_regions[i]->end)
				c='*' ;
			fprintf(stdout, "[perform_alignment] region %i: %i - %i (from_map=%i)\t%c\n", (int)i, current_regions[i]->start, current_regions[i]->end, current_regions[i]->from_map, c) ;
		}
	}
	
	if (myverbosity>=2)
		fprintf(stdout, "average alignment length %lint\n", _stats.qpalma_total_dna_length/(_stats.qpalma_total_alignments+1)) ;
	
	// TODO
	//current_regions[0]->strand = strand ;
	//current_regions[0]->chromosome = contig_idx ;

	if (myverbosity>=1)
		fprintf(stdout, "# [perform_alignment] performing alignment of read %s with %i regions and sequence of length %i, strand=%c, remapping=%i\n", 
				read.id(), (int)current_regions.size(), (int)dna.size(), strand, remapping) ;

	int nr_paths_p = 1;

	/*if (myverbosity >= 2) {
		for (size_t i = 0; i < current_regions.size(); i++) {
			fprintf(stdout, "# region %i: %i - %i (len=%i)\n", (int)i, current_regions[i]->start,
					current_regions[i]->end, current_regions[i]->end-current_regions[i]->start);
		}
		}*/

	/* setup read sequence and quality information */
	int est_len_p = read_string.length();
	//char* est = (char *) malloc(sizeof(char) * (est_len_p + 1));
	char est[est_len_p+1] ;

	strncpy(est, read_string.c_str(), est_len_p);
	est[est_len_p] = '\0';

	double prb[est_len_p] ;

	int offset_switched=0;
	for (size_t i = 0; i < (size_t) est_len_p; i++)
	{
		prb[i] = (read_quality[i] - alignment_parameters->quality_offset);
		if (prb[i]<-10 || prb[i]>70)
			fprintf(stderr, "prb[%i]=%f (offset=%i, %s, %s)\n", (int)i, prb[i], alignment_parameters->quality_offset, read_quality.c_str(), read.data()) ;

		if (_config.QPALMA_PRB_OFFSET_FIX)
		{
			if (alignment_parameters->quality_offset==33 && prb[i]>70)
			{
				fprintf(stderr, "setting prb offset from %i to %i\n", alignment_parameters->quality_offset, 64) ;
				_config.QPALMA_PRB_OFFSET_FIX=false ;
				alignment_parameters->quality_offset=64 ;
				read.set_quality_offset(64) ;
				offset_switched=1;
				
			}
			if (alignment_parameters->quality_offset==64 && prb[i]<-10)
			{
				fprintf(stderr, "setting prb offset from %i to %i\n", alignment_parameters->quality_offset, 33) ;
				_config.QPALMA_PRB_OFFSET_FIX=false ;
				alignment_parameters->quality_offset=33 ;
				read.set_quality_offset(33) ;
				offset_switched=1;
			}
		}
	}
	
	if (offset_switched==1){
		for (int i = 0; i < est_len_p; i++)
			prb[i] = (read_quality[i] - alignment_parameters->quality_offset);	
	}
	
	int num_qual_support = 0;
	if (alignment_parameters->num_qualityPlifs > 0)
		num_qual_support = alignment_parameters->qualityPlifs[0].len;
	//Alignment alignment(alignment_parameters->num_qualityPlifs, num_qual_support, true, myverbosity);
	Alignment alignment(alignment_parameters->num_qualityPlifs, num_qual_support, !(_config.NO_QPALMA), myverbosity);

	int a_len, d_len ;
	double * acceptor=NULL, *donor=NULL ;
	int hit_dna_converted ;
	{
		std::vector<size_t> cum_length ;
		int ret=get_splice_predictions<myverbosity>(current_regions, contig_idx, read, annotatedjunctions, prb, est,
													cum_length, a_len, acceptor, d_len, donor, dna, non_consensus_search, strand, ori, remapping) ;
		if (ret<0)
			return ret ;
	
		if (strand=='-')
		{
			dna = reverse(complement(dna)) ;
			positions = reverse(positions);
		}
		if (myverbosity>=4)
		{
			fprintf(stdout, "positions (%i):\n",(int)dna.length()) ;
			for (unsigned int i=0; i<positions.size(); i++)
				fprintf(stdout, "%i:%i ", i, positions[i]) ;
			fprintf(stdout, "\n") ;
		}
		
		if (myverbosity >= 3)
			fprintf(stdout, "# dna: %s\n", dna.c_str());
		if (myverbosity >= 3)
			fprintf(stdout, "# read: %s\n", read_string.c_str());
				
		// Convert real dna position of the hit into relative dna position 
		hit_dna_converted = convert_dna_position(hit_dna, cum_length, current_regions);
		if (strand=='-')
			hit_dna_converted=(int)dna.length()-1-hit_dna_converted;
		if (myverbosity >= 3)
			fprintf(stdout, "hit_dna_converted=%i (dna.length()=%ld, positions.size()=%ld)\n", hit_dna_converted, dna.length(), positions.size()) ;
	}
    
	TIME_CODE(clock_t start_time ;
	if (do_timing)
		start_time = clock();)

	bool remove_duplicate_scores = false;

    // create super-sequence and deletion list from variant list
	int seed_i;
	int seed_j;
	double best_match=alignment.init_seed_position (hit_read, hit_dna_converted, hit_length, seed_i, seed_j, est,
													est_len_p, (char*)dna.c_str(), (int) dna.length(), alignment_parameters->qualityPlifs, 
													alignment_parameters->matchmatrix, alignment_parameters->matchmatrix_dim[0]
													* alignment_parameters->matchmatrix_dim[1], prb);

	if (myverbosity >= 3)
		fprintf(stdout, "best_match=%lf\nhit_dna=%i, seed_j=%i\nhit_read=%i,seed_i=%i\n", best_match, hit_dna_converted, seed_j, hit_read, seed_i) ;
	
	std::vector<variant_cache_t *> variant_cache;
	
	std::vector<bool> ref_map;
	int pre_dna_size=dna.size() ;
	if (_config.USE_VARIANTS) // && variant_list.size()>0)
	{
		TIME_CODE(clock_t my_start_time ; )

		if (strand == '-')
		{
			//reverse super variants coordinates and sequences
			reverse_variants(variant_list, (int)dna.length());
		}
		try 
		{
			TIME_CODE(if (do_timing)
						  my_start_time = clock() ;)
			variant_cache = create_super_sequence_from_variants<myverbosity,false>(variant_list, dna, acceptor, a_len, donor, d_len, seed_j, ref_map) ;
		}
		catch (std::bad_alloc)
		{
			fprintf(stderr,	"[capture_hits] allocating memory for supersequence failed (dna.length()=%ld, variant_list.size()=%ld)\n", dna.length(), variant_list.size());
			//result.delete_regions();
			delete[] acceptor ;
			delete[] donor ;
			aln=NULL ;
			return -1;
		}
		
		TIME_CODE(if (do_timing)
					  _stats.variant_create_super_sequence_from_variants_time += clock() - my_start_time;)
	}
	if (myverbosity >= 3)
		fprintf(stdout, "after supersequence seed_j=%i\n", seed_j) ;
	if (dna.size()>100000)
		fprintf(stdout, "Warning: DNA sequence for alignment of read %s is long: dna.size()=%i (pre_dna_size=%i)\n", read.id(), (int)dna.size(), pre_dna_size) ;

	/* check whether we have scores for all donor and acceptor positions (first 10% of reads)*/
	if (!remapping)
	{
		postprocess_splice_predictions(contig_idx, a_len, acceptor, d_len, donor, dna, non_consensus_search, strand, ori, 
									   !_config.USE_VARIANTS, current_regions[0]->start, current_regions[current_regions.size()-1]->end) ;
		transform_splice_predictions(a_len, acceptor, d_len, donor);
	}

	if (myverbosity>=4)
	{
		fprintf(stdout, "positions (%i):\n",(int)dna.length()) ;
		for (unsigned int i=0; i<positions.size(); i++)
			fprintf(stdout, "%i:%i ", i, positions[i]) ;
		fprintf(stdout, "\n") ;
	}
	if (myverbosity >= 3)
		fprintf(stdout, "# dna: %s\n", dna.c_str());
	if (myverbosity >= 3)
		fprintf(stdout, "# read: %s\n", read_string.c_str());
	if (hit_read<0)
	{
		fprintf(stderr, "ERROR: hit_read=%i setting to 0 to recover\n", hit_read) ;
		//assert(0) ; BUG-TODO
		hit_read = 0 ;
	}

	if (perform_extra_checks && (hit_dna_converted<0 || hit_dna_converted>=(int)dna.size()))
	{
		fprintf(stderr, "ERROR: hit_dna_converted=%i, dnalen=%ld, strand=%c, ori=%i, remapping=%i\n", hit_dna_converted, dna.size(), strand, ori, remapping) ; // BUG-TODO
		assert(0) ;
		//result.delete_regions();
		delete[] acceptor ;
		delete[] donor ;
		delete[] aln ;
		aln=NULL ;
		return -1 ;
	}
	if (hit_dna_converted<0)
		hit_dna_converted=0 ;
	if (hit_dna_converted>=(int)dna.size())
		hit_dna_converted=dna.size()-1 ;
	
	if (perform_extra_checks && (seed_j<0 || seed_j>=(int)dna.size())) 
	{
		fprintf(stderr, "ERROR: seed_j=%i not in [0,%ld], strand=%c, ori=%i, remapping=%i\n", seed_j, dna.size(), strand, ori, remapping) ; // BUG-TODO
		//result.delete_regions();
		delete[] acceptor ;
		delete[] donor ;
		delete[] aln ;
		aln=NULL ;
		return -1 ;
	}
	if (seed_j<0)
		seed_j=0 ;
	if (seed_j>=(int)dna.size())
		seed_j=dna.size()-1 ;
	
	assert (hit_read >= 0);
	assert (hit_dna_converted >= 0);
	assert (seed_j >= 0);
	assert (seed_j >= 0 && seed_j<(int)dna.size());
	
	int qmm_value=MIN_NUM_MATCHES;
	if (non_consensus_search)
		qmm_value+= _config.MIN_NUM_MATCHES_PEN;

	{
		TIME_CODE(clock_t my_start_time ;
		if (do_timing)
			my_start_time = clock() ;)

		if ((_config.USE_VARIANTS && _config.IUPAC_SNPS) || _config.IUPAC_GENOME) // && (int)variant_list.size()>0)
            alignment.myalign_fast<true,true>(strand, contig_idx, nr_paths_p, (char*) dna.c_str(), (int) dna.length(), est,
                                              est_len_p, prb, alignment_parameters->h,
                                              alignment_parameters->matchmatrix,
                                              alignment_parameters->matchmatrix_dim[0]
                                              * alignment_parameters->matchmatrix_dim[1], donor, d_len,
                                              acceptor, a_len, alignment_parameters->qualityPlifs,
                                              remove_duplicate_scores, seed_i, seed_j, best_match, _config.SPLICED_MAX_INTRONS, _config.DP_MAX_DELETIONS,
                                              readMappings.topAlignments().get_max_gaps(), readMappings.topAlignments().get_max_mismatches(), readMappings.topAlignments().get_max_editops(),
                                              /* _config.NUM_GAPS, _config.NUM_MISMATCHES, readMappings.get_num_edit_ops(), */
                                              qmm_value, remapping,
                                              _config.USE_VARIANTS, _config.NO_GAP_END,_config.SPLICED_MIN_SEGMENT_LENGTH,_config.SPLICED_SHORTEST_INTRON_LENGTH,variant_cache);
        else if (_config.USE_VARIANTS) 
            alignment.myalign_fast<true,false>(strand, contig_idx, nr_paths_p, (char*) dna.c_str(), (int) dna.length(), est,
                                               est_len_p, prb, alignment_parameters->h,
                                               alignment_parameters->matchmatrix,
                                               alignment_parameters->matchmatrix_dim[0]
                                               * alignment_parameters->matchmatrix_dim[1], donor, d_len,
                                               acceptor, a_len, alignment_parameters->qualityPlifs,
                                               remove_duplicate_scores, seed_i, seed_j, best_match, _config.SPLICED_MAX_INTRONS, _config.DP_MAX_DELETIONS,
                                               readMappings.topAlignments().get_max_gaps(), readMappings.topAlignments().get_max_mismatches(), readMappings.topAlignments().get_max_editops(),
                                               /* _config.NUM_GAPS, _config.NUM_MISMATCHES, readMappings.get_num_edit_ops(), */
                                               qmm_value, remapping,
                                               _config.USE_VARIANTS, _config.NO_GAP_END,_config.SPLICED_MIN_SEGMENT_LENGTH,_config.SPLICED_SHORTEST_INTRON_LENGTH,variant_cache);
			
		else
			alignment.myalign_fast<false,false>(strand, contig_idx, nr_paths_p, (char*) dna.c_str(), (int) dna.length(), est,
												est_len_p, prb, alignment_parameters->h,
												alignment_parameters->matchmatrix,
												alignment_parameters->matchmatrix_dim[0]
												* alignment_parameters->matchmatrix_dim[1], donor, d_len,
												acceptor, a_len, alignment_parameters->qualityPlifs,
												remove_duplicate_scores, seed_i, seed_j, best_match,_config.SPLICED_MAX_INTRONS, _config.DP_MAX_DELETIONS,
												readMappings.topAlignments().get_max_gaps(), readMappings.topAlignments().get_max_mismatches(), readMappings.topAlignments().get_max_editops(),
												/* _config.NUM_GAPS, _config.NUM_MISMATCHES, readMappings.get_num_edit_ops(), */
												qmm_value, remapping,
												_config.USE_VARIANTS, _config.NO_GAP_END,_config.SPLICED_MIN_SEGMENT_LENGTH,_config.SPLICED_SHORTEST_INTRON_LENGTH,variant_cache);
		
		TIME_CODE(if (do_timing)
					  _stats.variant_myalign_fast_time += clock() - my_start_time;)
	}

	TIME_CODE(static pthread_mutex_t clock_mutex = PTHREAD_MUTEX_INITIALIZER;
	pthread_mutex_lock( &clock_mutex) ;
	if (do_timing)
		_stats.qpalma_align_time += clock() - start_time;
			  pthread_mutex_unlock( &clock_mutex) ;)

	int s_len = dna.length() ;
	int *s_align=NULL;
	int *e_align=NULL ;
	s_align=new int[s_len]; 
	e_align=new int[est_len_p] ;

	int mmatrix_p[alignment_parameters->matchmatrix_dim[0]
				  * alignment_parameters->matchmatrix_dim[1]];
	double alignscore;
	double *qScores = NULL;

	alignment.getAlignmentResults(s_align, e_align, mmatrix_p, &alignscore, qScores);
	int result_length = alignment.getResultLength();
	int *dna_align=NULL;
	int *est_align=NULL;
	dna_align=new int[result_length+1]; 
	for (int i=0; i<result_length+1; i++) dna_align[i]=0 ;
	est_align=new int[result_length+1];
	for (int i=0; i<result_length+1; i++) est_align[i]=0 ;

	alignment.getAlignmentArrays(dna_align, est_align);

	std::vector<FoundVariant> found_variants =alignment.getVariants();
	
	std::vector<int> exons;
	exons.clear();

	int min_exon_len = read.length();
	int max_intron_len = 0;
	int min_intron_len = 10000000;

	bool alignment_valid=true ;

	int alignment_matches_ref = 0;
	int alignment_gaps_ref = 0;
	int alignment_mismatches_ref = 0 ;
	int alignment_qual_mismatches_ref = 0 ;
	std::string read_anno_ref = std::string("");
	std::vector<int> aligned_positions_ref ;
	bool alignment_passed_filters_ref = false ;
	std::vector<Variant> align_variants_ref ;

	int alignment_gaps_var = 0;
	int alignment_mismatches_var = 0 ;
	bool alignment_passed_filters_var ;

	int used_variants=0;
	//bool alignment_variants_valid=true;
	bool determine_exons_successful=false ;
	bool considered_variants = false ;
	
	//Alignment not found with less than _config.NUM_EDIT_OPS or _config.NUM_GAPS or _config.NUM_MISMATCHES
	if (result_length<est_len_p)
	{
	    alignment_valid=false;
	    
	    if (myverbosity>=2)
		{
			fprintf(stdout, "result_length=%i, est_len_p=%i\n", result_length, est_len_p) ;
			fprintf(stdout, "No alignment found\n") ;
		}
	}
	else
	{  	
		if (_config.USE_VARIANTS) // && variant_list.size()>0)
		{
			considered_variants=true ;

			TIME_CODE(clock_t my_start_time ;
			if (do_timing)
				my_start_time = clock() ;)
			
			used_variants = reconstruct_reference_alignment<myverbosity>(variant_list,found_variants, dna, ref_map, s_align, s_len, e_align, est_len_p, 
																		 dna_align, est_align, result_length, remapping, 
																		 alignment_passed_filters_var, alignment_gaps_var, alignment_mismatches_var,
																		 variant_cache, _config.REPORT_USED_VARIANTS, read.id()) ;

			TIME_CODE(if (do_timing)
						  _stats.variant_reconstruct_reference_alignment_time += clock() - my_start_time;) 
			
			for (unsigned int i=0; i< variant_list.size();i++)
			{
				recover_variants_on_ref(variant_list[i], positions, strand, est_len_p,contig_idx);
				
				if (variant_list[i].used_count >=1)
					variant_list[i].read_id = read.id();
				else
					variant_list[i].non_used_count += 1;
			}
			
			if (myverbosity>=1)
			{
				fprintf(stdout,"Number of variants used for the alignment of %s: %i\n",read.id(),used_variants);				
				for (unsigned int v=0; v<variant_list.size();v++){
					if (variant_list[v].used_count >=1)
						fprintf(stdout,"variant to report: type=%i, position=%i, read pos=%i read_id=%s used_count=%i non_used_count=%i\n",variant_list[v].type,variant_list[v].position,variant_list[v].read_pos, (char *)variant_list[v].read_id.c_str(),variant_list[v].used_count,variant_list[v].non_used_count);
				}
			}

			alignment_valid = alignment_valid && alignment_passed_filters_var ;
		}

		if (alignment_valid)
			determine_exons_successful = determine_exons<myverbosity>(exons, dna, positions, remapping, strand, s_align, e_align, min_exon_len, max_intron_len, min_intron_len)  ;
		alignment_valid =  alignment_valid && determine_exons_successful ;
	}

	clear_variant_cache(variant_cache) ;

	bool non_consensus_alignment=false ;

	if (alignment_valid) //Should align consecutive positions and have at least one exon
	{

	    if (strand=='-')
			exons=reverse(exons) ;
	    for (size_t i=0; i<exons.size(); i+=2)
			exons[i+1]++ ;
	    
	    if (myverbosity >= 2) 
		{
			fprintf(stdout, "# exons:\n");
			for (size_t i = 0; i < exons.size(); i += 2)
				fprintf(stdout, "# %i. %i - %i\n", (int)i / 2, exons[i], exons[i + 1]);
	    }

		{
			int ret = determine_read_variants<myverbosity, discover_variants>(contig_idx, s_align, e_align, dna_align, est_align, positions, 
																			  variants, align_variants_ref, aligned_positions_ref,
																			  read, read_string, read_quality, read_anno_ref, est_len_p, result_length, strand, ori,
																			  alignment_matches_ref, alignment_gaps_ref, alignment_mismatches_ref, alignment_qual_mismatches_ref) ;
			
			alignment_passed_filters_ref = (ret==0) && alignment_pass_filters(min_intron_len,max_intron_len,alignment_mismatches_ref,alignment_gaps_ref,exons.size()/2,min_exon_len,remapping) ;
		}
		
		if (myverbosity >= 2)
			fprintf(stdout,
					"# alignment: score=%1.3f  matches=%i  gaps=%i  anno=%s\n",
					alignscore, alignment_matches_ref, alignment_gaps_ref,
					read_anno_ref.c_str());
		
		if (myverbosity >= 1)
		{
			
 			fprintf(stdout, "# alignment with %i exons on reference (%i, %i, %i, %i, %c, %i)\t%s\n", (int)exons.size()/2, 
					(int)alignment_mismatches_ref, (int)alignment_gaps_ref, (int)alignment_mismatches_ref+alignment_gaps_ref, (int)alignment_qual_mismatches_ref, strand, ori, read_anno_ref.c_str()) ;
			if (considered_variants)
				fprintf(stdout, "# alignment with %i exons on superseq. (%i, %i, %i, %c, %i)\n", (int)exons.size()/2, 
						(int)alignment_mismatches_var, (int)alignment_gaps_var, (int)alignment_mismatches_var+alignment_gaps_var, strand, ori) ;

			if (false && perform_extra_checks) // BUG-TODO (debug this without compiler optimization, gdb seems to have a bug and shows wrong value for alignment_gaps_var)
			{
				assert(alignment_gaps_var<=alignment_gaps_ref) ;
			}
		}

		if (!considered_variants)
		{ // if these variables where not previously assigned, then use the ref-version instead
			assert(alignment_mismatches_var==0) ;
			alignment_gaps_var = alignment_gaps_ref ;
			alignment_mismatches_var = alignment_mismatches_ref ;
		}
		//assert(alignment_mismatches_var<=alignment_mismatches_ref || alignment_gaps_ref>alignment_gaps_var) ;
		if (alignment_gaps_var>alignment_gaps_ref)
		{
			alignment_gaps_var=alignment_gaps_ref ; // BUG-TODO: reconstruct_reference_alignment counts gaps wrongly
			if (myverbosity>=1 || perform_extra_checks)
				fprintf(stderr, "Warning: BUG-TODO: reconstruct_reference_alignment counts gaps wrongly, using number of gaps in ref alignment instead \n") ;
		}
		if (alignment_mismatches_var>alignment_mismatches_ref)
		{
			alignment_mismatches_var=alignment_mismatches_ref ; // BUG-TODO: reconstruct_reference_alignment may count mismatches wrongly
			if (myverbosity>=1  || perform_extra_checks)
				fprintf(stderr, "Warning: BUG-TODO: reconstruct_reference_alignment may count mismatches wrongly, using number of mismatches in ref alignment instead (%s)\n", read.id()) ;
		}

		aln = new ALIGNMENT();
		{
			int ret=-1 ;
			if (aln)
				ret = construct_intron_strings(aln, contig_idx, exons, strand, ori, non_consensus_alignment, remapping, non_consensus_search) ;
			if (ret<0)
			{
				result.delete_regions();
				delete[] acceptor ;
				delete[] donor ;
				delete[] s_align;
				delete[] e_align;
				delete[] dna_align;
				delete[] est_align;
				delete aln ;
				aln=NULL ;
				return ret ;
			}
		}

		if (remapping && exons.size()>=4)
		{
			REAL score = 1.0 ; // high probability for donor and acceptor -> large ss score
			alignscore += lookup_penalty(&alignment_parameters->d, 0, &score) ; 
			alignscore += lookup_penalty(&alignment_parameters->a, 0, &score) ; 
		}
		
		aln->qpalma_score = alignscore;

		aln->num_matches_ref = alignment_matches_ref;
		aln->num_gaps_ref = alignment_gaps_ref;
        aln->num_mismatches_ref = alignment_mismatches_ref ;
        aln->qual_mismatches_ref = alignment_qual_mismatches_ref ;

		aln->considered_variants=considered_variants ;
		aln->num_gaps_var = alignment_gaps_var;
        aln->num_mismatches_var = alignment_mismatches_var ;

		//assert(read_anno_ref.length()<4*Config::MAX_READ_LENGTH) ;
		aln->read_anno.assign(read_anno_ref) ;
		aln->exons = exons;
		aln->chromosome = &contig_idx;
		aln->strand = strand;
		if (!_config.USE_VARIANTS_EDITOP_FILTER || !considered_variants)
			aln->passed_filters=alignment_passed_filters_ref ;
		else
			aln->passed_filters=alignment_passed_filters_var ;
		aln->non_consensus_alignment = non_consensus_alignment ; 
		aln->remapped = remapping ;
		aln->found_variants = variant_list ;
		aln->align_variants = align_variants_ref ;
		aln->aligned_positions = aligned_positions_ref ;
		aln->variant_positions = variant_pos;
		aln->from_gm = 3;

		if (ori == 0)
			aln->orientation = '+';
		else
			aln->orientation = '-' ;
		aln->min_exon_len = min_exon_len ;
		aln->max_intron_len = max_intron_len ;
		aln->spliced = (exons.size()>=4) ;
		if (!aln->spliced)
			isunspliced = true;

		strcpy(aln->read_id, read.id());
		
		if (myverbosity >= 2)
		{
			fprintf(stdout,
					"# alignment (%i/%i) with %i exons found for %s (score=%1.3f  ref_matches=%i  ref_mismatches=%i  ref_gaps=%i strand=%c): %s\npassed filters:%i\n",
					non_consensus_search,non_consensus_alignment, (int) exons.size() / 2, read.id(), alignscore,
					alignment_matches_ref, alignment_mismatches_ref, alignment_gaps_ref, strand, read_anno_ref.c_str(), aln->passed_filters?1:0);
			for (size_t i = 0; i < exons.size(); i += 2)
				fprintf(stdout, "# exon %i: %i - %i\n", (int)i / 2, exons[i], exons[i+ 1]);
		}
	} 
	else
		if (myverbosity>=1)
		{
			fprintf(stdout, "# dropped alignment with %i exons (%i, %i, %i)\n", (int)exons.size()/2, 
					(int)alignment_mismatches_ref, (int)alignment_gaps_ref, (int)alignment_mismatches_ref+alignment_gaps_ref) ;
			fprintf(stdout, "alignment_valid=%i\n", alignment_valid) ;
			fprintf(stdout, "alignment_passed_filters_var=%i\n", alignment_passed_filters_var) ;
			fprintf(stdout, "determine_exons_successful=%i\n", determine_exons_successful) ;
			fprintf(stdout, "max_intron_len<_config.SPLICED_LONGEST_INTRON_LENGTH=%i\n", max_intron_len<_config.SPLICED_LONGEST_INTRON_LENGTH) ;
			fprintf(stdout, "alignment_mismatches(%i) <= _config.NUM_MISMATCHES(%i)=%i\n", alignment_mismatches_ref, _config.NUM_MISMATCHES, alignment_mismatches_ref <= _config.NUM_MISMATCHES) ;
			fprintf(stdout, "alignment_gaps <= _config.NUM_GAPS=%i\n", alignment_gaps_ref <= _config.NUM_GAPS) ;
			fprintf(stdout, "alignment_mismatches(%i)+alignment_gaps(%i) <= _config.NUM_EDIT_OPS=%i\n",alignment_mismatches_ref, alignment_gaps_ref, alignment_mismatches_ref+alignment_gaps_ref <= _config.NUM_EDIT_OPS) ;
			fprintf(stdout, "exons.size()=%i\n", (int)exons.size()) ;
			fprintf(stdout, "((int)exons.size() <= (_config.SPLICED_MAX_INTRONS+1)*2)=%i\n", ((int)exons.size() <= (_config.SPLICED_MAX_INTRONS+1)*2)) ;
			fprintf(stdout, "min_exon_len >= _config.SPLICED_MIN_SEGMENT_LENGTH=%i\n", min_exon_len >= _config.SPLICED_MIN_SEGMENT_LENGTH) ;
		}

	if (acceptor != NULL){
		delete[] acceptor ;
		acceptor=NULL;
	}
	if (donor != NULL){
		delete[] donor ;
		donor = NULL;
	}
	if (s_align != NULL){
		delete[] s_align;
		s_align=NULL;
	}
	if (e_align != NULL){
		delete[] e_align;
		e_align =NULL;
	}
	if (dna_align != NULL){
		delete[] dna_align;
		dna_align = NULL;
	}
	
	if (est_align != NULL){
		delete[] est_align;
		est_align = NULL;
	}
	

	return (int) isunspliced;
}

double QPalma::score_unspliced(Read const &read, const char * read_anno, const char strand, const char ori) const
{
	if (alignment_parameters==NULL)
	{
		int num_matches = read.length() ;
		
		for (size_t i=0; i<strlen(read_anno); i++)
		{
			if (read_anno[i]=='[')
			{
				num_matches-- ;
				i+=3 ;
			}
		}
		return num_matches ;
	}

	int num_qual_support = 0;
	if (alignment_parameters->num_qualityPlifs > 0)
		num_qual_support = alignment_parameters->qualityPlifs[0].len;
//	Alignment alignment(alignment_parameters->num_qualityPlifs, num_qual_support, true, verbosity);
	Alignment alignment(alignment_parameters->num_qualityPlifs, num_qual_support, !_config.NO_QPALMA, verbosity);

	double prb[read.length()] ;
	int offset_switched=0;
	
	for (size_t i = 0; i < read.length(); i++)
	{
		prb[i] = (read.quality(0)[i] - alignment_parameters->quality_offset);
		if (prb[i]<-10 || prb[i]>70)
			fprintf(stderr, "prb[%i]=%f (offset=%i, %s, %s)\n", (int)i, prb[i], alignment_parameters->quality_offset, read.quality(0), read.data()) ;

		if (_config.QPALMA_PRB_OFFSET_FIX)
		{
			if (alignment_parameters->quality_offset==33 && prb[i]>70)
			{
				fprintf(stderr, "setting prb offset from %i to %i\n", alignment_parameters->quality_offset, 64) ;
				_config.QPALMA_PRB_OFFSET_FIX=false ;
				alignment_parameters->quality_offset=64 ;
				read.set_quality_offset(64) ;
				offset_switched=1;
				
			}
			if (alignment_parameters->quality_offset==64 && prb[i]<-10)
			{
				fprintf(stderr, "setting prb offset from %i to %i\n", alignment_parameters->quality_offset, 33) ;
				_config.QPALMA_PRB_OFFSET_FIX=false ;
				alignment_parameters->quality_offset=33 ;
				read.set_quality_offset(33) ;
				offset_switched=1;
			}
		}

		//assert(prb[i]>=-10 && prb[i]<=70) ;
	}

	// If offset has been switched: recompute qualities for all positions
	if (offset_switched==1){
		for (size_t i = 0; i < read.length(); i++)
			prb[i] = (read.quality(0)[i] - alignment_parameters->quality_offset);
	}
	

	if (ori=='-')
		reverse(prb,read.length());

	double score1 = alignment.scoreUnsplicedAlignment(read_anno, prb, read.length(), alignment_parameters->qualityPlifs, alignment_parameters->matchmatrix,strand) ;
	//float score2 = alignment.scoreUnsplicedAlignment(read_anno, prb, read.length(), alignment_parameters->qualityPlifs, alignment_parameters->matchmatrix, '-') ;
	
	//if (score1>score2)
	return score1 ;
		//return score2 ;
}

