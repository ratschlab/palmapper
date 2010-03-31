#include <assert.h>
#include <string>
#include <sys/time.h>
#include <time.h>

#include "genomemapper.h"
#include "IntervalQuery.h"
#include "dyn_prog/qpalma_dp.h"

int init_spliced_align(const char *fname, struct penalty_struct &h,
					   struct penalty_struct &a, struct penalty_struct &d,
					   struct penalty_struct *&qualityPlifs, int &num_qualityPlifs,
					   double*&matchmatrix, int dims[2], int &quality_offset)  ;
int init_alignment_parameters(std::string qpalma_file) ;
int clean_alignment_parameters() ;
int check_splice_files(std::string file_template) ;


int init_qpalma()
{
	if (_config.SPLICED_HITS)
	{
		if (_config.QPALMA_FILE.length()>0)
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
			if (!_config.NO_SPLICE_PREDICTIONS)
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
			if (!_config.NO_SPLICE_PREDICTIONS)
			{
				fprintf(stderr, "donor splice site predictions file not provided (use -acc option)\n") ;
				exit(-1) ;
			}
		}
	}
	return 0 ;
}

int compare_double(const void *a, const void *b) 
{
	double *ia = (double *) a;
	double *ib = (double *) b;

	return (int) (*ia - *ib);
}

int map_splice_sites(std::string file_template, char type, float &splice_site_threshold, bool estimate_thresh, bool do_report)
{
	char basename[1000] ;
	
	for (int chr=0; chr<(int)_genome.nrChromosomes() && (chr==0 || do_report); chr++)
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
				interval_matrix[1] = _genome.chromosome(chr).length();
				
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

				if (strand=='+' && chr==0 && estimate_thresh)
				{
					std::vector<double> acc_score2 ;
					
					try
						{
							for (int i=0; i<acc_size; i++)
								acc_score2.push_back(acc_score[i]) ;
						}
					catch (std::bad_alloc&)
						{
							fprintf(stderr, "[map_splice_sites] Could not allocate memory\n");
							delete[] acc_pos ;
							delete[] acc_index ;
							delete[] acc_score ;
							return -1;
						}
					
					std::sort(acc_score2.begin(), acc_score2.end()) ; 
					
					if (_config.VERBOSE)
						fprintf(stdout, "estimate_threshold: min=%1.3f  max=%1.3f %1.2f%% percentile=%1.3f\n", acc_score2[0], acc_score2[acc_size-1], splice_site_threshold*100, acc_score2[acc_size-1-(int)(acc_size*splice_site_threshold)]) ;
					
					splice_site_threshold=acc_score2[acc_size-1-(int)(acc_size*splice_site_threshold)] ;
				}

				/*int num_reported = 0 ;
				if (do_report)
					for (int i=0; i<acc_size; i++)
					{
						if (acc_score[i]>splice_site_threshold)
						{
							num_reported++ ;
							report_splice_site(chr, acc_pos[i], strand, type) ;
						}
						}*/
				
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

	return 0 ;
	
}

int clean_qpalma()
{
	return clean_alignment_parameters() ;
}

int init_alignment_parameters(std::string qpalma_file)
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

int clean_alignment_parameters()
{
	if (alignment_parameters==NULL)
		return -1 ;
	
	free(alignment_parameters->h.limits) ;
	free(alignment_parameters->h.penalties) ;
	free(alignment_parameters->a.limits) ;
	free(alignment_parameters->a.penalties) ;
	free(alignment_parameters->d.limits) ;
	free(alignment_parameters->d.penalties) ;
	for (int i=0; i<alignment_parameters->num_qualityPlifs; i++)
	{
		free(alignment_parameters->qualityPlifs[i].limits) ;
		free(alignment_parameters->qualityPlifs[i].penalties) ;
	}
	free(alignment_parameters->qualityPlifs) ;
	free(alignment_parameters) ;
	alignment_parameters = NULL ;

	return 0 ;
}


	
int check_splice_files(std::string file_template)
{
	char basename[1000] ;
	
	for (int chr=1; chr<(int)_genome.nrChromosomes(); chr++)
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
			fclose(fd) ;

			fd=fopen(confcumname, "r") ;
			if (fd==NULL)
			{
				fprintf(stderr, "%s does not exist\n", confcumname) ;
				return -1 ;
			}
			fclose(fd) ;
			
			if (strand=='+')
				strand='-' ;
			else
				break ;
		}
	}
	return 0 ;
}



int read_plif(FILE *fd, struct penalty_struct &plif) {
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
		fprintf(stderr,
				"[read_plif] Could not allocate memory (READ_ID = %s)\n",
				_read.id());
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

int read_matrix(FILE* fd, double *& matrix, char*& name, int dims[2]) {
	char buf[100];

	int narg = fscanf(fd, "%100s:\t", buf);
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
				"[read_matrix] Could not allocate memory (READ_ID = %s)\n",
				_read.id());
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

void skip_comment_lines(FILE* fd)
{
  const int buffer_len=1000 ;
  char buffer[buffer_len] ;
  
  while (!feof(fd))
    {
      long pos = ftell(fd) ;
      fgets(buffer, buffer_len, fd) ;
      if (buffer[0]!='#' && buffer[0]!=0)
	{
	  fseek(fd, pos, SEEK_SET) ;
	  break ;
	}
      //fprintf(stdout, "skipped comment line: %s\n", buffer) ;
    }
}

int init_spliced_align(const char *fname, struct penalty_struct &h,
		struct penalty_struct &a, struct penalty_struct &d,
		struct penalty_struct *&qualityPlifs, int &num_qualityPlifs,
		double*&matchmatrix, int dims[2], int &quality_offset) 
{
	num_qualityPlifs = 30;
	qualityPlifs = (struct penalty_struct*) malloc(
		sizeof(struct penalty_struct) * num_qualityPlifs);
	
	if (qualityPlifs == NULL) {
		fprintf(stderr, "[init_spliced_align] Could not allocate memory\n");
		return -1;
	}
	char *matrix_name = NULL;
	
	FILE * fd = fopen(fname, "r") ;
	skip_comment_lines(fd) ;
	int ret = read_plif(fd, h);
	if (ret < 0)
		return ret;

	assert(strcmp(h.name, "h")==0);
	h.max_len = _config.SPLICED_LONGEST_INTRON_LENGTH ;
	for (int i = 0; i < h.len; i++)
		if (h.limits[i] > 2)
			h.limits[i] = 2;

	skip_comment_lines(fd) ;
	ret = read_plif(fd, d);
	if (ret < 0)
		return ret;
	assert(strcmp(d.name, "d")==0);
	d.use_svm = 1;

	skip_comment_lines(fd) ;
	ret = read_plif(fd, a);
	assert(strcmp(a.name, "a")==0);
	if (ret < 0)
		return ret;
	a.use_svm = 1;
	for (int i = 0; i < num_qualityPlifs; i++)
	  {
	    skip_comment_lines(fd) ;
	    read_plif(fd, qualityPlifs[i]);
	  }

	skip_comment_lines(fd) ;
	ret = read_matrix(fd, matchmatrix, matrix_name, dims);
	if (ret < 0)
	{
		fprintf(stderr, "init_spliced_align: reading mmatrix failed\n") ;
		return ret;
	}
	assert(strcmp(matrix_name, "mmatrix")==0);
	free(matrix_name);
	matrix_name=NULL ;
	
	double *quality_offset_matrix ;
	int quality_offset_dims[2] ;
	skip_comment_lines(fd) ;
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

	return 0;
}
