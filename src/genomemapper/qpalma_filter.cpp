#include <assert.h>
#include <string>
#include <sys/time.h>
#include <time.h>

#include "genomemapper.h"
#include "IntervalQuery.h"
#include "dyn_prog/qpalma_dp.h"


static clock_t last_filter_report=0 ;

static const int verbosity = 0 ;

int qpalma_filter_reason = -1 ;
const int num_filter_reasons = 4 ;
int qpalma_filter_stat_spliced[num_filter_reasons] = {0,0,0,0} ;
int qpalma_filter_stat_unspliced[num_filter_reasons] = {0,0,0,0} ;


int get_num_splicesites(std::string file_template, const char* type, Chromosome const &chr, char strand, int start, int end, float thresh)
{
	std::vector<int> positions ;
	return get_splicesite_positions(file_template, type, chr, strand, start, end, thresh, false, positions) ;
}

void qpalma_filter_stat_report()
{
	for (int i=0; i<num_filter_reasons; i++)
	{
		fprintf(stdout, "[filter] reason %i:\t%i spliced\t%i unspliced\t%1.2f%%\n", i, qpalma_filter_stat_spliced[i], qpalma_filter_stat_unspliced[i], 100*((float)qpalma_filter_stat_spliced[i])/((float)qpalma_filter_stat_unspliced[i]+qpalma_filter_stat_spliced[i])) ;
	}
}

void qpalma_filter_stat(bool spliced)
{
	if (qpalma_filter_reason<0) 
		return ;
	
	if (spliced)
		qpalma_filter_stat_spliced[qpalma_filter_reason]++ ;
	else
		qpalma_filter_stat_unspliced[qpalma_filter_reason]++ ;
	
	if (((clock()-last_filter_report)/CLOCKS_PER_SEC>=10))
	{
		last_filter_report = clock() ;
		qpalma_filter_stat_report() ;
	}

	qpalma_filter_reason=-1 ;
}

int qpalma_filter(struct alignment_t *ali, int num_N)
{
	assert(ali->exons.size()<=2) ;
	static bool use_ss_for_filtering = true ;
	
	Chromosome const &chr = *ali -> chromosome ;
	int start = ali -> exons[0] ;
	int end = ali -> exons[1] ;
	
	unsigned int num_gaps = ali -> num_gaps ;
	unsigned int num_matches = ali -> num_matches ;

	if (num_gaps > _config.FILTER_BY_MAX_GAPS || _read.length()-num_matches > _config.FILTER_BY_MAX_MISMATCHES+num_N)
	{
		if (verbosity>=1)
			fprintf(stdout, "filter decides YES: num_gaps=%i, num_mismatches=%i, num_N=%i\n", num_gaps, _read.length()-num_matches, num_N) ;
		
		qpalma_filter_reason=2 ; // cheap positive filter
			
		return 1 ;
	}

	float thresh_acc = _config.FILTER_BY_SPLICE_SITES_THRESH_ACC ;
	float thresh_don = _config.FILTER_BY_SPLICE_SITES_THRESH_DON ;
	int region = _config.FILTER_BY_SPLICE_SITES_REGION ;

	if (region==-1 || !_config.FILTER_BY_SPLICE_SITES)
	{
		qpalma_filter_reason=0 ; // cheap negative filter
		return 0 ;
	}

	if ( _read.length()-num_matches < _config.FILTER_BY_SPLICE_SITES_EDIT_MIN+num_N )
	{
		qpalma_filter_reason=0 ; // cheap negative filter
		return 0 ;
	}
	
	if (use_ss_for_filtering)
		try
		{
			int num_ss = get_num_splicesites(_config.ACC_FILES, "acc", chr, '+', start-region, end+region, thresh_acc) + 
				get_num_splicesites(_config.ACC_FILES, "acc", chr, '-', start-region, end+region, thresh_acc) +
				get_num_splicesites(_config.DON_FILES, "don", chr, '+', start-region, end+region, thresh_don) +
				get_num_splicesites(_config.DON_FILES, "don", chr, '-', start-region, end+region, thresh_don) ;
			
			if ( num_ss>0 )
			{
				if (verbosity>=1)
					fprintf(stdout, "filter decides YES: num_ss=%i, gaps=%i, num_mismatches=%i, num_N=%i\n", num_ss, num_gaps, _read.length()-num_matches, num_N) ;
				
				qpalma_filter_reason=3 ; // positive splice site filter
				
				return 1 ;
			}
		}
		catch (IntervalQueryException & e)
		{
			fprintf(stdout, "Warning: do not use splice site files for triggering qpalma alignments\n") ;
			use_ss_for_filtering=false ;
		}
	
	if (verbosity>=1)
		fprintf(stdout, "filter decides NO: num_gaps=%i, num_mismatches=%i, num_N=%i\n", num_gaps, _read.length()-num_matches, num_N) ;
	
	qpalma_filter_reason=1 ; // negative splice site filter 

	return 0 ;
}

