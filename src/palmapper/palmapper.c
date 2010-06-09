// authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

// Authors: Gunnar R\"atsch and Lisa Thalheim
// Copyright (C) 2008 by Friedrich Miescher Laboratory of the Max Planck Society

#include "palmapper.h"
#include "print.h"
#include <palmapper/Mapper.h>

FILE *OUT_FP;
FILE *SP_OUT_FP;
FILE *TRIGGERED_LOG_FP; // #A#

Config _config;
Statistics _stats;

int main(int argc, char *argv[]) 
{
	_config.VersionHeader() ;

	init_shogun() ;
	
	// timing //////////////
	time_t timer, timer_mid, timer2;
  	timer=time(NULL);
  	//printf("The current time is %s",asctime(localtime(&timer)));
	////////////////////////

	// initialize variables
	_config.parseCommandLine(argc, argv);
	
	init_output_file(&_config);
	init_spliced_output_file(&_config);

	Genome genome;
	GenomeMaps genomemaps(&genome) ;

	_config.applyDefaults(&genome) ;
	_config.checkConfig() ;

	QueryFile queryFile(_config.QUERY_FILE_NAME);
	QPalma qpalma(&genome, &genomemaps, 0);
	Mapper hits(genome, genomemaps, queryFile, qpalma);


	// initialize GenomeMaps
	if (_config.REPORT_REPETITIVE_SEEDS || _config.REPORT_MAPPED_REGIONS || _config.REPORT_MAPPED_READS || _config.REPORT_FILE!=NULL || _config.FILTER_BY_SPLICE_SITES || _config.QPALMA_USE_SPLICE_SITES)
	{
		genomemaps.init_reporting() ;

		if (!_config.REPORT_RESET)
		{
			genomemaps.read_reporting() ;
			genomemaps.do_reporting(1) ;
		}
	}

	if (_config.SPLICED_HITS && _config.FILTER_BY_SPLICE_SITES && !_config.NO_SPLICE_PREDICTIONS)
	{
		if (_config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			fprintf(stdout, "Using splice sites with confidence in top %1.2f%% percentile for read filtering\n", 100*_config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC) ;
		else
			fprintf(stdout, "Using splice sites with confidence => %1.2f%%  for read filtering\n", 100*_config.FILTER_BY_SPLICE_SITES_THRESH_ACC) ;

		if (_config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			_config.FILTER_BY_SPLICE_SITES_THRESH_ACC = _config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC ;
		qpalma.map_splice_sites(_config.ACC_FILES, 'a', _config.FILTER_BY_SPLICE_SITES_THRESH_ACC, _config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC!=0.0, true) ;
		if (_config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			fprintf(stdout, " -> acceptor  splice sites with confidence >= %1.2f%% \n", 100*_config.FILTER_BY_SPLICE_SITES_THRESH_ACC) ;

		if (_config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			_config.FILTER_BY_SPLICE_SITES_THRESH_DON = _config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC ;
		qpalma.map_splice_sites(_config.DON_FILES, 'd', _config.FILTER_BY_SPLICE_SITES_THRESH_DON, _config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC!=0.0, true) ;
		if (_config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			fprintf(stdout, " -> donor  splice sites with confidence >= %1.2f%% \n", 100*_config.FILTER_BY_SPLICE_SITES_THRESH_DON) ;
	}

	if (_config.SPLICED_HITS && _config.QPALMA_USE_SPLICE_SITES && !_config.NO_SPLICE_PREDICTIONS)
	{
		if (_config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			fprintf(stdout, "Using splice sites with confidence in top %1.2f%% percentile to define QPALMA regions\n", 100*_config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC) ;
		else
			fprintf(stdout, "Using splice sites with confidence >= %1.2f to define QPALMA regions\n", _config.QPALMA_USE_SPLICE_SITES_THRESH_ACC) ;
		
		if (_config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			_config.QPALMA_USE_SPLICE_SITES_THRESH_ACC = _config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC ;
		qpalma.map_splice_sites(_config.ACC_FILES, 'a', _config.QPALMA_USE_SPLICE_SITES_THRESH_ACC, _config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC!=0.0, false) ;
		if (_config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			fprintf(stdout, "-> acceptor splice sites with confidence >= %1.2f%% \n", 100*_config.QPALMA_USE_SPLICE_SITES_THRESH_ACC) ;
			
		if (_config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			_config.QPALMA_USE_SPLICE_SITES_THRESH_DON = _config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC ;
		qpalma.map_splice_sites(_config.DON_FILES, 'd', _config.QPALMA_USE_SPLICE_SITES_THRESH_DON, _config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC!=0.0, false) ;
		if (_config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			fprintf(stdout, "-> donor splice sites with confidence >= %1.2f%% \n", 100*_config.QPALMA_USE_SPLICE_SITES_THRESH_DON) ;
	}

	if (_config.REPORT_GFF_FILE_NAME.size()>0)
		genomemaps.init_with_gff(_config.REPORT_GFF_FILE_NAME) ;
	
	// timing //////////////
	timer_mid=time(NULL);
	//printf("The current time is %s",asctime(localtime(&timer_mid)));
  	////////////////////////

 	if (_config.VERBOSE) { printf("Mapping reads\n"); }
	hits.map_reads(genome, genomemaps, &qpalma);

	if (_config.STATISTICS)	{
		print_stats(queryFile);
		printf("R E D U N D A N T : %d\n", hits.REDUNDANT);
		printf("Max used slots: %d\n", Mapper::MAX_USED_SLOTS);
		printf("\nList iterations: %d\n", _stats.listcount);
		printf("List iterations occurrences: %d\n",_stats.listocc);
	}
	else if (_config.VERBOSE) printf("Mapped Reads: %i of all %d reads\n", _stats.READS_MAPPED, _stats.NUM_READS);

	// timing //////////////
  	timer2=time(NULL);
  	//printf("The current time is %s",asctime(localtime(&timer2)));
  	int seconds = difftime(timer2, timer_mid);
  	int hours = seconds / 3600;
  	seconds %= 3600;
  	int minutes = seconds / 60;
  	seconds %= 60;
  	if (_config.STATISTICS) {
  		printf("Time needed to pre-process: %gs\n",difftime(timer_mid, timer));
  		printf("Time needed to map: %dh %dm %ds\n", hours, minutes, seconds);
  	}
  	seconds = difftime(timer2, timer);
  	hours = seconds / 3600;
  	seconds %= 3600;
  	minutes = seconds / 60;
  	seconds %= 60;
  	if (_config.STATISTICS) printf("Total time needed: %dh %dm %ds\n", hours, minutes, seconds);
  	////////////////////////

	if (_config.REPORT_REPETITIVE_SEEDS || _config.REPORT_MAPPED_REGIONS || _config.REPORT_MAPPED_READS || _config.REPORT_FILE!=NULL)
	{
		genomemaps.do_reporting(1) ;
		genomemaps.write_reporting() ;
		//_genomemaps.clean_reporting() ;
	}

	//if (_config.VERBOSE) { printf("La Fin.\n"); }

	return 0;
}
