// authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

// Authors: Gunnar R\"atsch and Lisa Thalheim
// Copyright (C) 2008 by Friedrich Miescher Laboratory of the Max Planck Society

#include "palmapper.h"
#include "print.h"
#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include <fcntl.h>

#include <lang/Thread.h>
#include <palmapper/FileReporter.h>
#include <palmapper/Mapper.h>
#include <palmapper/JunctionMap.h>
#include <palmapper/VariantMap.h>
#include <palmapper/shogun/Signal.h>
//#include <wait.h>


Config _config;
Statistics _stats;

using namespace lang;

void palmapper_cleanup() 
{
  if (_config.LOCK_FILE_NAME.size()>0)
    remove(_config.LOCK_FILE_NAME.c_str()) ;
}

 
class MapperThread : public Thread, public Mapper {
public:
MapperThread(Genome &genome,	GenomeMaps &genomemaps, QueryFile &queryFile, QPalma &qpalma, Reporter &reporter, JunctionMap &junctionmap, JunctionMap &annotated_junctions, VariantMap & variants)
	: 	Mapper(genome, genomemaps, queryFile, qpalma, reporter,junctionmap,annotated_junctions, variants) {}
	void run() {
		try 
		{
			map_reads();
		} 
		catch (ShogunException & e)
		{
			fprintf(stderr, "ERROR: caught ShogunException: %s in thread %lu.\n", e.get_exception_string(), pthread_self()) ;
			CSignal::do_show_read_ids() ;
			fprintf(stderr, "\nExiting\n") ;
			palmapper_cleanup() ;
			exit(-1) ;
		}
		catch (std::exception & e)
		{
			fprintf(stderr, "ERROR: caught std::exception: %s in thread %lu.\n", e.what(), pthread_self()) ;
			CSignal::do_show_read_ids() ;
			fprintf(stderr, "\nExiting\n") ;
			palmapper_cleanup() ;
			exit(-1) ;
		}
		catch (...)
		{
			fprintf(stderr, "ERROR: caught unknown exception in thread %lu.\n", pthread_self()) ;
			CSignal::do_show_read_ids() ;
			fprintf(stderr, "\nExiting\n") ;
			palmapper_cleanup() ;
			exit(-1) ;
		}
	}
};


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

	int lockfile=-232 ;
	if (_config.LOCK_FILE_NAME.size()>0)
	  {
	    lockfile=open(_config.LOCK_FILE_NAME.c_str(), O_WRONLY | O_CREAT | O_EXCL, 255 );
	    close(lockfile);
	    if (lockfile<0 && errno==EEXIST)
	      {
		//another process is holding the lock
		fprintf(stderr, "The lock file %s is locked; try again later\n", _config.LOCK_FILE_NAME.c_str());
		exit(-99) ;
	      }
	    else if (lockfile<0)
	      {
		fprintf(stderr, "Unable to lock the lock file %s (%i)\n", _config.LOCK_FILE_NAME.c_str(), lockfile);
		exit(-99) ;
	      }	      
	  }

	Genome genome;

	_config.applyDefaults(&genome) ;
	_config.checkConfig() ;

	if (_config.VERBOSE && !_config.MAP_REVERSE)
	  fprintf(stdout, "Warning: only trigger alignments on forward strand\n") ;

	FILE *OUT_FP = NULL ;
	if ( _config.OUTPUT_FORMAT != OUTPUT_FORMAT_BAM )
		OUT_FP = Util::openFile(_config.OUT_FILE_NAME, "w");
	else
	{
		OUT_FP = TopAlignments::open_bam_pipe(_config.OUT_FILE_NAME) ;
		
		if (OUT_FP==NULL)
		{
			fprintf(stderr, "popen failed: samtools subprocess could not be started (mapped reads)\n") ;
			exit(-1) ;
		}
	}
	
	FILE *SP_OUT_FP = OUT_FP ;
	if (_config.SPLICED_HITS && _config.SPLICED_OUT_FILE_NAME.length() > 0 && _config.SPLICED_OUT_FILE_NAME!=_config.OUT_FILE_NAME )
	{
		if ( _config.OUTPUT_FORMAT != OUTPUT_FORMAT_BAM )
			SP_OUT_FP = Util::openFile(_config.SPLICED_OUT_FILE_NAME, "w") ;
		else
		{
			SP_OUT_FP = TopAlignments::open_bam_pipe(_config.SPLICED_OUT_FILE_NAME) ;
			
			if (SP_OUT_FP==NULL)
			{
				fprintf(stderr, "popen failed: samtools subprocess could not be started (spliced reads)\n") ;
				exit(-1) ;
			}
		}
	}
    FILE *LEFTOVER_FP = _config.LEFTOVER_FILE_NAME.length() > 0 ? Util::openFile(_config.LEFTOVER_FILE_NAME, "w+") : NULL;

	gzFile USED_VARIANTS_FP = NULL;
	if (_config.REPORT_USED_VARIANTS)
	{
		USED_VARIANTS_FP= gzopen(_config.REPORT_USED_VARIANT_FILE_NAME.c_str(), "wb9") ; // Util::openFile(_config.REPORT_USED_VARIANT_FILE_NAME, "w") ;
	}

	QueryFile queryFile(_config.QUERY_FILE_NAMES,_config.QUERY_FILE_STRANDS);
	FileReporter reporter(OUT_FP, SP_OUT_FP, USED_VARIANTS_FP, LEFTOVER_FP);

	// junction maps
	JunctionMap junctionmap(genome, _config.MAP_JUNCTIONS_PSEUDO_ANNO_COV, _config.ACC_CONSENSUS,_config.DON_CONSENSUS,_config.ACC_CONSENSUS_REV,_config.DON_CONSENSUS_REV);
	JunctionMap annotated_junctions(genome, _config.MAP_JUNCTIONS_PSEUDO_ANNO_COV, _config.ACC_CONSENSUS,_config.DON_CONSENSUS,_config.ACC_CONSENSUS_REV,_config.DON_CONSENSUS_REV);	

	if ( _config.OUTPUT_FORMAT == OUTPUT_FORMAT_BAM)
	{
		TopAlignments::print_bam_header(genome, OUT_FP) ;
		if (SP_OUT_FP!=OUT_FP)
			TopAlignments::print_bam_header(genome, SP_OUT_FP) ;
	}

	// genomemaps
	GenomeMaps *genomemaps = NULL;
	if (!_config.NO_READ_MAPPING)
	  genomemaps = new GenomeMaps(genome);

	if (_config.MAP_JUNCTIONS){
		int ret=junctionmap.init_from_gffs(_config.MAP_JUNCTIONS_FILE);
		if (ret!=0)
			return -1;
		junctionmap.filter_junctions(_config.MAP_JUNCTIONS_COVERAGE, _config.MAP_JUNCTIONS_MIN_SEGMENT_LENGTH, _config.MAP_JUNCTIONS_MAP_WINDOW, *genomemaps, _config.VERBOSE);
	}
	
	if (_config.SCORE_ANNOTATED_SPLICE_SITES){
		int ret=annotated_junctions.init_from_gffs(_config.ANNOTATED_SPLICE_SITES_FILE);
		if (ret!=0)
			return -1;
		annotated_junctions.filter_junctions(_config.MAP_JUNCTIONS_COVERAGE, _config.MAP_JUNCTIONS_MIN_SEGMENT_LENGTH, _config.MAP_JUNCTIONS_MAP_WINDOW, *genomemaps, _config.VERBOSE);
	}


	// variant maps
	VariantMap variants(genome, _config.MERGE_VARIANT_SOURCE_IDS, _config.VALIDATE_VARIANTS);

	if (_config.USE_VARIANTS)
	{
		int ret=variants.init_from_files(_config.USE_VARIANT_FILE_NAME);
		if (ret!=0)
			return -1;

		if (_config.FILTER_VARIANT_JUNCTIONS)
			variants.filter_variants_junctions(junctionmap) ;
		if (_config.FILTER_VARIANTS)
			variants.filter_variants(_config.FILTER_VARIANT_MINSOURCECOUNT, _config.FILTER_VARIANT_MINCONFCOUNT, _config.FILTER_VARIANT_MAXNONCONFRATIO, 
									 0,
									 _config.FILTER_VARIANT_SOURCES, _config.FILTER_VARIANT_REQSOURCES,
									 _config.FILTER_VARIANT_MAXLEN, _config.FILTER_VARIANT_MAP_WINDOW, *genomemaps) ;
		if (_config.UNIQUE_VARIANT_SOURCE_IDS)
			variants.unique_variant_source_ids() ;
		if (_config.CONVERT_SUBSTITUTION_VARIANTS)
			variants.convert_substitutions() ;
	}

	QPalma *qpalma = NULL;
	if (!_config.NO_READ_MAPPING)
	  qpalma = new QPalma(&genome, genomemaps, _config.VERBOSE);

	if (_config.SPLICED_HITS && _config.FILTER_BY_SPLICE_SITES && !_config.NO_SPLICE_PREDICTIONS)
	{
		//genomemaps = new GenomeMaps(genome);
		if (_config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			fprintf(stdout, "Using splice sites with confidence in top %1.2f%% percentile for read filtering\n", 100*_config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC) ;
		else
			fprintf(stdout, "Using splice sites with confidence => %1.2f%%  for read filtering\n", 100*_config.FILTER_BY_SPLICE_SITES_THRESH_ACC) ;

		if (_config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			_config.FILTER_BY_SPLICE_SITES_THRESH_ACC = _config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC ;
		if (_config.ACC_FILES.length()>0)
			qpalma->map_splice_sites(_config.ACC_FILES, 'a', _config.FILTER_BY_SPLICE_SITES_THRESH_ACC, _config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC!=0.0, true) ;
		if (_config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			fprintf(stdout, " -> acceptor  splice sites with confidence >= %1.2f%% \n", 100*_config.FILTER_BY_SPLICE_SITES_THRESH_ACC) ;

		if (_config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			_config.FILTER_BY_SPLICE_SITES_THRESH_DON = _config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC ;
		if (_config.DON_FILES.length()>0)
			qpalma->map_splice_sites(_config.DON_FILES, 'd', _config.FILTER_BY_SPLICE_SITES_THRESH_DON, _config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC!=0.0, true) ;
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
		if (_config.ACC_FILES.length()>0)
			qpalma->map_splice_sites(_config.ACC_FILES, 'a', _config.QPALMA_USE_SPLICE_SITES_THRESH_ACC, _config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC!=0.0, false) ;
		if (_config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			fprintf(stdout, "-> acceptor splice sites with confidence >= %1.2f%% \n", 100*_config.QPALMA_USE_SPLICE_SITES_THRESH_ACC) ;
			
		if (_config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			_config.QPALMA_USE_SPLICE_SITES_THRESH_DON = _config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC ;
		if (_config.DON_FILES.length()>0)
			qpalma->map_splice_sites(_config.DON_FILES, 'd', _config.QPALMA_USE_SPLICE_SITES_THRESH_DON, _config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC!=0.0, false) ;
		if (_config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			fprintf(stdout, "-> donor splice sites with confidence >= %1.2f%% \n", 100*_config.QPALMA_USE_SPLICE_SITES_THRESH_DON) ;
	}

	// slightly hacky ... 
	if (qpalma)
		Read::set_quality_offset(qpalma->get_qpalma_quality_offset()) ;
	
//	if (_config.REPORT_GFF_FILE_NAME.size()>0)
//		genomemaps.init_with_gff(_config.REPORT_GFF_FILE_NAME) ;
	
	// timing //////////////
	timer_mid=time(NULL);
	//printf("The current time is %s",asctime(localtime(&timer_mid)));
  	////////////////////////

	if (!_config.NO_READ_MAPPING)
	  {
	    if (_config.VERBOSE) { printf("Mapping reads\n"); }
	    
	    CSignal::set_handler() ;
	    CSignal::toggle_show_read_ids(true) ;
	    
	    unsigned int numThreads = _config.NUM_THREADS;
	    MapperThread *threads[numThreads];
	    std::string threadIds(".+-:=!$'");
	    for (unsigned int i = 0; i < numThreads; ++i) {
	      threads[i] = new MapperThread(genome, *genomemaps, queryFile, *qpalma, reporter, junctionmap, annotated_junctions, variants);
	      threads[i]->setProgressChar(threadIds[i % threadIds.length()]);
	      //printf("Starting thread %d\n", i);
	      if (numThreads>1)
			  threads[i]->launch();
	      else
			  threads[i]->run();
	    }
	    for (unsigned int i = 0; i < numThreads; ++i) {
	      if (numThreads>1)
	        threads[i]->join();
	      delete threads[i];
	    }
	    reporter.done();
	    CSignal::toggle_show_read_ids(false) ;
	  }
	else
	  fprintf(stdout, "Not mapping reads ... \n") ;

	if (_config.OUT_FILE_NAME.length() > 0)
	  {
	    if ( _config.OUTPUT_FORMAT != OUTPUT_FORMAT_BAM )
	      fclose(OUT_FP);
	    else
	      {
		int ret=TopAlignments::close_bam_pipe(OUT_FP) ;
		if (ret!=0)
		  {
		    fprintf(stderr, "samtools pipe failed (ret=%i)\n", ret) ;
		    exit(-1) ;
		  }
	      }
	  }
	if (_config.SPLICED_OUT_FILE_NAME.length() > 0 && SP_OUT_FP!=OUT_FP)
	  {
	    if ( _config.OUTPUT_FORMAT != OUTPUT_FORMAT_BAM )
	      fclose(SP_OUT_FP);
	    else
	      {
		int ret=TopAlignments::close_bam_pipe(SP_OUT_FP) ;
		if (ret!=0)
		  {
		    fprintf(stderr, "samtools pipe failed (ret=%i)\n", ret) ;
		    exit(-1) ;
		  }
	      }
	  }
	if (_config.LEFTOVER_FILE_NAME.length() > 0) 
	  fclose(LEFTOVER_FP);
	
	if (_config.REPORT_USED_VARIANTS)
	  {
	    gzclose(USED_VARIANTS_FP);
	  }
	
	/*if (_config.TRANSCRIBE_GFF)
	  {
		variants.transcribe_gff(_config.TRANSCRIBE_GFF_FILE, _config.TRANSCRIBE_FASTA_FILE) ;
		}*/

	if (_config.STATISTICS)	{
		print_stats(queryFile);
//TODO:		printf("R E D U N D A N T : %d\n", mapper.REDUNDANT);
//		printf("Max used slots: %d\n", MAX_USED_SLOTS);
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

	if (genomemaps != NULL && (_config.REPORT_MAPPED_REGIONS || _config.REPORT_MAPPED_READS || _config.REPORT_FILE!=NULL))
	{
		genomemaps->do_reporting(1) ;
		genomemaps->write_reporting() ;
		//genomemaps.clean_reporting() ;
	}
	if (_config.REPORT_GENOME_COVERAGE)
		genomemaps->write_cov_reporting() ;

	if (_config.REPORT_JUNCTIONS)
		junctionmap.report_to_gff(_config.REPORT_JUNCTIONS_FILE);

	if (_config.REPORT_VARIANTS)
	{
		if (_config.FILTER_VARIANTS)
			variants.filter_variants(_config.FILTER_VARIANT_MINSOURCECOUNT, _config.FILTER_VARIANT_MINCONFCOUNT, _config.FILTER_VARIANT_MAXNONCONFRATIO, 
									 _config.FILTER_VARIANT_MINUSECOUNT,
									 _config.FILTER_VARIANT_SOURCES, _config.FILTER_VARIANT_REQSOURCES,
									 _config.FILTER_VARIANT_MAXLEN, 
									 _config.FILTER_VARIANT_MAP_WINDOW, *genomemaps) ;

		int ret=variants.report_to_file(_config.REPORT_VARIANTS_FILE_NAME); 
		if (ret!=0)
			return -1;

		ret=variants.stats_to_file("",10); 
		if (ret!=0)
			return -1;

		if (!_config.REPORT_VARIANTS_STATS_FILE_NAME.empty())
		{
			ret=variants.stats_to_file(_config.REPORT_VARIANTS_STATS_FILE_NAME, 1000); 
			if (ret!=0)
				return -1;
		}
	}
	
	IntervalQuery iq;
	iq.cleanup(true) ;

	if (qpalma != NULL)
		delete qpalma;
	if (genomemaps != NULL)
		delete genomemaps;

	fprintf(stdout, "\n") ;
	
	if (_config.VERBOSE) { printf("Mapping finished\n"); }
	
	CSignal::unset_handler() ;

	palmapper_cleanup() ;
	
	return 0;
}
