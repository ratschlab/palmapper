#include <palmapper/Mapper.h>
#include <palmapper/align.h>
#include <palmapper/print.h>

void Mapper::map_reads_timing(int count_reads, float this_read)
{
	fprintf(stdout, "# [map_reads] timing: %1.3f, %1.3f, %1.3f, %1.3f, %1.3f, avg/read: %1.3f (%i/%i spliced alignm., %2.1f%%)",
			((float)time1)/CLOCKS_PER_SEC, ((float)time2a)/CLOCKS_PER_SEC,
			((float)time2b)/CLOCKS_PER_SEC, ((float)time2c)/CLOCKS_PER_SEC, ((float)time3)/CLOCKS_PER_SEC, ((float)time3)/CLOCKS_PER_SEC/count_reads,
			num_spliced_alignments_triggered, count_reads, 100.0*num_spliced_alignments_triggered/count_reads) ;

	if (this_read>=0)
		fprintf(stdout, ", this read: %1.3f\n", this_read) ;
	else
		fprintf(stdout, "\n") ;
}

Mapper::Mapper(Genome const &genome_,	GenomeMaps &genomemaps_, QueryFile &queryFile, QPalma &qpalma, Reporter &reporter)
:	_genome(genome_), _genomeMaps(genomemaps_), _queryFile(queryFile), _qpalma(qpalma), _reporter(reporter)
{
	REDUNDANT = 0;
	GENOME = new CHROMOSOME_ENTRY *[_genome.LONGEST_CHROMOSOME];
	for (unsigned int i=0; i!=_genome.LONGEST_CHROMOSOME; ++i)
		GENOME[i] = NULL;

	// initialize with meta information
	if (_config.STATISTICS) {
		init_statistic_vars(); //updated
	}
	init_alignment_structures(&_config);
    _LEFTOVER_FP = _config.LEFTOVER_FILE_NAME.length() > 0 ? Util::openFile(_config.LEFTOVER_FILE_NAME, "w+") : NULL;
	_ADAPTERTRIM_LOG_FP = _config.ADAPTERTRIM_STRATEGY_LOG.length() > 0 ? Util::openFile(_config.ADAPTERTRIM_STRATEGY_LOG, "w+") : NULL;
	if (_ADAPTERTRIM_LOG_FP==NULL && _config.ADAPTERTRIM_STRATEGY_LOG.length() > 0)
	{
		fprintf(stderr, "ERROR : could not open log file %s\n", _config.ADAPTERTRIM_STRATEGY_LOG.c_str()) ;
		exit(1) ;
	}
	MAXHITS = 0;
	c_map_fast = 0;
	c_map_short_read = 0;
	time1 = time2a = time2c = time3 = 0;
	last_timing_report=0 ;
	num_spliced_alignments_triggered=0 ;
	_progressChar = '.';
}

Mapper::~Mapper() {
	delete[] GENOME;
}

int Mapper::init_alignment_structures(Config * config) {

	/////////////////////////////
	////// if you change config->NUM_EDIT_OPS here, then put alloc_hits_by_editops in init_hit_lists_operator() after this function!!!
	/////////////////////////////

	if (config->NUM_MISMATCHES > config->NUM_EDIT_OPS) {
		fprintf(
				stderr,
				"\n!!! WARNING: Nr. of allowed mismatches is set to nr. of allowed edit operations!\n\n");
		config->NUM_MISMATCHES = config->NUM_EDIT_OPS;
	}
	if (config->NUM_GAPS > config->NUM_EDIT_OPS) {
		fprintf(
				stderr,
				"\n!!! WARNING: Nr. of allowed gaps is set to nr. of allowed edit operations!\n\n");
		config->NUM_GAPS = config->NUM_EDIT_OPS;
	}

	if (config->M_SCORE == config->MM_SCORE) {
		fprintf(
				stderr,
				"ERROR: Sorry, until now, program cannot handle equal mismatch- and match-scores! Please restart with different values!\n");
		exit(1);
	}

	if (config->VERBOSE)
		printf(
				"allowed edit operations: %d / allowed mismatches: %d / allowed gaps: %d\n%s strategy\t%s overhang alignment\n------------------\n",
				config->NUM_EDIT_OPS, config->NUM_MISMATCHES, config->NUM_GAPS,
				config->OUTPUT_FILTER==OUTPUT_FILTER_ALL ? "all hit" : "top hit",
				config->OVERHANG_ALIGNMENT ? "with" : "without");

	if (config->GAP_SCORE > config->MM_SCORE)
		WORST_SCORE = config->NUM_GAPS * config->GAP_SCORE + (config->NUM_EDIT_OPS - config->NUM_GAPS)
				* config->MM_SCORE;
	else
		WORST_SCORE = config->NUM_MISMATCHES * config->MM_SCORE + (config->NUM_EDIT_OPS
				- config->NUM_MISMATCHES) * config->GAP_SCORE;
	WORST_MM_SCORE = config->NUM_MISMATCHES * config->MM_SCORE;


	//only for debugging purposes:
	//chrseq = (char *) malloc (config->MAX_READ_LENGTH * sizeof(char));
	//chrseq[0] = '\0';

	return (0);
}

int Mapper::map_reads()
{
	unsigned int count_reads = 0;
	//int first_slot = 0, first_pos = 0 ;
	if (_config.STATISTICS)
		MAX_USED_SLOTS = 0;

 	if (_config.VERBOSE) { printf("Start mapping: "); }

    if (_config.LOG_TRIGGERED && (_TRIGGERED_LOG_FP = fopen(_config.TRIGGERED_LOG_FILE.c_str(), "w")) == NULL) { // #A#
		fprintf(stderr, "ERROR : Couldn't open log file %s\n", _config.TRIGGERED_LOG_FILE.c_str());     // #A#
		exit(1);                                                                        // #A#
	}                                                                                   // #A#
	for (;;) {

		count_reads++;
		Read _read(_queryFile);
		Result result(_queryFile.read_count(), _read, *this);
		if (!_queryFile.next_read(_read))
			break;
		clock_t start_time = clock() ;
		if (_config.VERBOSE && (count_reads % 100 == 0))
			printf("%i..", count_reads) ;
		map_read(result, start_time);
		time3 += clock()-start_time ;
		if (_config.VERBOSE || ((clock()-last_timing_report)/CLOCKS_PER_SEC>=10))
		{
			last_timing_report = clock() ;
			map_reads_timing(result._nr, ((float)clock()-start_time)/CLOCKS_PER_SEC) ;
		}
		if (_config.READ_COUNT_LIMIT && count_reads >= _config.READ_COUNT_LIMIT)
			break ;
		// progress output, just for user convenience
		if ((count_reads % 1000 == 0)) {
			printf("%c", _progressChar);
			fflush(stdout);
		}
		if ((count_reads % 10000 == 0)) {
			printf("%i", count_reads);
			fflush(stdout);
		}
	}

	map_reads_timing(count_reads) ;
	_qpalma.capture_hits_timing();

	//TODO mt report leftovers in Reporter
	if (_config.LEFTOVER_FILE_NAME.length() > 0)
		fprintf(_LEFTOVER_FP, "#done\n");

	if (_config.LEFTOVER_FILE_NAME.length() > 0)
		fclose(_LEFTOVER_FP);

	if (_config.STATISTICS)
	{
		printf("\n\n    MAP_FAST  = %d\n", c_map_fast);
		printf("    MAP_COMPL = %d\n\n", c_map_short_read);
		printf("Maximal number of hits per read: %d\n",MAXHITS);
	}

	_stats.NUM_READS = count_reads - 1;

	if (_config.VERBOSE) printf("\n");

	return(0);
}

void Mapper::map_read(Result &result, clock_t start_time) {
	QPalma const *qpalma = &result._qpalma._qpalma;
	unsigned int rtrim_cut = 0 ;
	unsigned int polytrim_cut_start = 0 ;
	unsigned int polytrim_cut_end = 0 ;
	unsigned int polytrim_cut_start_curr = 0 ;
	unsigned int polytrim_cut_end_curr = 0 ;
	unsigned int poly_length_start=0 ;
	unsigned int poly_length_end=0 ;
	int read_mapped;
	int cancel = 0 ;
	bool FILTER_STAT = false ;

	Read* trim_orig_read = NULL ;

	Hits &hits(result._readMappings);
	Read &_read(result._read);


	unsigned int adapter_cut_start = 0 ;
	unsigned int adapter_cut_end = 0 ;
	if (_config.ADAPTERTRIM_STRATEGY)
	{
		int orig_len = _read.length() ;

		_read.find_adapter(adapter_cut_start, adapter_cut_end) ;
		if (_read.length()-(adapter_cut_start+adapter_cut_end) < _config.ADAPTERTRIM_STRATEGY_MIN_LEN)
		{
			if (_config.LEFTOVER_FILE_NAME.length() > 0)
				print_leftovers(_read, "(too short after trimming)", _LEFTOVER_FP);
			return;
		}

		if (adapter_cut_start!=0)
		{
			Read *trim_read=new Read(_read) ;
			_read.trim_read_start(trim_read, adapter_cut_start) ;
			delete trim_read ;
		}
		if (adapter_cut_end!=0)
		{
			Read *trim_read=new Read(_read) ;
			_read.trim_read_end(trim_read, adapter_cut_end) ;
			delete trim_read ;
		}
		if (_ADAPTERTRIM_LOG_FP)// && (adapter_cut_start!=0 || adapter_cut_end!=0))
			fprintf(_ADAPTERTRIM_LOG_FP, "%s\t%i\t%i\t%i\t%i\n", _read.id(), orig_len, adapter_cut_start, adapter_cut_end, _read.length()) ;
	}

restart:

	// make it somehow dependent on the read length, the index depth and the number of mismatches
	_genomeMaps.REPORT_REPETITIVE_SEED_DEPTH_EXTRA = _read.length() - _config.INDEX_DEPTH - _config.NUM_MISMATCHES  ;

	char const *READ = _read.data();
	if (_config.VERBOSE)
		printf("# _read.id()=%s READ=%s\n", _read.id(), READ) ;


	int num_N=0 ;
	for (unsigned int i=0; i < _read.length(); i++)
		if (READ[i]!='A' && READ[i]!='C' && READ[i]!='G' && READ[i]!='T')
			num_N++ ;
	if (num_N>_config.NUM_MISMATCHES)
	{
		if (_config.VERBOSE)
			fprintf(stdout, "read has %i non-ACGT characters, allowing only %i mismatches -> skip read\n", num_N, _config.NUM_MISMATCHES) ;
		return;
	}

	if (_read.length() < _config.HITLEN_LIMIT)
	{
		fprintf(stderr, "\n!!! WARNING! Read %d (%s) with length %d is shorter than the hitlength limit (=%d) and will not be processed!\n\n",
			result._nr, _read.id(), _read.length(), _config.HITLEN_LIMIT);
		return;
	}
	//printf("%d ", count_reads); fflush(stdout);

	if (_config.STATISTICS) _stats.HITS_PER_READ = 0;

	hits.HITS_IN_SCORE_LIST = 0;

	// map_fast IF 1) best hit strategy 2) only hits up to RL/ID mismatches without gaps should be found
	// READ_LENGTH / _config.INDEX_DEPTH is the number of seeds fitting in the current read
	int nr_seeds = (int) (_read.length() / _config.INDEX_DEPTH);
	if (!hits.ALL_HIT_STRATEGY || _config.OUTPUT_FILTER==OUTPUT_FILTER_TOP || (_config.NUM_MISMATCHES < nr_seeds && _config.NUM_GAPS == 0))
	{
		int ret	= hits.map_fast(_read);	// if no hits could have been found: _config.ALL_HIT_STRATEGY = -1, necessitating execution of normal mapping in the following
		if (ret<0)
			cancel = 1 ;
		else
			c_map_fast++;

		//if (first_slot<0) // in this case not a single seed could be found
		//	cancel = 1 ;
	}

	time1+= clock()-start_time ;
		// map_complete IF 1) all hit strategy 2) best hit strategy and no mappings found in map_fast BUT NOT IF MM < RL/ID AND gaps=0 (since map_fast has already found them)
	if (!cancel)
	{
		if (((hits.ALL_HIT_STRATEGY!=0) || //_config.NOT_MAXIMAL_HITS || // check again
			 (_config.OUTPUT_FILTER==OUTPUT_FILTER_TOP && hits.SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.size()==0)) &&
			(!(_config.NUM_MISMATCHES < nr_seeds && _config.NUM_GAPS == 0) ) )
		{
			c_map_short_read++;

			int ret = hits.map_short_read(_read, result._nr);

			if (ret<0)
				cancel=2 ;

			time2a+= clock()-start_time ;

			// removing duplicates:
			hits.dealloc_mapping_entries();

			hits.CHROMOSOME_ENTRY_OPERATOR.used = 0;

			time2b += clock()-start_time ;

			if (!_config.REPORT_REPETITIVE_SEEDS)
				ret = hits.browse_hits();
			if (ret<0)
				cancel = 3 ;

			if (hits.ALL_HIT_STRATEGY < 0)
				hits.ALL_HIT_STRATEGY = 0;		// resetting _config.ALL_HIT_STRATEGY
		} else
		{
			if (_config.VERBOSE)
				printf("skipping map_short_read\n") ;
		}
	}
	else
	{
		if (_config.VERBOSE)
			printf("canceled\n") ;
	}


	time2c += clock()-start_time ;

	if (hits.ALL_HIT_STRATEGY < 0)
		hits.ALL_HIT_STRATEGY = 0;         // resetting _config.ALL_HIT_STRATEGY

	if (_config.STATISTICS && _stats.HITS_PER_READ > MAXHITS)
		MAXHITS = _stats.HITS_PER_READ;

	if (_config.OUTPUT_FILTER==OUTPUT_FILTER_TOP)
		hits.reset_num_edit_ops() ;

	read_mapped = 0 ;

	if (!cancel && !_config.REPORT_REPETITIVE_SEEDS)
	{
		hits._topAlignments.start_top_alignment_record();

		read_mapped = hits.analyze_hits(qpalma);	// returns 1 if at least one hit is printed, 0 otherwise
		if (_config.VERBOSE)
			printf("%i unspliced alignment found\n", (int) hits._topAlignments.size());

		bool trigger = false ;
		if (_config.SPLICED_HITS || _config.LOG_TRIGGERED)
			trigger = hits._topAlignments.size()==0 ||
				qpalma->qpalma_filter(result._qpalma, hits._topAlignments.get_alignment(0), num_N)!=0 ;

		if ( trigger )
		{
			if (_config.LOG_TRIGGERED)
				_read.printOn(_TRIGGERED_LOG_FP);
		}

		if (_config.SPLICED_HITS && (trigger  || FILTER_STAT))
			{
				num_spliced_alignments_triggered++ ;

				//(top_alignments.size()==0 || top_alignments[0]->num_matches <= _read.lenght() - _config.NUM_EDIT_OPS/2) )
				try
					{
						int ret = qpalma->capture_hits(hits, result._qpalma);
						//fprintf(stderr, "capture_hits ret=%i\n", ret) ;
						if (ret<0)
							cancel=4 ;
						if (_config.VERBOSE)
							fprintf(stdout, "capture_hits generated %i alignments\n", ret) ;
						if (FILTER_STAT)
							qpalma->qpalma_filter_stat(result._qpalma, ret>0) ;
					}
				catch (std::bad_alloc&)
					{
						fprintf(stderr, "[map_reads] allocating memory in capture_hits failed\n") ;
						cancel=5 ;
					}
			}
		if (FILTER_STAT && num_spliced_alignments_triggered>=5000)
		{
			fprintf(stdout, "final filter stat report\n") ;
			qpalma->qpalma_filter_stat_report() ;
			FILTER_STAT=false ;
		}
	}

	if (!cancel)
	{
		if (hits._topAlignments.size()>0) {
			_stats.READS_MAPPED++ ;
			result._rtrim_cut = rtrim_cut;
			result._polytrim_cut_start = polytrim_cut_start_curr;
			result._polytrim_cut_end = polytrim_cut_end_curr;
			_reporter.report(result);
//			hits._topAlignments.end_top_alignment_record(_read, _OUT_FP, _SP_OUT_FP, rtrim_cut, polytrim_cut_start_curr, polytrim_cut_end_curr);
			_read.set_orig(NULL) ;
			delete trim_orig_read ;
			trim_orig_read=NULL ;
			return;
		}
			read_mapped = 1 ;
		//if (_config.VERBOSE && read_mapped)
		//	printf("unspliced or spliced alignment found\n");


		{
			if (_config.RTRIM_STRATEGY && (_read.length() > _config.RTRIM_STRATEGY_MIN_LEN))
			{
				if (rtrim_cut==0)
				{
					_read.set_orig(NULL) ;
					delete trim_orig_read ;
					trim_orig_read=new Read(_read) ;
					_read.set_orig(trim_orig_read) ;
				}

				for (int s=0; s<(int)_config.RTRIM_STRATEGY_STEP; s++)
				{
					_read.cutOffLast();
					rtrim_cut += 1 ;
				}
				hits.clear();
				goto restart ;
			}

			if (_config.POLYTRIM_STRATEGY && (_read.length() > _config.POLYTRIM_STRATEGY_MIN_LEN))
			{
				// intended logic: increase start and end alternatively
				// until the individual stopping conditions are reached

				if (polytrim_cut_start==0 && polytrim_cut_end==0)
				{
					// determine the number of T's at beginning or A's at end
					_read.find_poly(poly_length_start, poly_length_end) ;
					//fprintf(stdout, "poly_length_start=%i, poly_length_end=%i\n", poly_length_start, poly_length_end) ;

					// copy original read
					_read.set_orig(NULL) ;
					delete trim_orig_read ;
					trim_orig_read=new Read(_read) ;
					_read.set_orig(trim_orig_read) ;
					if (_read.is_full_poly())
						poly_length_start=poly_length_end=0 ;
				}
				assert(trim_orig_read!=NULL) ;

				if (poly_length_start <= _config.POLYTRIM_STRATEGY_POLY_MIN_LEN)
					poly_length_start=0 ;
				if (poly_length_end <= _config.POLYTRIM_STRATEGY_POLY_MIN_LEN)
					poly_length_end=0 ;

				bool restart=false ;
				if (poly_length_start>=_config.POLYTRIM_STRATEGY_POLY_MIN_LEN || poly_length_end>=_config.POLYTRIM_STRATEGY_POLY_MIN_LEN)
				{
					// determine which side to cut
					bool start_cond = (polytrim_cut_start < poly_length_start &&
									   _read.length() - polytrim_cut_start >= _config.POLYTRIM_STRATEGY_MIN_LEN) ;
					bool end_cond = (polytrim_cut_end < poly_length_end &&
									 _read.length() - polytrim_cut_end >= _config.POLYTRIM_STRATEGY_MIN_LEN) ;

					if (start_cond && (polytrim_cut_start<polytrim_cut_end || !end_cond))
					{
						polytrim_cut_start += _config.POLYTRIM_STRATEGY_STEP ;
						_read.trim_read_start(trim_orig_read, polytrim_cut_start) ;
						restart = true ;
						polytrim_cut_start_curr = polytrim_cut_start ;
						polytrim_cut_end_curr = 0 ;
					}
					if (end_cond && !restart)
					{
						polytrim_cut_end += _config.POLYTRIM_STRATEGY_STEP ;
						_read.trim_read_end(trim_orig_read, polytrim_cut_end) ;
						restart = true ;
						polytrim_cut_start_curr = 0 ;
						polytrim_cut_end_curr = polytrim_cut_end ;
					}
				}
				if (restart)
				{
					//fprintf(stdout, "polytrim_cut_start_curr=%i, polytrim_cut_end_curr=%i: %s\n", polytrim_cut_start_curr, polytrim_cut_end_curr, _read.data()) ;
					hits.clear();
					goto restart ;
				}
			}

			if (_config.LEFTOVER_FILE_NAME.length() > 0)
				print_leftovers(_read, "", _LEFTOVER_FP);
		}
	}
	else
	{
		if (read_mapped && _config.VERBOSE)
			fprintf(stderr, "lost unspliced alignments\n") ;
	}


	if (cancel)
	{
		fprintf(stderr, "read %s could not be mapped (cancel=%i): %s\n", _read.id(), cancel, READ) ;
		if (_config.LEFTOVER_FILE_NAME.length() > 0)
			print_leftovers(_read, " (read mapping failed)", _LEFTOVER_FP);
	}

	// forget about the original read
	_read.set_orig(NULL) ;
	delete trim_orig_read ;
	trim_orig_read=NULL ;
}
