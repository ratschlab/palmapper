// Authors: Korbinian Schneeberger, Joerg Hagmann, Gunnar Raetsch
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany
// Copyright (C) 2009-2010 by Friedrich Miescher Laboratory, Tuebingen, Germany

#include "genomemapper.h"

HIT::HIT() {
	readpos = 0;
	start = 0;
	end = 0;
	chromosome = 0;
	mismatches = 0;
	gaps = 0;
	start_offset = 0;
	end_offset = 0;
	aligned = 0;
	/*int i;
	//for (i=0; i!=MAX_MISMATCHES; ++i) mismatch[i] = READ_LENGTH + 1;
	for (i=0; i!=MAX_MISMATCHES; ++i) {
		edit_op[i].mm = 0;
	}*/
	orientation = ' ';
	same_eo_succ = NULL;
	next = NULL;
	last = NULL;
}

static std::vector<int> SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS ;

inline void reset_num_edit_ops(int num_edit_ops)
{
 	SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.resize(0) ;
}

inline void update_num_edit_ops(int num_edit_ops, char & all_hit_strategy, int & NUM_EDIT_OPS_)
{
	assert(num_edit_ops<Config::MAX_EDIT_OPS) ;

	if (!all_hit_strategy && num_edit_ops < NUM_EDIT_OPS_)
		NUM_EDIT_OPS_ = num_edit_ops ;

 	if (_config.SUMMARY_HIT_STRATEGY)
	{
		// keep a list of minimal edit operation (assuming that each hit is only reported once)
		// the list is twice as long as NUM_TOP_ALIGNMENTS, as will be filtered later according to
		// the qpalma-alignment score (this is a heuristic, which may work well enough; alternatively, one
		// would need to compute the qpalma score for each hit, which will be too expensive

		bool inserted = false ;

		std::vector<int>::iterator it = SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.begin();

		for (uint8_t i = 0; i < SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.size(); i++, it++)
		{
			if ( num_edit_ops < SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS[i] )
			{
				SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.insert(it, num_edit_ops);
				inserted = true;

				if (SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.size() > Config::NUM_TOP_ALIGNMENTS*2)
					SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.pop_back();

				break;
			}
		}
		if (!inserted && SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.size() < Config::NUM_TOP_ALIGNMENTS*2)
		{
			SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.push_back(num_edit_ops) ;
			inserted = true ;
		}
	}

	if (!all_hit_strategy && SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.size()>0)
		NUM_EDIT_OPS_ = SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS[SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.size()-1] ;
}

#define TIME_CODE(x) 

int size_hit(HIT *hit, unsigned int *oldlength, char num);
int browse_hits();
int duplicate(HIT* hit);
void printgenome();

clock_t time1 = 0, time2a=0, time2a_seek=0, time2a_part1=0, time2a_part2=0, time2a_part3=0, time2a_part4=0, time2a_part5=0, time2b=0, time2c, time3=0 ;
clock_t time2a_seed2genome = 0 ;
static clock_t last_timing_report=0 ;
int num_spliced_alignments_triggered=0 ;

void map_reads_timing(int count_reads, float this_read=-1)
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

int map_reads()
{
	char eof = 0, read_mapped;
	unsigned int count_reads = 0;
	int first_slot = 0, first_pos = 0 ;
	int num_edit_ops = _config.NUM_EDIT_OPS;
	if (_config.SUMMARY_HIT_STRATEGY)
		reset_num_edit_ops(num_edit_ops) ;
	if (_config.STATISTICS) MAX_USED_SLOTS = 0;
	unsigned int MAXHITS = 0;
	int c_map_fast = 0;
	int c_map_short_read = 0;
	bool FILTER_STAT = false ;

 	if (_config.VERBOSE) { printf("Start mapping: "); }
    
    if (_config.LOG_TRIGGERED && (TRIGGERED_LOG_FP = fopen(_config.TRIGGERED_LOG_FILE.c_str(), "w")) == NULL) { // #A#
		fprintf(stderr, "ERROR : Couldn't open log file %s\n", _config.TRIGGERED_LOG_FILE.c_str());     // #A#
		exit(1);                                                                        // #A#
	}                                                                                   // #A#
    FILE *QUERY_FP = Util::openFile(_config.QUERY_FILE_NAME, "r");
    FILE *LEFTOVER_FP = _config.LEFTOVER_FILE_NAME.length() > 0 ? Util::openFile(_config.LEFTOVER_FILE_NAME, "w") : NULL;

	while (!eof) {
		
		count_reads++;
		int cancel = 0 ;

		LONGEST_HIT = 0;
		int rtrim_cut = 0 ;
		int polytrim_cut_start = 0 ;
		int polytrim_cut_end = 0 ;
		int poly_length_start=0 ;
		int poly_length_end=0 ;
		Read* poly_orig_read = NULL ;

		clock_t start_time = clock() ;

		eof = _read.read_short_read(QUERY_FP);

	restart:

		// make it somehow dependent on the read length, the index depth and the number of mismatches 
		REPORT_REPETITIVE_SEED_DEPTH_EXTRA = _read.length() - _config.INDEX_DEPTH - _config.NUM_MISMATCHES  ;
	
		if (eof != 0)
			continue;

		char const *READ = _read.data();
		if (_config.VERBOSE)
			printf("# _read.id()=%s READ=%s\n", _read.id(), READ) ;

		if (_config.VERBOSE && (count_reads % 100 == 0)) 
			printf("%i..", count_reads) ;

		if (_config.READ_COUNT_LIMIT && count_reads>_config.READ_COUNT_LIMIT)
			break ;
		
		int num_N=0 ;
		for (int i=0; i<_read.length(); i++)
			if (READ[i]!='A' && READ[i]!='C' && READ[i]!='G' && READ[i]!='T')
				num_N++ ;
		if (num_N>_config.NUM_MISMATCHES)
		{
			if (_config.VERBOSE)
				fprintf(stdout, "read has %i non-ACGT characters, allowing only %i mismatches -> skip read\n", num_N, _config.NUM_MISMATCHES) ;
			continue ;
		}

		if (_read.length() < _config.HITLEN_LIMIT) {
			fprintf(stderr, "\n!!! WARNING! Read %d (%s) with length %d is shorter than the hitlength limit (=%d) and will not be processed!\n\n",
				count_reads, _read.id(), _read.length(), _config.HITLEN_LIMIT);
		}
		else {
			//printf("%d ", count_reads); fflush(stdout);

			// progress output, just for user convenience
			if ((count_reads % 1000 == 0)) {
				printf(".");
				fflush(stdout);
			}
			if ((count_reads % 10000 == 0)) {
				printf("%i", count_reads);
				fflush(stdout);
			}

			if (_config.STATISTICS) _stats.HITS_PER_READ = 0;

			HITS_IN_SCORE_LIST = 0;

			// map_fast IF 1) best hit strategy 2) only hits up to RL/ID mismatches without gaps should be found
			// READ_LENGTH / _config.INDEX_DEPTH is the number of seeds fitting in the current read
			int nr_seeds = (int) (_read.length() / _config.INDEX_DEPTH);
 			if (!_config.ALL_HIT_STRATEGY || _config.SUMMARY_HIT_STRATEGY || (_config.NUM_MISMATCHES < nr_seeds && _config.NUM_GAPS == 0)) 
			{
				int ret	= _read.map_fast(first_slot, first_pos);	// if no hits could have been found: _config.ALL_HIT_STRATEGY = -1, necessitating execution of normal mapping in the following
				if (ret<0)
					cancel = 1 ;
				else
					c_map_fast++;

				if (first_slot<0) // in this case not a single seed could be found
					cancel = 1 ;
 			}

			time1+= clock()-start_time ;
 			// map_complete IF 1) all hit strategy 2) best hit strategy and no mappings found in map_fast BUT NOT IF MM < RL/ID AND gaps=0 (since map_fast has already found them)
			if (!cancel)
			{
				if ((_config.ALL_HIT_STRATEGY!=0 ||
					 (_config.SUMMARY_HIT_STRATEGY && SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.size()==0)) &&
				    (!(_config.NUM_MISMATCHES < nr_seeds && _config.NUM_GAPS == 0) || _config.NOT_MAXIMAL_HITS) )
				{

					c_map_short_read++;
					int ret = _read.map_short_read(count_reads, first_slot, first_pos);
					
					if (ret<0)
						cancel=2 ;
					
					time2a+= clock()-start_time ;
					
					// removing duplicates:
					dealloc_mapping_entries();
					
					CHROMOSOME_ENTRY_OPERATOR->used = 0;
					
					time2b += clock()-start_time ;
					
					if (!_config.REPORT_REPETITIVE_SEEDS)
						ret = browse_hits();
					if (ret<0)
						cancel = 3 ;
					
					if (_config.ALL_HIT_STRATEGY < 0)
						_config.ALL_HIT_STRATEGY = 0;		// resetting _config.ALL_HIT_STRATEGY
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

			if (_config.ALL_HIT_STRATEGY < 0)
				_config.ALL_HIT_STRATEGY = 0;         // resetting _config.ALL_HIT_STRATEGY

			if (_config.STATISTICS && _stats.HITS_PER_READ > MAXHITS)
				MAXHITS = _stats.HITS_PER_READ;

			_config.NUM_EDIT_OPS = num_edit_ops;		// resetting _config.NUM_EDIT_OPS.. has been changed probably in alignments and map_fast
			if (_config.SUMMARY_HIT_STRATEGY)
				reset_num_edit_ops(num_edit_ops) ;
					
			read_mapped = 0 ;

			if (!cancel && !_config.REPORT_REPETITIVE_SEEDS)
			{
				_topalignments.start_best_alignment_record();

				read_mapped = print_hits();	// returns 1 if at least one hit is printed, 0 otherwise
				if (_config.VERBOSE)
					printf("%i unspliced alignment found\n", (int)_topalignments.size()); 

				bool trigger = false ;
				if (_config.SPLICED_HITS || _config.LOG_TRIGGERED)
					trigger = _topalignments.size()==0 || 
						qpalma_filter(_topalignments.get_alignment(0), num_N)!=0 ;

				if ( trigger )
				{
					if (_config.LOG_TRIGGERED) { // #A# begin
						if (_read.format() == 0)
							fprintf(TRIGGERED_LOG_FP, "@%s\n%s\n+\n%s\n", _read.id(), READ, _read.quality()[0]);
						else if (_read.format() == 1)
							fprintf(TRIGGERED_LOG_FP, ">%s\n%s\n", _read.id(), READ);
						else
							fprintf(TRIGGERED_LOG_FP, "%s\t%s\t%d\t%s\t%s\t%s\n", _read.id(), READ,
									_read.pe_flag(), _read.quality()[0], _read.quality()[1], _read.quality()[2]);
					}    // #A# end
					
				}
				
				if (_config.SPLICED_HITS && (trigger  || FILTER_STAT))
					{
						num_spliced_alignments_triggered++ ;

						//(top_alignments.size()==0 || top_alignments[0]->num_matches <= _read.lenght() - _config.NUM_EDIT_OPS/2) )
						try
							{
								int ret = capture_hits();
								//fprintf(stderr, "capture_hits ret=%i\n", ret) ;
								if (ret<0)
									cancel=4 ;
								if (_config.VERBOSE)
									fprintf(stdout, "capture_hits generated %i alignments\n", ret) ;
								if (FILTER_STAT)
									qpalma_filter_stat(ret>0) ;
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
					qpalma_filter_stat_report() ;
					FILTER_STAT=false ;
				}
			}

			if (!cancel)
			{
				if (_topalignments.size()>0)
					read_mapped = 1 ;
				//if (_config.VERBOSE && read_mapped)
				//	printf("unspliced or spliced alignment found\n"); 
				
				_topalignments.end_best_alignment_record(rtrim_cut, polytrim_cut_start, polytrim_cut_end);
				
				if (read_mapped)
					_stats.READS_MAPPED++ ;
				else 
				{
				    if (_config.RTRIM_STRATEGY && (_read.length() > _config.RTRIM_STRATEGY_MIN_LEN))
					{
				    	_read.cutOffLast();

						rtrim_cut += 1 ;
						
						dealloc_mapping_entries(); //muss eigentlich nur der container zaehler zurückgesetzt werden... optimization?
						dealloc_hits();
						dealloc_hits_by_score();
						CHROMOSOME_ENTRY_OPERATOR->used = 0;
						if (LONGEST_HIT != 0)
							dealloc_hit_lists_operator();

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
							// copy original read
							delete poly_orig_read ;
							poly_orig_read=new Read(_read) ;
						}
						assert(poly_orig_read!=NULL) ;

						// determine which side to cut
						bool restart=false ;
						bool start_cond = (polytrim_cut_start < poly_length_start &&
										   _read.length() - polytrim_cut_start >= _config.POLYTRIM_STRATEGY_MIN_LEN) ;
						bool end_cond = (polytrim_cut_end < poly_length_end &&
										 _read.length() - polytrim_cut_end >= _config.POLYTRIM_STRATEGY_MIN_LEN) ;

						if (start_cond && (polytrim_cut_start<polytrim_cut_end || !end_cond))
						{
							polytrim_cut_start += 1 ;
							_read.trim_read_start(poly_orig_read, polytrim_cut_start) ;
							restart = true ;
						}
						if (end_cond && !restart)
						{
							polytrim_cut_end += 1 ;
							_read.trim_read_end(poly_orig_read, polytrim_cut_end) ;
							restart = true ;
						}

						if (restart)
						{
							dealloc_mapping_entries();
							dealloc_hits();
							dealloc_hits_by_score();
							CHROMOSOME_ENTRY_OPERATOR->used = 0;
							if (LONGEST_HIT != 0)
								dealloc_hit_lists_operator();
							
							goto restart ;
						}
					}

					if (_config.LEFTOVER_FILE_NAME.length() > 0)
						print_leftovers("", LEFTOVER_FP);
				}
			}
			else
			{
				if (read_mapped && _config.VERBOSE)
					fprintf(stderr, "lost unspliced alignments\n") ;
			}

			time3 += clock()-start_time ;

			if (cancel)
			{
				fprintf(stderr, "read %s could not be mapped (cancel=%i): %s\n", _read.id(), cancel, READ) ;
				if (_config.LEFTOVER_FILE_NAME.length() > 0)
					print_leftovers(" (read mapping failed)", LEFTOVER_FP);
			}

			dealloc_mapping_entries(); //muss eigentlich nur der container zaehler zurückgesetzt werden... optimization?
			dealloc_hits();
			dealloc_hits_by_score();
			CHROMOSOME_ENTRY_OPERATOR->used = 0;
			if (LONGEST_HIT != 0)
				dealloc_hit_lists_operator();
			delete poly_orig_read ;

			if (_config.VERBOSE || ((clock()-last_timing_report)/CLOCKS_PER_SEC>=10))
			{
				last_timing_report = clock() ;
				map_reads_timing(count_reads, ((float)clock()-start_time)/CLOCKS_PER_SEC) ;
			}
		}
		
	}

	map_reads_timing(count_reads) ;
	capture_hits_timing() ;

	//TODO dd use otehr criteria
	fprintf(stdout, "\n#done\n") ;
	if (OUT_FP!=stdout)
		fprintf(OUT_FP, "#done\n") ;
	if (SP_OUT_FP!=stdout)
		fprintf(SP_OUT_FP, "#done\n") ;
	if (_config.LEFTOVER_FILE_NAME.length() > 0)
		fprintf(LEFTOVER_FP, "#done\n");

	if (_config.OUT_FILE_NAME.length() > 0)
		fclose(OUT_FP);
	if (_config.SPLICED_OUT_FILE_NAME.length() > 0)
		fclose(SP_OUT_FP);
	if (_config.LEFTOVER_FILE_NAME.length() > 0)
		fclose(LEFTOVER_FP);

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


int extend_seed(int direction, int seed_depth_extra, Chromosome &chr, int genome_pos, int readpos)
{
	int e = 0 ;
	if (direction==-1)	// forward strand
	{
		for (e=0; e<seed_depth_extra; e++)
		{
			int readpos_ = readpos+_config.INDEX_DEPTH+e-1 ;
			if (readpos_>=_read.length())
				break ;
			int genome_pos_ = genome_pos+_config.INDEX_DEPTH+e ;
			if (genome_pos_>=chr.length())
				break ;
			
//fprintf(stdout,"%c - %c +\n", READ[readpos_], CHR_SEQ(genome_chr,genome_pos_)) ;
			if (_read.data()[readpos_]!=chr[genome_pos_])
				break ;
		}
	}
	if (direction!=-1)	// reverse strand
	{
		for (e=0; e<seed_depth_extra; e++)
		{
			int readpos_ = readpos+_config.INDEX_DEPTH+e-1;
			if (readpos_>=_read.length())
				break ;
			int genome_pos_ = genome_pos-e-1;
			if (genome_pos_ < 0)
				break ;
			
//fprintf(stdout,"%c - %c -\n", get_compl_base(READ[readpos_]), CHR_SEQ(genome_chr,genome_pos_)) ;
			if (get_compl_base(_read.data()[readpos_])!=chr[genome_pos_])
				break ;
		}
	}

	return e ;
}

std::vector<bool> seed_covered ;
clock_t seed_covered_reporting_time = 0 ;

struct seedlist
{
// dd	int chr ;
	Chromosome *chr;
	int pos ;
} ;

int seed2genome(unsigned int num, unsigned int index_slot, unsigned int readpos, char reverse)
{
#ifndef BinaryStream_MAP
	STORAGE_ENTRY *index_mmap ;
#else
	CBinaryStream<STORAGE_ENTRY>* index_mmap;
#endif
	INDEX_ENTRY index_entry;

	unsigned int block;
	unsigned char pos;
	char strand;
	unsigned int direction;
	unsigned int mmoffset; //Mismatch-Offset
	int read_num = -1 ;
	char flag = 0;

	CHROMOSOME_ENTRY *chromosome_director, *chromosome_director_neighbor, *chromosome_director_new;
	MAPPING_ENTRY *mapping_entry, *existing_entry, *neighbor;

	HIT *hit = NULL;

	unsigned int i;
	unsigned int oldlength = _config.INDEX_DEPTH-1;

	// reverse = 1: only index, 2: only index_rev, 4: both
	while (reverse > 0) {

		if (reverse != 2) {
			index_mmap = INDEX_FWD_MMAP;
			index_entry = *(INDEX+index_slot);
			direction = -1;
		}
		else {
			index_mmap = INDEX_REV_MMAP;
			index_entry = *(INDEX_REV+index_slot);
			direction = 1;
		}

		if (read_num == num) 
		{
			printf("###############################################\n");
			printf("Add seed to genomepositions from slot # %d (%s) containing %d genomepositions\n", index_slot, get_seq(index_slot), index_entry.num);
		}

		
		if (index_entry.num) {

			STORAGE_ENTRY* se_buffer=NULL;

			try
			{
				se_buffer=new STORAGE_ENTRY[index_entry.num];
			}
			catch (std::bad_alloc)
			{
				return -1;
			}

			// make sure that every seed is only processed once
			int index_entry_num=index_entry.num ;
			int report_repetitive_seeds = _config.REPORT_REPETITIVE_SEEDS ;
			if (seed_covered.size()<index_slot)
				seed_covered.resize(index_slot+1, false) ;
			if (report_repetitive_seeds)
			{
				if (seed_covered[index_slot])
					report_repetitive_seeds=0 ;
				else
					seed_covered[index_slot] = true ;

				if ((clock()-seed_covered_reporting_time)/CLOCKS_PER_SEC>10)
				{
					seed_covered_reporting_time=clock() ;
					size_t num_covered = 0 ;
					for (size_t ii=0; ii<seed_covered.size(); ii++)
						if (seed_covered[ii])
							num_covered++ ;
					fprintf(stdout, "# seed coverage: %ld/%ld (%2.1f%%)\n", num_covered, seed_covered.size(), 100*((float)num_covered)/seed_covered.size()) ;
				}
			}
			
			if (index_entry.num>_config.SEED_HIT_CANCEL_THRESHOLD && !report_repetitive_seeds)
				index_entry_num=0 ;
			{
				TIME_CODE(clock_t start_time = clock()) ;
#ifndef BinaryStream_MAP
				index_pre_buffer(index_mmap, se_buffer, index_entry.offset-index_entry_num, index_entry_num);
#else
				index_mmap->pre_buffer(se_buffer, index_entry.offset-index_entry_num, index_entry_num);
#endif
				TIME_CODE(time2a_seek += clock()-start_time) ;
			}

			std::vector<struct seedlist> extended_seedlist ;
			
			for (i=0; i<index_entry_num; i++) { // foreach seed...

				TIME_CODE(clock_t start_time = clock()) ;
				
				if (read_num == num) {
					printf("############################\n");
					printf("Now adding seed # %d/%d of read %i (%s), slot %i, ori %d\n", i+1, index_entry.num, num, _read.id(), index_slot, reverse);
				}

				// Every mapping gets an entry
				mapping_entry = alloc_mapping_entry();
				if (!mapping_entry)
				{
					delete[] se_buffer;
					TIME_CODE(time2a_part1 += clock()-start_time) ;
					return -1 ;
				}

				mapping_entry->readpos = readpos;
				if (read_num == num) {
					printf("Readpos %d\n", readpos);
				}

				block = 0;
				pos = 0;
				unsigned char* p_block = (unsigned char*) &block;
				STORAGE_ENTRY se;
				unsigned char* p_id;

				// Get current position in the chromosome
				if (reverse != 2) {
					//memcpy(&block, &((index_mmap+(index_entry.offset-(i+1)))->id[0]), 3 * sizeof(char));
					//memcpy(&pos, &((index_mmap+(index_entry.offset-(i+1)))->id[3]), sizeof(unsigned char));
					//se = (*index_mmap)[index_entry.offset-(i+1)];
					se = se_buffer[index_entry.num-(i+1)];
					//se2 = (*index_mmap)[index_entry.offset-(i+1)];
					strand = 1;
				}
				else {
					//memcpy(&block, &((index_mmap+(index_entry.offset-(index_entry.num-i)))->id[0]), 3 * sizeof(char));
					//memcpy(&pos, &((index_mmap+(index_entry.offset-(index_entry.num-i)))->id[3]), sizeof(unsigned char));
					//se = (*index_mmap)[index_entry.offset-(index_entry.num-i)];
					se = se_buffer[i];
					strand = -1;
				}
				p_id=se.id;
				p_block[0]=p_id[0];
				p_block[1]=p_id[1];
				p_block[2]=p_id[2];
				pos = p_id[3];

				unsigned int genome_pos = pos + BLOCK_TABLE[block].pos;
				//unsigned int genome_chr = BLOCK_TABLE[block].chr;
				Chromosome &genome_chr = _genome.chromosome(BLOCK_TABLE[block].chr);

				// Check that read doesn't cross chrom borders
				if ((reverse != 2 && (genome_pos < readpos-1 || genome_pos+_read.length()-readpos >= genome_chr.length())) ||
				    (reverse == 2 && (genome_pos < _read.length()-(_config.INDEX_DEPTH+readpos-1) || genome_pos+_config.INDEX_DEPTH+readpos-2 >= genome_chr.length()))
					)
				{
					continue;
				}


				if (report_repetitive_seeds)
				{   // check every seed, whether it is extendable by REPORT_REPETITIVE_SEED_DEPTH_EXTRA nucleotides 
					/// and report it
					int e = extend_seed(direction, REPORT_REPETITIVE_SEED_DEPTH_EXTRA, genome_chr, genome_pos, readpos) ;

					if (e==REPORT_REPETITIVE_SEED_DEPTH_EXTRA)
					{
						struct seedlist seed;
						seed.chr= &genome_chr ;
						seed.pos=genome_pos ;
						extended_seedlist.push_back(seed) ;
					}					
				}

				if (index_entry.num>_config.SEED_HIT_CANCEL_THRESHOLD)
					continue ;

				if (index_entry.num>=_config.INDEX_DEPTH_EXTRA_THRESHOLD)
				{
					int ee = extend_seed(direction, _config.INDEX_DEPTH_EXTRA, genome_chr, genome_pos, readpos) ;
					if (ee!=_config.INDEX_DEPTH_EXTRA)
						continue ;
				}
				
				//if (i<10 && index_entry.num>10000)
				//	printf("block %d pos %d [block].pos %d genome_pos %d genome_chr %d\n", block, pos, BLOCK_TABLE[block].pos, genome_pos, genome_chr);

				if (read_num == num) {
					printf("Genome Entry: chr %d   pos %d   \n", genome_chr.nr()+1, genome_pos);
				}

				TIME_CODE(time2a_part1 += clock()-start_time) ;
				TIME_CODE(start_time = clock()) ;

				// Check if there is already a chromosome director and get this or a new one
				if (*(GENOME+genome_pos) == NULL) {
					flag = 1;

					if (read_num == num) {
						printf("Alloc new chromosome director genome_chr %d\n", genome_chr.nr()+1);
					}

					chromosome_director = alloc_chromosome_entry(genome_pos, genome_chr, strand);
					if (!chromosome_director)
					{
						delete[] se_buffer;
						TIME_CODE(time2a_part2 += clock()-start_time) ;
						return -1;    // chrom_container_size too small -> cancel this read
					}

					*(GENOME + genome_pos) = chromosome_director;

				}
				else {
					if (read_num == num) {
						printf("Found chromosome director\n");
					}

					chromosome_director = *(GENOME + genome_pos);

					if (read_num == num) {
						printf("ChrEntryOp:       %p\n", (CHROMOSOME_ENTRY_OPERATOR->entries));
						printf("ChrEntryOp[used]: %p\n", ((CHROMOSOME_ENTRY_OPERATOR->entries+CHROMOSOME_ENTRY_OPERATOR->used)));
						printf("chrom_director:   %p\n", chromosome_director);
						printf("used = %d\n",CHROMOSOME_ENTRY_OPERATOR->used);
						printf("chrom_dir->gen_pos: %d\n", chromosome_director->genome_pos);
						printf("genome_pos: %d\n", genome_pos);
					}

					// is chrom_director from the actual read or from a former one?
					if (chromosome_director >= (CHROMOSOME_ENTRY_OPERATOR->entries+CHROMOSOME_ENTRY_OPERATOR->used)
							|| (chromosome_director < (CHROMOSOME_ENTRY_OPERATOR->entries+CHROMOSOME_ENTRY_OPERATOR->used) &&
								(chromosome_director->genome_pos != genome_pos))) {

						// it's from a former read, thus we need a new chrom_director:
						chromosome_director = alloc_chromosome_entry(genome_pos, genome_chr, strand);
						if (!chromosome_director)
						{
							delete[] se_buffer;
							TIME_CODE(time2a_part2 += clock()-start_time) ;
							return -1;    // chrom_container_size too small -> cancel this read
						}

						*(GENOME + genome_pos) = chromosome_director;

						if (read_num == num) {
							printf("Overwrite chromosome director %p\n", chromosome_director);
						}
					}
				}

				// Parse the list of chromosome directors
				//printf("chrom_dir %d, genome_chr %d\n", chromosome_director->chromosome, genome_chr);
				int c=0;
				while (chromosome_director->next != NULL && (chromosome_director->chromosome != genome_chr.nr() ||
							(chromosome_director->chromosome == genome_chr.nr() && chromosome_director->strand != strand)))
				{
					if (_config.STATISTICS && c == 0) _stats.listocc++;
					if (_config.STATISTICS) _stats.listcount++;
					chromosome_director = chromosome_director->next;
				}

				// Chromosome director is set, but still it could be the wrong chromosome, if the right chromosome is not in there so far.
				if (chromosome_director->chromosome != genome_chr.nr() || chromosome_director->strand != strand) {
					if (read_num == num) printf("Expanding list with new chrom.director\n");
					chromosome_director_new = alloc_chromosome_entry(genome_pos, genome_chr, strand);
					if (!chromosome_director_new)
					{
						delete[] se_buffer;
						TIME_CODE(time2a_part2 += clock()-start_time) ;
						return -1;
					}

					chromosome_director_new->next = chromosome_director->next;
					chromosome_director->next = chromosome_director_new;

					chromosome_director = chromosome_director_new;
				}

				// Paste MAPPING_ENTRY in list of chromosome director slot
				if (read_num == num) {
					printf("Mapping entry, genome_chr %d\n", genome_chr.nr()+1);
					printf("Mapping entry, chromsome director %p\n", chromosome_director);
				}
				if (chromosome_director->mapping_entries == NULL) {

					if (read_num == num) {
						printf("List of mapping entries in the chromosome director was empty\n");
					}
					chromosome_director->mapping_entries = mapping_entry;
				}
				else {
					if (flag == 1) {
						printf("!!!!!!!!!!!Found entry in chr dir at genome pos %d  --  %p\n", genome_pos, chromosome_director->mapping_entries	);
						exit(1);
					}

					if (read_num == num) {
						printf("Already entries in the chromosome director\n");
					}
					existing_entry = chromosome_director->mapping_entries;
					mapping_entry->pred = existing_entry;
					existing_entry->succ = mapping_entry;
					chromosome_director->mapping_entries = mapping_entry;
				}


				TIME_CODE(time2a_part2 += clock()-start_time) ;
				TIME_CODE(start_time = clock()) ;
				
				// HIT EXTENSION

				//Check left (for plus strand) and right (for minus strand) neighbor at the genomeposition of the hit if there is a hit to join.
				if (genome_pos > 0 && readpos > 1) {
					if (read_num == num) {
						printf("Now checking if left neighbor exists and is willing to join %i\t%i\t%i\n",genome_pos, genome_pos+direction,reverse);
					}

					if (*(GENOME + (genome_pos + direction)) != NULL) { // Is there a chromosome director?
						if (read_num == num) {
							printf("  Found a neigboured chromosome director\n");
						}
						chromosome_director_neighbor = *(GENOME + (genome_pos + direction));

						// Is the chrom director from actual read?
						if (chromosome_director_neighbor < (CHROMOSOME_ENTRY_OPERATOR->entries+CHROMOSOME_ENTRY_OPERATOR->used) &&
								(chromosome_director_neighbor->genome_pos == genome_pos + direction)) {

							// is there a mapping entry in the right chromosome solt, mr. director?
							// Search for chromosome_director_neighbor for the correct chromosome and strand
							c=0;
							while (chromosome_director_neighbor->next != NULL && (chromosome_director_neighbor->chromosome != genome_chr.nr() ||
										(chromosome_director_neighbor->chromosome == genome_chr.nr() && chromosome_director_neighbor->strand != strand)))
							{
								if (_config.STATISTICS) _stats.listcount++;
								if (_config.STATISTICS && c == 0) _stats.listocc++;
								chromosome_director_neighbor = chromosome_director_neighbor->next;
							}

							if 	(chromosome_director_neighbor->chromosome == genome_chr.nr() && chromosome_director_neighbor->strand == strand)
							{

								if (read_num == num) {
									printf("  Found neighboured mapping entry list\n");
								}
								neighbor = chromosome_director_neighbor->mapping_entries;
								if (read_num == num) {
									printf("  Neighbour info: readpos %d \n", neighbor->readpos);
								}
								if (read_num == num) { printf("neighbor hit: "); if (neighbor->hit!=NULL) printhit(neighbor->hit); else printf("null\n"); }
								if (neighbor->readpos == mapping_entry->readpos-1) { // is the neighbored mapping entry also neighbored in the read?
									hit = neighbor->hit;
									if (hit != NULL) {
										oldlength = hit->end - hit->start + 1;
										mapping_entry->hit = neighbor->hit;
										if (reverse != 2) hit->end++;
										else hit->start--;
										if (read_num == num) printhit(hit);
									}
								}

							}
						}
					}

					if (read_num == num) {
						printf("Will the hit be combined with adjacent neighbor? %d\n", (hit!=NULL));
					}


					TIME_CODE(time2a_part3 += clock()-start_time) ;
					TIME_CODE(start_time = clock()) ;

					// MISMATCH extension

					//combine with possible hit at position seedlength+1 to the left(+) or right(-) to span hit over mismatch
					if (!_config.NOT_MAXIMAL_HITS && hit == NULL && _config.NUM_MISMATCHES != 0) {
						if (read_num == num) printf("Now checking if hit can be extended over mismatch\n");
						mmoffset = (reverse != 2)? -_config.INDEX_DEPTH - 1: _config.INDEX_DEPTH + 1;

						if ( (genome_pos + mmoffset > 0) && (genome_pos + mmoffset < genome_chr.length()) && (*(GENOME + (genome_pos + mmoffset)) != NULL) ) {
							chromosome_director_neighbor = *(GENOME + (genome_pos + mmoffset));



							// Is the chrom director from actual read?
							if  (chromosome_director_neighbor < (CHROMOSOME_ENTRY_OPERATOR->entries+CHROMOSOME_ENTRY_OPERATOR->used) && (chromosome_director_neighbor->genome_pos == genome_pos + mmoffset)) {

								// is there a mapping entry in the right chromosome solt, mr. director?
								c = 0;
								// Search for chromosome_director for the correct chromosome
								while (chromosome_director_neighbor->next != NULL && (chromosome_director_neighbor->chromosome != genome_chr.nr() ||
											(chromosome_director_neighbor->chromosome == genome_chr.nr() && chromosome_director_neighbor->strand != strand)))
								{
									if (_config.STATISTICS && c == 0) _stats.listocc++;
									if (_config.STATISTICS) _stats.listcount++;
									chromosome_director_neighbor = chromosome_director_neighbor->next;
								}

								if (chromosome_director_neighbor->chromosome == genome_chr.nr() && chromosome_director_neighbor->strand == strand) {
									neighbor = chromosome_director_neighbor->mapping_entries;
									if (read_num == num) printf("MM neighbor hit: ");
									if (read_num == num) printhit(neighbor->hit);
									if (read_num == num) printf("  Found a potential entry, readops(neighbor)=%d, (actual)=%d\n",neighbor->readpos,readpos);
									if (neighbor->readpos == mapping_entry->readpos - _config.INDEX_DEPTH - 1) {
										if (read_num == num) printf("  Readpos matches\n");

										if ((neighbor->hit)->mismatches < _config.NUM_EDIT_OPS && (neighbor->hit)->mismatches-(neighbor->hit)->gaps < _config.NUM_MISMATCHES) {
											hit = neighbor->hit;
											oldlength = hit->end - hit->start + 1;
											mapping_entry->hit = neighbor->hit;
											if (read_num == num) printf("  Fancy! Mismatches < 4\n");

											if (reverse != 2) {
												hit->end = hit->end + _config.INDEX_DEPTH + 1;
												hit->edit_op[hit->mismatches].pos = readpos - 1;
												hit->edit_op[hit->mismatches].mm = 1;
											}
											else {
												hit->start = genome_pos+1;	// +1 weg
												hit->edit_op[hit->mismatches].pos = ((INT)_read.length()) - readpos + 2;
												hit->edit_op[hit->mismatches].mm = 1;
											}
											assert(hit->edit_op[hit->mismatches].pos>=0 && hit->edit_op[hit->mismatches].pos<=_read.length()) ;
											hit->mismatches++;
											assert(hit->mismatches<Config::MAX_EDIT_OPS) ;
											if (read_num == num) printf("  Mismatch at pos %d, #mm=%d\n",hit->edit_op[hit->mismatches-1].pos, hit->mismatches);
											if (read_num == num) printhit(hit);
										}
									}
								}
							}
						}
					}
				}

				TIME_CODE(time2a_part4 += clock()-start_time) ;
				TIME_CODE(start_time = clock()) ;

				// for MM=0: if potential hit doesn't start at readpos 1, it cannot become perfect, thus it is not even allocated:
				if ( !(_config.NUM_MISMATCHES == 0 && readpos != 1 && !_config.NOT_MAXIMAL_HITS) ) {

					// create new hit:
					if (hit == NULL) {

						hit = alloc_hit();
						if (!hit)
						{
							return -1 ;
							TIME_CODE(time2a_part5 += clock()-start_time) ;
							delete[] se_buffer;
						}

						mapping_entry->hit = hit;

						hit->chromosome = &genome_chr;
						hit->end = genome_pos + _config.INDEX_DEPTH;	// -1
						hit->orientation = (reverse != 2)? '+': '-';

						oldlength = _config.INDEX_DEPTH - 1;

						hit->readpos = readpos;
						hit->start = genome_pos+1;	// +1 weg

					}

				}

				// for MM=0: if hit doesn't start at readpos=1, then do not report hit (due to: all perfect hits are extended from a hit starting at readpos 1!!!)
				if ( hit != NULL && !(_config.NUM_MISMATCHES == 0 && hit->readpos != 1 && !_config.NOT_MAXIMAL_HITS) )
				{
					size_hit(hit, &oldlength, (read_num==num));
				}

				TIME_CODE(time2a_part5 += clock()-start_time) ;

				//printgenome();

				flag = 0;
				hit = NULL;

			} //end of for each seed on read

			delete[] se_buffer;

			if (report_repetitive_seeds)
			{
				//fprintf(stdout, "report %i/%i repetitive seeds\n", (int)extended_seedlist.size(), (int)index_entry.num) ;
				for (int ii=0; ii<extended_seedlist.size(); ii++)
					report_repetitive_seed(*extended_seedlist[ii].chr, extended_seedlist[ii].pos, extended_seedlist.size())  ;
			}
			
			//if (index_entry.num>_config.INDEX_DEPTH_EXTRA_THRESHOLD)
			//fprintf(stdout, "dropped %i/%i entries\n", dropped_entries, index_entry.num) ;
		}

		reverse -= 2; // 1->-1, 2->0, 4->2

	} //end of while (for each strand)


	return(1);
}

void printgenome()
{
	printf("G E N O M E: \n");
	unsigned int i,c;
	HIT *hit;
	CHROMOSOME_ENTRY *ce;
	for (i=0; i!=LONGEST_CHROMOSOME; ++i) {
		c=0;
		ce = *(GENOME+i);
		if (ce != NULL) {
			printf("%d: ",i);
			while (ce != NULL) {
				++c;
				printf("(%d, %d, %d, %d) ", ce->chromosome+1, ce->strand, ce->genome_pos, (ce->next!=NULL));
				hit = ce->mapping_entries->hit;
				printhit(hit);
				ce = ce->next;
			}
			printf("((%d))\n",c);
		}
	}
}

//@TODO remove num
int size_hit(HIT *hit, unsigned int *oldlength, char num)
{
	HIT *last, *next;

	// close the gap where the hit is taken out
	// shortest possible hits are not under control of the operator so far
	if (*oldlength > _config.INDEX_DEPTH - 1) {

		if (hit->last == NULL) { //hit is the first in the list, the operator must be re-linked
			*(HIT_LISTS_OPERATOR+*oldlength) = hit->next;
			next = hit->next;
			if (next != NULL) {
				next->last = NULL;
			}
		}
		else {
			//sweet
			last = hit->last;
			next = hit->next;
			last->next = next;
			if (next != NULL) {
				next->last = last;
			}
		}
	}

	else {
		_stats.NUM_HITS++;
	}


	unsigned int length = hit->end - hit->start + 1;
	// add to new list
	if (*(HIT_LISTS_OPERATOR+length) != NULL) {
		next = *(HIT_LISTS_OPERATOR+length);
		next->last = hit;
		hit->next = next;
		hit->last = NULL;
		*(HIT_LISTS_OPERATOR+length) = hit;
	}
	else {
		hit->last = NULL;
		hit->next = NULL;
		*(HIT_LISTS_OPERATOR+length) = hit;
	}

	if (length > LONGEST_HIT) {
		LONGEST_HIT = length;
	}

	return(1);
}


int browse_hits()
{
	HIT* hit;
	int i;
	char perfect = 0;

	// browse hit_list foreach hitlength:
	for (i=_read.length(); i!=_config.INDEX_DEPTH - 1; --i)
	{
		
		int hitlength = 0 ;

		// if only perfect reads should be reported, break earlier:
		if ((_config.NUM_EDIT_OPS == 0) && (i < _read.length()))
			break;

		// if hitlength limit is reached, break earlier:
		if (i == _config.HITLEN_LIMIT - 1) break;

		if ((*(HIT_LISTS_OPERATOR + i)) != NULL) 
		{
			hit = *(HIT_LISTS_OPERATOR + i);
			
			// foreach hit with hitlength i:
			while (hit != NULL) 
			{
				hitlength = hit->end - hit->start + 1;
				
				if (_config.STATISTICS) 
					_stats.HITS_PER_READ++;

				// ##########################################################
				// Mismatch at first or last pos of read:
				// ##########################################################
				
				// if hit.readpos == 2, then spare alignment since if first base is mm, it's cheaper than a gap
				if (hit->readpos == 2) {
					if (_config.NOT_MAXIMAL_HITS || hit->mismatches < _config.NUM_MISMATCHES) {
						if (hit->orientation == '+' && hit->start != 1) {
							if (!_config.NOT_MAXIMAL_HITS || check_mm(*hit->chromosome, hit->start-2, 0, 1)) {
								hit->edit_op[hit->mismatches].pos = 1;
								hit->edit_op[hit->mismatches].mm = 1;
								hit->mismatches++;
							}
							hit->start--;
							hit->readpos--;
						}
						else if (hit->orientation == '-' && hit->end != hit->chromosome->length()) {
							if (!_config.NOT_MAXIMAL_HITS || check_mm(*hit->chromosome, hit->end, 0, -1)) {
								hit->edit_op[hit->mismatches].pos = _read.length();
								hit->edit_op[hit->mismatches].mm = 1;
								hit->mismatches++;
							}
							hit->end++;
							hit->readpos--;
						}
					}
					else {
						hit = hit->next;
						continue;
					}

					// update length:
					hitlength = hit->end - hit->start + 1;
				}

				// if hit ends at pos |read|-1, then spare alignment since if last base is mm, it's cheaper than a gap
				if (hit->readpos + hitlength == _read.length()) {
					if (_config.NOT_MAXIMAL_HITS || hit->mismatches < _config.NUM_MISMATCHES) {
						if (hit->orientation == '+' && hit->end != hit->chromosome->length()) {
							if (!_config.NOT_MAXIMAL_HITS || check_mm(*hit->chromosome, hit->end, _read.length()-1, 1)) {
								hit->edit_op[hit->mismatches].pos = _read.length();
								hit->edit_op[hit->mismatches].mm = 1;
								hit->mismatches++;
							}
							hit->end++;
						}
						else if (hit->orientation == '-' && hit->start != 1) {
							if (!_config.NOT_MAXIMAL_HITS || check_mm(*hit->chromosome, hit->start-2, _read.length()-1, -1)) {
								hit->edit_op[hit->mismatches].pos = 1;
								hit->edit_op[hit->mismatches].mm = 1;
								hit->mismatches++;
							}
							hit->start--;
						}
					}
					else {
						hit = hit->next;
						continue;
					}

					// update length:
					hitlength = hit->end - hit->start + 1;
				}
				// ##########################################################


				// 1) Hit spans the whole read:
				if (hitlength == _read.length())
				{
					if (_config.STATISTICS && hit->mismatches == 0) {
						
						// reporting perfect matching reads (only one count per read)
						if (!perfect) 
						{
							_stats.PERFECT_READS++;
							perfect = 1;
						}
						
						// reporting perfect hits
						if (hit->orientation == '+') 
							_stats.PERFECT_HITS++;
						else 
							_stats.PERFECT_HITS_REV++;
					}
					
					// report match:
					if (hit->mismatches <= _config.NUM_EDIT_OPS) 
					{
						if (_config.REPORT_MAPPED_REGIONS && hitlength >= Config::REPORT_MAPPED_REGIONS_MIN_LENGTH)
						{
							report_mapped_region(*hit->chromosome, hit->start, hit->end, hitlength - hit->mismatches) ;
						}
						
						// insert hit into HITS_BY_SCORE
						int ret = insert_into_scorelist(hit, 1) ;
						if (ret<0)
							return ret ;
						
						//printed = 1;
						if (_config.STATISTICS) _stats.NOT_ALIGNED[0]++;
						
						//if (!_config.ALL_HIT_STRATEGY && hit->mismatches < _config.NUM_EDIT_OPS)
						//	_config.NUM_EDIT_OPS = hit->mismatches;
						update_num_edit_ops(hit->mismatches, _config.ALL_HIT_STRATEGY, _config.NUM_EDIT_OPS) ;
					}
				}
				// 2) Hit has to be aligned:
				else 
				{
					unsigned int readstart;
					if (hit->orientation == '+') {
						readstart = hit->start - hit->readpos;	// start pos of read in genome	0-initialized
					}
					else {
						readstart = hit->start - (((int)_read.length()) - hit->readpos - hitlength + 2); 	// 0-initialized
					}

					// Alignment
					if (hit->mismatches < _config.NUM_EDIT_OPS || _config.NOT_MAXIMAL_HITS)
					{
						//if (_config.STATISTICS) (*(*(_stats.HITS_READPOS+(hitlength-_config.INDEX_DEPTH))+(hit->readpos-1)))++;

						if (_config.REPORT_MAPPED_REGIONS && hitlength>=Config::REPORT_MAPPED_REGIONS_MIN_LENGTH)
						{
							report_mapped_region(*hit->chromosome, hit->start, hit->end, hitlength - hit->mismatches) ;
						}

						if (_config.NUM_GAPS != 0) 
						{
							// KBOUND:
							int ret = prepare_kbound_alignment(hit, hit->start, hit->end, hit->readpos, *hit->chromosome, hit->orientation, hit->mismatches);
							if (ret<0)
								return ret ;
							hit->aligned = ret ;
						}
						else {
							// SIMPLE:
							int ret = align_hit_simple(hit, hit->start, hit->end, hit->readpos, *hit->chromosome, hit->orientation, hit->mismatches);
							if (ret<0)
								return ret ;
							hit->aligned = ret ;
						}
							
					}
					
				} // else has mismatches
				
				
				if (_config.STATISTICS) _stats.HITS_LEN[hitlength]++;

				hit = hit->next;

			} // while hitlist not empty

		}

	} //for each hitlength


	// store successfully aligned hits in HITS_BY_SCORE list for printout
	// by iterating a second time over HIT_LIST:
	for (i=_read.length(); i!=_config.INDEX_DEPTH - 1; --i) {

		// if only perfect reads should be reported, break earlier:
		if ((_config.NUM_EDIT_OPS == 0) && (i < _read.length()))
			return 1;

		// if hitlength limit is reached, break earlier:
		if (i == _config.HITLEN_LIMIT - 1)
			return 1;

		if ((*(HIT_LISTS_OPERATOR + i)) != NULL) {

			hit = *(HIT_LISTS_OPERATOR + i);

			// foreach hit with hitlength i:
			while (hit != NULL) {

				if (hit->aligned)
				{
					int ret = insert_into_scorelist(hit, 1);
					if (ret<0)
						return ret ;
				}
				hit = hit->next;
			}

		}

	}

	return 1;
}

/** Inserts a hit into the HITS_BY_SCORE list, which contains bins of hits that
 *  have the same score.
 *
 */

int insert_into_scorelist(HIT* hit, char d)
{
	if (d)
	{
		int ret = duplicate(hit) ;
		if (ret>0)
			return 0;
		else if (ret<0)
			return ret ;
	}


	int interval = (hit->mismatches-hit->gaps) * _config.MM_SCORE + hit->gaps * _config.GAP_SCORE - (((int)_read.length())-hit->mismatches) * _config.M_SCORE;
	if (HITS_BY_SCORE[interval].num == 0) {
		// first entry in list
		HITS_BY_SCORE[interval].hitpointer = hit;
	}
	else
	{
		// list has already some entries, insert to the front
		HIT *tmp_hit = HITS_BY_SCORE[interval].hitpointer;
		HITS_BY_SCORE[interval].hitpointer = hit;
		hit->same_eo_succ = tmp_hit;
	}
	HITS_BY_SCORE[interval].num++;
	HITS_IN_SCORE_LIST++;

	//if (!_config.ALL_HIT_STRATEGY && hit->mismatches < _config.NUM_EDIT_OPS)
	//_config.NUM_EDIT_OPS = hit->mismatches;
	update_num_edit_ops(hit->mismatches, _config.ALL_HIT_STRATEGY, _config.NUM_EDIT_OPS) ;

	return 1;
}


int duplicate(HIT* hit)
{
	int hitlength = hit->end - hit->start + 1;
	unsigned int readstart, readend;
	char strand;

	if (hit->orientation == '+') {
		readstart = hit->start - hit->readpos + hit->start_offset;	// start pos of read in genome	   0-initialized
		readend = hit->end + (((int)_read.length()) - hit->readpos - hitlength) + hit->end_offset;			// 0-initialized
		strand = 1;
	}
	else {
		readstart = hit->start - (((int)_read.length()) - hit->readpos - hitlength + 2) + hit->start_offset; 	// 0-initialized
		readend = hit->end + hit->readpos - 2 + hit->end_offset;									// 0-initialized
		strand = -1;
	}

	CHROMOSOME_ENTRY *chromosome_director, *chromosome_director_new;
	MAPPING_ENTRY *mapping_entry, *existing_entry;
	char flag = 0;

	if (*(GENOME+readstart) == NULL) {
				flag = 1;
				chromosome_director = alloc_chromosome_entry(readstart, *hit->chromosome, strand);
				if (!chromosome_director)
					return(-1);	// chrom_container_size too small -> cancel this read
				*(GENOME + readstart) = chromosome_director;
	}
	else {

		chromosome_director = *(GENOME + readstart);

		// is chrom_director from the actual read or from a former one?
		if (chromosome_director >= (CHROMOSOME_ENTRY_OPERATOR->entries+CHROMOSOME_ENTRY_OPERATOR->used)
			|| 	(chromosome_director < (CHROMOSOME_ENTRY_OPERATOR->entries+CHROMOSOME_ENTRY_OPERATOR->used) &&
		   		(chromosome_director->genome_pos != readstart))) {

		   			// it's from a former read, thus we need a new chrom_director:
					chromosome_director = alloc_chromosome_entry(readstart, *hit->chromosome, strand);

					if (!chromosome_director)
						return(-1);	// chrom_container_size too small -> cancel this read
					*(GENOME + readstart) = chromosome_director;

		}
	}


	int c=0;
 	// Search for chromosome_director for the correct chromosome
 	//printf("%d %d %d %d %d\n",(chromosome_director->next != NULL), chromosome_director->chromosome, hit->chromosome, chromosome_director->strand, strand);
	while (chromosome_director->next != NULL
			&& (chromosome_director->chromosome != hit->chromosome->nr()
					|| (chromosome_director->chromosome == hit->chromosome->nr()
							&& chromosome_director->strand != strand))) {
		if (_config.STATISTICS && c == 0)
			_stats.listocc++;
		if (_config.STATISTICS)
			_stats.listcount++;
		chromosome_director = chromosome_director->next;
	}

	if (chromosome_director->chromosome != hit->chromosome->nr() || chromosome_director->strand != strand)
	{
	    chromosome_director_new = alloc_chromosome_entry(readstart, *hit->chromosome, strand);
        if (!chromosome_director_new)
			return(-1);
		
     	chromosome_director_new->next = chromosome_director->next;
        chromosome_director->next = chromosome_director_new;
        chromosome_director = chromosome_director_new;
	}

	if (chromosome_director->mapping_entries == NULL)
	{
		mapping_entry = alloc_mapping_entry();
		if (!mapping_entry)
			return(-1) ;
		chromosome_director->mapping_entries = mapping_entry;

		mapping_entry->readpos = readend;
		mapping_entry->hit = hit;
	}
	else {

		if (flag == 1) {
			printf("!!!!!!!!!!!Found entry in chr dir at genome pos %d  --  %p\n", readstart, chromosome_director->mapping_entries	);
			exit(1);
		}

		double score2, score1 = (hit->mismatches-hit->gaps) * _config.MM_SCORE + hit->gaps * _config.GAP_SCORE - (((int)_read.length()) - hit->mismatches) * _config.M_SCORE;

		existing_entry = chromosome_director->mapping_entries;

		while (existing_entry != NULL) {
			score2 = (existing_entry->hit->mismatches-existing_entry->hit->gaps) * _config.MM_SCORE
					 + existing_entry->hit->gaps * _config.GAP_SCORE - (((int)_read.length()) - existing_entry->hit->mismatches) * _config.M_SCORE;

			if (existing_entry->readpos == readend && score1 == score2) {
				if (_config.STATISTICS) REDUNDANT++;
				return 1;
			}

			existing_entry = existing_entry->succ;
		}

		// no hit with same end pos of read and same score could have been found -> create another entry:
		existing_entry = chromosome_director->mapping_entries;
		mapping_entry = alloc_mapping_entry();
		if (!mapping_entry)
			return -1 ;

		mapping_entry->succ = existing_entry;
		existing_entry->pred = mapping_entry;
		//insert to the front of the mapping list
		chromosome_director->mapping_entries = mapping_entry;
		mapping_entry->readpos = readend;
		mapping_entry->hit = hit;

	}

	return 0;
}

// for debugging
char *get_seq(unsigned int n)
{
	char *seq = (char *) malloc ((_config.INDEX_DEPTH+1)*sizeof(char));
	if (seq == NULL) {
		fprintf(stderr, "[get_seq] Could not allocate memory (_read.id() = %s)\n", _read.id());
		exit(1);
	}
	int i, c;
	for (i=_config.INDEX_DEPTH-1; i>=0; --i) {
		c = (int) (n / Util::POWER[i]);
		switch (c)
		{
			case 0: seq[i] = 'A';
					break;
			case 1: seq[i] = 'C';
					break;
			case 2: seq[i] = 'G';
					break;
			case 3: seq[i] = 'T';
					break;
		}
		n -= (int) (c * Util::POWER[i]);
	}
	seq[_config.INDEX_DEPTH] = '\0';
	return seq;
}

