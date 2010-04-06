// Authors: Korbinian Schneeberger, Joerg Hagmann, Gunnar Raetsch
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany
// Copyright (C) 2009-2010 by Friedrich Miescher Laboratory, Tuebingen, Germany

#include <limits.h>
#include "genomemapper.h"
#include "Hits.h"
#include "align.h"
#include "print.h"

//config

Hits::Hits(Genome &genome_,	GenomeMaps &genomemaps_) {
/*
	readpos = 0;
	start = 0;
	end = 0;
	chromosome = 0;
	mismatches = 0;
	gaps = 0;
	start_offset = 0;
	end_offset = 0;
	aligned = 0;*/
	/*int i;
	//for (i=0; i!=MAX_MISMATCHES; ++i) mismatch[i] = READ_LENGTH + 1;
	for (i=0; i!=MAX_MISMATCHES; ++i) {
		edit_op[i].mm = 0;
	}*/
/*	orientation = ' ';
	same_eo_succ = NULL;
	next = NULL;
	last = NULL;
*/
	genome = &genome_;
	genomemaps = &genomemaps_;
	REDUNDANT = 0;
	SLOT = 0;
	HAS_SLOT = 0;
	// initialize with meta information
	init_from_meta_index(); // updated
	init_alignment_structures(&_config);
	init_hit_lists();
}

// for debugging purposes only

void printhit(HIT* hit) {
	printf("[%i: %i-%i/%c/chr%i/rp%i/%imm", hit->end - hit->start + 1,
			hit->start, hit->end, hit->orientation, hit->chromosome->nr() + 1,
			hit->readpos, hit->mismatches);
	if (hit->mismatches != 0) {
		printf(":");
		int i;
		for (i = 0; i != hit->mismatches; ++i) 
			printf(" %i%c", hit->edit_op[i].pos, (hit->edit_op[i].mm) ? 'M' : 'G');
	}
	printf("/ID: %s]", _read.id());
}

int Hits::init_alignment_structures(Config * config) {

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
				config->OUTPUT_FILTER==OUTPUT_FILTER_ALL ? "all hit" : ((config->OUTPUT_FILTER==OUTPUT_FILTER_ALL) ? "top hit" : "all hit"),
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


void Hits::printhits() 
{
	//printf("list:\n");
	//printf("print hitlist with readlength %i, read %s, last[rl-2]=%c\n",((int)_read.lenght()), READ, READ[((int)_read.lenght())-2]);
	int i;
	int c = 0;
	HIT* hit;
	for (i = _config.INDEX_DEPTH; i != ((int)_read.length()) + 1; ++i) {

		if (*(HIT_LISTS_OPERATOR + i) != NULL) {
			printf("%i: ", i);
			hit = *(HIT_LISTS_OPERATOR + i);
			do {
				//if (hit->orientation == '-') {
				printf("[%i: %i-%i/%c/chr%i/rp%i/%imm", hit->end - hit->start
						+ 1, hit->start, hit->end, hit->orientation,
						hit->chromosome->nr() + 1, hit->readpos, hit->mismatches);
				if (hit->mismatches != 0) {
					printf(":");
					int j;
					for (j = 0; j != hit->mismatches; ++j)
						printf(" %i%c", hit->edit_op[j].pos,
								(hit->edit_op[j].mm) ? 'M' : 'G');
				}
				printf("/%s]", _read.id());
				++c;
				//}
				hit = hit->next;
				if (c > 4)
					hit = NULL;
			} while (hit != NULL);
			printf("\n");
		}
		c = 0;
	}
	//printf("done\n");
}



#define TIME_CODE(x) 

int size_hit(HIT *hit, unsigned int *oldlength, char num);
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

int Hits::map_reads(TopAlignments * topalignments, QPalma* qpalma)
{
	char eof = 0, read_mapped;
	unsigned int count_reads = 0;
	int first_slot = 0, first_pos = 0 ;
	int num_edit_ops = _config.NUM_EDIT_OPS;
	if (_config.OUTPUT_FILTER==OUTPUT_FILTER_TOP)
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
		unsigned int rtrim_cut = 0 ;
		unsigned int polytrim_cut_start = 0 ;
		unsigned int polytrim_cut_end = 0 ;
		unsigned int polytrim_cut_start_curr = 0 ;
		unsigned int polytrim_cut_end_curr = 0 ;
		unsigned int poly_length_start=0 ;
		unsigned int poly_length_end=0 ;
		Read* poly_orig_read = NULL ;

		clock_t start_time = clock() ;

		eof = _read.read_short_read(QUERY_FP);

	restart:

		// make it somehow dependent on the read length, the index depth and the number of mismatches 
		genomemaps->REPORT_REPETITIVE_SEED_DEPTH_EXTRA = _read.length() - _config.INDEX_DEPTH - _config.NUM_MISMATCHES  ;
	
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
		for (unsigned int i=0; i<_read.length(); i++)
			if (READ[i]!='A' && READ[i]!='C' && READ[i]!='G' && READ[i]!='T')
				num_N++ ;
		if (num_N>_config.NUM_MISMATCHES)
		{
			if (_config.VERBOSE)
				fprintf(stdout, "read has %i non-ACGT characters, allowing only %i mismatches -> skip read\n", num_N, _config.NUM_MISMATCHES) ;
			continue ;
		}

		if (_read.length() < _config.HITLEN_LIMIT) 
		{
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
 			if (!_config.ALL_HIT_STRATEGY || _config.OUTPUT_FILTER==OUTPUT_FILTER_TOP || (_config.NUM_MISMATCHES < nr_seeds && _config.NUM_GAPS == 0)) 
			{
				int ret	= map_fast(_read, first_slot, first_pos);	// if no hits could have been found: _config.ALL_HIT_STRATEGY = -1, necessitating execution of normal mapping in the following
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
				if (((_config.ALL_HIT_STRATEGY!=0) || _config.NOT_MAXIMAL_HITS ||
					 (_config.OUTPUT_FILTER==OUTPUT_FILTER_TOP && SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.size()==0)) &&
				    (!(_config.NUM_MISMATCHES < nr_seeds && _config.NUM_GAPS == 0) ) )
				{

					c_map_short_read++;
					int ret = map_short_read(_read, count_reads, first_slot, first_pos);
					
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
			if (_config.OUTPUT_FILTER==OUTPUT_FILTER_TOP)
				reset_num_edit_ops(num_edit_ops) ;
					
			read_mapped = 0 ;
			
			if (!cancel && !_config.REPORT_REPETITIVE_SEEDS)
			{
				topalignments->start_top_alignment_record();

				read_mapped = analyze_hits(topalignments, qpalma);	// returns 1 if at least one hit is printed, 0 otherwise
				if (_config.VERBOSE)
					printf("%i unspliced alignment found\n", (int)topalignments->size()); 

				bool trigger = false ;
				if (_config.SPLICED_HITS || _config.LOG_TRIGGERED)
					trigger = topalignments->size()==0 || 
						qpalma->qpalma_filter(topalignments->get_alignment(0), num_N)!=0 ;

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
								int ret = qpalma->capture_hits();
								//fprintf(stderr, "capture_hits ret=%i\n", ret) ;
								if (ret<0)
									cancel=4 ;
								if (_config.VERBOSE)
									fprintf(stdout, "capture_hits generated %i alignments\n", ret) ;
								if (FILTER_STAT)
									qpalma->qpalma_filter_stat(ret>0) ;
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
				if (topalignments->size()>0)
					read_mapped = 1 ;
				//if (_config.VERBOSE && read_mapped)
				//	printf("unspliced or spliced alignment found\n"); 
				
				topalignments->end_top_alignment_record(rtrim_cut, polytrim_cut_start_curr, polytrim_cut_end_curr);
				
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
							_read.set_orig(NULL) ;
							delete poly_orig_read ;
							poly_orig_read=new Read(_read) ;
							_read.set_orig(poly_orig_read) ;
							if (_read.is_full_poly())
								poly_length_start=poly_length_end=0 ;
						}
						assert(poly_orig_read!=NULL) ;

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
								_read.trim_read_start(poly_orig_read, polytrim_cut_start) ;
								restart = true ;
								polytrim_cut_start_curr = polytrim_cut_start ;
								polytrim_cut_end_curr = 0 ;
							}
							if (end_cond && !restart)
							{
								polytrim_cut_end += _config.POLYTRIM_STRATEGY_STEP ;
								_read.trim_read_end(poly_orig_read, polytrim_cut_end) ;
								restart = true ;
								polytrim_cut_start_curr = 0 ;
								polytrim_cut_end_curr = polytrim_cut_end ;
							}
						}
						if (restart)
						{
							//fprintf(stdout, "polytrim_cut_start_curr=%i, polytrim_cut_end_curr=%i: %s\n", polytrim_cut_start_curr, polytrim_cut_end_curr, _read.data()) ;
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
			_read.set_orig(NULL) ;
			delete poly_orig_read ;

			if (_config.VERBOSE || ((clock()-last_timing_report)/CLOCKS_PER_SEC>=10))
			{
				last_timing_report = clock() ;
				map_reads_timing(count_reads, ((float)clock()-start_time)/CLOCKS_PER_SEC) ;
			}
		}
		
	}

	map_reads_timing(count_reads) ;
	qpalma->capture_hits_timing() ;

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


unsigned int extend_seed(int direction, unsigned int seed_depth_extra, Chromosome &chr, int genome_pos, int readpos)
{
	unsigned int e = 0 ;
	if (direction==-1)	// forward strand
	{
		for (e=0; e<seed_depth_extra; e++)
		{
			int readpos_ = readpos+_config.INDEX_DEPTH+e-1 ;
			if (readpos_>=(int)_read.length())
				break ;
			int genome_pos_ = genome_pos+_config.INDEX_DEPTH+e ;
			if (genome_pos_>=(int)chr.length())
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
			if (readpos_>=(int)_read.length())
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

int Hits::seed2genome(unsigned int num, unsigned int index_slot, unsigned int readpos, char reverse)
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
	unsigned int read_num = INT_MAX ;
	char flag = 0;

	CHROMOSOME_ENTRY *chromosome_director, *chromosome_director_neighbor, *chromosome_director_new;
	MAPPING_ENTRY *mapping_entry, *existing_entry, *neighbor;

	HIT *hit = NULL;

	unsigned int i;
	unsigned int oldlength = _config.INDEX_DEPTH-1;

	// reverse = 1: only index, 2: only index_rev, 4: both
	while (reverse > 0) {

		if (reverse != 2) {
			index_mmap = genome->INDEX_FWD_MMAP;
			index_entry = *(genome->INDEX+index_slot);
			direction = -1;
		}
		else {
			index_mmap = genome->INDEX_REV_MMAP;
			index_entry = *(genome->INDEX_REV+index_slot);
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
			unsigned int index_entry_num=index_entry.num ;
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
			
			if (index_entry.num > _config.SEED_HIT_CANCEL_THRESHOLD && !report_repetitive_seeds)
				index_entry_num=0 ;
			{
				TIME_CODE(clock_t start_time = clock()) ;
#ifndef BinaryStream_MAP
				genome->index_pre_buffer(index_mmap, se_buffer, index_entry.offset-index_entry_num, index_entry_num);
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

				unsigned int genome_pos = pos + genome->BLOCK_TABLE[block].pos;
				//unsigned int genome_chr = BLOCK_TABLE[block].chr;
				Chromosome &genome_chr = genome->chromosome(genome->BLOCK_TABLE[block].chr);

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
					int e = extend_seed(direction, genomemaps->REPORT_REPETITIVE_SEED_DEPTH_EXTRA, genome_chr, genome_pos, readpos) ;

					if (e==genomemaps->REPORT_REPETITIVE_SEED_DEPTH_EXTRA)
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
					unsigned int ee = extend_seed(direction, _config.INDEX_DEPTH_EXTRA, genome_chr, genome_pos, readpos) ;
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
				while (chromosome_director->next != NULL && (chromosome_director->chromosome != (int)genome_chr.nr() ||
															 (chromosome_director->chromosome == (int)genome_chr.nr() && chromosome_director->strand != strand)))
				{
					if (_config.STATISTICS && c == 0) _stats.listocc++;
					if (_config.STATISTICS) _stats.listcount++;
					chromosome_director = chromosome_director->next;
				}

				// Chromosome director is set, but still it could be the wrong chromosome, if the right chromosome is not in there so far.
				if (chromosome_director->chromosome != (int)genome_chr.nr() || chromosome_director->strand != strand) {
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
							while (chromosome_director_neighbor->next != NULL && (chromosome_director_neighbor->chromosome != (int)genome_chr.nr() ||
																				  (chromosome_director_neighbor->chromosome == (int)genome_chr.nr() && chromosome_director_neighbor->strand != strand)))
							{
								if (_config.STATISTICS) _stats.listcount++;
								if (_config.STATISTICS && c == 0) _stats.listocc++;
								chromosome_director_neighbor = chromosome_director_neighbor->next;
							}

							if 	(chromosome_director_neighbor->chromosome == (int)genome_chr.nr() && 
								 chromosome_director_neighbor->strand == strand)
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
					if (/*(!_config.NOT_MAXIMAL_HITS) &&*/ hit == NULL && _config.NUM_MISMATCHES != 0) {
						if (read_num == num) printf("Now checking if hit can be extended over mismatch\n");
						mmoffset = (reverse != 2)? -_config.INDEX_DEPTH - 1: _config.INDEX_DEPTH + 1;

						if ( (genome_pos + mmoffset > 0) && (genome_pos + mmoffset < genome_chr.length()) && (*(GENOME + (genome_pos + mmoffset)) != NULL) ) {
							chromosome_director_neighbor = *(GENOME + (genome_pos + mmoffset));


							// Is the chrom director from actual read?
							if  (chromosome_director_neighbor < (CHROMOSOME_ENTRY_OPERATOR->entries+CHROMOSOME_ENTRY_OPERATOR->used) && (chromosome_director_neighbor->genome_pos == genome_pos + mmoffset)) {

								// is there a mapping entry in the right chromosome solt, mr. director?
								c = 0;
								// Search for chromosome_director for the correct chromosome
								while (chromosome_director_neighbor->next != NULL && 
									   (chromosome_director_neighbor->chromosome != (int)genome_chr.nr() ||
										(chromosome_director_neighbor->chromosome == (int)genome_chr.nr() && 
										 chromosome_director_neighbor->strand != strand)))
								{
									if (_config.STATISTICS && c == 0) _stats.listocc++;
									if (_config.STATISTICS) _stats.listcount++;
									chromosome_director_neighbor = chromosome_director_neighbor->next;
								}

								if (chromosome_director_neighbor->chromosome == (int)genome_chr.nr() && 
									chromosome_director_neighbor->strand == strand) {
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
											assert(hit->edit_op[hit->mismatches].pos>=0 && hit->edit_op[hit->mismatches].pos<=(int)_read.length()) ;
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
				if ( !(_config.NUM_MISMATCHES == 0 && readpos != 1 && (!_config.NOT_MAXIMAL_HITS)) ) {

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
				for (unsigned int ii=0; ii<extended_seedlist.size(); ii++)
					genomemaps->report_repetitive_seed(*extended_seedlist[ii].chr, extended_seedlist[ii].pos, extended_seedlist.size())  ;
			}
			
			//if (index_entry.num>_config.INDEX_DEPTH_EXTRA_THRESHOLD)
			//fprintf(stdout, "dropped %i/%i entries\n", dropped_entries, index_entry.num) ;
		}

		reverse -= 2; // 1->-1, 2->0, 4->2

	} //end of while (for each strand)


	return(1);
}

void Hits::printgenome()
{
	printf("G E N O M E: \n");
	unsigned int i,c;
	HIT *hit;
	CHROMOSOME_ENTRY *ce;
	for (i=0; i!=genome->LONGEST_CHROMOSOME; ++i) {
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
int Hits::size_hit(HIT *hit, unsigned int *oldlength, char num)
{
	HIT *last, *next;

	// close the gap where the hit is taken out
	// shortest possible hits are not under control of the operator so far
	if ((int)*oldlength > (int)(_config.INDEX_DEPTH) - 1) {

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


int Hits::browse_hits()
{
	HIT* hit;
	int i;
	char perfect = 0;

	// browse hit_list foreach hitlength:
	for (i=_read.length(); i!=(int)_config.INDEX_DEPTH - 1; --i)
	{
		
		int hitlength = 0 ;

		// if only perfect reads should be reported, break earlier:
		if ((_config.NUM_EDIT_OPS == 0) && (i < (int)_read.length()))
			break;

		// if hitlength limit is reached, break earlier:
		if (i == (int)(_config.HITLEN_LIMIT) - 1) break;

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
							if (!_config.NOT_MAXIMAL_HITS || check_mm(_read, *hit->chromosome, hit->start-2, 0, 1)) 
							{
								hit->edit_op[hit->mismatches].pos = 1;
								hit->edit_op[hit->mismatches].mm = 1;
								hit->mismatches++;
							}
							hit->start--;
							hit->readpos--;
						}
						else if (hit->orientation == '-' && hit->end != hit->chromosome->length()) {
							if (!_config.NOT_MAXIMAL_HITS || check_mm(_read, *hit->chromosome, hit->end, 0, -1)) {
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
				if (hit->readpos + hitlength == (int)_read.length()) {
					if (_config.NOT_MAXIMAL_HITS || hit->mismatches < _config.NUM_MISMATCHES) {
						if (hit->orientation == '+' && hit->end != hit->chromosome->length()) {
							if (!_config.NOT_MAXIMAL_HITS || check_mm(_read, *hit->chromosome, hit->end, _read.length()-1, 1)) {
								hit->edit_op[hit->mismatches].pos = _read.length();
								hit->edit_op[hit->mismatches].mm = 1;
								hit->mismatches++;
							}
							hit->end++;
						}
						else if (hit->orientation == '-' && hit->start != 1) {
							if (!_config.NOT_MAXIMAL_HITS || check_mm(_read, *hit->chromosome, hit->start-2, _read.length()-1, -1)) {
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
				if (hitlength == (int)_read.length())
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
							genomemaps->report_mapped_region(*hit->chromosome, hit->start, hit->end, hitlength - hit->mismatches) ;
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
							genomemaps->report_mapped_region(*hit->chromosome, hit->start, hit->end, hitlength - hit->mismatches) ;
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
	for (i=_read.length(); i!=(int)(_config.INDEX_DEPTH) - 1; --i) {

		// if only perfect reads should be reported, break earlier:
		if ((_config.NUM_EDIT_OPS == 0) && (i < (int)_read.length()))
			return 1;

		// if hitlength limit is reached, break earlier:
		if (i == (int)(_config.HITLEN_LIMIT) - 1)
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

int Hits::insert_into_scorelist(HIT* hit, char d)
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


int Hits::duplicate(HIT* hit)
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
		   && (chromosome_director->chromosome != (int)hit->chromosome->nr()
			   || (chromosome_director->chromosome == (int)hit->chromosome->nr()
							&& chromosome_director->strand != strand))) {
		if (_config.STATISTICS && c == 0)
			_stats.listocc++;
		if (_config.STATISTICS)
			_stats.listcount++;
		chromosome_director = chromosome_director->next;
	}

	if (chromosome_director->chromosome != (int)hit->chromosome->nr() || chromosome_director->strand != strand)
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

// alloc

int Hits::alloc_genome_memory()
{
	if ((GENOME = (CHROMOSOME_ENTRY **) calloc (genome->LONGEST_CHROMOSOME, sizeof(CHROMOSOME_ENTRY**))) == NULL) {
		fprintf(stderr, "ERROR : not enough memory for genome memory\n");
		exit(1);
	}
		
	//@TODO is this really necessary? why isn't it already NULL
	unsigned int i;
	for (i=0; i!=genome->LONGEST_CHROMOSOME; ++i) {
		*(GENOME+i) = NULL;
	}

	return(0);
}

int Hits::init_from_meta_index() 
{
	genome->init_constants(); // updated

	alloc_genome_memory(); // updated

	init_operators(); //updated

	if (_config.STATISTICS) {
		init_statistic_vars(); //updated
	}
	if (!_config.HITLEN_LIMIT)
		_config.HITLEN_LIMIT = _config.INDEX_DEPTH;

	return (0);
}

int Hits::init_operators() {
	MAPPING_ENTRY_OPERATOR = alloc_mapping_entry_container();
	assert(MAPPING_ENTRY_OPERATOR!=NULL) ;

	MAPPING_ENTRY_OPERATOR_FIRST = MAPPING_ENTRY_OPERATOR;

	HIT_OPERATOR = alloc_hit_container();
	assert(HIT_OPERATOR!=NULL) ;
	HIT_OPERATOR_FIRST = HIT_OPERATOR;

	alloc_chromosome_entry_container();

	return (0);
}

int Hits::init_statistic_vars() {
	new (&_stats) Statistics();
	return (0);
}


int Hits::init_hit_lists() {
	alloc_hit_lists_operator();
	alloc_hits_by_score(); // A T T E N T I O N ! ! !   needs correct _config.NUM_EDIT_OPS which can be changed in init_alignment_structures() !!!

	return (0);
}



MAPPING_ENTRY* Hits::alloc_mapping_entry()
{
	MAPPING_ENTRY_CONTAINER *container;
	MAPPING_ENTRY *entry;

	if (MAPPING_ENTRY_OPERATOR->used >= CONTAINER_SIZE - 1) {
		container = alloc_mapping_entry_container();
		if (!container)
			return NULL ;
		MAPPING_ENTRY_OPERATOR->next = container;
		MAPPING_ENTRY_OPERATOR = container;
	}
	
	entry = &(MAPPING_ENTRY_OPERATOR->entries[MAPPING_ENTRY_OPERATOR->used]);
	MAPPING_ENTRY_OPERATOR->used++;

	entry->hit = NULL;
	entry->pred = NULL;
	entry->succ = NULL;
	entry->readpos = -1;
	
	return(entry);
}

MAPPING_ENTRY_CONTAINER* Hits::alloc_mapping_entry_container()
{
	MAPPING_ENTRY_CONTAINER *container;

	if ((container =  (MAPPING_ENTRY_CONTAINER *) calloc(1, sizeof(MAPPING_ENTRY_CONTAINER))) == NULL) {
		//fprintf(stderr, "ERROR : not enough memory for mapping entry container memory (alloc_mapping_entry_container)\n");
		fprintf(stderr, "WARNING : not enough memory for mapping entry container memory (alloc_mapping_entry_container)\n");
		return NULL ;
	}

	container->used = 0;
	container->next = NULL;

	return(container);
}


CHROMOSOME_ENTRY* Hits::alloc_chromosome_entry(unsigned int pos, Chromosome const &chr, char strand)
{
	if (CHROMOSOME_ENTRY_OPERATOR->used > _config.CHROM_CONTAINER_SIZE - 1) {
		fprintf(stderr, "\n!!! WARNING: Chromosome container size of %d is too small! Hits for read %s cannot be reported any more!\n\n", _config.CHROM_CONTAINER_SIZE, _read.id());
		return(NULL);
	}
	
	CHROMOSOME_ENTRY *entry;

	entry = &(CHROMOSOME_ENTRY_OPERATOR->entries[CHROMOSOME_ENTRY_OPERATOR->used]);
	entry->chromosome = chr.nr();
	entry->genome_pos = pos;
	entry->strand = strand;
	entry->next = NULL;
    entry->mapping_entries = NULL;
    
    CHROMOSOME_ENTRY_OPERATOR->used++;

    if (_config.STATISTICS && CHROMOSOME_ENTRY_OPERATOR->used > MAX_USED_SLOTS)
    	MAX_USED_SLOTS = CHROMOSOME_ENTRY_OPERATOR->used;

	return(entry);
}

CHROMOSOME_ENTRY_CONTAINER* Hits::alloc_chromosome_entry_container()
{

    if ((CHROMOSOME_ENTRY_OPERATOR = (CHROMOSOME_ENTRY_CONTAINER *) calloc (1, sizeof(CHROMOSOME_ENTRY_CONTAINER))) == NULL) {
    	fprintf(stderr, "ERROR : not enough memory for chromosome entry container memory\n");
		exit(1);
	}

	if ((CHROMOSOME_ENTRY_OPERATOR->entries = (CHROMOSOME_ENTRY *) calloc (_config.CHROM_CONTAINER_SIZE, sizeof(CHROMOSOME_ENTRY))) == NULL) {
		fprintf(stderr, "ERROR : not enough memory for chromosome entry container memory\n");
		exit(1);
	}
	
	CHROMOSOME_ENTRY_OPERATOR->used = 0;

	return(CHROMOSOME_ENTRY_OPERATOR);
}

HIT* Hits::alloc_hit()
{
	HIT_CONTAINER *container;
	HIT *hit;

	if (HIT_OPERATOR->used >= CONTAINER_SIZE - 1) {
		container = alloc_hit_container();
		if (!container)
			return NULL ;
		HIT_OPERATOR->next = container;
		HIT_OPERATOR = container;
	}
	
	hit = &(HIT_OPERATOR->entries[HIT_OPERATOR->used]);
	HIT_OPERATOR->used++;

	new (hit) HIT();
	return(hit);
}

HIT_CONTAINER* Hits::alloc_hit_container()
{
	HIT_CONTAINER *container;

	if ((container =  (HIT_CONTAINER *) calloc (1, sizeof(HIT_CONTAINER))) == NULL) {
		fprintf(stderr, "WARNING: not enough memory for mapping entry container memory (alloc_hit_container)\n");
		return NULL ;
	}

	container->used = 0;
	container->next = NULL;

	return(container);
}

int Hits::alloc_hit_lists_operator()
{
	if ((HIT_LISTS_OPERATOR = (HIT **) calloc (_config.MAX_READ_LENGTH, sizeof(HIT*))) == NULL) {
		fprintf(stderr, "ERROR : not enough memory for hitlist (alloc_hit_lists_operator)\n");
		exit(1);
	}

	return(0);
}

int Hits::alloc_hits_by_score()
{
	double max_score = _config.NUM_GAPS * _config.GAP_SCORE + (_config.NUM_EDIT_OPS - _config.NUM_GAPS) * _config.MM_SCORE;
	NUM_SCORE_INTERVALS = max_score / SCORE_INTERVAL;
	if (NUM_SCORE_INTERVALS * SCORE_INTERVAL != max_score) ++NUM_SCORE_INTERVALS;
	NUM_SCORE_INTERVALS++;
	
	if ((HITS_BY_SCORE = (HITS_BY_SCORE_STRUCT *) calloc (NUM_SCORE_INTERVALS, sizeof(HITS_BY_SCORE_STRUCT))) == NULL) {
		fprintf(stderr, "ERROR : not enough memory for hitlist by score (alloc_hits_by_score)\n");
		exit(1);
	}

	unsigned int i;
	for (i=0; i!=NUM_SCORE_INTERVALS; ++i) {
		HITS_BY_SCORE[i].hitpointer = NULL;
		HITS_BY_SCORE[i].num = 0;
	}

	return(0);
}

int Hits::alloc_readstart_bins()
{
	if ((READSTART_BINS = (HIT **) calloc ((genome->LONGEST_CHROMOSOME / ((_config.NUM_GAPS==0)? 1: _config.NUM_GAPS)), sizeof(HIT*) )) == NULL) {
		fprintf(stderr, "ERROR : not enough memory for readstart bin structure (alloc_readstart_bins)\n");
		exit(1);
	}

	return(0);
} 

int Hits::dealloc_mapping_entries()
{
	NUM_MAPPING_ENTRIES = 0;
	
	MAPPING_ENTRY_CONTAINER *container, *next;
	container = MAPPING_ENTRY_OPERATOR_FIRST->next;

	while (container != NULL) {
		next = container->next;
		NUM_MAPPING_ENTRIES += container->used;
		free(container);
		container = next;
	}

	NUM_MAPPING_ENTRIES += MAPPING_ENTRY_OPERATOR_FIRST->used;
	
	MAPPING_ENTRY_OPERATOR = MAPPING_ENTRY_OPERATOR_FIRST;
	MAPPING_ENTRY_OPERATOR->used = 0;
	MAPPING_ENTRY_OPERATOR->next = NULL;

	return(0);
}

int Hits::dealloc_hits()
{
	HIT_CONTAINER *container, *next;
	container = HIT_OPERATOR_FIRST->next;

	while (container != NULL) {
		next = container->next;
		free(container);
		container = next;
	}


	HIT_OPERATOR = HIT_OPERATOR_FIRST;
	HIT_OPERATOR->used = 0;
	HIT_OPERATOR->next = NULL;

	return(0);
}

int Hits::dealloc_chromosome_entries()
{
	unsigned int i;
	
	for (i=0; i < _config.CHROM_CONTAINER_SIZE; i++) {
		free(CHROMOSOME_ENTRY_OPERATOR->entries[i].mapping_entries);
		free(CHROMOSOME_ENTRY_OPERATOR->entries+i);
	}
	//free(CHROMOSOME_ENTRY_OPERATOR);

	return(0);
}

int Hits::dealloc_hit_lists_operator()
{
	unsigned int i;

	//for (i = _config.INDEX_DEPTH; i <= LONGEST_HIT; i++) {
	for (i = 0; i <= LONGEST_HIT; i++) 
	{
		*(HIT_LISTS_OPERATOR+i) = NULL;
	}
	
	return(0);
}

int Hits::dealloc_hits_by_score()
{
	unsigned int i;

	for (i = 0; i != NUM_SCORE_INTERVALS; i++) {
		HITS_BY_SCORE[i].hitpointer = NULL;
		HITS_BY_SCORE[i].num = 0;
	}

	return(0);
}



// read


int Hits::map_fast(Read & read, int & firstslot, int & firstpos)
{
#ifndef BinaryStream_MAP
	STORAGE_ENTRY *index_mmap ;
#else
	CBinaryStream<STORAGE_ENTRY>* index_mmap;
#endif
	INDEX_ENTRY index_entry;

	unsigned int pos, i, j, p, chars, block, chrom_overlap = 0, hits_reported = 0;
	int slot, run, nr_mms, readpos, chrpos, readstart, chrstart;
	firstslot = -1 ; firstpos=0 ;

	unsigned char position;
	char nr_runs, rev, nr_seeds = read.length() / _config.INDEX_DEPTH, mm, perfect = 0, cancel=0;

	if (/*_config.*/_config.NUM_GAPS != 0)
		nr_runs = (nr_seeds > 1)? 2: 1;
	else {
		if (nr_seeds > /*_config.*/_config.NUM_MISMATCHES) nr_runs = /*_config.*/_config.NUM_MISMATCHES + 1;
			else nr_runs = nr_seeds;
	}

	bool seed_already_inspected_fwd[nr_runs+1];
	bool seed_already_inspected_rev[nr_runs+1];

	int max_mms = nr_runs - 1;
	int mmpos[max_mms];
	//fprintf(stdout, "nr_runs=%i\n", nr_runs) ;

	for (run=1; run<=nr_runs; ++run) {

		if (nr_runs == 1) nr_runs = 0;	// a bit fishy, but nr_runs and run only have to be different, thats why this assignment is due
		if (run == nr_runs) readstart = ((int)read.length()) - /*_config.*/_config.INDEX_DEPTH;
			else		 	readstart = (run-1) * /*_config.*/_config.INDEX_DEPTH;

		slot = get_slot(read, readstart) ;

		//fprintf(stdout, "map_fast: slot=%i HAS_SLOT=%i read_start=%i seq=%s\n", slot, HAS_SLOT, readstart, get_seq(slot)) ;

		if (run==1)
		{
			if (slot<0)
			{
				for (int ii=1; ii<((int)read.length())-(int)_config.INDEX_DEPTH && firstslot<0; ii++)
				{
					firstpos = ii ;
					firstslot = get_slot(read, firstpos) ;
				}
			}
			else
			{
				firstslot = slot ;
				//firstpos = 0 ;
			}
		}

		// does not work when one of the first symbols are non-ACGT.
		//if (run == 1)
		//	firstslot = slot;


		if (slot >= 0)
		{	// tests if slot has an unallowed char!

			for (rev=0; rev <= /*_config.*/_config.MAP_REVERSE; ++rev) {
				if (!rev) {
					index_entry = genome->INDEX[slot];
					index_mmap = genome->INDEX_FWD_MMAP;
				}
				else  {
					index_entry = genome->INDEX_REV[slot];
					index_mmap = genome->INDEX_REV_MMAP;
				}

				// for each mapping position
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

					int index_entry_num=index_entry.num ;
					if (index_entry.num > _config.SEED_HIT_CANCEL_THRESHOLD) { // && !REPORT_REPETITIVE_SEEDS)
						index_entry_num=0 ;
						if (rev) seed_already_inspected_rev[run] = false;
						else seed_already_inspected_fwd[run] = false;
					}
					else {
						if (rev) seed_already_inspected_rev[run] = true;
						else seed_already_inspected_fwd[run] = true;
					}

					time_t start_time = clock() ;
#ifndef BinaryStream_MAP
					genome->index_pre_buffer(index_mmap, se_buffer, index_entry.offset-index_entry_num, index_entry_num);
#else
					index_mmap->pre_buffer(se_buffer, index_entry.offset-index_entry_num, index_entry_num);
#endif
					time2a_seek += clock()-start_time ;

					for (i=0; (int)i<index_entry_num; i++)
					{

						block = 0; position = 0;
						unsigned char* p_block = (unsigned char*) &block;
						STORAGE_ENTRY se;
						unsigned char* p_id;

						if (!rev) {
							//memcpy(&block, &((index_mmap+(index_entry.offset-(i+1)))->id[0]), 3 * sizeof(char));
							//memcpy(&position, &((index_mmap+(index_entry.offset-(i+1)))->id[3]), sizeof(unsigned char));
							se = se_buffer[index_entry.num-(i+1)];
						}
						else {
							//memcpy(&block, &((index_mmap+(index_entry.offset-(index_entry.num-i)))->id[0]), 3 * sizeof(char));
							//memcpy(&position, &((index_mmap+(index_entry.offset-(index_entry.num-i)))->id[3]), sizeof(unsigned char));
							se = se_buffer[i];//index_entry.offset-(index_entry.num-i)];
						}

						p_id=se.id;
						p_block[0]=p_id[0];
						p_block[1]=p_id[1];
						p_block[2]=p_id[2];
						position = p_id[3];
						pos = (unsigned int) position + genome->BLOCK_TABLE[block].pos;	// 0-initialized
						Chromosome &chr = genome->chromosome(genome->BLOCK_TABLE[block].chr);

						//if (_config.REPORT_REPETITIVE_SEEDS)
						//	report_repetitive_seed(chr, pos, index_entry.num)  ;

						//if (index_entry.num>_config.SEED_HIT_CANCEL_THRESHOLD)
						//	continue ;

						if (!rev) {
							chrstart = pos - (run!=nr_runs) * (run-1) * /*_config.*/_config.INDEX_DEPTH - (run==nr_runs) * (((int)read.length()) - /*_config.*/_config.INDEX_DEPTH);
						}	// 0-initialized
						else {
							chrstart = pos + (run!=nr_runs) *   run   * /*_config.*/_config.INDEX_DEPTH + (run==nr_runs) * ((int)read.length()) - 1;
						}

						// check if read can map on position in genome:
						if ( (!rev && chrstart < 0) ||
								(rev && chrstart < ((int)read.length()) - 1) ) {
							if (/*_config.*/_config.STATISTICS) chrom_overlap++;
						}
						else if ( (!rev && chrstart + ((int)read.length()) > (int)chr.length()) ||
								  ( rev && chrstart > (int)chr.length() - 1) ) {
							if (/*_config.*/_config.STATISTICS) chrom_overlap++;
						}
						else {

							nr_mms = 0;
							chars = 0;
							cancel = 0;

							readpos = 0;
							chrpos = chrstart;

							for (j=1; (int)j<run; ++j) {

								mm = 0;

								for (p=0; p!=/*_config.*/_config.INDEX_DEPTH; ++p) {

									char * read_data=read.data() ;
									read_data[readpos+p] = mytoupper(read_data[readpos+p]);
									if ( ( rev && get_compl_base(chr[chrpos-p]) != read_data[readpos+p]) ||
											(!rev &&                chr[chrpos+p]  != read_data[readpos+p]) ||
											!(unique_base(read_data[readpos+p])) ) {

										if (nr_mms < max_mms) {
											mmpos[nr_mms] = readpos + p + 1;
											++nr_mms;
										}
										else {
											cancel = 1;
											break;
										}
										++mm;

									}

									++chars;

								}

								if (!mm && ((rev && seed_already_inspected_rev[j]) || (!rev && seed_already_inspected_fwd[j]))) {
									cancel = 1;
									break;
								}

								if (cancel) break;

								chrpos += /*_config.*/_config.INDEX_DEPTH * (rev? -1: 1);
								readpos += /*_config.*/_config.INDEX_DEPTH;

							}


							chrpos  = chrpos  + (run!=nr_runs) * /*_config.*/_config.INDEX_DEPTH * (rev? -1: 1);
							readpos = readpos + (run!=nr_runs) * /*_config.*/_config.INDEX_DEPTH;
							while (!cancel && (int)chars != ((int)read.length()) - /*_config.*/(int)_config.INDEX_DEPTH) {

								read.data()[readpos] = mytoupper(read.data()[readpos]);
								if ( ( rev && get_compl_base(chr[chrpos]) != read.data()[readpos]) ||
										(!rev && 				 chr[chrpos]  != read.data()[readpos]) ||
										!(unique_base(read.data()[readpos])) )
								{
									if (nr_mms == max_mms) {
										cancel = 1;
										break;
									}
									mmpos[nr_mms] = readpos + 1;
									++nr_mms;
								}

								readpos++;
								chrpos += rev? -1: 1;
								chars++;

							}


							if ( !cancel && nr_mms <= max_mms ) {
								// create hit
								HIT* hit = new HIT();
								if (!hit)
								{
									delete[] se_buffer;
									return -1 ;
								}

								hit->chromosome = &chr;
								hit->readpos = 1;

								if (!rev) {
									hit->orientation = '+';
									hit->start = chrstart + 1;					// 1-initialized
									hit->end = chrstart + ((int)read.length());
								}
								else {
									hit->orientation = '-';
									hit->start = chrstart - ((int)read.length()) + 2;	// 1-initialized
									hit->end = chrstart + 1;
								}

								mm = 0;
								// create possible mismatches
								for (j=0; (int)j!=nr_mms; ++j) {
									hit->edit_op[j].mm = 1;
									if (hit->orientation == '+') hit->edit_op[j].pos = mmpos[j];
									else					 hit->edit_op[j].pos = ((int)read.length()) - mmpos[j] + 1;
									assert(hit->edit_op[j].pos >= -((int)read.length()) && hit->edit_op[j].pos<=(int)read.length()) ;
									hit->mismatches++;
									mm = 1;
								}

								//if (!_config.ALL_HIT_STRATEGY && nr_mms < max_mms)
								//	max_mms = nr_mms;

								update_num_edit_ops(nr_mms, _config.ALL_HIT_STRATEGY, max_mms) ;

								// perfect matching read
								if (/*_config.*/_config.STATISTICS) {
									if (!mm) {
										if (rev) /*_stats.*/_stats.PERFECT_HITS_REV++;
										else /*_stats.*/_stats.PERFECT_HITS++;
										if (!perfect) /*_stats.*/_stats.PERFECT_READS++;

										perfect = 1;
									}
									else /*_stats.*/_stats.NOT_ALIGNED[1]++;
								}

								int ret = insert_into_scorelist(hit, 0);
								assert(ret>=0) ;

								if (/*_config.*/_config.STATISTICS)	_stats.HITS_LEN[read.length()]++;

								hits_reported++;

							} // end of create hit

						} // end of no hit-overlap with chrom border

					} // end of for each mapping pos

					delete[] se_buffer;
				} // end of index entry num

			} // end of forward/reverse	rev

		} // end of slot != -1

	} // end of runs = different slots

	if (!/*_config.*/_config.ALL_HIT_STRATEGY && !hits_reported)
	{	//if best hit strategy, but no mappings found -> prepare for complete mapping!
		/*_config.*/_config.ALL_HIT_STRATEGY = -1;
	}
	else
	{
		if (/*_config.*/_config.STATISTICS)
		{
			/*_stats.*/_stats.NUM_HITS += hits_reported;
			/*_stats.*/_stats.HITS_PER_READ += hits_reported;
			/*_stats.*/_stats.ENDSTART_MAPPED[0] += chrom_overlap;
		}
	}
	//fprintf(stdout, "map_fast: firstslot=%i, firstpos=%i\n", firstslot, firstpos) ;

	return 0 ;
}


int Hits::map_short_read(Read& read, unsigned int num, int first_slot, int first_pos)
{
	char reverse;
	unsigned int readpos = first_pos ;
	unsigned int spacer = first_pos;
	unsigned int slot;

	//fprintf(stdout, "first_slot=%i\n", first_slot) ;

	if (first_slot < 0) {	// first slot has an unallowed char!
		return -1 ;
		/*assert(0) ; //obsoleted code
		spacer = -first_slot+1 +1 ;
		readpos = -first_slot+1 + 1 ;
		HAS_SLOT = 0;*/
	}
	else if (/*_config.*/_config.ALL_HIT_STRATEGY && !_config.OUTPUT_FILTER==OUTPUT_FILTER_TOP) {	// first slot hasn't been computed yet
		spacer = 0;
		slot = 0;
		HAS_SLOT = 0;
	}
	else {	// first slot has already been computed in map_fast and doesn't contain unallowed chars
		slot = first_slot;

		reverse = 0;
		if (genome->INDEX[slot].num != 0)
			reverse = 1;
		if (/*_config.*/_config.MAP_REVERSE && genome->INDEX_REV[slot].num != 0)
			reverse = (reverse + reverse) + 2;

		if (reverse > 0)
		{
			clock_t seed2genome_start=clock() ;
			//fprintf(stdout, "1: slot=%i\n", slot) ;
			int ret = seed2genome(num, slot, readpos + 1, reverse) ;
			time2a_seed2genome += clock() - seed2genome_start ;

			if (ret<0)
			{
				fprintf(stderr, "seed2genome<0 (1)\n")  ;
				return ret ; // add one to readpos for 1-initialization
			}
		}

		readpos++;
		spacer = /*_config.*/_config.INDEX_DEPTH;
		SLOT = first_slot;
		HAS_SLOT = 1;
	}

	while (spacer < read.length())
	{
		char * read_data = read.data() ;
		
		read_data[spacer] = mytoupper(read_data[spacer]);
		if (spacer < readpos + /*_config.*/_config.INDEX_DEPTH - 1)
		{
			if (read_data[spacer]=='A' || read_data[spacer]=='T' || read_data[spacer]=='C' || read_data[spacer]=='G')
			{
				spacer++;
			}
			else
			{
				spacer++;
				readpos = spacer;
				HAS_SLOT = 0;
			}
		}
		else {
			if (read_data[spacer]=='A' || read_data[spacer]=='T' || read_data[spacer]=='C' || read_data[spacer]=='G')
			{
				slot = get_slot(read, readpos);
				//fprintf(stdout, "map_short_read: slot=%i HAS_SLOT=%i readpos=%i seq=%s\n", slot, HAS_SLOT, readpos, get_seq(slot)) ;

				if (slot<0)
				{
					spacer++;
					readpos++;
					continue ;
				}

				// reverse: 0: slot doesnt match in either index or index_rev, 1: only index, 2: only index_rev, 4: both
				reverse = 0;
				//if (INDEX[slot].last_entry != NULL) {
				if (genome->INDEX[slot].num != 0) {
					reverse = 1;
				}
				if (/*_config.*/_config.MAP_REVERSE && genome->INDEX_REV[slot].num != 0) {
					reverse = (reverse + reverse) + 2;
				}

				if (reverse > 0) {
					clock_t seed2genome_start=clock() ;
					int ret = seed2genome(num, slot, readpos + 1, reverse) ;
					time2a_seed2genome += clock() - seed2genome_start ;
					if (ret<0)
					{
						fprintf(stderr, "seed2genome<0 (2)\n") ;
						return ret ; // add one to readpos for 1-initialization
					}
				}

				spacer++;
				readpos++;
				HAS_SLOT = 1;
			}
			else {
				spacer++;
				readpos=spacer;
				HAS_SLOT = 0;
			}
		}
	}

	HAS_SLOT = 0;

	return 1;
}


/**
 *
 *  \param pos Position in the current read for which to get a slot
 *  \return A slot
 */

int Hits::get_slot(Read & read, int pos)
{
    unsigned int slot = 0;
	unsigned int i;
	int c = 0;
	char * read_data = read.data() ;
	
	if (HAS_SLOT == 0)
	{

		for (i = 0; i < /*_config.*/_config.INDEX_DEPTH; i++)
		{
			read_data[pos + i] = mytoupper(read_data[pos + i]);

			switch (read_data[pos + i])
			{
			case 'A':
				c = 0;
				break;
			case 'C':
				c = 1;
				break;
			case 'G':
				c = 2;
				break;
			case 'T':
				c = 3;
				break;
			default:
				return -pos - i - 1 ; // subtract -1 to make sure the output is indeed negative (pos=i=0..?)
			}

			slot = slot + Util::POWER[i] * c;
		}

	} else {

		slot = SLOT;
		slot >>= 2;

		read_data[pos + /*_config.*/_config.INDEX_DEPTH - 1] = mytoupper(read_data[pos + /*_config.*/_config.INDEX_DEPTH - 1]);

		switch (read_data[pos + /*_config.*/_config.INDEX_DEPTH - 1]) {
		case 'A':
			slot = slot | genome->BINARY_CODE[0];
			break;
		case 'C':
			slot = slot | genome->BINARY_CODE[1];
			break;
		case 'G':
			slot = slot | genome->BINARY_CODE[2];
			break;
		case 'T':
			slot = slot | genome->BINARY_CODE[3];
			break;
		default:
			return -1;
		}
	}

	SLOT = slot;

	return (slot);
}

int Hits::align_hit_simple(HIT* hit, int start, int end, int readpos, Chromosome const &chromosome, int orientation, unsigned char mismatches)
{
	int Read_length = _read.length();
	char const * const READ = _read.data();
	//int Chr_length = CHR_LENGTH[chromosome];
	//int Num_edit_ops = _config.NUM_EDIT_OPS;
	int Num_mismatches = _config.NUM_MISMATCHES;
	// MAX_EDIT_OPS and _config.ALL_HIT_STRATEGY are used only once -> no global var

	EDIT_OPS edit_op[Config::MAX_EDIT_OPS];
	memcpy(edit_op, hit->edit_op, mismatches * sizeof(EDIT_OPS));

	int hitlength = end - start + 1;
	int afterhit_len = Read_length - readpos - hitlength + 1;

	// checking if read fits on the start or end of chromosome:
	//@TODO small overlaps should be mapped
	/*if (orientation == '+') {
		if (start < readpos) {
			//printf("Read cannot be mapped on plus strand at start of chromosome %d!\n", hit->chromosome+1);
			//if (_config.STATISTICS) _stats.ENDSTART_MAPPED[0]++;
			return 0;
		}
		if (end + afterhit_len > Chr_length) {
			//printf("Read cannot be mapped on plus strand at end of chromosome %d!\n", hit->chromosome+1);
			//if (_config.STATISTICS) _stats.ENDSTART_MAPPED[0]++;
			return 0;
		}
	}

	if (orientation == '-') {
		if (afterhit_len >= start) {
			//printf("Read cannot be mapped on minus strand at start of chromosome %d!\n", hit->chromosome+1);
			//if (_config.STATISTICS) _stats.ENDSTART_MAPPED[0]++;
			return 0;
		}
		if (end + readpos - 1 > Chr_length) {
			//printf("Read cannot be mapped on minus strand at end of chromosome %d!\n", hit->chromosome+1);
			//if (_config.STATISTICS) _stats.ENDSTART_MAPPED[0]++;
			return 0;
		}
	}*/

	//if (_config.STATISTICS) _stats.NUM_ALIGNMENTS++;

	int j;
	//int val;

	int readstart;
	if (orientation == '+') {
		readstart = start - readpos;	// 0-initialized
	}
	else {
		readstart = start - afterhit_len - 1;	//changed		0-initialized
	}


	// from read[0] to read[hit->readpos]
	for (j=0; j!=readpos-1; ++j) {

		if (	(orientation == '+')
				&& (	(chromosome[start - readpos + j] != READ[j])
				 	 || !(unique_base(READ[j]))		// [XX] should also be a mismatch!
					// if read[j]=X and chr_seq not, then first or-condition is automatically true -> genome sequence doesn't have to be checked
				   )
		   )
		{
			// create mismatch:
			if (mismatches < _config.NUM_EDIT_OPS) {
				(edit_op[mismatches]).pos = j+1;
				(edit_op[mismatches]).mm = 1;
			}

			mismatches++;
			assert(mismatches<Config::MAX_EDIT_OPS) ;
		}

		if (	(orientation == '-')
			 	&& (    (get_compl_base(chromosome[end + readpos - 2 - j]) != READ[j])
					 || !(unique_base(READ[j]))
				   )
		   )
		{
			// create mismatch:
			if (mismatches < _config.NUM_EDIT_OPS) {
				edit_op[mismatches].pos = Read_length - j;
				edit_op[mismatches].mm = 1;
			}

			mismatches++;
			assert(mismatches<Config::MAX_EDIT_OPS) ;
		}

		if (mismatches > _config.NUM_EDIT_OPS) {
			/*if (hit->orientation == '+') val = readstart;
				else val = readstart * 2;
			if (READSTART[val]) REDUNDANT++;
				else READSTART[val] = 1;*/
			return 0;
		}
	}

	// from read[hit->readpos + hitlength] to read[_read.lenght() - 1]
	int i = 0;
	j = readpos + hitlength - 1;


	while ((mismatches <= _config.NUM_EDIT_OPS) && (j < Read_length)) {

		if (	(orientation == '+')
			 	&& (    (hit->chromosome->operator [](end + i) != READ[j])
					 || !(unique_base(READ[j]))		// [XX] should also be a mismatch!
					// if read[j]=X and chr_seq not, then first or-condition is automatically true -> genome sequence doesn't have to be checked
				   )
		   )
		{
			// create mismatch:
			if (mismatches < _config.NUM_EDIT_OPS) {
				(edit_op[mismatches]).pos = j+1;
				(edit_op[mismatches]).mm = 1;
			}

			mismatches++;
			assert(mismatches<Config::MAX_EDIT_OPS) ;
		}


		if (	(orientation == '-')
			 	&& (    (get_compl_base(chromosome[start - 2 - i]) != READ[j])
					 || !(unique_base(READ[j]))
				   )
		   )
		{
			// create mismatch:
			if (mismatches < _config.NUM_EDIT_OPS) {
				(edit_op[mismatches]).pos = Read_length - j;
				(edit_op[mismatches]).mm = 1;
			}

			mismatches++;
			assert(mismatches<Config::MAX_EDIT_OPS) ;
		}

		if (mismatches > _config.NUM_EDIT_OPS) {
			/*if (hit->orientation == '+') val = readstart;
				else val = readstart * 2;
			if (READSTART[val]) REDUNDANT++;
				else READSTART[val] = 1;*/
			return 0;
		}

		++i;
		++j;
	}

	if (mismatches <= Num_mismatches) {	// there can't be gaps

		assert(mismatches<Config::MAX_EDIT_OPS) ;

		// write in hit-structure:
		hit->mismatches = mismatches;

//		memcpy(hit->edit_op, edit_op, mismatches * sizeof(EDIT_OPS));

		// this version leads to seg faults later on ... quite unclear why
		for (int ii=0; ii<mismatches; ii++)
		{
			assert(edit_op[ii].pos>=-((int)_read.length()) && edit_op[ii].pos<=((int)_read.length())) ;
			hit->edit_op[ii]=edit_op[ii] ;
			assert(hit->edit_op[ii].pos>=-((int)_read.length()) && hit->edit_op[ii].pos<=((int)_read.length())) ;
		}
		
		update_num_edit_ops(mismatches, _config.ALL_HIT_STRATEGY, _config.NUM_EDIT_OPS) ;

		return 1;
	}

	return 0;
}


// returns if aligned hit fulfills MM and gap-criterias, thus is printed out (alignments are called in this method)
int Hits::prepare_kbound_alignment(HIT* hit, int start, int end, int readpos, Chromosome const &chromosome, char orientation, char mismatches)
{
	// global vars -> local vars
	int Read_length = _read.length();
	int Chr_length = chromosome.length();
	char All_hit_strategy = _config.ALL_HIT_STRATEGY;
	//int Num_edit_ops = _config.NUM_EDIT_OPS;
	int Num_gaps = _config.NUM_GAPS;
	// variable _config.OVERHANG_ALIGNMENT is used only once -> no local variable

	int hitlength = end - start + 1;
	int afterhit_len = Read_length - readpos - hitlength + 1;

	// checking if read fits on the start or end of chromosome:
	//@TODO small overlaps should be mapped
	/*if (orientation == '+') {
		if (start < readpos) {
			//if (_config.STATISTICS) _stats.ENDSTART_MAPPED[0]++;
			return 0;
		}
		if (end + afterhit_len > Chr_length) {
			//if (_config.STATISTICS) _stats.ENDSTART_MAPPED[0]++;
			return 0;
		}
	}

	if (orientation == '-') {
		if (afterhit_len >= start) {
			//if (_config.STATISTICS) _stats.ENDSTART_MAPPED[0]++;
			return 0;
		}
		if (end + readpos - 1 > Chr_length) {
			//if (_config.STATISTICS) _stats.ENDSTART_MAPPED[0]++;
			return 0;
		}
	}*/

	// just perform global alignment if gap heuristic/speedup was disabled:
	if (!_config.OVERHANG_ALIGNMENT) 
	{
		//if (_config.STATISTICS) _stats.NUM_WHOLE_ALIGNMENTS++;
		/*if (kbound_global_alignment(hit, readpos, start, end, chromosome, orientation)) {
			mismatches = hit->mismatches;
			if (!All_hit_strategy && mismatches < _config.NUM_EDIT_OPS) _config.NUM_EDIT_OPS = mismatches;
			return insert_into_scorelist(hit, 1);
		}
		else return 0;*/
		mismatches = kbound_global_alignment(_read, hit, readpos, start, end, chromosome, orientation);
		if (mismatches < 0) return 0;
		assert(mismatches<Config::MAX_EDIT_OPS) ;

		//if (!All_hit_strategy && mismatches < _config.NUM_EDIT_OPS)
		//	_config.NUM_EDIT_OPS = mismatches;
		update_num_edit_ops(mismatches, All_hit_strategy, _config.NUM_EDIT_OPS) ;

		return 1;
	}


	// perform whole alignment pipeline:

	char k1_aligned;
	int offset;

	int readstart;	// start pos of read on genome
	if (orientation == '+') {
		readstart = start - readpos;		//	0-initialized
	}
	else {
		readstart = start - afterhit_len - 1;	//changed		0-initialized
	}
	int readend = readstart + Read_length;	// end pos of read on genome	1-initialized


	if (readpos != 1) {

			if (orientation == '+') {
				if (readstart - Num_gaps < 0) offset = readstart;
					else offset = Num_gaps;
			}
			else {
				if (readend + Num_gaps > Chr_length) offset = Chr_length - readend;
					else offset = Num_gaps;
			}

			// perform alignment
			k1_aligned = kbound_overhang_alignment(_read, hit, offset, 0, start, end, readpos, chromosome, orientation, mismatches);
			mismatches = hit->mismatches;
			assert(mismatches<Config::MAX_EDIT_OPS) ;

			// there are gaps on best path in alignment -> perform whole global alignment
			if (k1_aligned == 0) {

				//if (_config.STATISTICS) _stats.NUM_WHOLE_ALIGNMENTS++;
				/*if (kbound_global_alignment(hit, readpos, start, end, chromosome, orientation)) {
					mismatches = hit->mismatches;
					if (!All_hit_strategy && mismatches < _config.NUM_EDIT_OPS) _config.NUM_EDIT_OPS = mismatches;
					return insert_into_scorelist(hit, 1);
				}
				else return 0;*/
				mismatches = kbound_global_alignment(_read, hit, readpos, start, end, chromosome, orientation);
				if (mismatches < 0) return 0;
				assert(mismatches<Config::MAX_EDIT_OPS) ;
				//if (!All_hit_strategy && mismatches < _config.NUM_EDIT_OPS)
				//_config.NUM_EDIT_OPS = mismatches;
				update_num_edit_ops(mismatches, All_hit_strategy, _config.NUM_EDIT_OPS) ;
				return 1;
			}

			// too many mismatches in aligned read already:
			if (k1_aligned == -1) {
				return 0;
			}

	}

	if (readpos + hitlength != Read_length + 1) {

			if (orientation == '+') {
				if (readend + Num_gaps > Chr_length) offset = Chr_length - readend;
					else offset = Num_gaps;
			}
			else {
				if (readstart - Num_gaps < 0) offset = readstart;
					else offset = Num_gaps;
			}

			// perform alignment if at least one edit op can still be afforded:
			if (mismatches < _config.NUM_EDIT_OPS || _config.NOT_MAXIMAL_HITS) {
				k1_aligned = kbound_overhang_alignment(_read, hit, offset, readpos+hitlength-1, start, end, readpos, chromosome, orientation, mismatches);
				mismatches = hit->mismatches;
				assert(mismatches<Config::MAX_EDIT_OPS) ;
			}
			else {
				return 0;
			}

			// there are gaps on best path in alignment -> perform whole global alignment
			if (k1_aligned == 0) {

				//if (_config.STATISTICS) _stats.NUM_WHOLE_ALIGNMENTS++;
				/*if (kbound_global_alignment(hit, readpos, start, end, chromosome, orientation)) {
					mismatches = hit->mismatches;
					if (!All_hit_strategy && mismatches < _config.NUM_EDIT_OPS) _config.NUM_EDIT_OPS = mismatches;
					return insert_into_scorelist(hit, 1);
				}
				else return 0;*/
				mismatches = kbound_global_alignment(_read, hit, readpos, start, end, chromosome, orientation);
				if (mismatches < 0) return 0;
				assert(mismatches<Config::MAX_EDIT_OPS) ;
				//if (!All_hit_strategy && mismatches < _config.NUM_EDIT_OPS)
				//_config.NUM_EDIT_OPS = mismatches;
				update_num_edit_ops(mismatches, All_hit_strategy, _config.NUM_EDIT_OPS) ;

				return 1;
			}

			// too many mismatches?
			if (k1_aligned == -1) {
				return 0;
			}
	}

	// gapless alignment was successful -> insert hit into HITS_BY_SCORE:
	//insert_into_scorelist(hit, 1);

	//if (!All_hit_strategy && mismatches < _config.NUM_EDIT_OPS)
	//_config.NUM_EDIT_OPS = mismatches;
	update_num_edit_ops(mismatches, All_hit_strategy, _config.NUM_EDIT_OPS) ;

	// successful alignment:
	return 1;
}



// prints out all hits which have been inserted into HITS_BY_EDITOPS
// called once for each read (?)
int Hits::analyze_hits(TopAlignments* topalignments, QPalma* qpalma) 
{
	int i, printed = 0, nr;
	HIT *hit;

	//HIT *best_hit = NULL;
	//u_int32_t num_best_hits = 0;
	//vector<HIT *> found_hits ;
	//u_int32_t num_2nd_best_hits = 0;
	//u_int32_t num_other_hits = 0;
	int reported_reads = 0 ;

	for (i = 0; i != (int)NUM_SCORE_INTERVALS; ++i) {

		if (printed && !(_config.OUTPUT_FILTER==OUTPUT_FILTER_ALL) && !(_config.OUTPUT_FILTER==OUTPUT_FILTER_TOP))
			break; // best hit strategy

		if (HITS_BY_SCORE[i].hitpointer != NULL) 
		{
			// only _config.OUTPUT_FILTER_NUM_LIMIT numbers of alignment will be chosen randomly:
			if ((_config.OUTPUT_FILTER==OUTPUT_FILTER_RANDOM) && 
				(_config.OUTPUT_FILTER_NUM_LIMIT < 0) && (HITS_BY_SCORE[i].num > -_config.OUTPUT_FILTER_NUM_LIMIT)) 
			{
				srand((unsigned) time(NULL));
				
				int j, k, n;
				int lhits[-_config.OUTPUT_FILTER_NUM_LIMIT];
				for (j = 0; j != -_config.OUTPUT_FILTER_NUM_LIMIT; ++j) {
					n = 1;
					while (n != 0) {
						n = 0;
						lhits[j] = rand() % HITS_BY_SCORE[i].num;
						for (k = 0; k != j; ++k) {
							if (lhits[j] == lhits[k])
								++n;
						}
					}
				}

				qsort(lhits, -_config.OUTPUT_FILTER_NUM_LIMIT, sizeof(int), compare_int);

				hit = HITS_BY_SCORE[i].hitpointer;

				nr = 0;
				for (j = 0; j != HITS_BY_SCORE[i].num; ++j) {

					if (lhits[nr] == j) {
						printed += topalignments->report_unspliced_hit(hit, -_config.OUTPUT_FILTER_NUM_LIMIT, qpalma);
						nr++;
					}

					if (nr == -_config.OUTPUT_FILTER_NUM_LIMIT)
						break;

					hit = hit->same_eo_succ;
				}

			} else if (_config.OUTPUT_FILTER==OUTPUT_FILTER_TOP) {

				hit = HITS_BY_SCORE[i].hitpointer;

				// iterate over all hits for this score and collect data to generate a summary
				// This somewhat counter-intuitive code works because of the way the pre-existing
				// code was set up: we see the hits in the order of their score here - better hits first.
				while (hit != NULL) 
				{
					printed += topalignments->report_unspliced_hit(hit, 0, qpalma) ;
					hit = hit->same_eo_succ;
				}

			} else { // no random selection of output alignments:

				hit = HITS_BY_SCORE[i].hitpointer;

				while (hit != NULL) {

					if (!(_config.OUTPUT_FILTER==OUTPUT_FILTER_ALL))
						nr = HITS_BY_SCORE[i].num;
					else
						nr = HITS_IN_SCORE_LIST;

					if (_config.OUTPUT_FILTER_NUM_LIMIT == 0) { // no max nr of hits per read was specified, print all
						printed += topalignments->report_unspliced_hit(hit, nr, qpalma);
					} else if (_config.OUTPUT_FILTER_NUM_LIMIT > 0 && printed < _config.OUTPUT_FILTER_NUM_LIMIT) {
						printed += topalignments->report_unspliced_hit(hit, (nr < _config.OUTPUT_FILTER_NUM_LIMIT) ? nr : _config.OUTPUT_FILTER_NUM_LIMIT, qpalma);
					} else if (_config.OUTPUT_FILTER_NUM_LIMIT == printed) { // repeatmap many alignments already printed out -> stop printing -> next read
						return 1;
					}

					hit = hit->same_eo_succ;
				}
			}

			if (_config.REPORT_MAPPED_READS)
			{
				hit = HITS_BY_SCORE[i].hitpointer;

				while (hit != NULL) {

					nr = HITS_IN_SCORE_LIST;

					printed += report_read_alignment(hit, reported_reads);
					reported_reads++ ;
					hit = hit->same_eo_succ;
				}
			}
		}
	}

	if (printed != 0)
		return 1; // read could have been mapped
	else
		return 0; // read couldn't be mapped
}



int Hits::report_read_alignment(HIT* hit, int nbest) 
{
	int hitlength = hit->end - hit->start + 1;
	unsigned int readstart;
	if (hit->orientation == '+') 
	{
		readstart = hit->start - hit->readpos + hit->start_offset; // start pos of read in genome	0-initialized
	} else 
	{
		readstart = hit->start - (((int)_read.length()) - hit->readpos - hitlength + 2) + hit->start_offset; // 0-initialized
	}

	// PERFECT HITS:
	if (hit->mismatches == 0) 
	{
		genomemaps->report_mapped_read(*hit->chromosome, readstart, readstart+1+((int)_read.length()), ((int)_read.length()) - hit->mismatches, nbest) ;
	}
	// HITS WITH MISMATCHES:
	else 
	{
		char gap_offset = 0;
		char gap_in_read = 0;
		char gap_in_chr = 0;
		
		// sort mismatches in ascending order according to their abs positions and 'gap before mm'-strategy if equal
		qsort(hit->edit_op, hit->mismatches, sizeof(EDIT_OPS), compare_editops);
		
		int j ;
		for (j = 0; j != hit->mismatches; ++j) 
		{
			if (hit->edit_op[j].pos < 0) 
			{
				gap_in_chr = 1;
			}
			gap_in_read = 0;
			if (!hit->edit_op[j].mm) 
			{
				if (gap_in_chr) 
				{
					gap_offset--;
					gap_in_chr = 0;
				}
				else 
				{
					gap_offset++;
					gap_in_read = 1;
				}
			}
		}
		
		genomemaps->report_mapped_read(*hit->chromosome, readstart, readstart+1+((int)_read.length())+gap_offset, ((int)_read.length()) - hit->mismatches, nbest) ;
	}

	return 1;
}
