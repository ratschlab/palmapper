#pragma once

#include <stdio.h>
#include <palmapper/Config.h>
#include "IntervalQuery.h"
#include <time.h>

const int num_filter_reasons=4  ;

class Statistics {
public:
	Statistics();

	unsigned int PERFECT_READS;
	unsigned int PERFECT_HITS;
	unsigned int PERFECT_HITS_REV;
	unsigned long int NUM_HITS;
	unsigned int HITS_LEN[Config::MAX_READ_LENGTH]; // TODO Refactoring needed. Should be a parameter.
	unsigned int HITS_MM[Config::MAX_EDIT_OPS + 1];
	unsigned int READS_MAPPED;
	unsigned long int NUM_ALIGNMENTS;
	unsigned long int NUM_WHOLE_ALIGNMENTS;
	unsigned int ENDSTART_MAPPED[2];
	unsigned int NOT_ALIGNED[2];
	unsigned int NUM_READS;
	unsigned int HITS_PER_READ;
	unsigned long int GAPS_ENCOUNTERED[3];
	unsigned long int TOO_MANY_MMS[2];
	unsigned long int BREAK_GLOBAL_ALIGNMENT[2];
	unsigned long int BREAK_TB_IN_GLOBAL_ALIGNMENT;
	unsigned long int CELLS_GLOBAL;
	unsigned long int CELLS_OVERHANG;
	unsigned long int W;
	unsigned int listcount, listocc;

    // Alignments
	pthread_mutex_t alignment_num_mutex ;
	int alignment_num_unspliced_best ;
	int alignment_num_unspliced_suboptimal ;
	int alignment_num_spliced_best ;
	int alignment_num_spliced_suboptimal ;
	int alignment_num_unmapped;
	clock_t alignment_last_spliced_report ;

	// read stats
	pthread_mutex_t read_stats_mutex ;
	int reads_processed ;
	int reads_with_aligment ;
	int reads_with_unique_alignment ;
	int reads_with_best_spliced_alignment ;
	int reads_unmapped_editops ;
	int reads_unmapped_dropped ;
	int reads_unmapped_other ;	
	time_t read_stats_last_report ;
    time_t read_stats_start ;	
	int reads_time_step_size ;
	clock_t read_stats_seed_time, read_stats_UA_time, read_stats_SA_time, read_stats_other_time, read_stats_JA_time ;
	std::string read_stats_note ;
	

	void print_alignment_stats(bool force=false)
		{
			print_read_stats(force) ;
			return ;
			
			if (force || (clock()-alignment_last_spliced_report)/CLOCKS_PER_SEC>=100)
			{
				alignment_last_spliced_report = clock() ;
				fprintf(stdout, "\n# %i (%i) unspliced, %i (%i) spliced alignments, %i unmapped (spliced %2.1f%%, unmapped %2.1f%%)\n", 
						alignment_num_unspliced_best, alignment_num_unspliced_suboptimal, alignment_num_spliced_best, alignment_num_spliced_suboptimal, alignment_num_unmapped,
						100.0*alignment_num_spliced_best/(alignment_num_spliced_best+alignment_num_unspliced_best), 100.0*alignment_num_unmapped/(alignment_num_spliced_best+alignment_num_unspliced_best+alignment_num_unmapped)) ;

				float sum1=(variant_lower_bound_time+ variant_identification_time+ variant_create_super_sequence_from_variants_time+variant_myalign_fast_time+ variant_reconstruct_reference_alignment_time+1)/100.0 ;
				
				if (sum1>1)
					fprintf(stdout, "variant_timing:  lower_bound=%2.1f id=%2.1f create_super_sequence=%2.1f myalign_fast=%2.1f reconstruct_reference_alignment=%2.1f\n", 
							variant_lower_bound_time/sum1, variant_identification_time/sum1, variant_create_super_sequence_from_variants_time/sum1, 
							variant_myalign_fast_time/sum1, variant_reconstruct_reference_alignment_time/sum1) ;
				
				float sum2=(variant_create_super_sequence_from_variants_init_time+variant_create_super_sequence_from_variants_insertion_init_time+
							variant_create_super_sequence_from_variants_table_init_time+variant_create_super_sequence_from_variants_table_fill1_time+
							variant_create_super_sequence_from_variants_table_fill2_time + 1)/100 ;
				if (sum2>1)
				  fprintf(stdout, "create_super_sequence_timing:  init=%2.1f insert_init=%2.1f table_init=%2.1f fill1=%2.1f fill2=%2.1f\n", 
				  variant_create_super_sequence_from_variants_init_time/sum2, variant_create_super_sequence_from_variants_insertion_init_time/sum2,
				  variant_create_super_sequence_from_variants_table_init_time/sum2, variant_create_super_sequence_from_variants_table_fill1_time/sum2,
				  variant_create_super_sequence_from_variants_table_fill2_time/sum2) ;
				
			}
		}
	void reset_read_stats()
	{
		fprintf(stdout, "\nDate & time\t\t\tReads\t\t\tAligned Reads\t\t\tUnaligned Reads\t\t\t\tTime") ;
		fprintf(stdout, "\n\t\t\t\t  M/h\t       #\ttotal\tunique\tspliced\t\ttotal\t[editop\tMM\tother]\t\t[seed\tUA\tSA\tJA\tother]") ;
		read_stats_start=time(NULL)  ;
		reads_time_step_size=5 ;
		read_stats_seed_time=1;
		read_stats_UA_time=1 ;
		read_stats_SA_time=1 ;
		read_stats_JA_time=0 ;
		read_stats_other_time=0 ;
	}

	void print_read_stats(bool force=false)
		{
			if (force || difftime(time(NULL),read_stats_last_report)>=reads_time_step_size)
			{
				time_t tim=time(NULL);
				char *s=ctime(&tim);
				s[strlen(s)-1]=0;
				
				read_stats_last_report = time(NULL) ;
				if (difftime(time(NULL),read_stats_start)>10)
					reads_time_step_size=10 ;
				if (difftime(time(NULL),read_stats_start)>100)
					reads_time_step_size=30 ;
				if (difftime(time(NULL),read_stats_start)>1000)
					reads_time_step_size=120 ;
				if (difftime(time(NULL),read_stats_start)>10000)
					reads_time_step_size=600 ;
				
				fprintf(stdout, "\n%s\t%1.3f\t%8i\t%2.3f%%\t%2.1f%%\t%2.1f%%\t\t%2.3f%%\t%2.1f%%\t%2.1f%%\t%2.1f%%\t\t%2.1f%%\t%2.1f%%\t%2.1f%%\t%2.1f%%\t%2.1f%%", s, 
						reads_processed*3600.0/difftime(time(NULL),read_stats_start)/1000000, reads_processed, 
						100.0*reads_with_aligment/reads_processed, 100.0*reads_with_unique_alignment/reads_with_aligment,
						100.0*reads_with_best_spliced_alignment/reads_with_aligment,
						100.0*(reads_processed-reads_with_aligment)/reads_processed,
						100.0*reads_unmapped_editops/(reads_processed-reads_with_aligment),
						100.0*reads_unmapped_dropped/(reads_processed-reads_with_aligment),
						100.0*reads_unmapped_other/(reads_processed-reads_with_aligment),
						100.0*read_stats_seed_time/(read_stats_seed_time+read_stats_UA_time+read_stats_SA_time+read_stats_JA_time+read_stats_other_time),
						100.0*read_stats_UA_time/(read_stats_seed_time+read_stats_UA_time+read_stats_SA_time+read_stats_JA_time+read_stats_other_time),
						100.0*read_stats_SA_time/(read_stats_seed_time+read_stats_UA_time+read_stats_SA_time+read_stats_JA_time+read_stats_other_time),
						100.0*read_stats_JA_time/(read_stats_seed_time+read_stats_UA_time+read_stats_SA_time+read_stats_JA_time+read_stats_other_time),
						100.0*read_stats_other_time/(read_stats_seed_time+read_stats_UA_time+read_stats_SA_time+read_stats_JA_time+read_stats_other_time)) ;
				if (read_stats_note.size()>0)
					fprintf(stdout, "\t\t%s", read_stats_note.c_str()) ;
				read_stats_note="" ;
			}
		}


    // QPALMA Filter
	clock_t qpalma_last_filter_report;
	int qpalma_filter_stat_spliced[num_filter_reasons] ;
	int qpalma_filter_stat_unspliced[num_filter_reasons] ;

	void qpalma_filter_stat_report() const
		{
			for (int i=0; i<num_filter_reasons; i++)
			{
				fprintf(stdout, "[filter] reason %i:\t%i spliced\t%i unspliced\t%1.2f%%\n", i, qpalma_filter_stat_spliced[i], qpalma_filter_stat_unspliced[i], 100*((float)qpalma_filter_stat_spliced[i])/((float)qpalma_filter_stat_unspliced[i]+qpalma_filter_stat_spliced[i])) ;
			}
		}

	void qpalma_filter_stat(int &qpalma_filter_reason, bool spliced) 
		{
			if (qpalma_filter_reason<0)
				return ;
			
			if (spliced)
				qpalma_filter_stat_spliced[qpalma_filter_reason]++ ;
			else
				qpalma_filter_stat_unspliced[qpalma_filter_reason]++ ;
			
			if (((clock()-qpalma_last_filter_report)/CLOCKS_PER_SEC>=10))
			{
				qpalma_last_filter_report = clock() ;
				qpalma_filter_stat_report() ;
			}
			
			qpalma_filter_reason=-1 ;
		}

	// QPALMA
	int qpalma_read_count ;
	float qpalma_region1_time, qpalma_region_align_time, qpalma_align_time, intervalquery_total_time ;
	float variant_lower_bound_time, variant_identification_time, variant_create_super_sequence_from_variants_time, 
		variant_myalign_fast_time, variant_reconstruct_reference_alignment_time ;
	float variant_create_super_sequence_from_variants_init_time,	variant_create_super_sequence_from_variants_insertion_init_time, variant_create_super_sequence_from_variants_table_init_time, variant_create_super_sequence_from_variants_table_fill1_time, variant_create_super_sequence_from_variants_table_fill2_time;

	long int qpalma_total_dna_length, qpalma_total_alignments ;
	long int qpalma_total_num_threads, qpalma_total_num_thread_tasks ;

	clock_t qpalma_last_timing_report ;

	void qpalma_timing(float this_read=0) const
		{
			fprintf(stdout, "# [capture_hits] timing: %1.4f, %1.4f, %1.4f (%1.4f ss access; %lint and %1.1f threads per alignment)",
					((float) qpalma_region1_time) / CLOCKS_PER_SEC,
					((float) qpalma_region_align_time) / CLOCKS_PER_SEC, 
					((float) qpalma_align_time) / CLOCKS_PER_SEC, 
					((float) IntervalQuery::total_time) / CLOCKS_PER_SEC, 
					(qpalma_total_dna_length+1)/(qpalma_total_alignments+1),
					(float) qpalma_total_num_threads/(1e-6+qpalma_total_num_thread_tasks));
			
			if (this_read >= 0)
				fprintf(stdout, " this read: %1.4f", this_read);
			
			if (qpalma_read_count >= 0)
				fprintf(stdout, " average per read: %1.4f",
						((float) qpalma_region_align_time) / CLOCKS_PER_SEC / qpalma_read_count);
			
			fprintf(stdout, "\n");
		}

	// seed2genome
	clock_t hits_seed2genome, hits_last_report_total;
	int hits_seed2genome_cnt; 

	clock_t hits_seek, hits_part1, hits_part2, hits_part3, hits_part4, hits_part5,hits_last_report ;
	int hits_seek_cnt, hits_part1_cnt, hits_part2_cnt, hits_part3_cnt, hits_part4_cnt, hits_part5_cnt ;

	pthread_mutex_t hit_seed_mutex;

	void hits_timing_total()
		{
			if ((clock()-hits_last_report_total)/CLOCKS_PER_SEC>=1)
			{
				hits_last_report_total = clock() ;
				fprintf(stdout, "hits: total=%1.3f (avg. %1.6f)\n", 
						((double)hits_seed2genome)/CLOCKS_PER_SEC, ((double)hits_seed2genome)/(hits_seed2genome_cnt*CLOCKS_PER_SEC)) ;
			}
		}

	void hits_timing()
		{
			if ((clock()-hits_last_report)/CLOCKS_PER_SEC>=1)
			{
				hits_last_report = clock() ;
				fprintf(stdout, "hits: total=%1.3f (avg. %1.6f)\t seek=%1.3f (avg. %1.6f)\tpart1=%1.3f (avg. %1.6f)\tpart2=%1.3f (avg. %1.6f)\tpart3=%1.3f (avg. %1.6f)\tpart4=%1.4f (avg. %1.6f)\tpart5=%1.5f (avg. %1.6f)\n", 
						((double)hits_seed2genome)/CLOCKS_PER_SEC, ((double)hits_seed2genome)/(hits_seed2genome_cnt*CLOCKS_PER_SEC), 
						((double)hits_seek)/CLOCKS_PER_SEC, ((double)hits_seek)/(hits_seek_cnt*CLOCKS_PER_SEC), 
						((double)hits_part1)/CLOCKS_PER_SEC, ((double)hits_part1)/(CLOCKS_PER_SEC*hits_part1_cnt), 
						((double)hits_part2)/CLOCKS_PER_SEC, ((double)hits_part1)/(CLOCKS_PER_SEC*hits_part2_cnt), 
						((double)hits_part3)/CLOCKS_PER_SEC, ((double)hits_part1)/(CLOCKS_PER_SEC*hits_part3_cnt), 
						((double)hits_part4)/CLOCKS_PER_SEC, ((double)hits_part1)/(CLOCKS_PER_SEC*hits_part4_cnt), 
						((double)hits_part5)/CLOCKS_PER_SEC, ((double)hits_part1)/(CLOCKS_PER_SEC*hits_part5_cnt)) ;
			}
		}

};
