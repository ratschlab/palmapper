#pragma once

extern int compare_editops(const void *a, const void *b) ;
extern char get_compl_base(char c);
extern void print_stats();
extern int print_perfect_hits(unsigned int num);
extern int print_largest_hit();
extern void print_leftovers(const char *tag, FILE *LEFTOVER_FP);
extern void print_alignment_matrix(int chrstart, int readstart, int length, int offset_front, int offset_end, Chromosome const &chr, char ori, int K);
extern void print_alignment_stats(int num_unspliced_best, int num_unspliced_suboptimal, int num_spliced_best, int num_spliced_suboptimal) ;
extern int compare_int(const void *a, const void *b) ;
