#pragma once

#define DIAGONAL 'D'
#define LEFT 'L'
#define UP 'U'

extern double WORST_SCORE;
extern double WORST_MM_SCORE;

int kbound_overhang_alignment(HIT* hit, int offset, int readstart, int start, int end, unsigned short int hitreadpos, Chromosome const &chromosome, char orientation, unsigned char mismatches);
int kbound_global_alignment(HIT* hit, unsigned short int hitreadpos, unsigned int start, unsigned int end, Chromosome const &chromosome, char orientation);


//align.c
int check_mm(Chromosome const &chr, int genome_pos, int readpos, int ori);
extern int align_hit_simple(HIT* hit, int start, int end, int readpos, Chromosome const &chromosome, int orientation, unsigned char mismatches);
extern int prepare_kbound_alignment(HIT* hit, int start, int end, int readpos, Chromosome const &chromosome, char orientation, char mismatches);
