#pragma once

#define DIAGONAL 'D'
#define LEFT 'L'
#define UP 'U'

extern double WORST_SCORE;
extern double WORST_MM_SCORE;

int kbound_overhang_alignment(Read & read, HIT* hit, int offset, int readstart, int start, int end, unsigned short int hitreadpos, Chromosome const &chromosome, char orientation, unsigned char mismatches);
int kbound_global_alignment(Read & read, HIT* hit, unsigned short int hitreadpos, unsigned int start, unsigned int end, Chromosome const &chromosome, char orientation);


//align.c
int check_mm(Read & read, Chromosome const &chr, int genome_pos, int readpos, int ori);
