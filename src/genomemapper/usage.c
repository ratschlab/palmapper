// Authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include "genomemapper.h"

int usage() {
	printf("\nGenomeMapper/QPALMA v%s\n", VERSION);
	printf("written by Korbinian Schneeberger, Stephan Ossowski, Joerg Hagmann, Gunnar Raetsch, Lisa Thalheim, Fabio De Bona\n");
	printf("Max Planck Institute for Developmental Biology and Friedrich Miescher Laboratory, Tuebingen, Germany, 2008-2010\n\n");
	printf("USAGE: genomemapper [options]\n");
	printf("\n");
	printf("mandatory:\n");
	printf(" -i STRING      reference sequence (fasta)\n");
	printf(" -q STRING      query filename (fasta, fastq, SHORE flat file)\n");
	printf(" -cfg FILENAME  path to configuration file\n");
	printf("\n\n");
	printf("optional:\n");
        printf(" -a         report all hits (best hits only)\n");
        printf(" -z         report a summary of all hits\n");
        printf(" -S         report spliced hits\n\n");

	printf(" -f STRING  output format (\"shore\" or \"bed\")\n");
	printf(" -o STRING  output filename (stdout)\n");
	printf(" -H STRING  output filename for spliced hits (stdout)\n");
	printf(" -u STRING  unmapped reads filename\n\n");

	printf(" -r         disable reverse alignment\n");
	printf(" -h         always perform alignment on entire read\n");
	printf(" -d         align gaps most right (most left)\n");
	printf(" -w         allow more gaps for best hit\n");
	printf(" -e         report edit operations (alignment scores)\n");
	printf(" -l INT     seed length (index size)\n");
	printf(" -n INT     max number of best alignments (all)\n");
	printf(" -c INT     seed container size (15.000.000)\n\n");

	printf(" -seed-hit-cancel-threshold INT     number of hits of a seed that lead to its ignoration\n");
	printf(" -index-extend-threshold INT        number of hits of a seed that lead to an seed-length extension\n");
	printf(" -index-extend INT                  length of seed-length extension\n");
	printf(" -index-precache                    linearly read index file to fill caches\n\n");

	printf(" -rtrim INT                         shortens the read until a hit is found or the minimal length is reached (INT)\n\n");

	printf(" -rlim INT  limit the number of reads for alignment\n\n");

	printf(" -report STRING                  file for map reporting\n");
	printf(" -report-ro STRING               file for map reporting (read only)\n");
	printf(" -report-rep-seed                switch on reporting of repetitive seeds\n");
	printf(" -report-map-region              switch on reporting of mapped regions\n");
	printf(" -report-map-read                switch on reporting of mapped reads\n");
	printf(" -report-spliced-read            switch on reporting of spliced reads\n");
	printf(" -report-splice-sites FLOAT      report splice sites with confidence not less that threshold\n");
	printf(" -report-splice-sites-top-perc FLOAT   report splice sites with confidence in top percentile (between 0 and 1)\n");
	printf(" -qpalma-use-map                 use map for qpalma alignments\n");
	printf(" -qpalma-use-map-max-len         limit the map extension up- and downstream to the given length (100.000)\n");

	printf(" -acc STRING                           path name to acceptor splice site predictions\n");
	printf(" -don STRING                           path name to donor splice site predictions\n");
	printf(" -no-ss-pred                           indicates that no splice site predictions should be used\n");

	printf(" -filter-splice-sites-top-perc FLOAT   trigger spliced alignments, if read covers top percentile splice site (between 0 and 1)\n");
	printf(" -filter-max-mismatches INT            trigger spliced alignment, if unspliced alignment has at least this many mismatches\n");
	printf(" -filter-max-gaps INT                  trigger spliced alignment, if unspliced alignment has at least this many mismatches\n");

	printf(" -M INT     max number of mismatches (3)\n");
	printf(" -G INT     max number of gaps (1)\n");
	printf(" -E INT     max edit operations(3)\n");
	printf(" -m DOUBLE  mismatch penalty (4)\n");
	printf(" -g DOUBLE  gap penalty (5)\n");

	printf(" -v         verbose (silent)\n\n");

	printf("spliced hits definitions: (-S required)\n");
	printf(" -C INT     min combined length (25)\n");
	printf(" -L INT     min length of long hit (17)\n");
	printf(" -K INT     min length of short hit (12)\n");
	printf(" -I INT     longest intron length  (50000)\n");
	printf(" -SA INT    maximum number of spliced alignments per read (10)\n");	
	printf(" -NI INT    maximum number of introns in spliced alignments (2)\n");	
	printf(" -CT POSINT distance to tolerate between hit and existing hit cluster\n\n");

	return 0;
}

// int usage() {
// 	printf("genomemapper v%1.2f\n\n", VERSION);
// 	printf("USAGE: genomemapper [options] -g <genome fasta-file>\n");
// 	printf("                              -i <index-filename>\n");
// 	printf("                              -m <metaindex-filename>\n");
// 	printf("                              -q <query-filename>\n");
// 	printf("\n");
// 	printf("options:\n");
// 	printf(" -o <STRING>\t\toutput file\t[default: stdout]\n");
// 	printf(" -f <STRING>\t\toutput file format (\"shore\" or \"bed\")\t[default: SHORE flat file]\n");
// 	printf(" -u <STRING>\t\tleft_over file: file of unmapped reads\n\t\t\t(works only with -o option)\t[default: none]\n");
// 	printf(" -r \t\t\tdisables mapping against reverse complementary strand\n\t\t\t[default: enabled]\n");
// 	//printf(" -n <INT>\t\tmaximal allowed number of non-base symbols in reads\t[default: 4]\n");
// 	printf(" -A \t\t\tall hits strategy (report all hits)\t[default: best hits strategy]\n");
// 	printf(" -a <INT>\t\tmaximal number of best alignments printed out per read\n\t\t\t(doesn't work in combination with -ar and -A option)\n\t\t\t[default: print all]\n");
// 	printf(" -ar <INT>\t\tmaximal number of best alignments printed out per read, randomly chosen\n\t\t\t(doesn't work in combination with -a and -A option)\n\t\t\t[default: print all]\n");
// 	printf(" -h \t\t\talways perform alignments on the complete read\n\t\t\t(slightly more sensitive)\n\t\t\t[default: only gapped alignments on whole read]\n");
// 	printf(" -l <INT>\t\tseedlength increase: only extended seeds with lengths greater\n\t\t\tthan this value will be aligned (might miss better alignments!)\n\t\t\t[default: seedlength]\n");
// 	printf(" -M <INT>\t\tmaximal number of mismatches\t  [default: 3]\n");
// 	printf(" -G <INT>\t\tmaximal number of gaps\t\t  [default: 1]\n");
// 	printf(" -E <INT>\t\tmaximal number of edit operations [default: 3]\n");
// 	printf(" -gap_score <DOUBLE>\t\t(must be > 0)\t[default: 5]\n");
// 	printf(" -mismatch_score <DOUBLE>\t(must be > 0)\t[default: 4]\n");
// 	//printf(" -match_score <DOUBLE>\t\t[default: 0]\n");
// 	printf(" -e \t\t\tprints number of edit operations instead of scores in output file\n\t\t\t[default: scores]\n");
// 	printf(" -d \t\t\tgaps will be positioned most right\t[default: most left]\n");
// 	printf(" -w \t\t\tallows more gaps than specified if best alignment necessitates so\n\t\t\t[default: drop that hit]\n");
// 	printf(" -c <INT>\t\tchromosome container size\t[default: 1.000.000]\n");
// 	//printf(" -s \t\t\tcalculate and print statistics\t[default: disabled]\n");
// 	printf(" -v \t\t\tverbose mode\t[default: silent mode]\n");
// 	//printf(" -p \t\t\tdebug mode, long output\t[default: off]\n");
// 	printf("\n");
//
// 	return 0;
// }
