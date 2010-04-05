// Authors: Korbinian Schneeberger, Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include <genomemapper/Util.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <iostream>

#include <genomemapper/Hits.h>
#include <genomemapper/Read.h>
#include <genomemapper/genomemapper.h>

using namespace std;

extern clock_t time2a_seek ;
extern clock_t time2a_seed2genome  ;

int get_slot(int pos);

Read::Read() {
	READ_ID = new char[_config.MAX_READ_ID_LENGTH];
	for (int i = 0; i < 3; ++i) {
		READ_QUALITY[i] = (char*) calloc(_config.MAX_READ_LENGTH + 1, sizeof(char));
		if (READ_QUALITY[i] == NULL) {
			fprintf(stderr, "[init_defaults] Could not allocate memory\n");
			exit(1);
		}
	}
	linenr = 0;
	READ_LENGTH = 0;
}

Read::Read(const Read& read) 
{
	READ_ID = new char[_config.MAX_READ_ID_LENGTH];
	for (int i = 0; i < 3; ++i) {
		READ_QUALITY[i] = (char*) calloc(_config.MAX_READ_LENGTH + 1, sizeof(char));
		if (READ_QUALITY[i] == NULL) {
			fprintf(stderr, "[init_defaults] Could not allocate memory\n");
			exit(1);
		}
	}
	linenr = read.linenr ;
	READ_LENGTH = read.READ_LENGTH ;
	strncpy(READ_ID, read.READ_ID, _config.MAX_READ_ID_LENGTH) ;
	strncpy(READ, read.READ, _config.MAX_READ_LENGTH+1) ;
	for (int i = 0; i < 3; ++i) {
		strncpy(READ_QUALITY[i], read.READ_QUALITY[i], _config.MAX_READ_LENGTH + 1) ;
	}
	READ_FORMAT=read.READ_FORMAT ;
	READ_PE_FLAG = read.READ_PE_FLAG ;
	
}

/** Parses one line of read descriptions in either FASTA, FASTQ or FLATFILE
 *  format and sets a number of global variables describing this read.
 *
 *	The rest of the GenomeMapper code always operates on the current read.
 *	Information on this read is available through a number of global variables
 *	that are set in this routine.
 *
 *	\return A status code; 0 on success
 */
int Read::read_short_read(FILE *QUERY_FP)
{
	char line[10000];
	char *tmp;
	int linelen;

	if (fgets(line, 10000, QUERY_FP) == NULL) {
		if (READ_LENGTH == 0)
			cerr << "\n!!! WARNING: Input read file '" << /*_config.*/_config.QUERY_FILE_NAME << "' is empty!\n\n";
		return 1;
	}
	++linenr;

	if (strcspn(line, " \n\t") == 0) {
		do {
			if (fgets(line, 10000, QUERY_FP) == NULL) {
				if (READ_LENGTH == 0)
					cerr << "\n!!! WARNING: Input read file '" << /*_config.*/_config.QUERY_FILE_NAME << "' is empty!\n\n";
				return 1;
			}
			++linenr;
		} while (strcspn(line, " \n\t") == 0);
	}

	linelen = strlen(line);
	if (linelen < 3) {
		cerr << "ERROR: Unknown read input format! Do all the reads have an identifier?\n";
		exit(0);
	}

	if (line[0] == '@') {
		/////// FastQ input ///////

		// R E A D _ I D
		memset(READ, 0, /*_config.*/_config.MAX_READ_LENGTH) ;
		memset(READ_QUALITY[0], 0, /*_config.*/_config.MAX_READ_LENGTH) ;
		memset(READ_ID, 0, /*_config.*/_config.MAX_READ_ID_LENGTH) ;
		strncpy(READ_ID, line+1, strcspn(line, " \t\n")-1);

		do {
			if (fgets(line, 10000, QUERY_FP) == NULL) {
				fprintf(stderr, "ERROR: Read '%s' in line %lu is not complete in input query file '%s'! Missing read sequence and quality!\n", READ_ID, linenr, /*_config.*/_config.QUERY_FILE_NAME.c_str());
				exit(0);
			}
			++linenr;
		} while (strcspn(line, " \t\n") == 0);

		// R E A D
		strncpy(READ, line, strcspn(line, " \t\n"));
		//READ[36]=0 ;
		if (strlen(READ) > /*_config.*/_config.MAX_READ_LENGTH) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu is longer than the max read length (=%zu)! It will be omitted!\n\n", READ_ID, linenr, /*_config.*/_config.MAX_READ_LENGTH);
			return -1;
		}
		else if (strlen(READ) == 0) {
			fprintf(stderr, "ERROR: Cannot find read sequence of read '%s' in line %lu in input query file '%s'!\n", READ_ID, linenr, /*_config.*/_config.QUERY_FILE_NAME.c_str());
			exit(0);
		}
		if (strcspn(READ, "aAcCgGtTnNrRyYmMkKwWsSbBdDhHvV") != 0) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu contains non-IUPAC characters! It will be omitted!\n\n", READ_ID, linenr);
			return -1;
		}

		if (strlen(READ) < /*_config.*/_config.INDEX_DEPTH) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu is shorter than the specified seedlength! It will be omitted!\n\n", READ_ID, linenr);
			return -1;
		}

		do {
			if (fgets(line, 10000, QUERY_FP) == NULL) {
				fprintf(stderr, "ERROR: Read '%s' in line %lu is not complete in input query file '%s'! Missing quality!\n", READ_ID, linenr, /*_config.*/_config.QUERY_FILE_NAME.c_str());
				exit(0);
			}
			++linenr;
		} while (strcspn(line, " \t\n") == 0);

		// +
		if (strlen(line) < 1 || line[0] != '+') {
			fprintf(stderr, "ERROR: Read '%s' in line %lu is not in fastq format!\n", READ_ID, linenr);
			exit(0);
		}

		do {
			if (fgets(line, 10000, QUERY_FP) == NULL) {
				fprintf(stderr, "ERROR: Read '%s' in line %lu is not complete in input query file '%s'! Missing quality!\n", READ_ID, linenr, /*_config.*/_config.QUERY_FILE_NAME.c_str());
				exit(0);
			}
			++linenr;
		} while (strcspn(line, " \t\n") == 0);

		// Q U A L I T Y
		if (strlen(line) > 0)
			strncpy(READ_QUALITY[0], line, strcspn(line, " \t\n"));
		else {
			fprintf(stderr, "ERROR: Cannot find read quality of read '%s' in line %lu in input query file '%s'!\n", READ_ID, linenr, /*_config.*/_config.QUERY_FILE_NAME.c_str());
			exit(0);
		}

		/*if (strlen(READ_QUALITY[0]) != READ_LENGTH) {
			fprintf(stderr, "ERROR: Read quality 1 of read '%s' in line %lu hasn't length of read!\n", READ_ID, linenr);
			exit(0);
		}*/

		READ_LENGTH = strlen(READ);

		// O T H E R
		READ_PE_FLAG = 0;
		READ_QUALITY[1]=(char*)"" ;
		READ_QUALITY[2]=(char*)"";

		READ_FORMAT = 0;

	}
	else if (line[0] == '>') {
		/////// Fasta input ///////
		memset(READ, 0, /*_config.*/_config.MAX_READ_LENGTH) ;
		memset(READ_ID, 0, /*_config.*/_config.MAX_READ_LENGTH) ;

		strncpy(READ_ID, line+1, strcspn(line, " \t\n")-1);

		do {
			if (fgets(line, 10000, QUERY_FP) == NULL) {
				fprintf(stderr, "ERROR: Read '%s' in line %lu is not complete in input query file '%s'! Missing read sequence!\n", READ_ID, linenr, /*_config.*/_config.QUERY_FILE_NAME.c_str());
				exit(0);
			}
			++linenr;
		} while (strcspn(line, " \t\n") == 0);

		// R E A D
		strncpy(READ, line, strcspn(line, " \t\n"));
		if (strlen(READ) > /*_config.*/_config.MAX_READ_LENGTH) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu is longer than the max read length (=%zu)! It will be omitted!\n\n", READ_ID, linenr, /*_config.*/_config.MAX_READ_LENGTH);
			return -1;
		}
		else if (strlen(READ) == 0) {
			fprintf(stderr, "ERROR: Cannot find read sequence of read '%s' in line %lu in input query file '%s'!\n", READ_ID, linenr, /*_config.*/_config.QUERY_FILE_NAME.c_str());
			exit(0);
		}
		for (int i=0; i<(int)strlen(READ); i++)
			if (READ[i]=='-')
			{
				//fprintf(stderr, "replaced '-' with 'N'\n")  ;
				READ[i]='N' ;
			}

		if (strcspn(READ, "aAcCgGtTnNrRyYmMkKwWsSbBdDhHvV") != 0) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu contains non-IUPAC characters! It will be omitted!\n\n", READ_ID, linenr);
			return -1;
		}

		if (strlen(READ) < /*_config.*/_config.INDEX_DEPTH) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu is shorter than the specified seedlength! It will be omitted!\n\n", READ_ID, linenr);
			return -1;
		}

		READ_LENGTH = strlen(READ);

		READ_PE_FLAG = 0;
		strcpy(READ_QUALITY[0], READ) ;
		memset(READ_QUALITY[0], 'h', strlen(READ)) ;

		//READ_QUALITY[0] = (char*)"";
		READ_QUALITY[1] = (char*)"";
		READ_QUALITY[2] = (char*)"";

		READ_FORMAT = 1;

	}
	else {
		/////// Flatfile input ///////
		READ_ID = strtok(line, "\t");

		if ((int)strlen(READ_ID) == linelen) {
			fprintf(stderr, "ERROR: wrong read input data format, line %lu!\n", linenr);
			exit(0);
		}
		if (READ_ID == NULL) {
			fprintf(stderr, "ERROR: Read ID is empty, line %lu!\n", linenr);
			exit(0);
		}

		char *tok = strtok(NULL, "\t");
		if (tok == NULL) {
			fprintf(stderr, "ERROR: Read sequence is empty, line %lu!\n", linenr);
			exit(0);
		}

		if (strlen(tok) > /*_config.*/_config.MAX_READ_LENGTH) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu is longer than the max read length (=%zu)! It will be omitted!\n\n", READ_ID, linenr, /*_config.*/_config.MAX_READ_LENGTH);
			return -1;
		}
		//printf("%s sp: %d\n",READ,(int) strcspn(READ, "A"));
		if (strcspn(tok, "aAcCgGtTnNrRyYmMkKwWsSbBdDhHvV") != 0) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu contains non-IUPAC characters! It will be omitted!\n\n", READ_ID, linenr);
			return -1;
		}

		strcpy(READ, tok);
		READ_LENGTH = strlen(tok);
		if (READ_LENGTH < /*_config.*/_config.INDEX_DEPTH) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu is shorter than the specified seedlength! It will be omitted!\n\n", READ_ID, linenr);
			return -1;
		}

		tmp = strtok(NULL, "\t");
		if (tmp == NULL) {
			fprintf(stderr, "ERROR: Paired-end flag is empty, line %lu!\n", linenr);
			exit(0);
		}
		READ_PE_FLAG = atoi(tmp);

		READ_QUALITY[0] = strtok('\0', "\t");

		//fprintf(stderr, "hack!!!\n") ;
		//for (int i=0; i<strlen(READ_QUALITY[0]); i++)
		//	READ_QUALITY[0][i]='h' ;

		if (READ_QUALITY[0] == NULL) {
			fprintf(stderr, "ERROR: Read Quality 1 is empty, line %lu!\n", linenr);
			exit(0);
		}
		/*if (strlen(READ_QUALITY[0]) != READ_LENGTH) {
			fprintf(stderr, "ERROR: Read quality 1 hasn't length of read, line %lu!\n", linenr);
			exit(0);
		}*/

		// TODO it is very questionable to reassign the READ_QUALITY pointer since other branches of
		// this method function are memset()ting on it...
		READ_QUALITY[1] = strtok(NULL, "\t");
		if (READ_QUALITY[1] == NULL) {
			fprintf(stderr, "ERROR: Read Quality 2 is empty, line %lu!\n", linenr);
			exit(0);
		}
		/*if (strlen(READ_QUALITY[1]) != READ_LENGTH) {
			fprintf(stderr, "ERROR: Read quality 2 hasn't length of read, line %lu!\n", linenr);
			exit(0);
		}*/

		READ_QUALITY[2] = strtok(NULL, "\n");
		if (READ_QUALITY[2] == NULL) {
			fprintf(stderr, "ERROR: Read Quality 3 is empty, line %lu!\n", linenr);
			exit(0);
		}
		/*if (strlen(READ_QUALITY[2]) != READ_LENGTH) {
			fprintf(stderr, "ERROR: Read quality 3 hasn't length of read, line %lu!\n", linenr);
			exit(0);
		}*/

		READ_FORMAT = 2;
	}
	
	if (READ_LENGTH>ASSUMED_READ_LENGTH)
		ASSUMED_READ_LENGTH=READ_LENGTH ;

	return 0;
}

int Read::determine_read_length(const std::string & query_file)
{
    FILE *QUERY_FP = Util::openFile(query_file, "r") ;
	const unsigned sample_size = 10000 ;

	int sum_read_length=0 ;
	for (unsigned int i=0; i<sample_size; i++)
	{
		read_short_read(QUERY_FP) ;
		sum_read_length += READ_LENGTH ;
	}

	return sum_read_length/sample_size ;
}

