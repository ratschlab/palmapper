// Authors: Korbinian Schneeberger, Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include <palmapper/Util.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <iostream>

#include <palmapper/Hits.h>
#include <palmapper/QueryFile.h>
#include <palmapper/Read.h>
#include <palmapper/palmapper.h>

using namespace std;

int Read::PRB_QUALITY_OFFSET=33 ;

static std::string undefinedString = std::string("<no file>");
Read::Read(QueryFile &queryFile)
	:	_queryFile(queryFile), _location(QueryFile::Location(undefinedString, -1))
{
	READ_QUALITY[0][0] = READ_QUALITY[1][0] = READ_QUALITY[2][0] = '\0';
	READ_LENGTH = 0;
	READ_ID[0] = '\0';
	orig_read = NULL ;
	_nr = -1;
	prealigned_info=NULL ;
}

Read::Read(Read const &src, unsigned cutStart, unsigned cutEnd)
	:	_queryFile(src._queryFile), _location(src._location)
{
	copyFrom(src, cutStart, cutEnd);
	READ_FORMAT = src.READ_FORMAT ;
	READ_PE_FLAG = src.READ_PE_FLAG ;
	::strcpy(READ_ID, src.READ_ID);
	orig_read = &src;
	_nr = src._nr;
	if (src.prealigned_info!=NULL)
	{
		prealigned_info = new struct t_prealigned ;
		prealigned_info->exons=src.prealigned_info->exons ;
		prealigned_info->aligned_positions=src.prealigned_info->aligned_positions ;
		prealigned_info->contig=src.prealigned_info->contig ;
		prealigned_info->cigar=src.prealigned_info->cigar ;
		prealigned_info->paired_info=src.prealigned_info->paired_info ;
		prealigned_info->strand=src.prealigned_info->strand ;
		prealigned_info->quality=src.prealigned_info->quality ;
		prealigned_info->sam_flags=src.prealigned_info->sam_flags ;
		prealigned_info->start_position=src.prealigned_info->start_position ;
	}
}

Read::~Read() {
	if (prealigned_info)
		delete prealigned_info ;
}

void Read::copyFrom(Read const &src, unsigned cutStart, unsigned cutEnd) {
	assert(cutStart <= src.READ_LENGTH);
	assert(cutEnd <= src.READ_LENGTH);
	assert(cutStart + cutEnd <= src.READ_LENGTH);
	size_t len = src.READ_LENGTH - (cutStart + cutEnd);

	::memmove(READ, src.READ + cutStart, len);
	READ[len] = '\0';

	for (int i = 0; i < 3; ++i) {
		if (src.READ_QUALITY[i][0] != 0) {
			::memmove(READ_QUALITY[i], src.READ_QUALITY[i] + cutStart, len);
			READ_QUALITY[i][len] = '\0';
		} else {
			READ_QUALITY[i][0] = '\0';
		}
	}

	READ_LENGTH = len;
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
int Read::read_short_read()
{
	char line[10000];
	char *tmp;
	int linelen;

	char * saveptr=NULL ;
	
	do {
		if (!_queryFile.next_line(line, sizeof(line)))
			return 1;
	} while (strcspn(line, " \n\t") == 0);

	linelen = strlen(line);
	if (linelen < 3) {
		cerr << "ERROR: at " << _queryFile.getLocation() << " Unknown read input format! Do all the reads have an identifier?\n";
		exit(0);
	}

	_location = _queryFile.getLocation();
	if (prealigned_info!=NULL)
		delete prealigned_info ;
	prealigned_info=NULL ;
	
	if (line[0] == '@') {
		/////// FastQ input ///////

		// R E A D _ I D
		memset(READ, 0, _config.MAX_READ_LENGTH) ;
		memset(READ_QUALITY[0], 0, _config.MAX_READ_LENGTH) ;
		memset(READ_ID, 0, _config.MAX_READ_ID_LENGTH);

		strncpy(READ_ID, line+1, strcspn(line, " \t\n")-1);
		
		{
			char READ_ID_[_config.MAX_READ_ID_LENGTH] ;
			strcpy(READ_ID_, _config.READ_ID_PREFIX.c_str()) ;
			strcpy(&(READ_ID_[strlen(READ_ID_)]), READ_ID) ;
			strcpy(READ_ID, READ_ID_) ;
		}

		do {
			if (!_queryFile.next_line(line, sizeof(line))) {
				cerr << "ERROR: Read " << *this << " is not complete in input query file '%s'! Missing read sequence and quality!\n";
				exit(0);
			}
		} while (strcspn(line, " \t\n") == 0);

		// R E A D
		strncpy(READ, line, strcspn(line, " \t\n"));
		//READ[36]=0 ;
		if (strlen(READ) > _config.MAX_READ_LENGTH) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu is longer than the max read length (=%zu)! It will be omitted!\n\n", READ_ID, _queryFile.line_nr(), _config.MAX_READ_LENGTH);
			return -1;
		}
		else if (strlen(READ) == 0) {
			cerr << "ERROR: Cannot find read sequence of read " << *this;
			exit(0);
		}
		if (strcspn(READ, "aAcCgGtTnNrRyYmMkKwWsSbBdDhHvV") != 0) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu contains non-IUPAC characters! It will be omitted!\n\n", READ_ID, _queryFile.line_nr());
			return -1;
		}

		if (strlen(READ) < _config.INDEX_DEPTH) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu is shorter than the specified seedlength! It will be omitted!\n\n", READ_ID, _queryFile.line_nr());
			return -1;
		}

		do {
			if (!_queryFile.next_line(line, sizeof(line))) {
				cerr << "ERROR: Read " << *this << " is not complete'! Missing quality!\n";
				exit(0);
			}
		} while (strcspn(line, " \t\n") == 0);

		// +
		if (strlen(line) < 1 || line[0] != '+') {
			fprintf(stderr, "ERROR: Read '%s' in line %lu is not in fastq format!\n", READ_ID, _queryFile.line_nr());
			exit(0);
		}

		do {
			if (!_queryFile.next_line(line, sizeof(line))) {
				cerr << "ERROR: Read " << *this << " is not complete'! Missing quality!\n";
				exit(0);
			}
		} while (strcspn(line, " \t\n") == 0);

		// Q U A L I T Y
		if (strlen(line) > 0)
			strncpy(READ_QUALITY[0], line, strcspn(line, " \t\n"));
		else {
			cerr <<  "ERROR: Cannot find read quality of read " << *this << "'!\n";
			exit(0);
		}

		/*if (strlen(READ_QUALITY[0]) != READ_LENGTH) {
			fprintf(stderr, "ERROR: Read quality 1 of read '%s' in line %lu hasn't length of read!\n", READ_ID, linenr);
			exit(0);
		}*/

		READ_LENGTH = strlen(READ);
		
		//Fixtrim global strategy
		if (_config.FIXTRIM_STRATEGY_LEN < READ_LENGTH) {
			READ[_config.FIXTRIM_STRATEGY_LEN] = '\0';
			READ_QUALITY[0][_config.FIXTRIM_STRATEGY_LEN] = '\0';
			READ_LENGTH = _config.FIXTRIM_STRATEGY_LEN;
		}


		//Fixtrim left and right strategy
		int trimborders= _config.FIXTRIMRIGHT_STRATEGY_LEN + _config.FIXTRIMLEFT_STRATEGY_LEN;
		if (trimborders >0){
			if (READ_LENGTH - trimborders < _config.INDEX_DEPTH){
				cerr <<	"ERROR: fixtrim alignment length too short\n";
				exit(0) ;
			}
			copyFrom(*this, _config.FIXTRIMLEFT_STRATEGY_LEN,_config.FIXTRIMRIGHT_STRATEGY_LEN);
			
		}
		
		

		// O T H E R
		READ_PE_FLAG = 0;
		READ_FORMAT = 0;
	}
	else if (line[0] == '>') {
		/////// Fasta input ///////
		memset(READ, 0, _config.MAX_READ_LENGTH) ;
		memset(READ_ID, 0, _config.MAX_READ_LENGTH) ;

		strncpy(READ_ID, line+1, strcspn(line, " \t\n")-1);

		do {
			if (!_queryFile.next_line(line, sizeof(line))) {
				cerr << "ERROR: Read " << *this << " is not complete'! Missing read sequence!\n";
				exit(0);
			}
		} while (strcspn(line, " \t\n") == 0);

		// R E A D
		strncpy(READ, line, strcspn(line, " \t\n"));
		if (strlen(READ) > _config.MAX_READ_LENGTH) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu is longer than the max read length (=%zu)! It will be omitted!\n\n", READ_ID, _queryFile.line_nr(), _config.MAX_READ_LENGTH);
			return -1;
		}
		else if (strlen(READ) == 0) {
			cerr << "ERROR: Cannot find read sequence of read " << *this << "'!\n";
			exit(0);
		}
		for (int i=0; i<(int)strlen(READ); i++)
			if (READ[i]=='-')
			{
				//fprintf(stderr, "replaced '-' with 'N'\n")  ;
				READ[i]='N' ;
			}

		if (strcspn(READ, "aAcCgGtTnNrRyYmMkKwWsSbBdDhHvV") != 0) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu contains non-IUPAC characters! It will be omitted!\n\n", READ_ID, _queryFile.line_nr());
			return -1;
		}

		if (strlen(READ) < _config.INDEX_DEPTH) {
			fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu is shorter than the specified seedlength! It will be omitted!\n\n", READ_ID, _queryFile.line_nr());
			return -1;
		}

		READ_LENGTH = strlen(READ);

		if (_config.FIXTRIM_STRATEGY_LEN < READ_LENGTH) {
			READ[_config.FIXTRIM_STRATEGY_LEN] = '\0';
			READ_LENGTH = _config.FIXTRIM_STRATEGY_LEN;
		}

		//Fixtrim left and right strategy
		int trimborders= _config.FIXTRIMRIGHT_STRATEGY_LEN + _config.FIXTRIMLEFT_STRATEGY_LEN;
		if (trimborders >0){
			if (READ_LENGTH - trimborders < _config.INDEX_DEPTH){
				cerr <<	"ERROR: fixtrim alignment length too short\n";
				exit(0) ;
			}
			copyFrom(*this, _config.FIXTRIMLEFT_STRATEGY_LEN, _config.FIXTRIMRIGHT_STRATEGY_LEN);
			
		}

		READ_PE_FLAG = 0;
		strcpy(READ_QUALITY[0], READ) ;
		memset(READ_QUALITY[0], 'g', strlen(READ)) ;

		//READ_QUALITY[0] = (char*)"";
		READ_FORMAT = 1;
	}
	else 
		if (_location.has_extension(".sam") || _config.REALIGN_READS)
		{
			prealigned_info = new struct t_prealigned ;
			
			/////// SAM input file ///////
			char const *rid = strtok_r(line, "\t", &saveptr);
			if (rid == NULL) {
				fprintf(stderr, "ERROR: Read ID is empty, line %lu!\n", _queryFile.line_nr());
				exit(0);
			}
			if ((int)strlen(rid) == linelen) {
				fprintf(stderr, "ERROR: wrong read input data format, line %lu!\n", _queryFile.line_nr());
				exit(0);
			}
			strcpy(READ_ID, rid);

			tmp = strtok_r(NULL, "\t", &saveptr);
			if (tmp == NULL) {
				fprintf(stderr, "ERROR: SAM flags are empty, line %lu!\n", _queryFile.line_nr());
				exit(0);
			}
			prealigned_info->sam_flags = atol(tmp);

			char *tok = strtok_r(NULL, "\t", &saveptr); 
			if (tok == NULL) {
				fprintf(stderr, "ERROR: sam contig name field is empty, line %lu!\n", _queryFile.line_nr());
				exit(0);
			}
			prealigned_info->contig=std::string(tok) ;

						
			tmp = strtok_r(NULL, "\t", &saveptr);
			if (tmp == NULL) {
				fprintf(stderr, "ERROR: sam position field is empty, line %lu!\n", _queryFile.line_nr());
				exit(0);
			}
			prealigned_info->start_position = atoi(tmp);


			tmp = strtok_r(NULL, "\t", &saveptr);
			if (tmp == NULL) {
				fprintf(stderr, "ERROR: sam quality field is empty, line %lu!\n", _queryFile.line_nr());
				exit(0);
			}
			prealigned_info->quality = atoi(tmp);

			tok = strtok_r(NULL, "\t", &saveptr);
			if (tok == NULL) {
				fprintf(stderr, "ERROR: sam cigar field is empty, line %lu!\n", _queryFile.line_nr());
				exit(0);
			}
			prealigned_info->cigar=std::string(tok) ;

			for (int k=0; k<3; k++)
			{
				tok = strtok_r(NULL, "\t", &saveptr);
				if (tok == NULL) {
					fprintf(stderr, "ERROR: paired field %i is empty, line %lu!\n", k, _queryFile.line_nr());
					exit(0);
				}
				prealigned_info->paired_info.push_back(std::string(tok)) ;
			}

			tok = strtok_r(NULL, "\t", &saveptr);
			if (tok == NULL) {
				fprintf(stderr, "ERROR: Read sequence is empty, line %lu!\n", _queryFile.line_nr());
				exit(0);
			}
			
			if (strlen(tok) > _config.MAX_READ_LENGTH) {
				fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu is longer than the max read length (=%zu)! It will be omitted!\n\n", READ_ID, _queryFile.line_nr(), _config.MAX_READ_LENGTH);
				return -1;
			}
			if (strcspn(tok, "aAcCgGtTnNrRyYmMkKwWsSbBdDhHvV") != 0) {
				fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu contains non-IUPAC characters! It will be omitted!\n\n", READ_ID, _queryFile.line_nr());
				return -1;
			}
			strcpy(READ, tok);
			READ_LENGTH = strlen(tok);
			if (READ_LENGTH < _config.INDEX_DEPTH) {
				fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu is shorter than the specified seedlength! It will be omitted!\n\n", READ_ID, _queryFile.line_nr());
				return -1;
			}

			// optional read qualities:
			for (int i = 0; i < _maxNrQualities; ++i) {
				char const *qual = strtok_r('\0', "\t\n", &saveptr);
				if (qual == NULL)
					break;
				//TODO Jörg removed the length test: intentionally?
				if (::strlen(qual) !=  READ_LENGTH)
					fprintf(stderr, "WARNING: quality length != read length for read %s in line %lu (going ahead with possibly indeterministic values)\n", READ_ID, _queryFile.line_nr());
				::strncpy(READ_QUALITY[i], qual, READ_LENGTH);
			}
			
			reconstruct_exons_from_cigar() ;

			READ_FORMAT = 3;
		}
		else {
			/////// Flatfile input ///////
			char const *rid = strtok_r(line, "\t", &saveptr);
			if (rid == NULL) {
				fprintf(stderr, "ERROR: Read ID is empty, line %lu!\n", _queryFile.line_nr());
				exit(0);
			}
			if ((int)strlen(rid) == linelen) {
				fprintf(stderr, "ERROR: wrong read input data format, line %lu!\n", _queryFile.line_nr());
				exit(0);
			}
			strcpy(READ_ID, rid);
			
			char *tok = strtok_r(NULL, "\t", &saveptr);
			if (tok == NULL) {
				fprintf(stderr, "ERROR: Read sequence is empty, line %lu!\n", _queryFile.line_nr());
				exit(0);
			}
			
			if (strlen(tok) > _config.MAX_READ_LENGTH) {
				fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu is longer than the max read length (=%zu)! It will be omitted!\n\n", READ_ID, _queryFile.line_nr(), _config.MAX_READ_LENGTH);
				return -1;
			}
			//printf("%s sp: %d\n",READ,(int) strcspn(READ, "A"));
			if (strcspn(tok, "aAcCgGtTnNrRyYmMkKwWsSbBdDhHvV") != 0) {
				fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu contains non-IUPAC characters! It will be omitted!\n\n", READ_ID, _queryFile.line_nr());
				return -1;
			}
			
			strcpy(READ, tok);
			READ_LENGTH = strlen(tok);
			if (READ_LENGTH < _config.INDEX_DEPTH) {
				fprintf(stderr, "\n!!! WARNING: Read '%s' in line %lu is shorter than the specified seedlength! It will be omitted!\n\n", READ_ID, _queryFile.line_nr());
				return -1;
			}
			
			tmp = strtok_r(NULL, "\t", &saveptr);
			if (tmp == NULL) {
				fprintf(stderr, "ERROR: Paired-end flag is empty, line %lu!\n", _queryFile.line_nr());
				exit(0);
			}
			READ_PE_FLAG = atoi(tmp);
			
			// optional read qualities:
			for (int i = 0; i < _maxNrQualities; ++i) {
				char const *qual = strtok_r('\0', "\t\n", &saveptr);
				if (qual == NULL)
					break;
				//TODO Jörg removed the length test: intentionally?
				if (::strlen(qual) !=  READ_LENGTH)
					fprintf(stderr, "WARNING: quality length != read length for read %s in line %lu (going ahead with possibly indeterministic values)\n", READ_ID, _queryFile.line_nr());
				::strncpy(READ_QUALITY[i], qual, READ_LENGTH);
			}
			
			if (_config.FIXTRIM_STRATEGY_LEN < READ_LENGTH) {
				READ[_config.FIXTRIM_STRATEGY_LEN] = '\0';
				READ_QUALITY[0][_config.FIXTRIM_STRATEGY_LEN] = '\0';
				READ_QUALITY[1][_config.FIXTRIM_STRATEGY_LEN] = '\0';
				READ_QUALITY[2][_config.FIXTRIM_STRATEGY_LEN] = '\0';
				READ_LENGTH = _config.FIXTRIM_STRATEGY_LEN;
			}
			
			//Fixtrim left and right strategy
			int trimborders= _config.FIXTRIMRIGHT_STRATEGY_LEN + _config.FIXTRIMLEFT_STRATEGY_LEN;
			if (trimborders >0){
				if (READ_LENGTH - trimborders < _config.INDEX_DEPTH){
					cerr <<	"ERROR: fixtrim alignment length too short\n";
					exit(0) ;
				}
				copyFrom(*this, _config.FIXTRIMLEFT_STRATEGY_LEN,_config.FIXTRIMRIGHT_STRATEGY_LEN);
				
			}
			
			READ_FORMAT = 2;
		}
	
	return 0;
}

void Read::printOn(std::ostream &out) const {
	out << READ_ID << " in " << _location;
}

void Read::printOn(FILE *file) const 
{
	if (READ_FORMAT == 0)
		fprintf(file, "@%s\n%s\n+\n%s\n", READ_ID, READ, READ_QUALITY[0]);
	else if (READ_FORMAT == 1)
		fprintf(file, ">%s\n%s\n", READ_ID, READ);
	else if (READ_FORMAT == 3)
	{
		assert(prealigned_info!=NULL) ;
		
		fprintf(file, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n", READ_ID, prealigned_info->sam_flags, prealigned_info->contig.c_str(), 
				prealigned_info->start_position, prealigned_info->quality, prealigned_info->cigar.c_str(), 
				prealigned_info->paired_info[0].c_str(), prealigned_info->paired_info[1].c_str(), prealigned_info->paired_info[2].c_str(), 
				READ, READ_QUALITY[0]);
	}
	else
		fprintf(file, "%s\t%s\t%d\t%s\t%s\t%s\n", READ_ID, READ,
				READ_PE_FLAG, READ_QUALITY[0], READ_QUALITY[1], READ_QUALITY[2]);

}

void Read::reconstruct_exons_from_cigar()
{
	// use prealigned_info->cigar and prealigned_info->start_position to fill
	// prealigned_info->exons and prealigned_info->aligned_positions

    size_t cig_pos = prealigned_info->cigar.find_first_of("HSMIDN");
	int genome_pos = prealigned_info->start_position;

	prealigned_info->exon.push_back(-1) ;
	prealigned_info->exon.push_back(-1) ;
    while (cig_pos != std::string::npos) 
	{
        unsigned int length = atoi(prealigned_info->cigar.substr(0, op_pos).c_str());
        switch (prealigned_info->cigar.at(cig_pos)) 
        {
        case 'M': 
            if (prealigned_info->exons[prealigned_info->exons.size()-2]==-1)
				prealigned_info->exons[prealigned_info->exons.size()-2].push_back(genome_pos) ;
            for (unsigned int idx = 0; idx < length; idx++)
                prealigned_info->aligned_positions.push_back(genome_pos++);
			prealigned_info->exons[prealigned_info->exons.size()-1].push_back(genome_pos-1) ;
			break;
		case 'D': 
			genome_pos += length; 
			prealigned_info->exons[prealigned_info->exons.size()-1] = genome_pos-1 ; // extend last exon
			break;
        case 'N': 
			genome_pos += length; 
			prealigned_info->exon.push_back(-1) ;
			prealigned_info->exon.push_back(-1) ;
			break;
        case 'I': 
            for (unsigned int idx = 0; idx < length; idx++)
                prealigned_info->aligned_positions.push_back(-1);
			break;
        case 'S': case 'H': break;
        }
        cig_pos = prealigned_info->cigar.find_first_of("HSMIDN", cig_pos + 1);
    }
	for (unsigned int i=0; i<prealigned_info->exon.size(); i++)
		assert(prealigned_info->exon[i]>=0) ;
	for (unsigned int i=0; i<prealigned_info->exon.size()-1; i++)
		assert(prealigned_info->exon[i]<prealigned_info->exon[i+1]) ;

}

