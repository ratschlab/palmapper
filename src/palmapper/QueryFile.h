#pragma once

#include <stdio.h>

#include <palmapper/Read.h>

class QueryFile {
public:
	QueryFile(std::string filename);
	~QueryFile();

	Read *next_read();

	bool next_line(char *buf, int maxLen);

	unsigned long line_nr() const {
		return _lineNr;
	}

	unsigned int max_length() const {
		return _maxReadLen;
	}

	static int determine_read_length(std::string const &filename);
private:
	FILE *_file;
	int _lineNr;
	unsigned int _maxReadLen;
	std::string _filename;
};
