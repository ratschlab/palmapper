#pragma once

#include <stdio.h>

#include <palmapper/Read.h>
#include <lang/Thread.h>

using lang::Mutex;

class QueryFile {

public:
	QueryFile(std::string filename);
	~QueryFile();

	Read *next_read();

	bool next_read(Read &read);

	bool next_line(char *buf, int maxLen);

	unsigned long line_nr() const {
		return _lineNr;
	}

	unsigned int max_length() const {
		return _maxReadLen;
	}

	int read_count() const {
		return _readCount;
	}

	void reset_read_count() {
		_readCount=0;
	}

	static int determine_read_length(std::string const &filename);

private:
	FILE *_file;
	int _lineNr;
	unsigned int _maxReadLen;
	std::string _filename;
	int _readCount;
	Mutex _mutex;
};
