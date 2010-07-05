#include <iostream>

#include <palmapper/QueryFile.h>
#include <palmapper/Read.h>

using namespace std;

QueryFile::QueryFile(std::string filename) {
	_filename = filename;
	_file = Util::openFile(filename, "r");
	_lineNr = 0;
	_maxReadLen = 0;
	_readCount = 0;
}

QueryFile::~QueryFile() {
	::fclose(_file);
}

bool QueryFile::next_line(char *buf, int maxLen) {
	if (::fgets(buf, maxLen, _file) == NULL)
		return false;
	++_lineNr;
	return true;
}

Read *QueryFile::next_read() {
	Read *read = new Read(*this);
	if (!next_read(*read)) {
		delete read;
		return NULL;
	}
	return read;
}

bool QueryFile::next_read(Read &read) {
	Mutex::Locker locker(_mutex);
	if (read.read_short_read() > 0) {
		if (_readCount == 0)
			cerr << "\n!!! WARNING: Input read file '" << _filename << "' contains no usable read!\n\n";
		return false;
	}
	read._nr = _readCount++;
	if (read.length() > _maxReadLen)
		_maxReadLen = read.length();
	return true;
}

int QueryFile::determine_read_length(std::string const &filename) {
	QueryFile file(filename);
	int const sample_size = 10000;
	int sum_read_length = 0;
	int nr_read = 0;

	Read r(file);
	for (;  nr_read < sample_size; ++nr_read) {
		if (!file.next_read(r))
			break;
		sum_read_length += r.length();
	}
	return sum_read_length / nr_read;
}
