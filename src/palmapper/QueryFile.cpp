#include <palmapper/QueryFile.h>
#include <palmapper/Read.h>

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
	if (read.read_short_read() > 0)
		return false;
	++_readCount;
	if (read.length() > _maxReadLen)
		_maxReadLen = read.length();
	return true;
}

int QueryFile::determine_read_length(std::string const &filename) {
	QueryFile file(filename);
	int const sample_size = 10000;
	int sum_read_length = 0;
	int nr_read = 0;

	for (;  nr_read < sample_size; ++nr_read) {
		Read *r = file.next_read();
		if (r == NULL)
			break;
		sum_read_length += r->length();
	}
	return sum_read_length / nr_read;
}
