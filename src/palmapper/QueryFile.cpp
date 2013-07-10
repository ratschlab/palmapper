#include <iostream>
#include <palmapper/palmapper.h>

#include <palmapper/QueryFile.h>
#include <palmapper/Read.h>

using namespace std;

QueryFile::QueryFile(std::vector<std::string> const &filenames,std::vector<int> const &strands) : _filenames(filenames), _strands(strands) {
	_currentFile = -1;
	_current_strand=-1;
	_readCount = 0;
	_file = NULL;
	open_next_file();
	passed_first_read=false ;
	passed_last_read=false ;
}


QueryFile::~QueryFile() {
	::gzclose(_file);
}

bool QueryFile::next_line(char *buf, int maxLen) {
	if (::gzgets(_file, buf, maxLen) == NULL)
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

bool QueryFile::next_read(Read &read, int &strand) {
	Mutex::Locker locker(_mutex);
	int ret ;

	while ((ret=read.read_short_read(passed_first_read, passed_last_read)) > 0) {
		if (ret==2)
			continue ;
		//reset first/last read logic
		passed_first_read=false;
		passed_last_read=false ;
		if (!open_next_file()) {
			if (_readCount == 0 && _config.VERBOSE>0)
				cerr << "\n!!! WARNING: None of the given file(s) contain any usable read!\n\n";
			return false;
		}
	}
	strand=_current_strand;
	read._nr = _readCount++;
	if (read.length() > _maxReadLen)
		_maxReadLen = read.length();
	return true;
}

bool QueryFile::next_read(Read &read) {
	Mutex::Locker locker(_mutex);
	int ret = 0 ;
	
	while ((ret=read.read_short_read(passed_first_read, passed_last_read)) > 0) {
		if (ret==2)
			continue ;
		//reset first/last read logic
		passed_first_read=false ;
		passed_last_read=false ;
		if (!open_next_file()) {
			if (_readCount == 0  && _config.VERBOSE>0)
				cerr << "\n!!! WARNING: None of the given file(s) contain any usable read!\n\n";
			return false;
		}
	}
	read._nr = _readCount++;
	if (read.length() > _maxReadLen)
		_maxReadLen = read.length();
	return true ;
}


bool QueryFile::open_next_file() {
	if (++_currentFile >= _filenames.size())
		return false;
	if (_file != NULL)
		::gzclose(_file);
	_file=NULL ;

	_file = Util::gzopenFile(_filenames[_currentFile].c_str(), "r");
	_lineNr = 0;
	_maxReadLen = 0;
	_current_strand=_strands[_currentFile];
	return true;
}

int QueryFile::determine_read_length(std::vector<std::string> const &filenames,std::vector<int> const &strands) {
	QueryFile file(filenames,strands);
	int const sample_size = 10000;
	int sum_read_length = 0;
	int nr_read = 0;

	Read r(file);
	for (;  nr_read < sample_size; ++nr_read) {
		if (!file.next_read(r))
			break;
		sum_read_length += r.length();
	}
	if (nr_read==0)
		return -1 ;
	return sum_read_length / nr_read;
}

void QueryFile::Location::printOn(std::ostream &out) const {
	out << *_filename;
	out << '(';
	out << _lineNr;
	out << ')';
}
