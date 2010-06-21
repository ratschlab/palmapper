#include <assert.h>

#include <palmapper/FileReporter.h>

FileReporter::FileReporter(FILE *out, FILE *sp_out) {
	_out = out;
	_sp_out = sp_out;
	_lastResult = -1;
	memset(_results, 0, sizeof(Mapper::Result*) * _nrResults);
}

void FileReporter::report(Mapper::Result &result) {
	Mutex::Locker locker(_mutex);
	//printf("Delivering result %i\n", result._read.getNr());
	if (result._read.getNr() >= _lastResult + _nrResults) {
		// buffer overrun
		//TODO: the following is not really an option...
		printf("!!!!!!!! Bad thing: result lost due to buffer overrun\n");
		delete &result;
		return;
	}
	assert(_results[result._read.getNr() % _nrResults] == NULL);
	_results[result._read.getNr() % _nrResults] = &result;
	for (int i = _lastResult + 1;  ; ++i) {
		int pos = i % _nrResults;
		if (_results[pos] == NULL)
			return;
		//printf("Writing result %i\n", i);
		Mapper::Result &r(*_results[pos]);
		assert(&r != NULL);
		assert(r._read.getNr() == i);
		if (r._state == Mapper::HitFound)
			r._readMappings.topAlignments().end_top_alignment_record(r._read, _out, _sp_out, r._rtrim_cut, r._polytrim_cut_start, r._polytrim_cut_end);
		_lastResult = i;
		delete &r;
		_results[pos] = NULL;
	}
}
