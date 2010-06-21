#pragma once

#include <palmapper/Mapper.h>

class FileReporter : public Mapper::Reporter {
public:
	FileReporter(FILE *out, FILE *sp_out, FILE *left_overs);

	void done() {
		if (_out != stdout)
			fprintf(_out, "#done\n") ;
		if (_sp_out !=stdout)
			fprintf(_sp_out, "#done\n") ;
	}

	void report(Mapper::Result &result);

private:
	static int const _nrResults = 1024;
	Mutex _mutex;
	FILE *_out;
	FILE *_sp_out;
	FILE *_left_overs;
	Mapper::Result *_results[_nrResults];
	volatile int _lastResult;
};
