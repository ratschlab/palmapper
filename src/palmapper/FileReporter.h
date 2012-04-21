#pragma once

#include <palmapper/Mapper.h>
#include <palmapper/JunctionMap.h>
#include <palmapper/VariantMap.h>
#include <zlib.h>

class FileReporter : public Mapper::Reporter {
	struct Entry {
		Entry() {_used= false;}
		std::string _out;
		std::string _sp_out;
		std::string _variants_out;
		std::string _left_overs;
		bool _used;
	};

public:
	FileReporter(FILE *out, FILE *sp_out, 	gzFile _variants_out, FILE *left_overs);

	void done() {
		if (_out != stdout)
			fprintf(_out, "\n") ;
		//fprintf(_out, "#done\n") ;
		if (_sp_out != NULL && _sp_out != stdout && _sp_out!=_out)
			fprintf(_sp_out, "\n") ;
		if (_variants_out != NULL && _variants_out != stdout)
			gzprintf(_variants_out, "\n") ;

		//fprintf(_sp_out, "#done\n") ;
	}

	void report(Mapper::Result &result, JunctionMap &junctionmap, VariantMap& variants);

private:
	static int const _nrResults = 1024*100;
	static void print(FILE *file, std::string &str);
	Mutex _mutex;
	lang::Signal _roomLeft;
	FILE *_out;
	FILE *_sp_out;
	gzFile _variants_out;
	FILE *_left_overs;
	Entry _results[_nrResults];
	volatile int _lastResult;
};
