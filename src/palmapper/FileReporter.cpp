#include <assert.h>
#include <sstream>

#include <palmapper/FileReporter.h>
#include <palmapper/print.h>

static clock_t last_warning_time = clock();

FileReporter::FileReporter(FILE *out, FILE *sp_out, FILE *variants_out, FILE *left_overs) {
	_out = out;
	_sp_out = sp_out;
	_variants_out = variants_out;
	_left_overs = left_overs;
	_lastResult = -1;

	//static clock_t last_warning_time = clock();
}

void FileReporter::report(Mapper::Result &result, JunctionMap &junctionmap, VariantMap & variants) 
{
//	fprintf(stdout,"Delivering result %i\n", result._orig.getNr());

	int readNr = result._work.getNr();
	std::stringstream out;
	std::stringstream sp_out;
	std::stringstream variants_out;
	std::stringstream leftOvers;

	if (result._state == Mapper::ReadMapped || (result._state < Mapper::IgnoreResultBound && _config.INCLUDE_UNMAPPED_READS_SAM)){
		result._readMappings.topAlignments().end_top_alignment_record(result._work, &out, &sp_out, &variants_out, result._rtrim_cut, result._polytrim_cut_start, result._polytrim_cut_end,
																	  junctionmap, variants);
	} else {
		if (result._state < Mapper::IgnoreResultBound && _config.LEFTOVER_FILE_NAME.length() > 0) {
			char const *text = "";
			switch (result._state) {
				case Mapper::MappingFailed:
					text = " (read mapping failed)";
					break;
				case Mapper::TooShortAfterTrimming:
					text = "(too short after trimming)";
					break;
				default:
					;
			}

			if (!_config.INCLUDE_UNMAPPED_READS_SAM)
				print_leftovers(result._work, text, &leftOvers);
		}
	}
	delete &result;

	Mutex::Locker locker(_mutex);
	while (readNr >= _lastResult + _nrResults) 
	{
		if ((clock()-last_warning_time)/CLOCKS_PER_SEC>=5) // wait for at most 5 seconds before outputing a new warning
		{
			printf("Warning: small result buffer may degrade performance\n");
			last_warning_time=clock() ;
		}
		_roomLeft.wait(_mutex);
	}
	int index = readNr % _nrResults;
	Entry &e = _results[index];
	assert(!e._used);
	_results[index]._used = true;

	e._out = out.str();
	e._sp_out = sp_out.str();
	e._left_overs = leftOvers.str();
	e._variants_out = variants_out.str();
	

	for (int i = _lastResult + 1;  ; ++i) {
		int pos = i % _nrResults;
		if (!_results[pos]._used)
			break;
		//printf("Writing result %i\n", i);
		Entry &en(_results[pos]);
		assert(en._used);
		print(_out, en._out);
		print(_sp_out, en._sp_out);
		print(_variants_out, en._variants_out);
		print(_left_overs, en._left_overs);
		en._used = false;
		_lastResult = i;
	}
	_roomLeft.notifyAll();
}


void FileReporter::print(FILE *file, std::string &s){
	if (s.length() == 0 )
		return;
	if (::fwrite(s.c_str(), 1, s.length(), file) != s.length()) {
		fprintf(stderr, "Could not write results\n");
		exit(1);
	}
	s.assign("");
}
