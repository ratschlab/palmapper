#include <genomemapper/Util.h>

bool Util::doInit() {
	POWER[0] = 1;
	for (int i = 1; i <= MAX_INDEX_DEPTH; i++) {
		POWER[i] = POWER[i - 1] * 4;
	}
	return true;
}

bool Util::_isInitialized = doInit();

int Util::POWER[MAX_INDEX_DEPTH+1];

FILE *Util::openFile(std::string const &name, char const *mode) {
	FILE *ret = fopen(name.c_str(), mode);
	if (ret == NULL) {
		fprintf(stderr, "ERROR : Couldn't open input file %s\n", name.c_str());
		exit(1);
	}
	return ret;
}
