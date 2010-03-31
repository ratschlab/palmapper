#pragma once

#include <stdio.h>
#include <stdlib.h>

#include <string>

#include <genomemapper/Config.h>

class Util {
public:
	static int POWER[MAX_INDEX_DEPTH+1];
	static FILE *openFile(std::string const &name, char const *mode);
private:
	static bool _isInitialized;
	static bool doInit();
};
