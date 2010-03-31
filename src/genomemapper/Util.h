#pragma once

#include <stdlib.h>

#include <genomemapper/Config.h>

class Util {
public:
	static int POWER[MAX_INDEX_DEPTH+1];

private:
	static bool _isInitialized;
	static bool doInit();
};
