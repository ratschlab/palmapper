#include <new>
#include <stdarg.h>
#include <stdio.h>
#include <iostream>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <cxxabi.h>

#include <palmapper/Util.h>



bool Util::doInit() {
	POWER[0] = 1;
	for (int i = 1; i <= MAX_INDEX_DEPTH; i++) {
		POWER[i] = POWER[i - 1] * 4;
	}

	return true;
}

bool Util::_isInitialized = doInit();

int Util::POWER[MAX_INDEX_DEPTH+1];

FILE *Util::openFile(char const *name, char const *mode) {
	FILE *ret = fopen(name, mode);
	if (ret == NULL) {
		fprintf(stderr, "ERROR : Couldn't open file %s for %s\n", name, mode);
		exit(1);
	}
	return ret;
}

gzFile Util::gzopenFile(char const *name, char const *mode) {
	gzFile ret = gzopen(name, mode);
	if (ret == NULL) {
		fprintf(stderr, "ERROR : Couldn't open file %s for %s\n", name, mode);
		exit(1);
	}
	return ret;
}

#ifndef __APPLE__ 
// this goes wrong with Mac's gcc

void *operator new(size_t size) {
	static std::bad_alloc ex;
	void *ret = operator new(size, std::nothrow);
	if (ret == NULL)
	{
	  fprintf(stderr, "Warning: allocation of %ld bytes failed\n", size) ;
	  print_stacktrace(stderr, 10) ;

	  throw ex;
	}
	return ret;
	}

#endif

void fprintf(std::ostream *out, char const *format, ...) {
	va_list ap;
	va_start(ap, format);
	char *strp;
	int len = vasprintf(&strp, format, ap);
	va_end(ap);

	if (len == -1) {
		fprintf(stderr, "Unable to format string. Exiting!!\n");
		exit(1);
	}
	out->write(strp, len);
	free(strp);
}
