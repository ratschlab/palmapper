// Authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include "mkindex.h"

int usage() {
	printf("\ngmindex v%s\n", VERSION);
	printf("developed by Korbinian Schneeberger, Stephan Ossowski and Joerg Hagmann\n");
	printf("Max Planck Institute for Developmental Biology, Tübingen, Germany, 2008\n\n");
	printf("USAGE: mkindex [options]\n");
	printf("\n");
	printf("mandatory:\n");
	printf(" -i STRING  input fastafile\n");
	//printf(" -x STRING  index filename\n");
	//printf(" -m STRING  map forward index filename\n");	
	//printf(" -c STRING  map reverse index filename\n");	
	//printf(" -t STRING  meta index filename\n");
	printf("\n");
	printf("optional:\n");
	printf(" -s INT     seed length, range 5 to 13 (12)\n");
	printf(" -r         do not build reverse index\n");
	printf(" -v         verbose (silent)\n");
	printf("\n");

	return 0;
}
