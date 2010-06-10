// Authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include <palmapper/Util.h>
#include <palmapper/palmapper.h>

FILE *init_output_file(Config &config)
{
	if (config.OUT_FILE_NAME.length() > 0)
		return Util::openFile(config.OUT_FILE_NAME, "w");
	return stdout;
}

FILE *init_spliced_output_file(Config &config, FILE *OUT_FP)
{
	if (config.SPLICED_HITS && config.SPLICED_OUT_FILE_NAME.length() > 0)
		return Util::openFile(config.SPLICED_OUT_FILE_NAME, "w");
	return OUT_FP;
}
