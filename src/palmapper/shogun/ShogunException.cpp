/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Written (W) 1999-2009 Soeren Sonnenburg
 * Copyright (C) 1999-2009 Fraunhofer Institute FIRST and Max-Planck-Society
 */

#include "ShogunException.h"
#include "Signal.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

using namespace shogun;

ShogunException::ShogunException(const char* str)
{
#ifndef WIN32
	CSignal::unset_handler();
#endif
   val = (char*) malloc(sizeof(char)*4096);
   if (val)
       strncpy(val,str,4096);
   else
   {
       fprintf(stderr, "Could not even allocate memory for exception - dying.\n");
       exit(1);
   }
}

ShogunException::~ShogunException()
{
}
