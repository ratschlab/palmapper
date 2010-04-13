#pragma once

#include <stdio.h>
#include <stdlib.h>

#include <string>

#include <genomemapper/Config.h>

class Util {
public:
	static int POWER[MAX_INDEX_DEPTH+1];
	static FILE *openFile(std::string const &name, char const *mode) {
		return openFile(name.c_str(), mode);
	}

	static FILE *openFile(char const *name, char const *mode);

	static void skip_comment_lines(FILE* fd)
		{
			const int buffer_len=1000 ;
			char buffer[buffer_len] ;
			
			while (!feof(fd))
			{
				long pos = ftell(fd) ;
				fgets(buffer, buffer_len, fd) ;
				if (buffer[0]!='#' && buffer[0]!=0)
				{
					fseek(fd, pos, SEEK_SET) ;
					break ;
				}
				//fprintf(stdout, "skipped comment line: %s\n", buffer) ;
			}
		}

	static std::string read_line(FILE* fd)
		{
			char buf[10000] ;
			
			int narg = fscanf(fd, "%10000s:\t", buf);
			assert(narg==1);
			return std::string(buf) ;
		}

private:
	static bool _isInitialized;
	static bool doInit();
};
