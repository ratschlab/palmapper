#pragma once

#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <iostream>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <cxxabi.h>

#include <zlib.h>

#include <palmapper/Config.h>

class Util {
public:
	static int POWER[MAX_INDEX_DEPTH+1];
	static FILE *openFile(std::string const &name, char const *mode) {
		return openFile(name.c_str(), mode);
	}

	static FILE *openFile(char const *name, char const *mode);
	static gzFile gzopenFile(char const *name, char const *mode);

	static int skip_comment_lines(FILE* fd)
		{
			const int buffer_len=5000 ;
			char buffer[buffer_len] ;
			int line_num=0 ;
			
			while (!feof(fd))
			{
				long pos = ftell(fd) ;
				if (fgets(buffer, buffer_len, fd) == NULL)
					return line_num ;
				if (buffer[0]!='#' && buffer[0]!=0)
				{
					fseek(fd, pos, SEEK_SET) ;
					break ;
				}
				line_num++ ;
				//fprintf(stdout, "skipped comment line: %s\n", buffer) ;
			}
			return line_num ;
		}

	static std::string read_line(FILE* fd)
		{
			char buf[10000] ;
			int narg = fscanf(fd, "%10000s", buf);
			if (narg!=1)
				fprintf(stderr, "narg=%i\n", narg) ;
			assert(narg==1);
			return std::string(buf) ;
		}

private:
	static bool _isInitialized;
	static bool doInit();
};

void fprintf(std::ostream *out, char const *format, ...);


// stacktrace.h (c) 2008, Timo Bingmann from http://idlebox.net/
// published under the WTFPL v2.0

/** Print a demangled stack backtrace of the caller function to FILE* out. */
inline void print_stacktrace(FILE *out = stderr, unsigned int max_frames = 63)
{
	
    fprintf(out, "stack trace:\n");

    // storage array for stack trace address data
    void* addrlist[max_frames+1];

    // retrieve current stack addresses
    int addrlen = backtrace(addrlist, sizeof(addrlist) / sizeof(void*));

    if (addrlen == 0) {
        fprintf(out, "  <empty, possibly corrupt>\n");
        return;
    }
	

    // resolve addresses into strings containing "filename(function+address)",
    // this array must be free()-ed
    char** symbollist = backtrace_symbols(addrlist, addrlen);

    // allocate string which will be filled with the demangled function name
    size_t funcnamesize = 256;
	
    char* funcname = (char*)malloc(funcnamesize);
	

    // iterate over the returned symbol lines. skip the first, it is the
    // address of this function.
    for (int i = 1; i < addrlen; i++)
    {
        char *begin_name = 0, *begin_offset = 0, *end_offset = 0;

        // find parentheses and +address offset surrounding the mangled name:
        // ./module(function+0x15c) [0x8048a6d]
        for (char *p = symbollist[i]; *p; ++p)
        {
            if (*p == '(')
                begin_name = p;
			
            else if (*p == '+')
                begin_offset = p;
			
            else if (*p == ')' && begin_offset) {
                end_offset = p;
                break;
            }
        }
		

        if (begin_name && begin_offset && end_offset
            && begin_name < begin_offset)
        {
            *begin_name++ = '\0';
            *begin_offset++ = '\0';
            *end_offset = '\0';

            // mangled name is now in [begin_name, begin_offset) and caller
            // offset in [begin_offset, end_offset). now apply
            // __cxa_demangle():

            int status;
			
            char* ret = abi::__cxa_demangle(begin_name,
                                            funcname, &funcnamesize, &status);
			
            if (status == 0) {
                funcname = ret;
				// use possibly realloc()-ed string
                fprintf(out, "  %s : %s+%s\n",
                        symbollist[i], funcname, begin_offset);
            }
            else {
                // demangling failed. Output function name as a C function with
                // no arguments.
                fprintf(out, "  %s : %s()+%s\n",
                        symbollist[i], begin_name, begin_offset);
            }
        }
        else
        {
            // couldn't parse the line? print the whole line.
            fprintf(out, "  %s\n", symbollist[i]);
        }
    }
	

    free(funcname);
    free(symbollist);
}
