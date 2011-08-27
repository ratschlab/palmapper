#include <vector>
#include <map>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>

typedef void handler_t(int);

static std::map<pthread_t, std::string> current_read_ids ;
static bool show_read_ids=false ;
static bool handler_installed = false ;


void report_current_read_id(std::string read_id) 
{
	pthread_t id=pthread_self() ;
	current_read_ids[id]=read_id ;
}

void handler(int sig)
{
	fprintf(stderr, "\nERROR: SEGSEGV in thread %lu encountered\n\n", pthread_self()) ;
	if (show_read_ids)
	{
		fprintf(stderr, "Read IDs that were recently processed:\n") ;
		
		for (std::map<pthread_t, std::string>::iterator it=current_read_ids.begin(); it!=current_read_ids.end(); it++)
		{
			char c=' ' ;
			if (pthread_self()==it->first)
				c='*' ;
			fprintf(stderr, "\t%lu\t%s\t%c\n", it->first, it->second.c_str(), c) ;
		}
	}
	
	fprintf(stderr, "\n\nTerminating process\n\n") ;
	
    exit(-1) ;

	return ;
}

/* Signal: wrapper for sigaction */
handler_t * Signal(int sig, handler_t *h)
{
    struct sigaction action, old_action;
	
    action.sa_handler = h;
	
    sigemptyset(&action.sa_mask);
	
    action.sa_flags = SA_NODEFER| SA_RESTART ;

    if (sigaction(sig, &action, &old_action)<0)
        fprintf(stderr, "Error installing signal handler!\n");
	

    return old_action.sa_handler;
}

void init_signal_handler(bool show_read_ids_)
{
	if (!handler_installed)
		Signal(SIGSEGV, handler);

	show_read_ids=show_read_ids_ ;
}

