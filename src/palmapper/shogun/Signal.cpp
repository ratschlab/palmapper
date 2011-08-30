/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Written (W) 1999-2009 Soeren Sonnenburg
 * Copyright (C) 1999-2009 Fraunhofer Institute FIRST and Max-Planck-Society
 */

//#include "lib/config.h"

#ifndef WIN32

#include <stdlib.h>
#include <signal.h>
#include <string.h>

#include "io.h"
#include "Signal.h"

using namespace shogun;

int CSignal::signals[NUMTRAPPEDSIGS]={SIGINT, SIGURG, SIGSEGV, SIGBUS};
struct sigaction CSignal::oldsigaction[NUMTRAPPEDSIGS];
bool CSignal::active=false;
bool CSignal::cancel_computation=false;
bool CSignal::cancel_immediately=false;
bool CSignal::block_computations=false;

std::map<pthread_t, std::string> CSignal::current_read_ids ;
bool CSignal::show_read_ids = false ;

Mutex CSignal::_mutex;
lang::Signal CSignal::_roomLeft;

CSignal::CSignal()
: CSGObject()
{
}

CSignal::~CSignal()
{
	if (!unset_handler())
		fprintf(stderr, "error uninitalizing signal handler\n");
}

void CSignal::do_show_read_ids()
{
	if (!show_read_ids)
		return ;
	
	fprintf(stderr, "Read IDs that were recently processed:\n") ;
	
	for (std::map<pthread_t, std::string>::iterator it=current_read_ids.begin(); it!=current_read_ids.end(); it++)
	{
		char c=' ' ;
		if (pthread_self()==it->first)
			c='*' ;
		fprintf(stderr, "\t%lu\t%s\t%c\n", it->first, it->second.c_str(), c) ;
	}
}

void CSignal::handler(int signal)
{
	if (signal == SIGINT)
	{
		if (block_computations)
			return ;
		block_computations=true ;
		if (show_read_ids)
			fprintf(stderr, "\nJust exit / Immediately return to prompt / Prematurely finish computations / Show Read IDs / Do nothing (J/I/P/S/D)? ");
		else
			fprintf(stderr, "\nJust exit / Immediately return to prompt / Prematurely finish computations / Do nothing (J/I/P/D)? ");
		char answer=fgetc(stdin);

		if (answer == 'I')
		{
			unset_handler();
			set_cancel(true);
			if (sg_print_error)
				sg_print_error(stdout, "sg stopped by SIGINT\n");
		}
		else if (answer == 'J')
		{
			fprintf(stderr, "Exiting.") ;
			exit(-1) ;
		}
		else if (answer == 'P')
			set_cancel();
		else if (answer == 'S')
			do_show_read_ids();
		else
			fprintf(stderr, "Continuing...\n");
		block_computations=false ;
		_roomLeft.notifyAll();
	}
	else if (signal == SIGURG)
		set_cancel();
	else  if (signal == SIGSEGV || signal==SIGBUS)
	{
		if (signal == SIGSEGV)
			fprintf(stderr, "\nERROR: SEGSEGV in thread %lu encountered\n\n", pthread_self()) ;
		else
			fprintf(stderr, "\nERROR: SEGBUS in thread %lu encountered\n\n", pthread_self()) ;
		do_show_read_ids() ;
		fprintf(stderr, "\n\nTerminating process.\n\n") ;
		
		exit(-1) ;
	}
	else
		fprintf(stderr, "unknown signal %d received\n", signal);
}

bool CSignal::set_handler()
{
	if (!active)
	{
		struct sigaction act;
		sigset_t st;

		sigemptyset(&st);
		for (int32_t i=0; i<NUMTRAPPEDSIGS; i++)
			sigaddset(&st, signals[i]);

#ifndef __INTERIX
		act.sa_sigaction=NULL; //just in case
#endif
		act.sa_handler=CSignal::handler;
		act.sa_mask = st;
		act.sa_flags = 0;

		for (int32_t i=0; i<NUMTRAPPEDSIGS; i++)
		{
			if (sigaction(signals[i], &act, &oldsigaction[i]))
			{
				fprintf(stderr, "Error trapping signals!\n");
				for (int32_t j=i-1; j>=0; j--)
					sigaction(signals[i], &oldsigaction[i], NULL);

				clear();
				return false;
			}
		}

		active=true;
		return true;
	}
	else
		return false;
}

bool CSignal::unset_handler()
{
	if (active)
	{
		bool result=true;

		for (int32_t i=0; i<NUMTRAPPEDSIGS; i++)
		{
			if (sigaction(signals[i], &oldsigaction[i], NULL))
			{
				SG_SPRINT("error uninitalizing signal handler for signal %d\n", signals[i]);
				result=false;
			}
		}

		if (result)
			clear();

		return result;
	}
	else
		return false;
}

void CSignal::clear_cancel()
{
	cancel_computation=false;
	cancel_immediately=false;
}

void CSignal::set_cancel(bool immediately)
{
	cancel_computation=true;

	if (immediately)
		cancel_immediately=true;
}

void CSignal::clear()
{
	clear_cancel();
	active=false;
	memset(&CSignal::oldsigaction, 0, sizeof(CSignal::oldsigaction));
}

void CSignal::toggle_show_read_ids(bool show_ids)
{
	show_read_ids=show_ids ;
}

void CSignal::report_current_read_id(std::string read_id) 
{
	pthread_t id=pthread_self() ;
	current_read_ids[id]=read_id ;
}

#endif //WIN32
