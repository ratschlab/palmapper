#pragma once

#include <palmapper/Genome.h>
#include <list>
#include <string>
#include <stdlib.h> 

struct junction_str {
	int start; 
	int end; 
	int coverage;
	bool consensus ;
	char strand;
};
typedef struct junction_str Junction;


class JunctionMap
{

public:
	JunctionMap(Genome const &genome_, int min_coverage_) ;
	~JunctionMap() ;

	void insert_junction(char strand, int chr, int start, int end, bool consensus, int coverage);
	int init_from_gffs(std::string &gff_fname);
	int report_to_gff(std::string &gff_fname);
	void filter_junctions();
	

	std::list<Junction> * junctionlist ;

protected:

	int init_from_gff(std::string &gff_fname);

	Genome const *genome;
	int min_coverage;

	pthread_mutex_t junction_mutex;
	
};

