#pragma once

#include <palmapper/Genome.h>
#include <list>
#include <string>
#include <stdlib.h> 

struct junction_str {
	int start; 
	int end; 
	int coverage;
	char strand;
};
typedef struct junction_str Junction;


class JunctionMap
{

public:
	JunctionMap(Genome const &genome_) ;
	~JunctionMap() ;

	void insert_junction(char strand, int chr, int start, int end,int coverage);
	int init_from_gff(std::string &gff_fname);
	int report_to_gff(std::string &gff_fname);
	
	std::list<Junction> * junctionlist ;

protected:
	Genome const *genome;
};

