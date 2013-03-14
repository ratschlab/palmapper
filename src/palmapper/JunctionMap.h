#pragma once

#include <assert.h>
#include <palmapper/Genome.h>
#include <palmapper/GenomeMaps.h>
#include <deque>
#include <string>
#include <stdlib.h> 

struct junction_str {
	int start; 
	int end; 
	int coverage;
	bool consensus ;
	char strand;
	std::string read_id ;
	std::string intron_string ;
	int junction_qual ;
};
typedef struct junction_str Junction;

struct exon_str {
	int start; 
	int end; 
	int chr;
	char strand;
};
typedef struct exon_str Exon;

class JunctionMap
{

public:
	JunctionMap(Genome const &genome_, int anno_pseudo_coverage, std::vector<const char*> ACC_CONSENSUS_, std::vector<const char*> DON_CONSENSUS_, std::vector<const char*> ACC_CONSENSUS_REV_, std::vector<const char*> DON_CONSENSUS_REV_ ) ;
	~JunctionMap() ;

	void insert_junction(char strand, int chr, int start, int end, bool consensus, const char* intron_string, int junction_qual, const char* read_id, int coverage);
	int init_from_gffs(std::string &gff_fname);
	int report_to_gff(std::string &gff_fname);
	void filter_junctions(int min_coverage, int min_junction_qual, int filter_by_map, const GenomeMaps & genomemaps, int verbosity);

	std::deque<Junction> * junctionlist_by_start ;
	std::deque<Junction> * junctionlist_by_end ;

	void lock() 
	{ 
		pthread_mutex_lock( &junction_mutex) ; 
	}

	void unlock() 
	{ 
		pthread_mutex_unlock( &junction_mutex) ; 
	}
	

protected:
	bool is_consensus_intron(char strand, int chr, int start, int end);
	
	int init_from_gff(std::string &gff_fname);

	Genome const *genome;
	int anno_pseudo_coverage;

	pthread_mutex_t junction_mutex;
	std::vector<const char*> ACC_CONSENSUS, DON_CONSENSUS, ACC_CONSENSUS_REV, DON_CONSENSUS_REV ;	
   	


};

inline std::deque<Junction>::iterator  my_lower_bound_by_start ( std::deque<Junction>::iterator first, std::deque<Junction>::iterator  last, const int& value )
{
	std::deque<Junction>::iterator it;
	long int count, step;
	count = distance(first,last);

	while (count>0)
	{
		it = first;
		step = count/2;
		advance (it,step);
		
		if ( (*it).start < value) 
		{
			first=it; //++it;
			count-=step+1;
		}
		else count=step;
	}
	return first;
}

inline std::deque<Junction>::iterator  my_lower_bound_by_end ( std::deque<Junction>::iterator first, std::deque<Junction>::iterator  last, const int& value )
{
	std::deque<Junction>::iterator it;
	long int count, step;
	count = distance(first,last);

	while (count>0)
	{
		it = first;
		step = count/2;
		advance (it,step);
		
		if ( (*it).end < value) 
		{
			first=it; //++it;
			count-=step+1;
		}
		else count=step;
	}
	return first;
}


inline int get_combination_span ( std::vector<std::deque<Junction>::iterator> combination, bool ascending)
{
    if (combination.size() < 2)
        return 0;

    int span = 0;
    // this assumes that the vector is sorted by ascending /descending (parameter)
    if (ascending) {
        for (size_t i = 0; i < combination.size() - 1; i++) {
            assert(combination.at(i+1)->start > combination.at(i)->end);
            span += (combination.at(i+1)->start - combination.at(i)->end - 1);
        }
    } else {
        for (size_t i = combination.size() - 1; i > 0; i--) {
            assert(combination.at(i)->start > combination.at(i-1)->end);
            span += (combination.at(i)->start - combination.at(i-1)->end - 1);
        }
    }
    return span;
}
