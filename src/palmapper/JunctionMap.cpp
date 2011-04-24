#include <palmapper/JunctionMap.h>
#include <string>
#include <list>
#include <palmapper/Genome.h>
#include <stdlib.h> 
#include <palmapper/Util.h>
#include <pthread.h>

JunctionMap::JunctionMap(Genome const &genome_, int min_coverage_)
{
	genome = &genome_ ;
	unsigned int nbchr = genome->nrChromosomes();
	
	junctionlist = new std::deque<Junction>[nbchr];
	
	min_coverage=min_coverage_;

	int ret = pthread_mutex_init(&junction_mutex, NULL) ;// = PTHREAD_MUTEX_INITIALIZER ;
	assert(ret==0) ;
	
}


JunctionMap::~JunctionMap()
{
	delete[] junctionlist;	
}


void JunctionMap::filter_junctions()
{
	lock() ;

	int total=0, 
		used_nonconsensus=0, 
		used_consensus=0,
		filtered_consensus=0, 
		filtered_nonconsensus=0 ;
	
	for (unsigned int chr=0; chr < genome->nrChromosomes(); chr++)
	{
		if (junctionlist[chr].empty())
			continue;
		
		// create copy of list
		std::deque<Junction>::iterator it=junctionlist[chr].begin(); 
		std::deque<Junction> list  ;
		while (!junctionlist[chr].empty() and it!=junctionlist[chr].end())
		{
			list.push_back(*it) ;
			it++ ;
		}
		junctionlist[chr].clear() ;

		it=list.begin(); 
		
		while (!list.empty() and it!=list.end())
		{
			assert((*it).coverage>0);
			if ( ((*it).coverage<min_coverage) || (((*it).coverage < 2*min_coverage || ((*it).junction_qual<30)) && (!(*it).consensus) ) )
			{
				if ((*it).consensus)
					filtered_consensus++ ;
				else
					filtered_nonconsensus++ ;
			}
			else
			{
				if ((*it).consensus)
					used_consensus++ ;
				else
					used_nonconsensus++ ;
				
				junctionlist[chr].push_back(*it) ;
			}
			it++;
		}
		total+=junctionlist[chr].size();
	}
	unlock() ;
	
	fprintf(stdout,"Number of junctions in database (min support=%i): %i/%i consensus, %i/%i nonconsensus, %i total\n", 
			min_coverage, used_consensus, used_consensus+filtered_consensus, used_nonconsensus, used_nonconsensus+filtered_nonconsensus, total);
}
	
bool comp_junction(std::deque<Junction>::iterator &a, std::deque<Junction>::iterator & b)
{
	if ((*a).start<(*b).start)
		return true ;
	return false ;
}

void JunctionMap::insert_junction(char strand, int chr, int start, int end, bool consensus,  const char* intron_string,
								  int junction_qual, const char *read_id, int coverage = 1)
{
	lock() ;
	
	//Sorted list by donor positions first and then acceptor positions
	Junction j;

	//fprintf(stdout,"%c %i %i %i\n",strand, chr, start, end);

	if (junctionlist[chr].empty())
	{
		j.start=start;
		j.end=end;
		j.coverage=coverage;
		j.strand=strand;
		j.consensus=consensus ;
		j.intron_string = intron_string ;
		j.read_id = read_id ;
		j.junction_qual = junction_qual ;
		junctionlist[chr].push_back(j);

		unlock() ;
		return;
		
	}

	std::deque<Junction>::iterator it = my_lower_bound(junctionlist[chr].begin(), junctionlist[chr].end(), start) ;
	//std::deque<Junction>::iterator it = junctionlist[chr].begin() ;
	
	for (; it!=junctionlist[chr].end(); it++)
	{
		if (start <  (*it).start)
		{
			j.start=start;
			j.end=end;
			j.coverage=coverage;
			j.strand=strand;
			j.consensus = consensus ;
			j.intron_string = intron_string ;
			j.read_id = read_id ;
			j.junction_qual = junction_qual ;
			junctionlist[chr].insert(it,j);

			unlock() ;
			return;
		}
		if (start ==  (*it).start)
		{
			if (end < (*it).end)
			{
				j.start=start;
				j.end=end;
				j.coverage=coverage;
				j.strand=strand;
				j.consensus = consensus ;
				j.intron_string = intron_string ;
				j.read_id = read_id ;
				j.junction_qual = junction_qual ;
				junctionlist[chr].insert(it,j);

				unlock() ;
				return;
			}
			
			if (end == (*it).end)
			{
				if (strand == (*it).strand)
				{
					if ((*it).consensus!=consensus)
					{
						fprintf(stderr, "WARNING: consensus mismatch:\n%s:\t%i-%i %c %i %i %i\n%s:\t%i-%i %c %i %i %i\n", 
								(*it).read_id.c_str(), (*it).start, (*it).end, (*it).strand, (*it).coverage, (*it).consensus, (*it).junction_qual,
								read_id, start, end, strand, coverage, consensus, junction_qual) ;

						assert(0) ; // this should not happen -> please report this bug and try commenting out the assertion
						
						if (!consensus) // try to handle this case 
							(*it).consensus=false ;
					}
					if (junction_qual > (*it).junction_qual)
					{
						(*it).junction_qual=junction_qual ;
						(*it).read_id=read_id ;
					}
					
					(*it).coverage += coverage;

					unlock() ;
					return;
				}
				if (strand == '+')
				{
					j.start=start;
					j.end=end;
					j.coverage=coverage;
					j.consensus = consensus ;
					j.intron_string = intron_string ;
					j.strand = strand;
					j.read_id = read_id ;
					j.junction_qual = junction_qual ;
					junctionlist[chr].insert(it,j);

					unlock() ;
					return;
				}
			}
		}
		continue;
	}
	
	j.start=start;
	j.end=end;
	j.coverage=coverage;
	j.strand=strand;
	j.consensus = consensus ;
	j.intron_string = intron_string ;
	j.read_id = read_id ;
	j.junction_qual = junction_qual ;
	junctionlist[chr].push_back(j);
	
	unlock() ;
	return ;
}

int JunctionMap::init_from_gff(std::string &gff_fname)
{

	fprintf(stdout, "initializing splice site junction list with GFF file %s\n", gff_fname.c_str()) ;

	FILE * fd=Util::openFile(gff_fname.c_str(), "r") ;
	if (!fd)
		return -1 ;


	int exon_lines=0 ;
	int intron_lines=0 ;

	while (!feof(fd))
	{
		char chr_name[1000], source[1000], type[1000], properties[1000], strand, tmp1, tmp2 ;
		int start, end ;

		Util::skip_comment_lines(fd) ;
		
		int num = fscanf(fd, "%1000s\t%1000s\t%1000s\t%i\t%i\t%c\t%c\t%c\t%1000s\n", chr_name, source, type, &start, &end, &tmp1, &strand, &tmp2, properties) ;  
		if (num!=9)
		{
			if (feof(fd))
				break ;
			fprintf(stdout, "gff line only contained %i columns, aborting\n", num) ;
		}
		
		if (strcmp(type, "intron")==0)
		{

			int chr_idx = genome->find_desc(chr_name) ;
			if (chr_idx==-1)
			{
				fprintf(stderr, "chromosome %s not found. known chromosome names:\n", chr_name) ;
				genome->print_desc(stderr) ;
				return -1 ;
			}
			
			std::string tmp(properties);
			
			int pos_cov=tmp.find("Confirmed=");
			if (pos_cov>0)
				pos_cov += strlen("Confirmed=") ;
			else
			{
				pos_cov=tmp.find("Note=");
				if (pos_cov>0)
					pos_cov += strlen("Note=") ;
			}
			int coverage = 1;
			if (pos_cov>0)
				coverage= atoi(tmp.substr(pos_cov).c_str());

			bool nonconsensus=false ;
			char* intron_string = strdup("") ;
			int pos_cons=tmp.find("Nonconsensus=");

			if (pos_cons>0)
			{
				pos_cons += strlen("Nonconsensus=") ;
				nonconsensus = atoi(tmp.substr(pos_cons).c_str());

				int pos_intron=tmp.find("IntronSeq=");
				if (pos_intron>0)
				{
					pos_intron += strlen("IntronSeq=") ;
					intron_string = strdup(tmp.substr(pos_intron).c_str()) ;
				}
			}
			
			int junction_qual = 0 ;
			int pos_qual = tmp.find("BestSplit=") ;
			if (pos_qual>0)
			{
				pos_qual += strlen("BestSplit=") ;
				junction_qual = atoi(tmp.substr(pos_qual).c_str());
			}

			char * read_id = strdup("gff") ;
			int pos_id=tmp.find("ReadID=");
			if (pos_id>0)
			{
				pos_id += strlen("ReadID=") ;
				read_id = strdup(tmp.substr(pos_id).c_str()) ;
			}
			
			insert_junction(strand,chr_idx,start, end, !nonconsensus, intron_string, junction_qual, read_id, coverage);

			free(intron_string) ;
			free(read_id) ;
			
			intron_lines++;
		}
		

		if (strcmp(type, "exon")==0)
		{
			exon_lines++ ;
		}

		if (intron_lines%100000==0)
			fprintf(stdout, "read %i intron lines\n", intron_lines) ;
		
	}
	fclose(fd) ;

	fprintf(stdout, "read %i intron lines\n", intron_lines) ;
	
	return 0 ;

}

int JunctionMap::report_to_gff(std::string &gff_fname)
{
	lock() ;

	int nb_introns=0;
	
	fprintf(stdout, "report splice site junction list in GFF file %s\n", gff_fname.c_str()) ;
	
	FILE * fd=Util::openFile(gff_fname.c_str(), "w") ;
	if (!fd)
		return -1 ;	
	for (unsigned int i=0; i<genome->nrChromosomes(); i++){
		
		const char * chr= genome->get_desc(i);
		std::deque<Junction>::iterator it;
		
		for (it=junctionlist[i].begin(); it!=junctionlist[i].end(); it++){			
			fprintf(fd,"%s\tpalmapper\tintron\t%i\t%i\t.\t%c\t.\tID=intron_%i;Confirmed=%i;BestSplit=%i;ReadID=%s",
					chr,(*it).start,(*it).end,(*it).strand,nb_introns,(*it).coverage, (*it).junction_qual, (*it).read_id.c_str());
			if (!(*it).consensus)
				fprintf(fd,";Nonconsensus=1;IntronSeq=%s\n", (*it).intron_string.c_str());
			else
				fprintf(fd,"\n");
			nb_introns++;
		}
	
	}
	fclose(fd) ;
	fprintf(stdout, "report %i introns\n", nb_introns) ;	

	unlock() ;

	return 0;
	
}

int JunctionMap::init_from_gffs(std::string &gff_fname)
{

	int previousfound=0;
	int found=gff_fname.find(",");
	std::string filename;
	
	while (found >= 0)
	{
		
		filename = gff_fname.substr(previousfound, found-previousfound);
		int ret = init_from_gff(filename);
		if (ret!=0)
			return ret;
	   
		previousfound=found+1;
		found=gff_fname.find(",",found+1);
	}
	
	filename=gff_fname.substr(previousfound);
	int ret=init_from_gff(filename);
	if (ret!=0)
		return  ret;
	
	filter_junctions();
	
	return ret;
	
}
