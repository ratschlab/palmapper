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
	
	junctionlist = new std::list<Junction>[nbchr];
	
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
	pthread_mutex_lock( &junction_mutex) ;

	int total=0, total_nonconsensus=0 ;
	
	for (unsigned int chr=0; chr < genome->nrChromosomes(); chr++)
	{
		if (junctionlist[chr].empty())
			continue;
		
		std::list<Junction>::iterator it=junctionlist[chr].begin(); 

		while (!junctionlist[chr].empty() and it!=junctionlist[chr].end())
		{
			assert((*it).coverage>0);
			if ( ((*it).coverage<min_coverage) || (((*it).coverage < 2*min_coverage) && (!(*it).consensus) ) )
			{
				it=junctionlist[chr].erase(it);
			}
			else
			{
				it++;
				if (!(*it).consensus)
					total_nonconsensus++ ;
			}
		}
		
		total+=junctionlist[chr].size();
		
	}

	pthread_mutex_unlock( &junction_mutex) ;
	
	fprintf(stdout,"Number of junctions in database (min support=%i): %i consensus, %i nonconsensus\n", 
			min_coverage,total-total_nonconsensus,total_nonconsensus);
}
	

void JunctionMap::insert_junction(char strand, int chr, int start, int end, bool consensus, int coverage=1)
{
	pthread_mutex_lock( &junction_mutex) ;

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
		junctionlist[chr].push_back(j);

		pthread_mutex_unlock( &junction_mutex) ;
		return;
		
	}

	std::list<Junction>::iterator it;
	for (it=junctionlist[chr].begin(); it!=junctionlist[chr].end(); it++)
	{
		if (start <  (*it).start)
		{
			j.start=start;
			j.end=end;
			j.coverage=coverage;
			j.strand=strand;
			j.consensus = consensus ;
			junctionlist[chr].insert(it,j);

			pthread_mutex_unlock( &junction_mutex) ;
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
				junctionlist[chr].insert(it,j);

				pthread_mutex_unlock( &junction_mutex) ;
				return;
			}
			
			if (end == (*it).end)
			{
				if (strand == (*it).strand)
				{
					if ((*it).consensus!=consensus)
					{
						fprintf(stderr, "WARNING: consensus mismatch:\n%i-%i %c %i %i\n%i-%i %c %i %i\n", (*it).start, (*it).end, (*it).strand, (*it).coverage, (*it).consensus, start, end, strand, coverage, consensus) ;
						//assert(0) ;

						if (!consensus)
							(*it).consensus=false ;
					}

					(*it).coverage += coverage;

					pthread_mutex_unlock( &junction_mutex) ;
					return;
				}
				if (strand == '+')
				{
					j.start=start;
					j.end=end;
					j.coverage=coverage;
					j.consensus = consensus ;
					j.strand = strand;
					junctionlist[chr].insert(it,j);

					pthread_mutex_unlock( &junction_mutex) ;
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
	junctionlist[chr].push_back(j);
	
	pthread_mutex_unlock( &junction_mutex) ;
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
			int pos_cons=tmp.find("Nonconsensus=");
			if (pos_cons>0)
				pos_cons += strlen("Nonconsensus=") ;
			if (pos_cons>0)
				nonconsensus = atoi(tmp.substr(pos_cons).c_str());
			
			insert_junction(strand,chr_idx,start, end, !nonconsensus, coverage);

			intron_lines++;
		}
		

		if (strcmp(type, "exon")==0)
		{
			exon_lines++ ;
		}
		
	}
	fclose(fd) ;

	fprintf(stdout, "read %i intron lines\n", intron_lines) ;
	
	return 0 ;

}

int JunctionMap::report_to_gff(std::string &gff_fname)
{
	pthread_mutex_lock( &junction_mutex) ;

	int nb_introns=0;
	
	fprintf(stdout, "report splice site junction list in GFF file %s\n", gff_fname.c_str()) ;
	
	FILE * fd=Util::openFile(gff_fname.c_str(), "w") ;
	if (!fd)
		return -1 ;	
	for (unsigned int i=0; i<genome->nrChromosomes(); i++){
		
		const char * chr= genome->get_desc(i);
		std::list<Junction>::iterator it;
		
		for (it=junctionlist[i].begin(); it!=junctionlist[i].end(); it++){			
			fprintf(fd,"%s\tpalmapper\tintron\t%i\t%i\t.\t%c\t.\tID=intron_%i;Confirmed=%i",chr,(*it).start,(*it).end,(*it).strand,nb_introns,(*it).coverage);
			if (!(*it).consensus)
				fprintf(fd,";Nonconsensus=1\n");
			else
				fprintf(fd,"\n");
			nb_introns++;
		}
	
	}
	fclose(fd) ;
	fprintf(stdout, "report %i introns\n", nb_introns) ;	

	pthread_mutex_unlock( &junction_mutex) ;

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
