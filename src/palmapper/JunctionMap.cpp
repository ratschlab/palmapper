#include <palmapper/JunctionMap.h>
#include <string>
#include <list>
#include <palmapper/Genome.h>
#include <stdlib.h> 
#include <palmapper/Util.h>

JunctionMap::JunctionMap(Genome const &genome_, int min_coverage_)
{
	genome = &genome_ ;
	unsigned int nbchr = genome->nrChromosomes();
	
	junctionlist = new std::list<Junction>[nbchr];
	
	min_coverage=min_coverage_;
	
}


JunctionMap::~JunctionMap()
{
	delete[] junctionlist;	
}


void JunctionMap::filter_junctions()
{

	int total=0;
	
	for (unsigned int chr=0; chr < genome->nrChromosomes(); chr++){

		if (junctionlist[chr].empty())
			continue;
		
		std::list<Junction>::iterator it=junctionlist[chr].begin(); 

		while(!junctionlist[chr].empty() and it!=junctionlist[chr].end()){

			if ((*it).coverage>0 and (*it).coverage<min_coverage){
				it=junctionlist[chr].erase(it);
			}
			else
				it++;
		}
		
		total+=junctionlist[chr].size();
		
	}
	
	fprintf(stdout,"Number of junctions in database (min support=%i): %i\n",min_coverage,total);
	

}
	

void JunctionMap::insert_junction(char strand, int chr, int start, int end, int coverage=1)
{

	//Sorted list by donor positions first and then acceptor positions
	Junction j;

	//fprintf(stdout,"%c %i %i %i\n",strand, chr, start, end);

	if (junctionlist[chr].empty())
	{
		j.start=start;
		j.end=end;
		j.coverage=coverage;
		j.strand=strand;
		junctionlist[chr].push_back(j);
		return;
		
	}


	std::list<Junction>::iterator it;
	for (it=junctionlist[chr].begin(); it!=junctionlist[chr].end(); it++){
		
	  
		if (start <  (*it).start){
			j.start=start;
			j.end=end;
			j.coverage=coverage;
			j.strand=strand;
			junctionlist[chr].insert(it,j);
			return;
		}
	
		if (start==  (*it).start){
			if (end < (*it).end){
				j.start=start;
				j.end=end;
				j.coverage=coverage;
				j.strand=strand;
				junctionlist[chr].insert(it,j);
				return;
			}
			
			if (end == (*it).end){
				if (strand == (*it).strand){
					(*it).coverage+=coverage;
					return;
				}
				if (strand == '+'){
					j.start=start;
					j.end=end;
					j.coverage=coverage;
					j.strand=strand;
					junctionlist[chr].insert(it,j);
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
	junctionlist[chr].push_back(j);
	
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
			
			int pos_cov=tmp.find("Note=");
			
			int coverage;
			if (pos_cov>0){
				coverage= atoi(tmp.substr(pos_cov+5).c_str());
			}
			else
				coverage=1;
			
			insert_junction(strand,chr_idx,start, end, coverage);
			
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

	int nb_introns=0;
	
	fprintf(stdout, "report splice site junction list in GFF file %s\n", gff_fname.c_str()) ;
	
	FILE * fd=Util::openFile(gff_fname.c_str(), "w") ;
	if (!fd)
		return -1 ;	
	for (unsigned int i=0; i<genome->nrChromosomes(); i++){
		
		const char * chr= genome->get_desc(i);
		std::list<Junction>::iterator it;
		
		for (it=junctionlist[i].begin(); it!=junctionlist[i].end(); it++){			
			fprintf(fd,"%s\tpalmapper\tintron\t%i\t%i\t.\t%c\t.\tID=intron_%i;Note=%i\n",chr,(*it).start,(*it).end,(*it).strand,nb_introns,(*it).coverage);
				nb_introns++;
		}
	
	}
	fclose(fd) ;
	fprintf(stdout, "report %i introns\n", nb_introns) ;	
	return 0;
	
}

int JunctionMap::init_from_gffs(std::string &gff_fname)
{

	int previousfound=0;
	int found=gff_fname.find(",");
	std::string filename;
	
	while (found >= 0){
		
		filename=gff_fname.substr(previousfound,found-previousfound);
		int ret=init_from_gff(filename);
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
