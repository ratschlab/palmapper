#include <palmapper/JunctionMap.h>
#include <string>
#include <list>
#include <palmapper/Genome.h>
#include <stdlib.h> 
#include <palmapper/Util.h>

JunctionMap::JunctionMap(Genome const &genome_)
{
	genome = &genome_ ;
	unsigned int nbchr = genome->nrChromosomes();
	
	junctionlist = new std::list<Junction>[nbchr];

}


JunctionMap::~JunctionMap()
{
	delete[] junctionlist;	
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
			insert_junction(strand,chr_idx,start, end );
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
	for (int i=0; i<genome->nrChromosomes(); i++){
		
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
