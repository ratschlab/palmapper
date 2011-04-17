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
	
	junctionlist[0] = new std::list<Junction>[nbchr];
	junctionlist[1] = new std::list<Junction>[nbchr];

}


JunctionMap::~JunctionMap()
{
	delete[] junctionlist[1];
	delete[] junctionlist[0];

	
}


void JunctionMap::insert_junction(char strand, int chr, int start, int end, int coverage=1)
{

	//Sorted list by donor positions first and then acceptor positions
	Junction j;

	//fprintf(stdout,"%c %i %i %i\n",strand, chr, start, end);
	int num_strand=0;
	
	if (strand=='-')
		num_strand=1;
		


	if (junctionlist[num_strand][chr].empty())
	{
		j.start=start;
		j.end=end;
		j.coverage=coverage;
		junctionlist[num_strand][chr].push_back(j);
		return;
		
	}


	std::list<Junction>::iterator it;
	for (it=junctionlist[num_strand][chr].begin(); it!=junctionlist[num_strand][chr].end(); it++){
	
		if (start <  (*it).start){
			j.start=start;
			j.end=end;
			j.coverage=coverage;
			junctionlist[num_strand][chr].insert(it,j);
			return;
		}
	
		if (start==  (*it).start){
			if (end < (*it).end){
				j.start=start;
				j.end=end;
				j.coverage=coverage;
				junctionlist[num_strand][chr].insert(it,j);
				return;
			}
			
			if (end == (*it).end){
				(*it).coverage+=coverage;
				return;
			}
		}
		
		continue;
	}
	
	
	j.start=start;
	j.end=end;
	j.coverage=coverage;
	junctionlist[num_strand][chr].push_back(j);
	
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

	char strand;
	
	for (int i=0; i<2; i++){
		
		if ( i==0)
			strand='+';
		else
			strand='-';

		for (int j=0; j<genome->nrChromosomes(); j++){
			
			const char * chr= genome->get_desc(j);
			std::list<Junction>::iterator it;
			
			for (it=junctionlist[i][j].begin(); it!=junctionlist[i][j].end(); it++){			
				fprintf(fd,"%s\tpalmapper\tintron\t%i\t%i\t.\t%c\t.\tID=intron_%i;Note=%i\n",chr,(*it).start,(*it).end,strand,nb_introns,(*it).coverage);
				nb_introns++;
			}
		}		
	}
	fclose(fd) ;
	fprintf(stdout, "report %i introns\n", nb_introns) ;	
	return 0;
	
}
