#include <palmapper/JunctionMap.h> 
#include <string>
#include <list>
#include <map>
#include <palmapper/Genome.h>
#include <palmapper/Config.h>
#include <stdlib.h> 
#include <palmapper/Util.h>
#include <pthread.h>

bool JunctionMap::is_consensus_intron(char strand, int chr, int start, int end)
{
	const Chromosome & chrom = genome->chromosome(chr);
	//fprintf(stdout,"%c-%c %c-%c\n",chrom[start], chrom[start+1],chrom[end-1],chrom[end]);

	if (strand=='+'){
	  assert(start<=end) ;
		bool is_consensus_don = false ;
		for (unsigned int j=0; j < DON_CONSENSUS.size(); j++){
			
			if (chrom[start] == DON_CONSENSUS[j][0] && chrom[start+1] == DON_CONSENSUS[j][1])
			{
				is_consensus_don = true ;
				break ;
			}
		}
		
		if (! is_consensus_don)
			return false;
	
		bool is_consensus_acc = false ;
		for (unsigned int j=0; j < ACC_CONSENSUS.size(); j++)
			if (chrom[end-1] == ACC_CONSENSUS[j][0] && chrom[end] == ACC_CONSENSUS[j][1])
			{
				is_consensus_acc = true ;
				break ;
			}
		//fprintf(stdout, "+ is_consensus_don=%i, is_consensus_acc=%i [%c%c]\n", is_consensus_don, is_consensus_acc, chrom[end-1], chrom[end]) ;
		return is_consensus_acc;
	}
	else{
		bool is_consensus_acc = false ;
		for (unsigned int j=0; j < ACC_CONSENSUS_REV.size(); j++)
			if (chrom[start] == ACC_CONSENSUS_REV[j][0] && chrom[start+1] == ACC_CONSENSUS_REV[j][1])
			{
				is_consensus_acc = true ;
				break ;
			}
	
		if (! is_consensus_acc)
			return false;
	
		bool is_consensus_don = false ;
		for (unsigned int j=0; j < DON_CONSENSUS_REV.size(); j++)
			if (chrom[end-1] == DON_CONSENSUS_REV[j][0] && chrom[end] == DON_CONSENSUS_REV[j][1])
			{
				is_consensus_don = true ;
				break ;
			}
		//fprintf(stdout, "- is_consensus_don=%i, is_consensus_acc=%i\n", is_consensus_don, is_consensus_acc) ;
		return is_consensus_don;
	}
	
	
}

bool compare_exons(Exon ex1,Exon ex2)
{
	if (ex1.start < ex2.start)
		return true;
	if (ex1.start > ex2.start)
		return false;
	return (ex1.end < ex2.end);
}


JunctionMap::JunctionMap(Genome const &genome_, int anno_pseudo_coverage_, 
						 std::vector<const char*> ACC_CONSENSUS_, std::vector<const char*> DON_CONSENSUS_, 
						 std::vector<const char*> ACC_CONSENSUS_REV_, std::vector<const char*> DON_CONSENSUS_REV_ )
{
	genome = &genome_ ;
	unsigned int nbchr = genome->nrChromosomes();
	
	junctionlist = new std::deque<Junction>[nbchr];
	
	anno_pseudo_coverage = anno_pseudo_coverage_ ;
	
	ACC_CONSENSUS= ACC_CONSENSUS_;
	ACC_CONSENSUS_REV= ACC_CONSENSUS_REV_;
	DON_CONSENSUS= DON_CONSENSUS_;
	DON_CONSENSUS_REV= DON_CONSENSUS_REV_;
	int ret = pthread_mutex_init(&junction_mutex, NULL) ;// = PTHREAD_MUTEX_INITIALIZER ;
	assert(ret==0) ;
	
}


JunctionMap::~JunctionMap()
{
	for (unsigned int i=0; i< genome->nrChromosomes(); i++)
		junctionlist[i].clear();
	
	delete[] junctionlist;	
}


void JunctionMap::filter_junctions(int min_coverage, int min_junction_qual, int filter_by_map, const GenomeMaps & genomemaps, int verbosity)
{
	lock() ;
	if (verbosity>0)
	  {
	    fprintf(stdout, "Filtering junctions, requiring\n* %i as minimum confirmation count\n", min_coverage) ;
	    if (min_junction_qual>0)
	      fprintf(stdout, "* requiring minimum junction quality of %i (usually distance to border)\n", min_junction_qual) ;
	    if (filter_by_map>=0)
	      fprintf(stdout, "* requiring junction next to mapped read or annotated exon with distance at most %i bp\n", filter_by_map) ;
	  }
	
	int total=0, 
	  used_nonconsensus=0, 
	  used_consensus=0,
	  filtered_consensus=0, 
	  filtered_nonconsensus=0 ;

	int N=0, T=0 ;
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
			assert((*it).coverage>=0);

			bool take = true ;
			if ((*it).junction_qual<min_junction_qual)
				take = false ;
			if ((*it).coverage<min_coverage)
				take = false ;
			if (((*it).coverage < 2*min_coverage || ((*it).junction_qual<30)) && min_junction_qual!=0 && (!(*it).consensus))
				take = false ;

			if (take && filter_by_map>=0)
			{
				bool map=false ;
				for (int p=-filter_by_map; p<=filter_by_map && !map; p++)
				  {
				    if ((*it).start+p>=0 && (*it).start+p<(int)genome->chromosome(chr).length())
				      map |= genomemaps.CHR_MAP(genome->chromosome(chr), (*it).start+p) ;
				    if ((*it).end+p>=0 && (*it).end+p<(int)genome->chromosome(chr).length())
				      map |= genomemaps.CHR_MAP(genome->chromosome(chr), (*it).end+p) ;
				  }
				if (!map)
				  take=false ;
			}
			
			if (!take)
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
		int n=filtered_consensus+filtered_nonconsensus+used_consensus+used_nonconsensus ;
		int t=used_consensus+used_nonconsensus ;
		if (verbosity>0)
		  fprintf(stdout, "%s: analyzed %i junctions, accepted %i junctions (%2.1f%%)\n", genome->chromosome(chr).desc(), n, t, 100.0*t/n) ;
		total+=junctionlist[chr].size();
		N+=n ;
		T+=t ;
	}
	unlock() ;	
	if (verbosity>0)
	  fprintf(stdout, "All: analyzed %i junctions, accepted %i junctions (%2.1f%%)\n", N, T, 100.0*T/N) ;

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

	if (coverage<0) // annotation
		coverage = anno_pseudo_coverage ;
	if (junction_qual<0) // annotation
		junction_qual = anno_pseudo_coverage ;

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
						fprintf(stderr, "ERROR: consensus mismatch:\n%s:\t%i-%i %c %i %i %i\n%s:\t%i-%i %c %i %i %i\n", 
								(*it).read_id.c_str(), (*it).start, (*it).end, (*it).strand, (*it).coverage, (*it).consensus, (*it).junction_qual,
								read_id, start, end, strand, coverage, consensus, junction_qual) ;

						//assert(0) ; // this should not happen -> please report this bug and try commenting out the assertion
						
						if (!consensus) // try to handle this case 
							(*it).consensus=false ;
					}
					if (junction_qual > (*it).junction_qual)
					{
						(*it).junction_qual=junction_qual ;
						(*it).read_id=read_id ;
					}
					
					if ((*it).coverage!=0 && coverage!=0)
						(*it).coverage += coverage;
					else
						(*it).coverage = 0;
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
	int num_introns=0 ;
	std::string prev_parent;
	prev_parent.assign("");
	std::map<std::string,std::list<Exon> > transcript_map;
	int prev_chr=-1;
	
		
	while (!feof(fd))
	{
		char chr_name[1000], source[1000], type[1000], properties[1000], strand, tmp1[1000], tmp2[1000] ;
		int start, end ;

		Util::skip_comment_lines(fd) ;
		
		//Scan gff3 line
		int num = fscanf(fd, "%1000s\t%1000s\t%1000s\t%i\t%i\t%1000s\t%c\t%1000s\t%1000s\n", chr_name, source, type, &start, &end, tmp1, &strand, tmp2, properties) ;  
		if (num!=9)
		{
			if (feof(fd))
				break ;
			fprintf(stdout, "gff line in file %s only contained %i columns, aborting\n", gff_fname.c_str(), num) ;
			//exit(-1) ;
		}
		

		// Line comes from a gff3 file built by PALMapper
		if ((strcmp(source, "palmapper")==0  || strcmp(source, "TopHat")==0) and strcmp(type, "intron")==0)
		{

			int chr_idx = genome->find_desc(chr_name) ;
			if (chr_idx==-1)
			{
				fprintf(stderr, "chromosome %s not found. known chromosome names:\n", chr_name) ;
				genome->print_desc(stderr) ;
				fclose(fd) ;
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
					free(intron_string) ;
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
				free(read_id) ;
				read_id = strdup(tmp.substr(pos_id).c_str()) ;
			}

			// check intron consensus
			bool consensus_intron= is_consensus_intron(strand,chr_idx,start,end);
			nonconsensus=!consensus_intron ;
			
			//Attention: positions in this file start at 0! :S
			insert_junction(strand,chr_idx,start, end, !nonconsensus, intron_string, junction_qual, read_id, coverage);
			num_introns++ ;

			free(intron_string) ;
			free(read_id) ;

			
			intron_lines++;
		}
		
		//Line comes from annotation
		if (strcmp(type, "exon")==0)
		{
			exon_lines++ ;
			
			//Get Parent name
			std::string parent;
			std::string tmp(properties);
			int pos_parent=tmp.find("Parent=");
			if (pos_parent <0){
				fprintf(stderr, "No parent information for exon in gff3 file\n") ;
				fclose(fd) ;
				return -1 ;
			}
			else
				pos_parent+=7;
			
			int pos_parent_end=tmp.find(";",pos_parent);
			if (pos_parent<0)
				parent=tmp.substr(pos_parent);
			else
				parent=tmp.substr(pos_parent,pos_parent_end-pos_parent);

			//if (intron_lines%100000==0)
			//	fprintf(stdout, "read %i intron lines\n", intron_lines) ;
		
			//Get chromosome index
			int chr_idx = genome->find_desc(chr_name) ;
			if (chr_idx==-1)
			{
				fprintf(stderr, "chromosome %s not found. known chromosome names:\n", chr_name) ;
				genome->print_desc(stderr) ;
				fclose(fd) ;
				return -1 ;
			}
			

			//Build exon
			Exon ex;
			ex.start=start;
			ex.end=end;
			ex.strand=strand;
			ex.chr=chr_idx;

			if (prev_chr==-1)
				prev_chr=chr_idx;
			else{
				
				//Deal with transcripts of the current chromosome
				if (prev_chr!=chr_idx){
					//fprintf(stdout,"prev chr=%i curr chr=%i\n",prev_chr, chr_idx);	
					for (std::map<std::string,std::list<Exon> >::iterator it=transcript_map.begin();it!=transcript_map.end();it++){
						
						std::list<Exon> exons_list = (*it).second;
						exons_list.sort(compare_exons);
						std::list<Exon>::iterator it_prev=exons_list.begin();
						std::list<Exon>::iterator it_next=exons_list.begin();
						it_next++;

						//fprintf(stdout,"%s\n",(*it).first.c_str());	
						while( it_next!=exons_list.end() ){
							
							//fprintf(stdout,"Exon1: %i-%i %c %i\n", (*it_prev).start,(*it_prev).end,(*it_prev).strand,(*it_prev).chr);
							//fprintf(stdout,"Exon2: %i-%i %c %i\n", (*it_next).start,(*it_next).end,(*it_next).strand,(*it_next).chr);
							if ((*it_prev).strand != (*it_next).strand  || (*it_prev).chr != (*it_next).chr || (*it_prev).end >= (*it_next).start){
								fprintf(stderr, "No consistent information between exons from the same transcript in gff3 file\n") ;
								fclose(fd) ;
								return -1 ;		
							}
							
							bool consensus_intron= is_consensus_intron((*it_next).strand,(*it_next).chr,(*it_prev).end,(*it_next).start-2);

							//In annotation, positions on sequence starts at 1 (coverage set to 0 when from annotation)
							insert_junction((*it_next).strand,(*it_next).chr,(*it_prev).end,(*it_next).start-2,consensus_intron,"", -1, (*it).first.c_str(), -1);
							num_introns++ ;
							it_prev++;
							it_next++;
							
						}
						
					}
					
					transcript_map.clear();
					prev_chr=chr_idx;
				}
			}
			transcript_map[parent].push_back(ex);		
			
		}	
	}

	//Last transcripts
	for (std::map<std::string,std::list<Exon> >::iterator it=transcript_map.begin();it!=transcript_map.end();it++){
		
		std::list<Exon> exons_list = (*it).second;
		exons_list.sort(compare_exons);
		std::list<Exon>::iterator it_prev=exons_list.begin();
		std::list<Exon>::iterator it_next=exons_list.begin();
		it_next++;
		
		//fprintf(stdout,"%s\n",(*it).first.c_str());	
		while( it_next!=exons_list.end() ){
			
			//fprintf(stdout,"Exon1: %i-%i %c %i\n", (*it_prev).start,(*it_prev).end,(*it_prev).strand,(*it_prev).chr);
			//fprintf(stdout,"Exon2: %i-%i %c %i\n", (*it_next).start,(*it_next).end,(*it_next).strand,(*it_next).chr);
			if ((*it_prev).strand != (*it_next).strand  || (*it_prev).chr != (*it_next).chr || (*it_prev).end >= (*it_next).start){
				fprintf(stderr, "No consistent information between exons from the same transcript in gff3 file\n") ;
				fclose(fd) ;
				return -1 ;		
			}
			
			bool consensus_intron= is_consensus_intron((*it_next).strand,(*it_next).chr,(*it_prev).end,(*it_next).start-2);

			
			//In annotation, positions on sequence starts at 1 (coverage set to 0 when from annotation)
			insert_junction((*it_next).strand,(*it_next).chr,(*it_prev).end,(*it_next).start-2,consensus_intron,"", -1, (*it).first.c_str(), -1);
			num_introns++ ;
			it_prev++;
			it_next++;
			
		}
						
	}
					
	transcript_map.clear();

	fclose(fd) ;

	fprintf(stdout, "read %i GFF intron lines\n", intron_lines) ;
	fprintf(stdout, "read %i GFF exon lines\n", exon_lines) ;

	fprintf(stdout, "Read/inferred %i intron junctions\n", num_introns) ;
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
	
	return ret;
	
}
