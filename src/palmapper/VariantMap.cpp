#include <palmapper/VariantMap.h>
#include <string>
#include <list>
#include <map>
#include <palmapper/Genome.h>
#include <palmapper/Config.h>
#include <stdlib.h> 
#include <palmapper/Util.h>
#include <pthread.h>


VariantMap::VariantMap(Genome const &genome_)
{
	genome = &genome_ ;
	unsigned int nbchr = genome->nrChromosomes();
	
	variantlist = new std::deque<Variant>[nbchr];
	
	int ret = pthread_mutex_init(&variant_mutex, NULL) ;// = PTHREAD_MUTEX_INITIALIZER ;
	assert(ret==0) ;
	
}

VariantMap::~VariantMap()
{
	delete[] variantlist;	
}


void VariantMap::insert_variant(int chr, int pos, int ref_len, int variant_len, const std::string & ref_str, const std::string & variant_str)
{
	lock() ;

	enum polytype pt = pt_unknown ;
	if (ref_len==1 && variant_len==1)
		pt = pt_SNP ;
	else if (ref_len==0 && variant_len>0)
		pt = pt_insertion ;
	else if (variant_len==0 && ref_len>0)
		pt = pt_deletion ;
	else if (variant_len>0 && ref_len>0)
		pt = pt_substitution ;

	Variant j;
	int end = pos+ref_len ;

	if (variantlist[chr].empty())
	{
		j.type = pt ;
		j.position = pos ;
		j.end_position=end ;
		j.ref_len = ref_len ;
		j.variant_len = variant_len ;
		j.ref_str=ref_str ;
		j.variant_str = variant_str ;

		variantlist[chr].push_back(j);

		unlock() ;
		return;
	}

	std::deque<Variant>::iterator it = my_lower_bound(variantlist[chr].begin(), variantlist[chr].end(), pos) ;
	
	for (; it!=variantlist[chr].end(); it++)
	{
		if (pos <  (*it).position)
		{
			j.type = pt ;
			j.position = pos ;
			j.end_position=end ;
			j.ref_len = ref_len ;
			j.variant_len = variant_len ;
			j.ref_str=ref_str ;
			j.variant_str = variant_str ;

			variantlist[chr].insert(it,j);

			unlock() ;
			return;
		}
		if (pos ==  (*it).position)
		{
			if (end <= (*it).end_position)
			{
				j.type = pt ;
				j.position = pos ;
				j.end_position=end ;
				j.ref_len = ref_len ;
				j.variant_len = variant_len ;
				j.ref_str=ref_str ;
				j.variant_str = variant_str ;
				
				if (!variant_identical(j, *it))
					variantlist[chr].insert(it,j);

				unlock() ;
				return;
			}
		}
		continue;
	}

	j.type = pt ;
	j.position = pos ;
	j.end_position=end ;
	j.ref_len = ref_len ;
	j.variant_len = variant_len ;
	j.ref_str=ref_str ;
	j.variant_str = variant_str ;
	
	variantlist[chr].push_back(j);

	unlock() ;
	return ;
}

int VariantMap::init_from_sdi(std::string &sdi_fname)
{

	fprintf(stdout, "initializing genome variant list with SDI file %s\n", sdi_fname.c_str()) ;

	FILE * fd=Util::openFile(sdi_fname.c_str(), "r") ;
	if (!fd)
		return -1 ;
	int variant_lines = 0, variant_lines_checked = 0 ;

	while (!feof(fd))
	{
		char chr_name[1000], ref_str[100000], variant_str[100000], tmp[1000], buf[250000] ;
		int position, lendiff ;
		
		Util::skip_comment_lines(fd) ;
		
		fgets(buf, 250000, fd) ;

		//Scan sdi3 line
		int num = sscanf(buf, "%1000s\t%i\t%i\t%100000s\t%100000s\t%1000s\t%1000s\n", chr_name, &position, &lendiff, ref_str, variant_str, tmp, tmp) ;  
		if (num<5)
		{
			if (feof(fd))
				break ;
			fprintf(stdout, "sdi line only contained %i columns (5 expected), aborting (%s)\n", num, chr_name) ;
		}

		int chr_idx = genome->find_desc(chr_name) ;
		if (chr_idx==-1)
		{
			fprintf(stderr, "chromosome %s not found. known chromosome names:\n", chr_name) ;
			genome->print_desc(stderr) ;
			fclose(fd) ;
			return -1 ;
		}
		if (strcmp(ref_str, "-")==0)
			strcpy(ref_str, "") ;
		if (strcmp(variant_str, "-")==0)
			strcpy(variant_str, "") ;

		int ref_len = strlen(ref_str) ;
		int variant_len = strlen(variant_str) ;
		assert(lendiff==variant_len-ref_len) ;

		// validate variants on genome sequence
		if (ref_len>0)
		{
			//fprintf(stdout, "pos=%i\tref_len=%i\tvariant_len=%i\n", pos, ref_len, variant_len) ;
			for (int i=0; i<ref_len; i++)
			{
				if (genome->chromosome(chr_idx)[position+i-1]!=ref_str[i] && 
					ref_str[i]!='N' && ref_str[i]!='Y' && ref_str[i]!='W' && ref_str[i]!='K' && ref_str[i]!='S' && ref_str[i]!='M' && ref_str[i]!='R' && ref_str[i]!='D')
				{
					fprintf(stdout, "error: variant map disagrees with genome: %i\t%i\t%c\t%c\t%s\n", i, position+i, genome->chromosome(chr_idx)[position+i-1], ref_str[i], ref_str) ;
					exit(-1) ;
				}
			}
			variant_lines_checked++ ;
		}
			
		insert_variant(chr_idx, position-1, ref_len, variant_len, ref_str, variant_str);
		variant_lines++ ;
	}		

	fclose(fd) ;

	fprintf(stdout, "read %i variant lines (checked %i)\n", variant_lines, variant_lines_checked) ;

	return 0 ;

}

int VariantMap::report_to_sdi(std::string &sdi_fname)
{
	lock() ;

	int nb_variants=0;
	
	fprintf(stdout, "report genome list in SDI file %s\n", sdi_fname.c_str()) ;
	
	FILE * fd=Util::openFile(sdi_fname.c_str(), "w") ;
	if (!fd)
		return -1 ;	

	for (unsigned int i=0; i<genome->nrChromosomes(); i++){
		
		const char * chr= genome->get_desc(i);
		std::deque<Variant>::iterator it;
		
		for (it=variantlist[i].begin(); it!=variantlist[i].end(); it++)
		{			
			fprintf(fd,"%s\t%i\t%i\t%s\t%s\n",
					chr,(*it).position,(*it).variant_len-(*it).ref_len,(*it).ref_str.c_str(), (*it).variant_str.c_str());
			nb_variants++;
		}
	}
	fclose(fd) ;
	fprintf(stdout, "reported %i variants\n", nb_variants) ;	

	unlock() ;

	return 0;
	
}

int VariantMap::init_from_sdis(std::string &sdi_fname)
{

	int previousfound=0;
	int found=sdi_fname.find(",");
	std::string filename;
	
	while (found >= 0)
	{
		
		filename = sdi_fname.substr(previousfound, found-previousfound);
		int ret = init_from_sdi(filename);
		if (ret!=0)
			return ret;
	   
		previousfound=found+1;
		found=sdi_fname.find(",",found+1);
	}
	
	filename=sdi_fname.substr(previousfound);
	int ret=init_from_sdi(filename);
	if (ret!=0)
		return  ret;
	
	return ret;
	
}
