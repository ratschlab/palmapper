#! /usr/bin/python

"""
Runs Palmapper on single-end or paired-end data.
"""

import optparse, os, sys, tempfile, shutil

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()
 
def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option('', '--threads', dest='threads', help='The number of threads to run')
    parser.add_option('', '--input1', dest='input1', help='The (forward or single-end) reads file in Sanger FASTQ format')
    parser.add_option('', '--input2', dest='input2', help='The reverse reads file in Sanger FASTQ format')
    parser.add_option('', '--unspliced-output', dest='unspliced_output', help='The output file for unspliced reads')
    parser.add_option('', '--spliced-output', dest='spliced_output', help='The output file for spliced reads')
    parser.add_option('', '--paired', dest='paired', help='Whether the data is single- or paired-end')
    parser.add_option('', '--genomeSource', dest='genomeSource', help='The type of reference provided')
    parser.add_option('', '--ref', dest='ref', help='The reference genome to use or index')
    parser.add_option('', '--indexSettings', dest='index_settings', help='Whether or not indexing options are to be set')
    parser.add_option('', '--params', dest='params', help='Whether to use default or specified parameters for GenomeMapper')
    parser.add_option('', '--seedlength', dest='seedlength', help='Index Seed Length')
    parser.add_option('', '--alignseedlength', dest='alignseedlength', help='Alignment Seed Length')
    parser.add_option('', '--format', dest='format', help='Output format (bed or shore)')
    parser.add_option('', '--maxmismatches', dest='maxmismatches', help='Maximal number of mismatches')
    parser.add_option('', '--maxgaps', dest='maxgaps', help='Maximal number of gaps')
    parser.add_option('', '--maxedits', dest='maxedits', help='Maximal number of edit operations')
    parser.add_option('', '--reportall', dest='reportall', help='Report all hits')
    parser.add_option('', '--acc', dest='acc', help='Acceptor SS predictions')
    parser.add_option('', '--don', dest='don', help='Donor SS predictions')
    parser.add_option('', '--qpalma', dest='qpalma', help='QPALMA parameter file')
    parser.add_option('', '--qpalma-params', dest='qpalma_params', help='Whether to use default or specified parameters for QPALMA')
    parser.add_option('', '--aligntrigger', dest='aligntrigger', help='aligntrigger')
    parser.add_option('', '--alignshorttrigger', dest='alignshorttrigger', help='alignshorttrigger')
    parser.add_option('', '--aligncombinedtrigger', dest='aligncombinedtrigger', help='aligncombinedtrigger')
    parser.add_option('', '--maxintronlength', dest='maxintronlength', help='maxintronlength')
    parser.add_option('', '--maxintronnum', dest='maxintronnum', help='maxintronnum')
    parser.add_option('', '--clustertol', dest='clustertol', help='clustertol')
    parser.add_option('', '--maxalign', dest='maxalign', help='maxalign')
    #parser.add_option('', '--filter_ss_tophit', dest='filter_ss_tophit', help='filter_ss_tophit')
    parser.add_option('', '--report_ss', dest='report_ss', help='report_ss')
    parser.add_option('', '--qpalma-use-map-max-len', dest='mapmaxlen', help='qpalma-use-map-max-len')
    parser.add_option('', '--reportmappedread', dest='reportmappedread', help='reportmappedread')
    parser.add_option('', '--reportsplicedread', dest='reportsplicedread', help='reportsplicedread')
    (options, args) = parser.parse_args()
    
    # index if necessary
    if options.genomeSource == 'history':
        # set up commands
        if options.index_settings =='index_pre_set':
            indexing_cmds = ''
        else:
            try:
                indexing_cmds = '%s ' % \
                                (('','-s %s'%options.seedlength)[options.seedlength!='None' and options.seedlength>=1])
            except ValueError:
                indexing_cmds = ''
                
        # make temp directory for placement of indices and copy reference file there
        tmp_dir = tempfile.gettempdir()
        try:
            os.system('cp %s %s' % (options.ref, tmp_dir))
        except Exception, erf:
            stop_err('Error creating temp directory for indexing purposes\n' + str(erf))
        options.ref = os.path.join(tmp_dir, os.path.split(options.ref)[1])
        cmd1 = 'pmindex -v -i %s %s' % (options.ref, indexing_cmds)
        try:
            os.system(cmd1)
        except Exception, erf:
            stop_err('Error indexing reference sequence\n' + str(erf))
    
    if options.params == 'pre_set':
        aligning_cmds = ' -M 3 -G 1 -E 3 -l 18 '
    else:
        try:
            aligning_cmds = '%s %s %s %s %s %s' % \
                            (('','-f %s' % options.format)[options.format!='None'],
                             ('','-a')[options.reportall!='None'],
                             ('','-M %s' % options.maxmismatches)[options.maxmismatches!='None'],
                             ('','-G %s' % options.maxgaps)[options.maxgaps!='None'],
                             ('','-E %s' % options.maxedits)[options.maxedits!='None'],
                             ('','-l %s' % options.alignseedlength)[options.alignseedlength!='None'])
        except ValueError, erf:
            stop_err('Something is wrong with the alignment parameters and the alignment could not be run\n' + str(erf))

    aligning_cmds+=('','-report-map-read')[options.reportmappedread!='None' and options.reportmappedread!='false' and options.reportmappedread!='False']

    if options.qpalma_params == 'pre_set':
        qpalma_cmds = '-L 18 -K 13 -C 24 -I 50000 -NI 2 -SA 10 -CT 10 -z -S -seed-hit-cancel-threshold 10000  -report-map-read -report-spliced-read  -report-splice-sites 0.9 -report-map-region -qpalma-use-map'
    else:
        try:
            print options
            qpalma_cmds = '%s %s %s %s %s %s %s %s %s %s -S -qpalma-use-map ' % \
                          (('','-L %s' % options.aligntrigger)[options.aligntrigger!='None'],
                           ('','-K %s' % options.alignshorttrigger)[options.alignshorttrigger!='None'],
                           ('','-C %s' % options.aligncombinedtrigger)[options.aligncombinedtrigger!='None'],
                           ('','-I %s' % options.maxintronlength)[options.maxintronlength!='None'],
                           ('','-NI %s' % options.maxintronnum)[options.maxintronnum!='None'],
                           ('','-CT %s' % options.clustertol)[options.clustertol!='None'],
                           ('','-SA %s' % options.maxalign)[options.maxalign!='None'],
                           ('','-report-splice-sites %s' % options.report_ss)[options.report_ss!='None'],
                           ('','-qpalma-use-map-max-len %s' % options.mapmaxlen)[options.mapmaxlen!='None'],
                           ('','-report-spliced-read')[options.reportsplicedread!='None' and options.reportsplicedread!='false' and options.reportsplicedread!='False'])
        except ValueError, erf:
            stop_err('Something is wrong with the QPALMA alignment parameters and the alignment could not be run\n' + str(erf))

    (unmapped_tmp_file, unmapped_tmp_fname) = tempfile.mkstemp(suffix='', prefix='unmapped_', dir=None) ;
    os.close(unmapped_tmp_file) ;

    # creating copy of critical files on local tmp file system
    index_tmp_dname = tempfile.mkdtemp(suffix='', prefix='gmindex_tmp_', dir=None) ;
    if options.ref[0:5]!='/tmp/':
        try:
            print 'cp %s* %s' % (options.ref, index_tmp_dname)
            os.system('cp %s* %s' % (options.ref, index_tmp_dname))
        except Exception, erf:
            stop_err('Error creating temp directory for indexing purposes\n' + str(erf))
        options.ref = os.path.join(index_tmp_dname, os.path.split(options.ref)[1])

    acc_tmp_dname = tempfile.mkdtemp(suffix='', prefix='acc_', dir=None) ;
    try:
        print 'cp %s/pred/contig_* %s' % (options.acc, acc_tmp_dname)
        os.system('cp %s/pred/contig_* %s' % (options.acc, acc_tmp_dname))
    except Exception, erf:
        stop_err('Error creating temp directory for indexing purposes\n' + str(erf))
    options.acc = os.path.join(acc_tmp_dname, 'contig_%i%c') 

    don_tmp_dname = tempfile.mkdtemp(suffix='', prefix='don_', dir=None) ;
    try:
        print 'cp %s/pred/contig_* %s' % (options.don, don_tmp_dname)
        os.system('cp %s/pred/contig_* %s' % (options.don, don_tmp_dname))
    except Exception, erf:
        stop_err('Error creating temp directory for indexing purposes\n' + str(erf))
    options.don = os.path.join(don_tmp_dname, 'contig_%i%c') 
    
    # prepare actual aligning commands
    (report_file, report_fname) = tempfile.mkstemp(suffix='', prefix='report_', dir=None) 
    os.close(report_file)
    try:
        os.unlink(report_fname)
    except:
        pass
    
    if options.paired == 'paired':
        print "Sorry, paired end alignments are not supported yet"
        return
    else:
        cmd2 = 'genomemapper -z %s -i %s -q %s -o %s -u %s -report %s' % (aligning_cmds, options.ref, options.input1, options.unspliced_output, unmapped_tmp_fname, report_fname) 
	cmd2a = 'genomemapper -z %s %s -i %s -q %s -o %s -H %s -u %s -qpalma %s -acc %s -don %s -report %s' % (aligning_cmds, qpalma_cmds, options.ref, options.input1, options.unspliced_output, options.spliced_output, unmapped_tmp_fname, options.qpalma, options.acc, options.don, report_fname) 

    # align
    try:
        print cmd2
        #os.environ['LD_LIBRARY_PATH']=''
        os.system(cmd2)
        print cmd2a
        os.system(cmd2a)
    except Exception, erf:
        stop_err("Error aligning sequence\n" + str(erf))

    try:
        shutil.rmtree(index_tmp_dname)
        shutil.rmtree(acc_tmp_dname)
        shutil.rmtree(don_tmp_dname)
        os.unlink(report_fname)
    except:
        pass

if __name__=="__main__": __main__()
