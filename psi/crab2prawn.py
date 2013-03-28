#!/usr/bin/env python

import psitools
import csv
import optparse
import subprocess
import getpass

def c2p( args ):
    # need
    # csv file to get the dataset names
    # id range
    # njobs
    # destination path
    # cmssw dir? can we guess it?
    # template
    # createSession.py -s $nick -g $version -l $originaldataset -i $inFile -o $oFile -w $cmsswdir [-t analyze/Analysis] -j 10 $template
    
    req = ['ids','workingDir','outputDir','nJobs','template','csvFile']
    opt = psitools.ReqOptions(args,req)

    ids = psitools.strToNumbers(opt.ids)
    print ' - Ids considered for inclusion: ', opt.ids
    datasets =  psitools.getDatasetsFromCSV( opt.csvFile, ids )
    print [d.nick for d in datasets]

    for d in datasets:
        shortName = opt.prefix+'id'+d.id+'.'+d.nick
        inputPath = shortName+'.input' 
        inputFile = open(inputPath,'w')
        allFiles =  psitools.dbsGetFileList( d )
        inputFile.write('\n'.join(allFiles))
        inputFile.close()

        groups = d.ver+':'+opt.groups if opt.groups is not None else d.ver 
        cmd = ['pwn_CreateSession.py','-n',
               '-s',shortName,
               '-g',groups,
               '-l',d.name,
               '-i',inputPath,
               '-w',opt.workingDir,
               '-o',opt.outputDir,
               '-j',str(opt.nJobs),
               opt.template
              ]
        if opt.optArgs is not None:
            cmd.extend(['-a',opt.optArgs])
        print ' - Executing',' '.join(cmd)
        create = subprocess.Popen(cmd,stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (stdout,stderr) = create.communicate()
        hline = '-'*80
        print create.poll()
        if create.poll() != 0:
            print hline
            print stdout
            print hline
            print stderr
            print hline
            raise RuntimeError('Something went wrong while creating %s' % d.nick)
        print stdout




if __name__ == '__main__':
    usage = 'usage: %prog -p [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--ids', dest='ids',  help='Dataset ids to be grabbed from the spreadsheet', default='')
    parser.add_option('-w', '--workingDir', dest='workingDir', help='Working directory', default='.')
    parser.add_option('-o', '--outputDir', dest='outputDir', help='Output directory', default='.')
    parser.add_option('-g', '--groups', dest='groups', help='Column separated list of groups') 
    parser.add_option('-j', '--nJobs', type='int', dest='nJobs', help='Number of jobs')
    parser.add_option('-a', '--optArgs', dest='optArgs', help='optional arguments')
    parser.add_option('-t', '--template', dest='template', help='Script template')
    parser.add_option('-c', '--csv', dest='csvFile', help='CSV source file')
    parser.add_option('-p', '--prefix', dest='prefix', help='Prefix to the job name', default='')
    
    (opt,args) = parser.parse_args()
    

    try: 
        c2p( opt.__dict__ )
    except RuntimeError as rt:
        print 'Runtime error!',rt
