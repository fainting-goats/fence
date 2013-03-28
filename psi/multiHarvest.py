#!/usr/bin/env python

import os
import sys
import re
import subprocess
import optparse

def checkRes( files ):
    # regex: basename_job_retry_randoms.root
    #           g1     g2   g3    g4
    gridRe = re.compile('(.*)_([0-9]{1,})_([0-9]{1,})_([0-9a-zA-Z]{3}\.root)')
    
    #group by id
    res = {}
    firstId = 1

    for f in files:
        print f
        m = gridRe.search(f)
        if m is None:
            continue

        job = int(m.group(2))
        retry = int(m.group(3))
        if job not in res.iterkeys():
            res[job] = [retry]
        else:
            res[job].append(retry)

    if len(res.keys()) == 0:
        print 'No files found'
        return False

    lastId = 0
    lastId = sorted(res.keys())[-1]
    print 'Highest id found',lastId

    missingJobs = set(range(1,lastId+1))-set(res.keys())
    missingJobs = sorted(missingJobs)
    
    if len(missingJobs) == 0:
        missStr = 'None'
    else:
        missStr='['
        for i in missingJobs:
            missStr += str(i)+','
        missStr = missStr[:-1]+']'

    print 'Missing files:',missStr
    
    isOk = True
    for (job,retries) in res.iteritems():
        if len(retries) != 1:
            print 'Warning: job',job,'has duplicates:',retries
            isOk &= False

    return isOk

def multiHarvest( basepath ):
    datasets =  os.listdir(basepath)
    targetPath = 'ntuples'
    try:
        os.mkdir(targetPath)
    except:
        print targetPath,'already exists'


    for d in datasets:
        print 'Harvesting results in',d
        path = basepath+'/'+d+'/res/' 
        print path
        files = os.listdir(path)
        expr = re.compile(d+'.*\.root')
        target = 'ntuples/'+d+'.root'
        sources = []
        names = []
        for f in files:
            if expr.match(f) is not None:
                names.append(f)
                sources.append(path+f)

        if not checkRes(names):
            continue

#         sys.exit(0)
        hline = '-'*80

        cmd = ['hadd','-f',target]
        cmd.extend(sources)
        print sources
        hadd = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        (stdout, stderr) = hadd.communicate()
        print hline
        print '|  hadd return code:',hadd.returncode
        print hline
        if hadd.returncode is not 0:
            print '|  hadd stdout'
            print hline
            print stdout
            print hline
            print '|  hadd stderr'
            print hline
            print stderr
            print hline




if __name__ == '__main__':
    usage = 'usage: %prog [options] <file.csv>'
    parser = optparse.OptionParser(usage)
    parser.add_option('--site', '-s', dest='site', help='Destination. Can be [t3psi,t2cscs]')

    (opt,args) = parser.parse_args()

    if len(args) == 0:
        parser.error('Path not specified')

    path = args[0]

    multiHarvest( path )



