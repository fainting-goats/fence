#!/usr/bin/env python

import optparse
import csv
import subprocess
import re
import os
import threading
import time
import sys
import psitools

sema = threading.Semaphore(5)
DBS_SERVER = {"ph01":"http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_01/servlet/DBSServlet",
              "ph02":"http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_02/servlet/DBSServlet"}
#                      'http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_02/servlet/DBSServlet',
class Chdir:         
    def __init__( self, newPath ):  
        self.savedPath = os.getcwd()
        os.chdir(newPath)

    def __del__( self ):
        os.chdir( self.savedPath )

class Dataset:
    def __init__(self, id, nick, ver, name):
        self.id = id
        self.nick = nick
        self.ver  = ver
        self.name = name
        self.transfer = False

#     def __str__(self):
#         pass
#         return '[%(id)s] %(nick)s %(ver)s %(name)d' % self.__dict__

class ThreadTransfer( threading.Thread ):
    def __init__(self, ds, dbs, site, blacklist, whitelist):
        threading.Thread.__init__(self)
        self.dataset     = ds
        self.dbs         = dbs
        self.site        = site
        self.go          = False
        self.kill        = False
        self.code        = None
        self.blacklist   = blacklist
        self.whitelist   = whitelist
        self.t1          = 0.
        self.t2          = 0.
        self.tLastChange = 0.

    def run(self):
        sema.acquire()
        if self.kill:
            print ' * Transfer',self.dataset.nick,'aborted'
            sema.release()
            return
        
        self.go = True
        self.t1 = time.time()
        self.t2 = self.t1
        print ' - Starting transfer',self.dataset.nick
        wd = os.getcwd()+'/'+self.dataset.id
        try:
            if not os.path.exists(wd):
                os.mkdir(wd)
            cmd = ['dbs_transferRegister.py','--dbs=%s' % self.dbs,'--debug','--to-site=%s' % self.site,self.dataset.name]
            if self.blacklist is not None:
                cmd.append('--blacklist=%s' % self.blacklist) 
            if self.whitelist is not None:
                cmd.append('--whitelist=%s' % self.whitelist) 
            print '- Executing:',' '.join(cmd)
            logFile = open(wd+'/grab.log','w')
            logFile.write('#'*80+'\n')
            logFile.write(' '.join(cmd)+'\n')
            logFile.write('#'*80+'\n')
            dbsTransfer = subprocess.Popen(cmd, cwd=wd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            retcode = None
            self.tLastChange = time.time()
            while( not retcode ): 
                retcode = dbsTransfer.poll()
#                 print 'retcode',retcode
#                 print 'kill',self.kill
                if retcode is not None:
                    self.code = retcode
                    break
                if self.kill:
                    try:
                        dbsTransfer.kill()
                    except OSError:
                        pass
                    break
                else:
                    time.sleep(1.)

                out = psitools.readAllSoFar(dbsTransfer)
                if len(out) != 0:
                    logFile.write(out)
                    logFile.flush()
                    os.fsync(logFile)
                    self.tLastChange = time.time()

#                 if time.time()-self.tLastChange > 600: #10 mins?
#                     print ' ? Warning: no output since 600 seconds'




        finally:
            self.t2 = time.time()
            print ' + Transfer',self.dataset.nick,'completed' if not self.kill else 'killed'
            sema.release()


class Grabber:
    def __init__(self, csvFile, destination, ids):
        self.versions = set()
        self.datasets = []
        self.existingDatasets = []
        self.ids = ids
        self.csvFile = csvFile
        self.destination = destination
        self.threads = []
        self.blacklist = None
        self.whitelist = None
        self.dbs = None
        
    def configure(self, version=None, forceTransfer=False):
        self.datasets = psitools.getDatasetsFromCSV( self.csvFile, self.ids )
#         print '\n'.join([ d.name for d in self.datasets])
        for d in self.datasets:
            if d.ver != '':
                self.versions.add(d.ver)

#         print self.versions
        if not forceTransfer:
            for v in self.versions:
                self.existingDatasets.extend( getExistingDatasets(self.dbs, v, self.destination) )

#         print  self.existingDatasets

        for  ds in self.datasets:
            if ds.name in self.existingDatasets:
                continue
            numId = int(re.search('[0-9]+',ds.id).group(0))
            if not numId in self.ids:
                continue

            print 'Registering',ds.id, ds.nick,'for transfer'
            ds.transfer = True
        print '---'

    def start(self):
        for ds in self.datasets:
            if not ds.transfer:
                continue

            t = ThreadTransfer(ds,self.dbs,self.destination, self.blacklist, self.whitelist)
            self.threads.append(t)
            if threading.activeCount() > 1:
                pause = 1.
                print ' - Wait',pause,'sec before queuing',ds.id,ds.nick
                time.sleep(pause)

            t.start()
            print ' -',ds.nick,'Ready.'

        print ' - All datasets ready or running'
        print '---'

    def monitor(self):
        counter = 0
        reportEvery = 600
        while threading.activeCount() > 1:
            time.sleep(.1)
            if (counter % reportEvery) == 0:
                self.report()
            counter += 1
        # report one last time
        self.report()

    def kill(self):
        for t in self.threads:
            t.kill = True

        while threading.activeCount() > 1:
            time.sleep(.1)

    def report(self):
        timeFmt = '%a, %d %b %Y %H:%M:%S' 
        print '- report ('+time.strftime( timeFmt, time.gmtime())+')---'
        lmax = max([len(t.dataset.nick) for t in self.threads])
        for t in self.threads:
            if not t.go:
                statusstr = 'Not started'
            elif t.code is None:
                statusstr = 'Running'
            elif t.code == 0:
                statusstr = 'Finished successfully'
            else:
                statusstr = 'Terminated with error ('+str(t.code)+')'
            timeStr = time.strftime(timeFmt,time.gmtime(t.tLastChange)) 
            dt = '%.2f' % (time.time()-t.tLastChange,)
            print '| {0:^5} {1}'.format(t.dataset.id,t.dataset.nick.ljust(lmax)),
            if t.go:
                print statusstr,'(last change ',timeStr,dt,'sec ago)'
            else:
                print statusstr


    def summary(self):
        for t in self.threads:
            if t.kill:
                exitstr = 'Killed'
            elif not t.go:
                exitstr = 'Not started'
            elif t.code==0:
                exitstr = 'OK'
            else:
                exitstr = 'Error('+str(t.code)+')'

            timeStr = time.strftime('%a, %d %b %Y %H:%M:%S',time.gmtime(t.t1))
            print t.dataset.id,t.dataset.nick,t.site,timeStr,'[',t.t2-t.t1,']',exitstr


def getExistingDatasets(dbs, version, site):

    cmd=['dbs', 
         'search',
         '--url=%s' % DBS_SERVER[dbs],
         '--query=find dataset where dataset like /*%s* and site like %s' % (version,site)]
    print ' '.join(cmd)
    dbsQuery = subprocess.Popen(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout,stderr) = dbsQuery.communicate()
    outRows = stdout.splitlines()

    # no dataset found
    if len(outRows) == 1:
        return []

    for i in range(4):
        outRows.pop(0)

    print '---'
    print '\n'.join(outRows)
    print '---'

    return outRows
    


if __name__ == '__main__':
    usage = 'usage: %prog [options] <file.csv>'
    parser = optparse.OptionParser(usage)
    parser.add_option('--site', '-s', dest='site', help='Destination. Can be [t3psi,t2cscs]')
    parser.add_option('--dbs',  '-d', dest='dbs',  help='DBS database to use', default='ph01')
    parser.add_option('--ids',  '-i', dest='ids',  help='Dataset ids to be grabbed from the spreadsheet', default='')
    parser.add_option('--blacklist',  dest='blacklist', help='Blacklisted sites, comma separated')
    parser.add_option('--whitelist',  dest='whitelist', help='Whitelisted sites, comma separated')
    parser.add_option('--version',    dest='skimVersion', help='Grab samples of this version')
    parser.add_option('--forceTransfer', dest='forceTransfer', action='store_true', help='don\'t check if the dataset is already available at destination')

    (opt,args) = parser.parse_args()

    if len(args) == 0:
        parser.error('Input csv file not defined')

    csvFile = args[0]
    if opt.site == 't3psi' or opt.site == 'T3_CH_PSI':
        site = 'T3_CH_PSI'
    elif opt.site == 't2cscs' or opt.site == 'T2_CH_CSCS':
        site = 'T2_CH_CSCS'
    elif opt.site == 't2ifca' or opt.site == 'T2_ES_IFCA':
        site = 'T2_ES_IFCA'
    else:
        parser.error('site can be either t3psi or t2cscs')

    ids = psitools.strToNumbers(opt.ids)
    print ' - Ids considered for transfer: ', opt.ids
#     print '   ',','.join([str(id) for id in ids])
#     print ','.join(ids)

    try:
        g = Grabber(csvFile, site, ids)
        g.blacklist = opt.blacklist
        g.whitelist = opt.whitelist
        g.dbs = opt.dbs
        g.configure( opt.skimVersion, opt.forceTransfer )
        g.start()
        g.monitor()
    except ValueError as v:
        print v
        g.kill()
    except KeyboardInterrupt:
        print 'Ctrl-C detected. Terminate them all'
        g.kill()
    finally:
        print 'The End'
        print '--- Transfer Summary'
        g.summary()

