#!/bin/env python

from sys import exit, argv
import subprocess

import optparse

def main():
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)

    parser.add_option('--site', '-s', dest='site', help='Site where files are located. Can be [t3psi,t2cscs,t2pisa,t2ifca]')

    (opt, args) = parser.parse_args()

    if not opt.site:
        parser.error('No site selected')

    if len(args)!=1:
        parser.error('Wrong number of arguments')

    dir = args[0]
    sites = {}
    sites['t3psi'] =('srm://t3se01.psi.ch:8443/srm/managerv2?SFN=','/pnfs/psi.ch/cms/trivcat/')
    sites['t2cscs']=('srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=','/pnfs/lcg.cscs.ch/cms/trivcat/')
    sites['t2pisa']=('srm://cmsdcache.pi.infn.it:8443/srm/managerv2?SFN=','/pnfs/pi.infn.it/data/cms/')
    sites['t2ifca']=('srm://srm01.ifca.es:8444/srm/managerv2?SFN=','/cms/')


#     if opt.site not in sites:
#         parser.error('site can be either t3psi or t2cscs or t2pisa')

#     if opt.site == 't3psi':
#         srmpath = "srm://t3se01.psi.ch:8443/srm/managerv2?SFN="
#         rootpath = srmpath+'/pnfs/psi.ch/cms/trivcat/'
#     elif opt.site == 't2cscs':
#         srmpath = "srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN="
#         rootpath = srmpath+"/pnfs/lcg.cscs.ch/cms/trivcat/"
#     else:
#         parser.error('site can be either t3psi or t2cscs')

#     srmpath  = sites[opt.site][0]
    rootpath = sites[opt.site][0]+sites[opt.site][1]

    srmPath = rootpath+dir
    print "srmls "+srmPath
    srmls = subprocess.Popen(['srmls',srmPath])
    srmls.wait()

if __name__ == '__main__':
    main()
