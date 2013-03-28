#!/bin/env python

import os
import sys
import subprocess

import optparse
import string

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
    It must be "yes" (the default), "no" or None (meaning
    an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes":True,   "y":True,  "ye":True,
             "no":False,     "n":False}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while 1:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid.keys():
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                "(or 'y' or 'n').\n")

#-------------------------------------------------------------------------------
def srmList( siteRoot, srmPath ):

    srmlscmd =  ['srmls',siteRoot+srmPath]
#     print ' '.join(srmlscmd)
    srmls = subprocess.Popen(srmlscmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
#     srmls.wait()
    (out,err) = srmls.communicate()
    print out
    if len(out) is 0:
        return None

    outLines = out.splitlines()
    header = outLines.pop(0).split()[1]

    entries = []
    for line in outLines:
        if len(line) is 0:
            continue
        tokens = line.split()
        if tokens[0] == '0' and tokens[1][-1] != '/':
            tokens[1] = tokens[1]+'/'
        entries.append(tokens[1])
    return (header, entries)

#-------------------------------------------------------------------------------
def sortByType( entries ):
    dirs = []
    files = []
    for entry in entries:
        if entry[-1] == '/':
            dirs.append(entry)
        else:
            files.append(entry)

    return (dirs,files)
    
#-------------------------------------------------------------------------------
def recursiveList(siteRoot, relPath ):
    allFiles = []
    allDirs = []
    (header,entries) = srmList(siteRoot, relPath)
#     print entries


    (dirs,files) = sortByType(entries)
    if len(dirs) is not 0:
        print '-- directories\n',string.join(dirs,'\n'),'\n'
    if len(files) is not 0:
        print '-- files\n',string.join(files,'\n'),'\n'

    for d in dirs:
        print '- Listing',d
        (subDirs,subFiles) = recursiveList(siteRoot, d )
        allDirs.extend(subDirs)
        allFiles.extend(subFiles)

    # add the files and directories at the end
    allFiles.extend(files)
    allDirs.append(relPath)
#     allDirs.extend(dirs)

    return (allDirs,allFiles)

#-------------------------------------------------------------------------------
def performDelete(siteRoot,dirs,files):
    for file in files:
        print 'Removing',file
        srmrm = subprocess.Popen(['srmrm','-2',siteRoot+file],stdout=subprocess.PIPE)
#         srmrm.wait()
        (out,err) = srmrm.communicate()
        print out

    for dir in dirs:
        print 'Removing',dir
        srmrmdir = subprocess.Popen(['srmrmdir',siteRoot+dir],stdout=subprocess.PIPE)
#         srmrmdir.wait()
        (out,err) = srmrmdir.communicate()
        print out

#-------------------------------------------------------------------------------
def main():
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)

    parser.add_option('-s','--site', dest='site', help='Site where files are located. Can be [t3psi,t2cscs]')

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

#     srmpath  = sites[opt.site][0]
    siteRoot = sites[opt.site][0]
    rootpath = sites[opt.site][1]

    srmPath = rootpath+dir

#     if opt.site == 't3psi':
#         siteRoot = 'srm://t3se01.psi.ch:8443/srm/managerv2?SFN='
#         srmPath = '/pnfs/psi.ch/cms/trivcat'+dir
#     elif opt.site == 't2cscs':
#         siteRoot = 'srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN='
#         srmPath = "/pnfs/lcg.cscs.ch/cms/trivcat"+dir
#     else:
#         parser.error('site can be either t3psi or t2cscs')
    print 'siteRoot \''+siteRoot+'\''
    print 'srmPath  \''+srmPath+'\''
    (dirs,files) = recursiveList(siteRoot,srmPath)
    if len(dirs) is not 0:
        print 'Directories found:'
        print string.join(dirs,'\n')
    if len(files) is not 0:
        print 'File found:'
        print string.join(files,'\n')

    if query_yes_no('Do you want to proceed?','no'):
        performDelete(siteRoot,dirs,files)
    else:
        print 'Removal aborted'
#     (header,entries) = srmList(srmPath)

#     (dirs,files) = sortByType(entries)
#     print 'Dirs:\n',dirs
#     print 'Files:\n',files

#     for dir in dirs:
#         print dir
#         srmList(dir)
            
if __name__ == '__main__':
        main()
