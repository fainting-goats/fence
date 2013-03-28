#!/usr/bin/env python

import csv
import optparse
import re
import ConfigParser
import odict
import os
import psitools


# def strToNumbers( numString ):
#     numbers = []
#     if len(numString)==0:
#         return numbers;

#     items = numString.split(',')
#     for item in items:
#         nums = item.split('-')
#         if len(nums) == 1:
#             single number
#             numbers.append(int(item))
#         elif len(nums) == 2:
#             i = int(nums[0])
#             j = int(nums[1])
#             if i > j:
#                 raise ValueError('Invalid interval '+item)
#             numbers.extend(range(i,j+1))
#         else:
#             raise ValueError('Argument '+item+'is not valid')

#     return set(numbers)

# class Dataset:
#     def __init__(self, id, nick, ver, name):
#         self.id = id
#         self.nick = nick
#         self.ver  = ver
#         self.name = name
#         self.transfer = False
#         self.t0 = None
#         self.t1 = None


class MultiCrabConfigurator:
    def __init__(self, csvFile, datatype,crabcfg, ids):
        self.csvFile = csvFile
        self.crabcfg = crabcfg
        self.datatype = datatype
        self.ids = ids
        self.datasets = psitools.getDataSets( self.csvFile, self.ids )
        self.parser = ConfigParser.ConfigParser( dict_type=odict.OrderedDict)

#     def configure(self):
#         reader = csv.reader(open(self.csvFile,'rb'),delimiter=',')
#         header = reader.next()
#         l = len(header)
#                 print 'header len',l

#         header = [ h.lower() for h in header]

#         cols = dict(zip(header,range(len(header))))
#         print cols
#         for row in reader:
#             if len(row) != l:
#                 continue

#             d = Dataset(row[cols['id']],row[cols['nickname']],row[cols['skim version']],row[cols['output dataset']])

#             numId = int(re.search('[0-9]+',d.id).group(0))
#             if len(self.ids) != 0 and not numId in self.ids:
#                 continue

#             self.datasets.append(d);

#             if d.ver == '':
#                 print 'Dataset',d.nick,'doesn\'t have a version number'
#                 continue

    def build(self):

        # check that the cfg file exists
        if not os.path.exists(self.crabcfg):
            raise NameError('File '+self.crabcfg+' does not exists')

        crabparser = ConfigParser.ConfigParser()
        crabparser.read(self.crabcfg)

        crabdir = os.path.dirname(self.crabcfg)
        if crabdir == '':
            crabdir = '.'
        pypath = crabparser.get('CMSSW','pset')

        if pypath[0]!='/':
            print 'relative'
            pypath = crabdir+'/'+pypath

        print 'crabdir =',crabdir
        print 'pypath  =',pypath

        self.parser.add_section('MULTICRAB')
        self.parser.set('MULTICRAB','cfg',self.crabcfg)

        self.parser.add_section('COMMON')
        self.parser.set('COMMON','CMSSW.total_number_of_events','-1')
        if self.datatype=='data':
            self.parser.set('COMMON','CMSSW.total_number_of_lumis','-1')
        self.parser.set('COMMON','CMSSW.number_of_jobs','10')
        self.parser.set('COMMON','CMSSW.pset',pypath)
        self.parser.set('COMMON','USER.ui_working_dir','./jobs_'+self.datatype)

        for d in self.datasets:
            print '- Adding',d.nick
            self.parser.add_section(d.nick)
            self.parser.set(d.nick, 'CMSSW.datasetpath', d.name)
            self.parser.set(d.nick, 'CMSSW.pycfg_params', 'outputFile='+d.nick+'.root')


    def write(self, cfgfile):
        print cfgfile
        configfile = open(cfgfile, 'wb')
        self.parser.write(configfile)
        

if __name__ == '__main__':
    usage = 'usage: %prog -p [options] <file.csv> '
    parser = optparse.OptionParser(usage)
#     parser.add_option('--site', '-s', dest='site', help='Site where files are located. Can be [t3psi,t2cscs]')
    parser.add_option('--ids',  '-i', dest='ids',  help='Dataset ids to be grabbed from the spreadsheet', default='')
    parser.add_option('--cfg',  '-c', dest='crabcfg',  help='cmsRun python config file') 
    parser.add_option('--type', '-t', dest='datatype',  help='data type, can be [data,mc]') 
    (opt,args) = parser.parse_args()

    if len(args) == 0:
        parser.error('Input csv file not defined')

    csvFile = args[0]
    
    if opt.crabcfg is None:
        parser.error('Please specify the python config file')

    allowedDataypes = ['mc','data']
    if opt.datatype is None or opt.datatype not in allowedDataypes :
        parser.error('Please specify the datatype')

    ids = strToNumbers(opt.ids)
    print ' - Ids considered for inclusion: ', opt.ids

    c = MultiCrabConfigurator( csvFile, opt.datatype, opt.crabcfg, ids)
    c.build()
    c.write('multicrab.cfg')


