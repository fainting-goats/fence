#!/usr/bin/env python

import sys
import re

import HiggsAnalysis.CombinedLimit.DatacardParser as dcparser
import optparse

#
#
#

def interPollate( path0, pathout, srcmass, intmass):
    file0 = open(path0,'r')

    endheader = False

    bins=[]
    pids=[]
    channels=[]
    rates=[]

    names = []
    ilines = []

    alllines = file0.readlines()
    for i,line in enumerate(alllines):
        if line.lower().startswith('observation'):
            endheader=True
            continue
        if not endheader:
            continue
        tokens = line.split()
        if line.startswith('process'):
            if not channels:
                channels = tokens[1:]
            else:
                pids = tokens[1:]
            names.append(line[:line.find(tokens[1])])
            ilines.append(i)
            
        if line.startswith('bin'):
            bins = tokens[1:]
            names.append(line[:line.find(tokens[1])])
            ilines.append(i)

        if line.startswith('rate'):
            rates = map(float,tokens[1:])
            names.append(line[:line.find(tokens[1])])
            ilines.append(i)


#     print 'bins:',bins
#     print 'pids:',pids
#     print 'chan:',channels
#     print 'rate:',rates
#     print names
#     print ilines

    if not len(bins)==len(pids)==len(channels)==len(rates):
        raise RuntimeError('channels mismatch!!!')

    l = len(bins)

    # alls = []
    # for i in xrange(l):
    #     alls.append( (bins[],pids
    temp =  zip(bins,pids,channels,rates)
        
    signals = ['ZH', 'WH', 'qqH', 'ggH']

    '''
    For GGF WW2l, it's $2*$7*0.108*0.108*9
    For VBF WW2l, it's $3*$7*0.108*0.108*9
    For VH->WW, it's ($4+$5+$6)*$7
    '''

    srcX = xsecs[srcmass]
    intX = xsecs[intmass]


    for i,(bin,pid,chan,rate) in enumerate(temp):
        a = 1
        b = 0
        if chan not in signals:
            continue
    #         print i
            temp[i] = (bin,pid,chan,0.)
        if chan in ['ggH']:
            a = srcX[1]*srcX[6]*108.*108*9
            b = intX[1]*intX[6]*108.*108*9
        if chan in ['qqH']:
            a = srcX[2]*srcX[6]*108.*108*9
            b = intX[2]*intX[6]*108.*108*9
        if chan in ['WH','ZH']:
            a = srcX[6]*(srcX[3]+srcX[4]+srcX[5])
            b = intX[6]*(intX[3]+intX[4]+intX[5])
        irate = rate/a*b
        temp[i] = (bin,pid,chan,irate)

    # print temp
    short = min(map(len,names))
    replaced = []
    replaced.append( names[0][:short]+' '.join([ bin.rjust(8) for (bin,pid,chan,rate) in temp ])+'\n' )
    replaced.append( names[1][:short]+' '.join([ pid.rjust(8) for (bin,pid,chan,rate) in temp ])+'\n' )
    replaced.append( names[2][:short]+' '.join([ chan.rjust(8) for (bin,pid,chan,rate) in temp ])+'\n' )
    replaced.append( names[3][:short]+' '.join([ ('%.3f' % rate).rjust(8) for (bin,pid,chan,rate) in temp ])+'\n' )

#     print ' '.join([ str(x).ljust(8) for x in range(1,len(srcX))] )
#     print ' '.join([ ('%.3f' % x).ljust(8) for x in srcX] )
#     print ' '.join([ ('%.3f' % x).ljust(8) for x in intX] )

    out = open(pathout,'w') 

    print 'mass:',srcmass,'->',intmass,' | ',path0,'->',pathout
    for i,line in enumerate(alllines):
        if i not in ilines:
            out.write( line )
            continue
        out.write( replaced[ilines.index(i)] )

    out.close()
    file0.close()


#-------

parser = optparse.OptionParser()
dcparser.addDatacardParserOptions(parser)

(opts,args) = parser.parse_args()

pathx = '/shome/thea/HWW/work/dc/xs_br_hww_ecm8tev.txt'
path0 = '/shome/thea/HWW/work/dc/cebby/150/hwwof_0j_cut_8TeV.txt'
path0 = '/shome/thea/HWW/work/dc/maiko/{mass}/'
outdir = '/shome/thea/HWW/work/dc/p5'

#                                      
xfile = open(pathx,'r') 
xsecs = {}
for i,line in enumerate(xfile):
    if i==0: continue
    tokens = line.split()
    xsecs[float(tokens[0])]=map(float,tokens)

xfile.close()
#                                      

import numpy
import os.path

masses=[110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 170, 180, 190, 200, 250, 300, 350, 400, 450, 500, 550, 600]

r0 = set(numpy.arange(110,150,0.5))
r1 = set(numpy.arange(150,200,1))
r0 = set(numpy.arange(130,135,0.5))
r1 = set()
# print r0
newmasses = sorted( (r0 | r1) - set(masses))

refmasses = map( lambda nm: min(masses,key=lambda m: abs(nm-m)),newmasses)

newAndRef = zip(newmasses,refmasses)
# print newAndRef

for nm,mref in newAndRef:
    print nm, mref
    pathm = path0.format( mass=mref)
    if not os.path.isdir(pathm):
        raise RuntimeError(mref+' is NOT a directory')

    newdir = os.path.join(outdir,'%.1f' % nm)
    if not os.path.exists(newdir):
        os.system('mkdir -p '+newdir)
    print newdir
    print pathm
    for root, subFolders, files in os.walk(pathm):
         
        for f in files:
            srcdc = os.path.join(pathm,f)
            newdc = os.path.join(newdir,f)
            interPollate(srcdc, newdc, nm, mref)

        


# sys.exit(0)

# massin  = 150
# massout = 152



# interPollate(path0, pathout, massin, massout)











#     dc0 = dcparser.parseCard(file0,opts)


#     print 'bins              ',dc0.bins
#     print 'obs               ',dc0.obs
#     print 'processes         ',dc0.processes
#     print 'signals           ',dc0.signals
#     print 'isSignal          ',dc0.isSignal
#     print 'keyline           ',dc0.keyline
#     print 'exp               ',dc0.exp
#     print 'systs             ',dc0.systs
#                               
#     print 'shapeMap          ',dc0.shapeMap
#     print 'hasShape          ',dc0.hasShape
#     print 'flatParamNuisances',dc0.flatParamNuisances

#     cmax  = 5
#     pollo = dc0

#     keyline = []
#     pids = {}
#     pids.update( dict( [(s,-1-i) for (i,s) in enumerate(pollo.signals)]) )
#     pids.update( dict( [(s,i)    for (i,s) in enumerate([p for p in pollo.processes if p not in pollo.signals]) ] ) ) 


#     for b,p,s in pollo.keyline:
#         if not s: continue
#         keyline.append( (b,p,pids[p],pollo.exp[b][p]) )

#     for b,p,s in pollo.keyline:
#         if s: continue
#         keyline.append( (b,p,pids[p],pollo.exp[b][p]) )

#     print keyline

#     sys.exit(0)
#     print
#     print '#Interpolation of',path0
#     print '#             and',path0
#     print
#     print "imax %d number of bins" % len(pollo.bins)
#     print "jmax %d number of processes minus 1" % (len(pollo.processes)-1)
#     print "kmax %d number of nuisance parameters" % len(pollo.systs)
#     print "-" * 130

#     cmax = max([cmax]+[len(l) for l in pollo.obs.iterkeys()]+[len(str(x)) for x in pollo.obs.itervalues()])
#     cfmt = "%-"+str(cmax)+"s";
#     print "bin         ", "  ".join([cfmt % x for x in pollo.obs.iterkeys()])
#     print "observation ", "  ".join([cfmt % x for x in pollo.obs.itervalues()])
#     print "-" * 130
#     print 'bin',' '.join([ bin for bin,chan,pid,rate in keyline])
#     print 'process',' '.join([ chan for bin,chan,pid,rate in keyline])
#     print 'process',' '.join([ str(pid) for bin,chan,pid,rate in keyline])
#     print 'rate',' '.join([ str(rate) for bin,chan,pid,rate in keyline])

