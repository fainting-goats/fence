#!/usr/bin/env python


import HiggsAnalysis.CombinedLimit.DatacardParser as cardparser
import hwwinfo
import glob
import sys
import operator
import optparse
import os.path
import re

class Option: pass

def main():
    input = 'all_syst/datacards/hww-4.63fb.mH130.sf_0j_shape.txt'
#     input = 'all_syst/datacards/hww-4.63fb.mH130.comb_0j_shape.txt'

    # input directory
    # output directory
    # fake rate (factor)
    # dy extr (factor)

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('--prefix','-p',dest='prefix',help='prefix',default=None)
    parser.add_option('--output','-o',dest='output',help='output',default=None)
    parser.add_option('--FakeRate', dest='fakeRate',help='FakeRate ratio',default=None, type='float')
    parser.add_option('--DYLLextr', dest='dyllExtr',help='DYLL extr ratio',default=None, type='float')
    parser.add_option('--filter','-f',   dest='filter',  help='Datacard filter',default='*.txt')
#     hwwtools.addOptions(parser)
#     hwwtools.loadOptDefaults(parser)
    (opt, args) = parser.parse_args()

    if not opt.prefix:
        parser.print_help()
        parser.error('Prefix not defined')

    if not opt.output:
        parser.print_help()
        parser.error('Need to specify the output path')

    inputpath = os.path.join(opt.prefix,'datacards')

    datacards = glob.glob(inputpath+'/'+opt.filter)

    output = os.path.join(opt.output,'datacards')
    os.system('mkdir -p '+output)

    for dc in datacards:
        print dc
        indc = open(dc,'r')
        
        outdc = open(os.path.join(output,os.path.basename(dc)),'w')

        tarocCard(indc, outdc, opt)
    


def tarocCard(input, output, opt):
#     f = open(input,'r')
    
    dc_opt = Option()

    dc_opt.bin = True
    dc_opt.stat = False

    dc = cardparser.parseCard(input,dc_opt)

    # tarocching going on here

    f = operator.itemgetter(0)
    effects = map(f,dc.systs)

    # fake rate slimming
    if opt.fakeRate:
        try:
            fr = dc.systs[effects.index('FakeRate')][4]

            for bin,proc in fr.iteritems():
                for p,v in proc.iteritems():
                    if v != 0.:
                        v = (v-1.)*opt.fakeRate+1
                        proc[p] = v
        except ValueError:
            pass

    if opt.dyllExtr:
        expr = re.compile('CMS_hww_DYLL_\dj_extr')
        for syst in dc.systs:
#             extr = dc.systs[effects.index('CMS_hww_DYLL_0j_extr')][4]
            if not expr.match(syst[0]):
                continue
            extr = syst[4]

            for bin,proc in extr.iteritems():
                for p,v in proc.iteritems():
                    if v != 0.:
                        v = (v-1.)*opt.dyllExtr+1
                        proc[p] = v

    writeDatacard(dc,output)

def writeDatacard( dc, out):
    print >> out, 'Tarocched datacard' 
    print >> out, "imax %d number of bins" % len(dc.bins)
    print >> out, "jmax %d number of processes minus 1" % (len(dc.processes) - 1)
    print >> out, "kmax %d number of nuisance parameters" % (len(dc.systs))
    print >> out, "-" * 130



    clens = []
    for b,vals in dc.shapeMap.iteritems():
        for filter,args in vals.iteritems():
            clens.append(max([len(filter),len(b)]))

    cmax = max(clens)
    cfmt = "%-"+str(cmax)+"s ";

    for b,vals in dc.shapeMap.iteritems():
        for f,args in vals.iteritems():
            print >> out, 'shapes',cfmt % f,cfmt % b,' '.join(args)

#     print >> out, "-" * 130
    cmax = max([len(b) for b in dc.bins]+[len(p) for p in dc.processes])
    cfmt = "%-"+str(cmax)+"s ";
    
    print >> out, "bin         ", "  ".join([cfmt % x for x in dc.bins])
    print >> out, "observation ", "  ".join([cfmt % int(dc.obs[x]) for x in dc.bins])
    print >> out, "-" * 130

    pidline = []; signals = []; backgrounds = []
    tmpsignals = [];
#     for (b,p,s) in dc.keyline:
#         if s:
#             if p not in tmpsignals: tmpsignals.append(p)
    for (b,p,s) in dc.keyline:
        if s:
#             if p not in signals: signals.append(p)
            pidline.append(dc.signals.index(p)-len(dc.signals)+1)
        else:
            if p not in backgrounds: backgrounds.append(p)
            pidline.append(1+backgrounds.index(p))

    smax = max([len(l) for (l,nf,p,a,e) in dc.systs])

    hmax = max([10] + [len( ('%-'+str(smax)+'s[nofloat]  %s %s') % (l,p,' '.join([str(x) for x in a]))) for (l,nf,p,a,e) in dc.systs])
    hfmt = "%-"+str(hmax)+"s "


    systlen = []
    for (l,nf,pdf,a,e) in dc.systs:
        for b,p,s in dc.keyline:
            systlen.append(len(str(e[b][p])))
        
#         print >> out, [ i for i in x.itervalues() for x in e.itervalues() for (l,nf,pdf,a,e) in dc.systs]
    cmax = max([cmax]+[max(len(p),len(b)) for p,b,s in dc.keyline]+systlen)
 
    cfmt = "%-"+str(cmax)+"s ";

    print >> out, hfmt % "bin",     "  ".join([cfmt % b for b,p,s in dc.keyline])
    print >> out, hfmt % "process", "  ".join([cfmt % p for b,p,s in dc.keyline])
    print >> out, hfmt % "process", "  ".join([cfmt % x for x in pidline])
    print >> out, hfmt % "rate",    "  ".join([cfmt % dc.exp[b][p] for b,p,s in dc.keyline])
    
    print >> out, "-" * 130
    # (lsyst,nofloat,pdf,args,errline)
    for (l,nf,p,a,e) in dc.systs:
        print >> out, hfmt % ( ('%-'+str(smax)+'s %s %s') % (l,p,' '.join([str(x) for x in a]))) ,
        effects = []
        for b,p,s in dc.keyline:
            if e[b][p] == 0:
                effects.append('-')
            else:
                effects.append(e[b][p])
        print >> out, '  '.join([ cfmt % e for e in effects])
        


if __name__ == '__main__':
    main()
