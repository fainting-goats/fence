#!/usr/bin/env python

import HWWAnalysis.Misc.odict as odict
from HWWAnalysis.Misc.ROOTAndUtils import Tee

outdir = 'www/test/truth_xxx'
logpath = outdir+'/truth.log'
tee = Tee(logpath)

import sys
sys.argv.append('-b')

import hwwlatino

from hwwinfo2g import wwnamedcuts as wwcuts
from ginger.analysis import CutFlow,AnalysisView

orchard = '/scratch/thea/test_b'

samples = hwwlatino.samples(125,'8TeV','Data2012','SM','topplots')

wwflow = CutFlow(wwcuts.wwcommon)

del wwflow['bveto_mu']
del wwflow['bveto_ip']
wwflow['ptll'] = 'ptll>45'

print '------------------------------'
for n,c in wwflow.iteritems():
    print '%-30s: %s'% (n,c)

topflow = CutFlow()

topflow['base']   = wwflow.string()

lumi = 19.468

for mc in samples.keys():
    if mc in ['ttbar','tW']: continue
    if mc in ['tW']: continue

    samples.pop(mc)

analysers = hwwlatino.makeanalysers(samples,orchard,topflow,lumi)

from ginger.filters import UnderOverTucker
uo = UnderOverTucker()

for a in analysers.itervalues(): a.filters.append(uo)
import ROOT


flavours = odict.OrderedDict([
    ('bquarks',('abs(jetflv{0})==5' , ROOT.kOrange)),
    ('gluons' ,('abs(jetflv{0})==21', ROOT.kRed-7)),
    ('others' ,('abs(jetflv{0})!=5 && abs(jetflv{0})!=21',ROOT.kGreen-7)),
])

flavours = odict.OrderedDict([
    ('bquarks','abs(jetflv1)==5'),
    ('gluons' ,'abs(jetflv1)==21'),
    ('others' ,'abs(jetflv1)!=5 && abs(jetflv1)!=21'),
])

ptbins  = range(30, 70, 10) + range( 70, 150, 20) + range(150,250, 50) + [250,500,1000]
etabins = [0.,0.75,1.5,2.8,5]
etabins = [0.,1.4,2.8,5]
etabins = [0.,2.8,5]
#ptbins  = range(30,300,10)#+range(300,1001,100)
ptbins  = range(30,300,10)#+range(300,1001,100)
#ptbins  = [30,1000]
#etabins = [0.,2.8,5]

plt = analysers['ttbar'].splitplot('pippo','jetpt2',flavours,bins=(ptbins,))
ylds= analysers['ttbar'].splityields(flavours)

for n,p in plt.iteritems():
    print n,p.Integral()
    print n,ylds[n]

sys.exit(0)


abcd1j = {}
abcd1j['A']   = 'njet==2 && bveto_munj30 && softtche<=2.1 && jettche2>2.1'
abcd1j['B']   = 'jettche1>2.1'
abcd1j['C']   = 'njet==1 && bveto_munj30 && softtche<=2.1'
abcd1j['D']   = 'jettche1>2.1 '
abcd1j['C-D'] = 'nbjettche==0'
abcd1j['mu']  = 'bveto_mu'

abcd1j_mu = {}
abcd1j_mu['A']   = 'njet==2 && bveto_mu && softtche<=2.1 && jettche2>2.1'
abcd1j_mu['B']   = 'jettche1>2.1'
abcd1j_mu['C']   = 'njet==1 && bveto_mu && softtche<=2.1'
abcd1j_mu['D']   = 'jettche1>2.1 '
abcd1j_mu['C-D'] = 'nbjettche==0'
abcd1j_mu['mu']  = 'bveto_mu'

regions = abcd1j_nosoft

from ginger.filters import UnderOverTucker
from commons import AlienDict

# specific cuts
cutsAB = CutFlow()
cutsAB['A'] = regions['A']
cutsAB['B'] = regions['B']

cutsCD = CutFlow()
cutsCD['C'] = regions['C']
cutsCD['D'] = regions['D']

def plotregios( analysers, cuts):

    # print the cuts because it's nice
    print '--1j cuts--'.ljust(80,'-')
    print '\n'.join(['%-30s: %s'% i for i in cuts.iteritems()])

    # and we want to tuck the under and overflows in
    uo = UnderOverTucker()

    print '--Updating buffers--'.ljust(80,'-')
    for n,a in analysers.iteritems():
        old = a.selectedentries()
        a.update(cuts)
        a.bufferentries()
        # and tuck in underflow/overflows by default
        a.filters.append(uo)
        print '  ',n,':',old,'>>',a.selectedentries(),'...done'

    ptbins  = range(30, 70, 10) + range( 70, 150, 20) + range(150,250, 50) + [250,500,1000]
    etabins = [0.,0.75,1.5,2.8,5]
    etabins = [0.,1.4,2.8,5]
    etabins = [0.,2.8,5]
    #ptbins  = range(30,300,10)#+range(300,1001,100)
    ptbins  = range(30,300,10)#+range(300,1001,100)
    #ptbins  = [30,1000]
    #etabins = [0.,2.8,5]

    import ROOT


    flavours1 = {k:v[0].format(1) for k,v in flavours.iteritems()}
    flavours2 = {k:v[0].format(2) for k,v in flavours.iteritems()}

    print 'plots: [',
    plots = AlienDict()
    for n,a in analysers.iteritems():
        print '%s,' % n,
        # plot the 2d histograms for the efficiency
        for f,extra in flavours1.iteritems():
            pf = a.plotsflow(n+'_etapt_'+f,'fabs(jeteta1):jetpt1',extra=extra, bins=(ptbins,etabins))

            for c,h in pf.iteritems():
                # rearrange the plots as cut,process
                plots[c][n][f]['2D'] = h
                for i,eta in enumerate(etabins):
                    hp = h.ProjectionX('%s_%d' % (h.GetName(),i),i,i,'e')
                    hp.SetLineWidth(2)
                    hp.SetLineColor(flavours[f][1]+1)
                    hp.SetFillColor(flavours[f][1])
                    plots[c][n][f][i] = hp

        sys.stdout.flush()

    plots.lock()
    print ']'

    #print plots

    i = 1

    for reg in cuts.keys()+['base']:
        for n in analysers.iterkeys():
            print ('--- '+reg+':'+n+' ').ljust(80,'-')
            ys = {f:plots[reg][n][f][i].Integral() for f in flavours.iterkeys()}
            y_all = sum(ys.itervalues())
            #print 'Total yield:',y_all
            print '| '+' | '.join([' '*10       ,'y_all'.rjust(10)]+         ['{0:>10}'.format(k)       for k in ys.iterkeys()   ])+' |'
            print '| '+' | '.join([reg.ljust(10),'{0:>10.1f}'.format(y_all)]+['{0:>10.4f}'.format(y/y_all) for y in ys.itervalues() ])+' |'
            #for f,y in ys.iteritems():
                #print f,'%.3f %%' % (100*y/y_all)

        for n in analysers.iterkeys():
            stack = ROOT.THStack()
            for f in flavours:
                stack.Add(plots[reg][n][f][i])

            c = ROOT.TCanvas()
            stack.Draw('hist')
            stack.GetXaxis().SetTitle('jetpt1')
            leg = c.BuildLegend()
            leg.SetFillColor(ROOT.kWhite)
            leg.SetBorderSize(0)
            c.SaveAs('%s/%s_%s.pdf'% (outdir,reg,n))
            c.SaveAs('%s/%s_%s.png'% (outdir,reg,n))


import copy

#plotregios(copy.deepcopy(analysers),cutsAB)
#plotregios(copy.deepcopy(analysers),cutsCD)
