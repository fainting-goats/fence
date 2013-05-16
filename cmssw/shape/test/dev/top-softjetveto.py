#!/usr/bin/env python

import os
macrofile = './loadsamples.py'
if macrofile.startswith('./'):
    mydir = os.path.dirname(__file__)
    macrofile = mydir+macrofile[1:]
#execfile(mydir+'/loadsamples.py')
execfile(macrofile)

from ginger.analysis import CutFlow
from ginger.plotter import H1RatioPlotter
import ROOT
import copy

cf = CutFlow()
cf['bveto_mu'] = 'bveto_mu'
cf['tagpt']  = 'jetpt2>20 && (njet>=1 && njet<=2)'
cf['subtag'] = 'jettche2>2.1'


steps = [25, 20, 15, 10, 5]
for s in steps:
    cf[str(s)] = 'jetpt3 < %d' % s

for a in analysers.itervalues():
    a.extend( cf )

attbar = analysers['ttbar']
atW    = analysers['tW']

# here would be a nice breakpoint.

# plot the distributions of tche for ttbar
pttB = attbar.plotsflow('btag1','jettche1',bins=(500,-20,30))
#pttB = atW.plotsflow('btag1','jettche1',bins=(500,-20,30))

    #h.GetXaxis().SetTitle('jettche1')

for i,h in enumerate(pttB[4:].itervalues()):
    h.SetLineColor(i+2)

pttN = copy.deepcopy(pttB)
for h in pttN[1:].itervalues():
    h.Scale(1./h.Integral())



colors = [1,2,3,4,5,6,7,8,9,10]
markers = [ROOT.kFullCircle]*10
#h1 = H1RatioPlotter(colors=colors,markers=markers)
h1 = H1RatioPlotter()
h1.plotratio=False
h1.legalign = ('r','t')
h1.legtextsize = 20
h1.legboxsize = 35
h1.userrangex = (-5,15)
h1.set( *(pttB[3:].values()))
c1 = h1.plot('hist')


h2 = H1RatioPlotter()
h2.plotratio=False
h2.legalign = ('r','t')
h2.legtextsize = 20
h2.legboxsize = 35
h2.userrangex = (-5,15)
h2.set( *(pttN[3:].values()))
c2 = h2.plot('hist')

jetpt3_A  = attbar.views['subtag'].plot('subtag'         , 'jetpt3',bins=(30,0,30))
jetpt3_B  = attbar.views['subtag'].plot('jet1 btag'      , 'jetpt3',bins=(30,0,30), cut='jettche1>2.1')

h3 = H1RatioPlotter()
h3.colors    = [ROOT.kRed+1      , ROOT.kAzure-6    , ROOT.kRed+1      , ROOT.kAzure-7    , ROOT.kOrange-2  ]
h3.markers   = [ROOT.kFullCircle , ROOT.kFullCircle , ROOT.kOpenCircle , ROOT.kOpenCircle , ROOT.kFullCircle]
h3.set(jetpt3_A,jetpt3_B)#,jetpt3_A1, jetpt3_B1)
h3.legalign = ('r','b')
c3 = h3.plot()

jettche3_A  = attbar.views['subtag'].plot('subtag'         , 'jettche3 > jettche1',bins=(2,0,2))
jettche3_B  = attbar.views['subtag'].plot('jet1 btag'      , 'jettche3 > jettche1',bins=(0,0,2), cut='jettche1>2.1')

h4 = H1RatioPlotter()
h4.colors    = [ROOT.kRed+1      , ROOT.kAzure-6    , ROOT.kRed+1      , ROOT.kAzure-7    , ROOT.kOrange-2  ]
h4.markers   = [ROOT.kFullCircle , ROOT.kFullCircle , ROOT.kOpenCircle , ROOT.kOpenCircle , ROOT.kFullCircle]
h4.set(jettche3_A,jettche3_B)#,jetpt3_A1, jetpt3_B1)
h4.legalign = ('r','b')
c4 = h4.plot()
#jetpt4_A = attbar.views['subtag'].plot('stoca','jetpt4',bins=(30,0,30))
#jetpt4_B = attbar.views['subtag'].plot('stoca','jetpt4',cut='jettche1>2.1',bins=(30,0,30))

#h4 = H1RatioPlotter()
#h4.colors    = [ROOT.kRed+1      , ROOT.kAzure-6    , ROOT.kAzure+9    , ROOT.kOrange+7   , ROOT.kOrange-2  ]
#h4.markers   = [ROOT.kFullCircle , ROOT.kFullCircle , ROOT.kOpenCircle , ROOT.kFullCircle , ROOT.kFullCircle]
#h4.set(jetpt4_A,jetpt4_B)
#c4 = h4.plot()


#c2 = ROOT.TCanvas()
#pttN['subtag'].Draw('hist')
#for i,h in enumerate(pttN[4:].itervalues()): h.SetLineColor(i+2); h.Draw('hist same')
#c2.BuildLegend()
