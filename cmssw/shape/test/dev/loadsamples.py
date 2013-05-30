#!/usr/bin/env python

import sys
sys.argv.append('-b')

import hwwlatino

from hwwinfo2g import wwnamedcuts as wwcuts
from ginger.analysis import CutFlow,AnalysisView

orchard = '/shome/thea/HWW/work/dds/trees/top'

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

analysers = hwwlatino.makeanalysers(samples,orchard,topflow,lumi)
