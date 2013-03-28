#!/bin/env python

import ROOT
import sys
import optparse

files = sys.argv[1:]
n = max([len(f) for f in files])
for rootfile in files:
    treename = 'latino'

    f = ROOT.TFile.Open(rootfile)
    t = f.Get(treename)

    print '-'*80
    print '-->',rootfile,':',', '.join([ b.GetName() for b in t.GetListOfBranches()])



