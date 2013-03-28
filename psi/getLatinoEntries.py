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

    print rootfile.ljust(n+2),':',t.GetName(), t.GetTitle(), 'Entries:',t.GetEntries()



