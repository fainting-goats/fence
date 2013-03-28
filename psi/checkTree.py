#!/bin/env python

import ROOT
import sys
import optparse

rootfile = sys.argv[1]
treename = 'latino'

f = ROOT.TFile.Open(rootfile)
latino = f.Get('latino')
# probe_tree = f.Get('probe_tree')

print rootfile.ljust(50),':',latino.GetName(), 'Entries:',latino.GetEntries()
# print rootfile.ljust(50),':',probe_tree.GetName(), 'Entries:',probe_tree.GetEntries()



