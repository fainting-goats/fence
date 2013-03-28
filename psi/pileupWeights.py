#!/usr/bin/env python
import ROOT
import sys
import optparse

def main(puROOT ):
    f = ROOT.TFile.Open(puROOT)
#     f.ls()
    pu = f.Get('pileup')
    entries = pu.GetEntries()
    integral = pu.Integral()
    print 'Entries',entries
    print 'Integral',integral
    weights = []
    for c in range(1,pu.GetNbinsX()+1):
        print c,pu.GetXaxis().GetBinCenter(c),pu.GetBinContent(c)
        weights.append(pu.GetBinContent(c)/integral)

    print 'Number of weights computed',len(weights)
    for i in range(len(weights)):
        print i,weights[i]
    print weights


if __name__ == '__main__':
    sys.argv.extend(['-b','-n'])
    main(sys.argv[1])

