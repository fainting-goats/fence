#!/usr/bin/env python


import os.path
import numpy

import ROOT
treename = 'latino'
var = 'ptll'
njet = '0'
mstyle = ROOT.kFullCircle
msize  = .7
colors = [ROOT.kBlack,ROOT.kBlue, ROOT.kRed, ROOT.kGreen+1,ROOT.kViolet-1]
opt=''

postfix='_'+opt.replace(' ','-') if opt else ''

path = '/shome/thea/HWW/work/shape2012/trees/nominals/'
dyfiles = [
    path+'latino_036_DY10toLLMad.root',
    path+'latino_037_DY50toLLMad.root',
]

mypath = os.path.dirname(os.path.abspath(__file__))
ROOT.gInterpreter.ExecuteMacro(mypath+'/LatinoStyle2.C')
ROOT.TH1.SetDefaultSumw2()

from ROOTtree import sample, TreeWorker
from ROOTplotter import plot

samples = [
    sample(
        title     = 'lomet',
        tree      = 'latino',
        selection = 'mpmet > 20',
        files     = dyfiles
    ),

    sample(
        title     = 'dymva >0.6',
        tree      = 'latino',
        selection = 'dymva1 > 0.6',
        files     = dyfiles
    ),

    sample(
        title     = 'dymva >0.3',
        tree      = 'latino',
        selection = 'dymva1 > 0.3',
        files     = dyfiles
    ),

    sample(
        title     = 'dymva >0.0',
        tree      = 'latino',
        selection = 'dymva1 > 0.0',
        files     = dyfiles
    ),

    sample(
        title     = 'dymva >-0.3',
        tree      = 'latino',
        selection = 'dymva1 > -0.3',
        files     = dyfiles
    )
]

if njet=='1': 
    del samples[1]
    del colors[1]

cut = 'ptll > 20 && njet == '+njet


workers = []
for s in samples:
    print 'Making',s.title
    w = TreeWorker(s.tree,s.files)
    w.setselection( s.selection )
    w.setweight(s.weight)

    workers.append(w)
    print ' Entries',w.entries()

w0 = workers[0]
wmvas = workers[1:]

ROOT.TH1.AddDirectory(False)
edges = range(0,60,5)+range(60,100,10)+range(100,200,20)
bins = numpy.array(edges,dtype='d')
ht = ROOT.TH1F('h_ptll',';ptll',len(bins)-1,bins)

print 'Plottong'
# hs = [ w.plot('h_ptll',var,cut,bins=(100,0,200), options=opt) for w in workers]
hs = [ w.fill(ht.Clone(),var,cut, options=opt) for w in workers]

stack = ROOT.THStack('xx',';'+var)

for i,h in enumerate(hs): 
    h.SetTitle(samples[i].title)
    h.SetLineColor(colors[i])
    h.SetMarkerColor(colors[i])
    h.SetMarkerStyle(mstyle)
    h.SetMarkerSize(msize)
    stack.Add(h)

c = ROOT.TCanvas('check','check',5)
c.SetLogy()
stack.Draw('nostack')

legend = c.BuildLegend(0.7)
legend.SetFillStyle(0)
legend.SetBorderSize(1)
for e in legend.GetListOfPrimitives():
    e.SetOption('lp')

c.SaveAs('www/test/dyrw_'+njet+'j_ptll'+postfix+'.png')
c.SaveAs('www/test/dyrw_'+njet+'j_ptll'+postfix+'.pdf')
c.Clear()

h0 = hs[0]
hx = hs[1:]

rstack = ROOT.THStack('yy',';'+var)

for i,h in enumerate(hx):
    hr = h.Clone()
    hr.Divide(h0)

    hr.SetLineColor(colors[i+1])
    hr.SetMarkerColor(colors[i+1])
    hr.SetMarkerStyle(mstyle)
    hr.SetMarkerSize(msize)
    rstack.Add(hr)

rstack.Draw('nostack')

legend = c.BuildLegend(0.7)
legend.SetFillStyle(0)
legend.SetBorderSize(1)
for e in legend.GetListOfPrimitives():
    e.SetOption('lp')

c.SaveAs('www/test/dyrw_'+njet+'j_ratio'+postfix+'.png')
c.SaveAs('www/test/dyrw_'+njet+'j_ratio'+postfix+'.pdf')

import sys
sys.exit(0)
