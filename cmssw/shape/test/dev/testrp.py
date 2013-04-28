#!/bin/env python

import ROOT
import os
import numpy as np
import array

from ginger.plotter import H1RatioPlotter

def testplotter():

    ROOT.TH1.SetDefaultSumw2()
 
    bins = (20,-3,3)
    edges = list(np.linspace(-3,0,10))+list(np.linspace(1.,4.,4))
    bins = (len(edges)-1,array.array('d',edges))

    fs1 = ROOT.TF1('s1','5*exp(-0.5*((x-0)/1)**2)',-3,3)
    fb1 = ROOT.TF1('b1','x*x',-3,3)
    fsb1 =ROOT.TF1('sb1','s1+b1',-3,3)

    h = ROOT.TH1F('aa','aa',*bins)
    h.FillRandom('s1',100)

    g = ROOT.TH1F('bb','bb',*bins)
    g.FillRandom('b1',50)

    data = ROOT.TH1F('data_obs','data_obs',*bins)
    data.FillRandom('sb1',150)

    hr = H1RatioPlotter(legmargin=10,markersize=20,errsty=3001,errcol=ROOT.kRed)
    hr.set(h+g,data)
    hr.scalemax = 6
#     hr.userrange = (-2.5,3)
#     hr.ltitle   = 'stocazz'
#     hr.rtitle   = 'stameng'
#     hr.ytitle2  = '[cm]'
#     hr.colors   = [ROOT.kRed,ROOT.kAzure]
#     hr.textsize = 40
#     hr.legboxsize = 24
#     hr.legtextsize = 15

#     hr.markersize = 32
#     hr.legalign = ('r','b')
#     hr.legalign = ('r','t')
    

    c = hr.plot('E1')
    c[0,0].SetLogy()
    c.Print('testrp.png')
    c.Print('testrp.pdf')

if __name__ == '__main__':
    testplotter()

