#!/bin/env python

import ROOT
import os

import hwwplot
shape_path = os.path.join(os.getenv('CMSSW_BASE'),'src/HWWAnalysis/ShapeAnalysis/')
print 'Shape directory is',shape_path
ROOT.gInterpreter.ExecuteMacro(shape_path+'macros/LatinoStyle2.C');


defaults  = {
    'ggH':  { 'color':633, 'label':'ggH m_{H}=125', },
    'vbfH': { 'color':634, 'label':'qqH m_{H}=125', },
    'VH':   { 'color':635, 'label':'VH m_{H}=125' , },

    'VV':   { 'color':858, 'label':'WZ/ZZ' , },  
    'WJet': { 'color':921, 'label':'W+jets' , },
    'Vg':   { 'color':617, 'label':'V+#gamma' , },
    'VgS':  { 'color':618, 'label':'V+#gamma*' , },
    'Top':  { 'color':400, 'label':'top' , },
    'DYTT': { 'color':418, 'label':'DY+jets' , },
#     'DYLL': { 'color':419, 'label':'DY+jets' , },
    'WW':   { 'color':851, 'label':'WW' , },
    'ggWW': { 'color':853, 'label':'ggWW' , },

}

shapes = {}
names = ['data','ggH', 'vbfH', 'VH', 'VV', 'WJet', 'Vg', 'VgS', 'Top', 'DYTT', 'WW', 'ggWW']
l = len(names)
for i,n in enumerate(names):
    h = ROOT.TH1F('name_'+n,'title '+n,l,0,l)
    h.Fill(i+0.5,i+1)
    shapes[n] = h

d = shapes['data']
for i in xrange(l):
    d.Fill(i+0.5,i+2*(i%2) )


plot = hwwplot.HWWPlot( defaults )
# print shapes

plot.setdata(shapes['data'])

plot.addsig('ggH',  shapes['ggH'])
plot.addsig('vbfH', shapes['vbfH'])
plot.addsig('VH',   shapes['VH'], label='stocazz')

plot.addbkg('VV',   shapes['VV'])
plot.addbkg('WJet', shapes['WJet'])
plot.addbkg('Vg',   shapes['Vg'])
plot.addbkg('VgS',  shapes['VgS'])
plot.addbkg('Top',  shapes['Top'])
plot.addbkg('DYTT', shapes['DYTT'])
# plot.addbkg('DYLL', shapes['DYLL'])
plot.addbkg('WW',   shapes['WW'])
plot.addbkg('ggWW', shapes['ggWW'])

plot.prepare()

## 1 = signal over background , 0 = signal on its own
plot.set_addSignalOnBackground(0);

## 1 = merge signal in 1 bin, 0 = let different signals as it is
plot.set_mergeSignal(0);

plot.setMass(125); 

c1 = ROOT.TCanvas("mll","mll",500,600)

# hookDebugger()
plot.Draw(c1, 1, True)

c1.Print('testplot.pdf')

if False
    shape_path = os.path.join(os.getenv('CMSSW_BASE'),'src/HWWAnalysis/ShapeAnalysis/')
    print 'Shape directory is',shape_path

    ROOT.gROOT.LoadMacro(shape_path+'macros/PlotVHqqHggH.C+')
    ROOT.gInterpreter.ExecuteMacro(shape_path+'macros/LatinoStyle2.C');

    fs1 =ROOT.TF1('s1','5*exp(-0.5*((x-0)/1)**2)',-3,3)
    fb1 =ROOT.TF1('b1','x*x',-3,3)
    fsb1 =ROOT.TF1('sb1','s1+b1',-3,3)

    h = ROOT.TH1F('aa','aa',100,-3,3)
    h.FillRandom('s1',1000)

    g = ROOT.TH1F('bb','bb',100,-3,3)
    g.FillRandom('b1',500)

    data = ROOT.TH1F('data_obs','data_obs',100,-3,3)
    data.FillRandom('sb1',1500)


    hs = ROOT.PlotVHqqHggH()

    hs.setLumi(12.103);
    hs.setLabel("unrolled");
    hs.addLabel("    #sqrt{s} = 8 TeV");

    hs.setDataHist( data )


    vectColourBkg        = ROOT.vector('int')()
    vectSystBkg          = ROOT.vector('double')()
    vectScaleBkg         = ROOT.vector('double')()
    vectNameBkg          = ROOT.vector('std::string')()
    vectNormalizationBkg = ROOT.vector('double')()
    vectTHBkg            = ROOT.vector('TH1F*')()

    vectColourSig        = ROOT.vector('int')()
    vectSystSig          = ROOT.vector('double')()
    vectScaleSig         = ROOT.vector('double')()
    vectNameSig          = ROOT.vector('std::string')()
    vectNormalizationSig = ROOT.vector('double')()
    vectTHSig            = ROOT.vector('TH1F*')()

    #----

    vectColourSig.push_back(632)    
    vectNameSig.push_back('sig')     
    vectTHSig.push_back(h)

    vectColourBkg.push_back(ROOT.kBlue)    
    vectNameBkg.push_back('bkg')     
    vectTHBkg.push_back(g)

    #----
    hs.set_vectTHSig     (vectTHSig);      
    hs.set_vectNameSig   (vectNameSig);    
    hs.set_vectColourSig (vectColourSig);  

    hs.set_vectTHBkg     (vectTHBkg);      
    hs.set_vectNameBkg   (vectNameBkg);    
    hs.set_vectColourBkg (vectColourBkg);  

    hs.prepare();


    hs.set_addSignalOnBackground(1); ## 1 = signal over background , 0 = signal on its own

    hs.set_mergeSignal(0);    ## 1 = merge signal in 1 bin, 0 = let different signals as it is

    hs.setMass(125); 

    c1 = ROOT.TCanvas("mll","mll",500,600)

    hs.Draw(c1)

    c1.Print('testplot.pdf')
