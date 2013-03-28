#!/bin/env python



import sys
import os.path

from HiggsAnalysis.CombinedLimit.DatacardParser import *

dcpath = 'pippo/datacards/hww-5.06fb.mH125.of_0j_shape.txt'
file = open(dcpath,'r')
class dummy: pass

options = dummy()
options.stat = False
options.bin = True
options.noJMax = False
options.nuisancesToExclude = []
options.nuisancesToRescale = []

DC = parseCard(file, options)

if not DC.hasShapes:
    sys.exit(-1)


# print DC.shapeMap

try:
    #here one has to check what channel to use?
    # what if there are more channels?

    hpath, hname =  DC.shapeMap['*']['data_obs']
except KeyError as e:
    print 'Shape for '+str(e) +'not found!' 
    sys.exit(-1)

import ROOT
ROOT.TH1.AddDirectory(False)

print os.path.join(os.path.dirname(dcpath),hpath)
hfile = ROOT.TFile.Open(os.path.join(os.path.dirname(dcpath),hpath))


hdata = hfile.Get(hname)
print hdata

template = hdata.Clone('shape_template')
template.SetTitle('shape_template')
template.Reset()

template.Print()
print 'AAAAA',template.GetXaxis().GetTitle()

wspath = os.path.splitext(dcpath)[0]+'.root'
if not os.path.exists(wspath):
    # workspace + parameters = shapes
    print 'Making the workspace...',
    sys.stdout.flush()
    os.system('text2workspace.py '+dcpath)
    print 'done.'


ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
wsfile = ROOT.TFile.Open(wspath)
if wsfile.__nonzero__():
    print 'Workspace file found'

wsfile.ls()
w = wsfile.Get('w')
# w.allPdfs().printMultiline(ROOT.std.cout, 3)

#prendi il morph
#mettici i parametr
#estrai i valori

channel = DC.bins[0]#'of_0j'
process = 'Top'

# data = w.data('data_obs')
# fast = w.pdf('shapeSig_{0}_{1}_morph'.format(channel,process))
# obs = fast.getObservables(data)
# pars = fast.getParameters(data)

# for i in xrange(data.numEntries()):
#     obs.__assign__(data.get(i))
#     print i,obs.find('CMS_th1x').getVal(),fast.getVal()




mlfdir = '.'
print 'Fitting the workspace...',
sys.stdout.flush()
# os.system('combine -M MaxLikelihoodFit --out ' + mlfdir + 
#           ' --saveNormalizations '+wspath)
print 'done.'

mlfpath = os.path.join(mlfdir,'mlfit.root')
mlffile = ROOT.TFile.Open(mlfpath)

norm_fit_b = mlffile.Get('norm_fit_b')
fit_b = mlffile.Get('fit_b')
norm_fit_s = mlffile.Get('norm_fit_s')
fit_s = mlffile.Get('fit_s')


# norm_fit_s.printMultiline(ROOT.std.cout,7)

norm_b = norm_fit_b.find('n_exp_bin{0}_proc_{1}'.format(channel,process))
norm_s = norm_fit_s.find('n_exp_bin{0}_proc_{1}'.format(channel,process))


print 'shapeBgk_{0}_{1}_morph'.format(channel,process)

data = w.data('data_obs')
morph = w.pdf('shapeBkg_{0}_{1}_morph'.format(channel,process))
# obs = fast.getObservables(data)
# pars = fast.getParameters(data)
# pars.printMultiline(ROOT.std.cout,3)

def list2set( alist ):
#     for p in alist:
#         pass

    next =  alist.createIterator()
    aset = ROOT.RooArgSet()
    while True: 
        p = next.Next()
        if not p: break
        aset.add( p )
    return aset

def pdf2TH1(h, pdf, data, norm=None, pars=None):
    print pdf
    pdf_obs  = pdf.getObservables(data)
    pdf_pars = pdf.getParameters(data)

    if h.GetNbinsX() != data.numEntries():
        raise ValueError('bins mismatch!')
    
    if pars:
        if isinstance(pars,ROOT.RooArgList):
            pars = list2set(pars)

        pdf_pars.__assign__(pars)

    for i in xrange(data.numEntries()):
        pdf_obs.__assign__(data.get(i))
        h.SetBinContent(i+1,pdf.getVal())

    if norm:
        if isinstance(norm,ROOT.RooAbsReal):
            norm = norm.getVal()

        h.Scale(norm/h.Integral())

    pass

print 'Fit Init'

# for i in xrange(data.numEntries()):
#     obs.__assign__(data.get(i))
#     print i,obs.find('CMS_th1x').getVal(),fast.getVal()

h = template.Clone('histo_'+process)
h.SetTitle(process)
h_s = template.Clone('s_histo_'+process)
h_s.SetTitle(process)
h_b = template.Clone('b_histo_'+process)
h_b.SetTitle(process)

pdf2TH1(h,morph,data, DC.exp[channel][process])
pdf2TH1(h_s,morph,data,norm_s,fit_s.floatParsFinal())
pdf2TH1(h_b,morph,data,norm_b,fit_b.floatParsFinal())

h_b.SetLineColor(ROOT.kBlue)
h_s.SetLineColor(ROOT.kRed)

c1 = ROOT.TCanvas()
h.Draw()
h_b.Draw('same')
h_s.Draw('same')
c1.GetListOfPrimitives().ls()
c1.SaveAs('test_'+process+'.pdf')

# pars.__assign__( list2set(fit_s.floatParsFinal() )

# print 'Fit_s'
# for i in xrange(data.numEntries()):
#     obs.__assign__(data.get(i))
#     print i,obs.find('CMS_th1x').getVal(),fast.getVal()



# TFile *_file0 = TFile::Open("pippo/datacards/hww-5.06fb.mH125.of_0j_shape.root")
# gSystem->Load("libHiggsAnalysisCombinedLimit")
# hpdf = w->pdf("shapeSig_ggH_of_0jPdfshapeSig_ggH_of_0jPdf")
# data = w->data("data_obs")
# RooArgSet * obs
# obs = hpdf->getObservables(data)
# *obs = *data.get(14)
# hpdf->getVal(obs)

