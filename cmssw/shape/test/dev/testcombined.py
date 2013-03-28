#!/bin/env python

from HiggsAnalysis.CombinedLimit.DatacardParser import *
from HiggsAnalysis.CombinedLimit.ShapeTools import *
import ROOT
import os
import logging

class dcoptions: pass

class roofiter:
    def __init__(self,collection):
        self._iter = collection.fwdIterator()

    def __iter__(self):
        return self

    def next(self):
        o = self._iter.next()
        if not o: raise StopIteration
        return o
#---
def getnorms(pdf, obs, norms = None ):
    '''helper function to exctact the normalisation factors'''

    out = norms if norms!=None else {}

    logging.debug('searching norms in class: %s' % pdf.__class__.__name__ )

    if isinstance(pdf,ROOT.RooSimultaneous):
        cat = pdf.indexCat()
        idx = cat.getIndex()
        for i in xrange(cat.numBins('')):
            cat.setBin(i)
            pdfi = pdf.getPdf(cat.getLabel());
            if pdfi.__nonzero__(): getnorms(pdfi, obs, out);
        # restore the old index
        cat.setIndex(idx)
        #pass


    if isinstance(pdf,ROOT.RooProdPdf):
        pdfs = ROOT.RooArgList(pdf.pdfList())
        for pdfi in roofiter(pdfs):
            if pdfi.dependsOn(obs): getnorms(pdfi,obs,out)

    if isinstance(pdf,ROOT.RooAddPdf):
        coefs = ROOT.RooArgList(pdf.coefList())
        for c in roofiter(coefs):
            out[c.GetName()] =  c.getVal(obs)

    return out

def printstats( dcpath ):
    wspath = os.path.splitext(dcpath)[0]+'.root'


    options = dcoptions()
    options.stat = False
    options.bin = True
    options.noJMax = False
    options.nuisancesToExclude = []
    options.nuisancesToRescale = []

    options.fileName = dcpath
    options.out = None
    options.cexpr = False
    options.fixpars = False
    options.libs = []
    options.verbose = 0
    options.poisson = 0
    options.mass = 250

    dcfile = open(dcpath,'r')
    DC = parseCard(dcfile, options)
    print DC.bins
    
    f = ROOT.TFile.Open(wspath)
    w = f.Get('w')
    print w
    
    print '-'*80
    c = w.cat("CMS_channel")
    c.Print("V")


    print '-'*80
    model_s = w.pdf('model_s')
    model_s.Print('V')
    print '-'*80
    data = w.data('data_obs')

    return DC,w,model_s,data,c


if __name__ == '__main__':
    print 'Test!'
    ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')

    
    
    dc,ws,model_s,data,cat = printstats('datacards/hww-24.41fb.mH250.comb_0j1j_shape.txt')
#     dc,ws,model_s,data,cat = printstats('datacards/hww-19.47fb.mH125.comb_of_shape.txt')
#     printstats('datacards/hww-19.47fb.mH250.of_0j_shape.txt')
