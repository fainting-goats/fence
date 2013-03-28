#!/bin/env python



import sys
import os.path

from HiggsAnalysis.CombinedLimit.DatacardParser import *
from HWWAnalysis.Misc.ROOTAndUtils import TH1AddDirSentry

import ROOT


mlfdir = '.'

def RooList2Set( alist ):

    next =  alist.createIterator()
    aset = ROOT.RooArgSet()
    while True: 
        p = next.Next()
        if not p: break
        aset.add( p )
    return aset

class ShapeGluer:
    def __init__(self, bin, DC, ws, fit=None): 
        self._DC = DC
        self._bin = bin
        self._ws = ws
        self._fit = fit
        
        self._build()

    def _build(self):
        
        if self._bin in self._DC.shapeMap:  bintag = self._bin
        elif '*' in self._DC.shapeMap:     bintag = '*'
        else:
            raise ValueError('Couldn\'t find '+bin+' or * in shapeMap')

        try:
            hpath, hname =  self._DC.shapeMap[bintag]['data_obs']
        except KeyError as e:
            raise KeyError('Shape for '+str(e) +'not found!') 

        sentry = TH1AddDirSentry()
        hpath = os.path.join(os.path.dirname(dcpath),hpath)
        hfile = ROOT.TFile.Open(os.path.join(hpath))
        if not hfile.__nonzero__():
            raise IOError('Could not open '+wspath)

        hdata = hfile.Get(hname)

        # the Xaxis label has to be cheked
        self._template = hdata.Clone('shape_template')
        self._template.SetTitle('shape_template')
        self._template.Reset()

    
    def glue(self):
        exp = self._DC.exp[self._bin]
        shapes = dict([ (p,self._glueprocess(p)) for p in self._DC.processes if exp[p] != 0])
        shapes['Data'] = self._gluedata()

        return shapes

    def _makeHisto(self, name, title):

        sentry = TH1AddDirSentry()
        h = self._template.Clone(name)
        h.SetTitle(title)

        return h
        


    def _gluedata(self):
        from math import sqrt
        h = self._makeHisto('histo_Data','Data')

        data = self._ws.data('data_obs')
        for i in xrange(data.numEntries()):
            data.get(i)
            h.SetBinContent(i+1,data.weight())
            h.SetBinError(i+1,sqrt(data.weight()))
        
        return h


    def _glueprocess(self, process):

        data = self._ws.data('data_obs')
        tag = 'Sig' if process in self._DC.signals else 'Bkg'

#         self._ws.allPdfs().printMultiline(ROOT.std.cout,3)
        mname = 'shape{0}_{1}_{2}_morph'.format(tag,self._bin,process)
        sname = 'shape{0}_{2}_{1}Pdf'.format(tag,self._bin,process)
        morph  = self._ws.pdf('shape{0}_{1}_{2}_morph'.format(tag,self._bin,process))
        static = self._ws.pdf('shape{0}_{2}_{1}Pdf'.format(tag,self._bin,process))
        errs = self._ws.pdf('shape{0}_{2}_{1}_CMS_hww_{2}_{1}_stat_shapeUpPdf'.format(tag,self.bin,process))
#         for the errors one needs to take the 2 original histograms
#  'shape{0}_{2}_{1}.format(tag,self.bin,process)
#  'shape{0}_{2}_{1}_CMS_hww_{2}_{1}_stat_shapeUpPdf'.format(tag,self.bin,process)

        if morph.__nonzero__():
            shape = morph
        elif static.__nonzero__():
            shape = static
        else:
            self._ws.allPdfs().printMultiline(ROOT.std.cout,3)
            print morph.__nonzero__(),morph, mname
            print static.__nonzero__(),static, sname
            raise ValueError('Can\'t find the nether the morph nor the shape!!! '+process)

        h = self._makeHisto('histo_'+process, process)

        if self._fit:
            norms, res = self._fit
            norm = norms.find('n_exp_bin{0}_proc_{1}'.format(self._bin, process))
            ShapeGluer._rooPdf2TH1(h,shape,data, norm, self._fit[1].floatParsFinal())
        else:
            ShapeGluer._rooPdf2TH1(h,shape,data, self._DC.exp[self._bin][process])

        return h


    @staticmethod
    def _rooPdf2TH1(h, pdf, data, norm=None, pars=None):

        pdf_obs  = pdf.getObservables(data)
        pdf_pars = pdf.getParameters(data)

        if h.GetNbinsX() != data.numEntries():
            raise ValueError('bins mismatch!')
        
        if pars:
            if isinstance(pars,ROOT.RooArgList):
                pars = RooList2Set(pars)

            pdf_pars.__assign__(pars)

        for i in xrange(data.numEntries()):
            pdf_obs.__assign__(data.get(i))
            h.SetBinContent(i+1,pdf.getVal())

        if norm:
            if isinstance(norm,ROOT.RooAbsReal):
                norm = norm.getVal()

            h.Scale(norm/h.Integral())


def THSum( shapes, labels, name=None, title=None):

    labels = [ p for p in shapes.iterkeys() if p in labels ]
    if not labels: return None

    sentry = TH1AddDirSentry()
    h = shapes[labels[0]].Clone()

    for l in labels[1:]:
        h.Add(shapes[l])

    if name: h.SetName(name)
    if title: h.SetName(title)

    return h


def main( dcpath, opts ):
    '''
1. read the datacard
2. convert to ws
3. run combine
4. open get the mlfit rootfile

1-4 don't need to know the content of the card, only to check that there are shapes inside

'''

    mypath = os.path.dirname(os.path.abspath(__file__))
    ROOT.gInterpreter.ExecuteMacro(mypath+'/LatinoStyle2.C')
    try:
        ROOT.gROOT.LoadMacro(mypath+'/MWLPlot.C+g')
    except RuntimeError:
        ROOT.gROOT.LoadMacro(mypath+'/MWLPlot.C++g')


    # 1. load the datacard
    dcfile = open(dcpath,'r')
    class dummy: pass

    options = dummy()
    options.stat = False
    options.bin = True
    options.noJMax = False
    options.nuisancesToExclude = []
    options.nuisancesToRescale = []

    DC = parseCard(dcfile, options)

    if not DC.hasShapes:
        sys.exit(-1)

    if len(DC.bins) != 1:
        raise ValueError('Only 1 bin datacards supported at the moment: '+', '.join(DC.bins))
        

    # 2. convert to ws
    wspath = os.path.splitext(dcpath)[0]+'.root'
    if not os.path.exists(wspath):
        # workspace + parameters = shapes
        print 'Making the workspace...',
        sys.stdout.flush()
        os.system('text2workspace.py '+dcpath)
        print 'done.'

    ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
    wsfile = ROOT.TFile.Open(wspath)
    if not wsfile.__nonzero__():
        raise IOError('Could not open '+wspath)
    
    w = wsfile.Get('w')

    # 3. run combine
    import tempfile
    mlfdir = tempfile.mkdtemp(prefix='mlfit_')
    print 'Fitting the workspace...',
    sys.stdout.flush()
    os.system('combine -M MaxLikelihoodFit --out ' + mlfdir + 
              ' --saveNormalizations '+wspath)
    print 'done.'

    # open the output and get the normalizations
    mlfpath = os.path.join(mlfdir,'mlfit.root')
    mlffile = ROOT.TFile.Open(mlfpath)
    if not mlffile.__nonzero__():
        raise IOError('Could not open '+wspath)

    sig_fit = ( mlffile.Get('norm_fit_s'),  mlffile.Get('fit_s') )
    bkg_fit = ( mlffile.Get('norm_fit_b'),  mlffile.Get('fit_b') )

    print DC.bins
    mode = 'input' 

    bin = DC.bins[0]

    if mode=='sig':
        fit = sig_fit
    elif mode == 'bkg':
        fit = bkg_fit
    elif mode == 'input':
        fit = None
    else:
        raise ValueError('mode can be only sig, bkg or input')

    gluer = ShapeGluer(bin, DC, w, fit)

    shapes = gluer.glue()

    shapes['Hsum']  = THSum(shapes,['ggH','vbfH','wzttH'],'histo_higgs','higgs')
    shapes['WWsum'] = THSum(shapes,['WW','qqWW'],'histo_WWsum','WWsum')
    shapes['VVsum'] = THSum(shapes,['VV','Vg'],'histo_VVsum','VVsum')
    shapes['DYsum'] = THSum(shapes,['DYLL','DYTT'],'histo_DYsum','DYsum')

    plot = ROOT.MWLPlot()
    plot.setDataHist(shapes['Data'])
    if mode != 'bkg':
        plot.setHWWHist(shapes['Hsum'])

    plot.setWWHist(shapes['WWsum'])  
    plot.setZJetsHist(shapes['DYsum'])
    plot.setTopHist(shapes['Top'])
    plot.setVVHist(shapes['VVsum'])
    plot.setWJetsHist(shapes['WJet'])

    cName = 'c_test'
    ratio = False
    c = ROOT.TCanvas(cName,cName) if ratio else ROOT.TCanvas(cName,cName,2)
    plot.setMass(opt.mass)
    plot.setLumi(opt.lumi)
    plot.setLabel(opt.xlabel)
    plot.Draw(c,1,ratio)

    c.Print(cName+'_'+'_'.join([bin,mode])+'.pdf')


if __name__ == '__main__':
    import optparse
    import hwwtools
    ## option parser
    usage = 'usage: %prog [options] datacard'
    parser = optparse.OptionParser(usage)
    
    parser.add_option('-i','--input',dest='inputdir',help='Input dir')
    parser.add_option('-m','--mass',dest='mass',help='Mass',default=-1)
#     parser.add_option('-o','--output',dest='outputdir',help='Output dir',default='.')
#     parser.add_option('-v','--variations',dest='variations',help='make the scale up/down stacks',action='store_true',default=False)
    parser.add_option('-x','--xlabel',dest='xlabel',help='X-axis label',default='')
    parser.add_option('-r','--ratio',dest='ratio',help='Plot the data/mc ration', action='store_false',default=True)
    parser.add_option('-l','--lumi', dest='lumi', type='float', help='Luminosity', default=None)
#     parser.add_option('-f','--filter',dest='filter', help='Filter on the variations', default='*')
 

    hwwtools.loadOptDefaults(parser)
    (opt, args) = parser.parse_args()
    sys.argv.append('-b')

    try:
        dcpath = args[0]
    except IndexError:
        parser.print_usage()
    
    try:
        main(dcpath, opt)
    except:
        import traceback
        exc_type, exc_value, exc_traceback = sys.exc_info()
#         print "*** print_tb:"
#         traceback.print_tb(exc_traceback, limit=1, file=sys.stdout)
#         print "*** print_exception:"
#         traceback.print_exception(exc_type, exc_value, exc_traceback,
#                                   limit=2, file=sys.stdout)
        print "*** print_exc:"
        traceback.print_exc()
        print "*** format_exc, first and last line:"
        formatted_lines = traceback.format_exc().splitlines()
        print formatted_lines[0]
        print formatted_lines[-1]
        print "*** format_exception:"
        print repr(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))
#         print "*** extract_tb:"
#         print repr(traceback.extract_tb(exc_traceback))
#         print "*** format_tb:"
#         print repr(traceback.format_tb(exc_traceback))
#         print "*** tb_lineno:", exc_traceback.tb_lineno

