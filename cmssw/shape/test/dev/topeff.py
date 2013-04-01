#!/usr/bin/env python
from ginger.analysis import TreeAnalyser,Cut,CutFlow,Latino,Sample
import HWWAnalysis.Misc.odict as odict
from hwwinfo2g import wwnamedcuts as wwcuts
from ginger.painter import Canvas,Pad,Legend

# orchard = '/shome/mtakahashi/HWW/Tree/ShapeAna/53x_195fb/tree_skim_wwmin/'
# orchard = '/shome/thea/HWW/work/shapeMoriond/trees/dileptons'
orchard = '/shome/thea/HWW/work/dds/trees/top'

class AlienDict(dict):
    """Implementation of perl's autovivification feature."""
    def __init__(self,*args, **kwargs):
        # init the dict
        super(self.__class__,self).__init__(self, *args, **kwargs)
        self._lock = False
    
    def lock(self):
        self._lock = True
        for a in self.itervalues():
            if type(a) == type(self):
                a.lock()

    def unlock(self):
        self._lock = False
        for a in self.itervalues():
            if type(a) == type(self):
                a.unlock()

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            if self._lock:
                raise
            else:
                value = self[item] = type(self)()
                return value

# ---
def main( opt ):
    import hwwlatino
    samples = hwwlatino.samples(125,'8TeV','Data2012','SM','topestimate')
#     samples = hwwlatino.samples(125,'8TeV','Data2012','SM','topplots')

    wwflow = CutFlow(wwcuts.wwcommon)

    del wwflow['bveto_mu']
    del wwflow['bveto_ip']
    wwflow['ptll'] = 'ptll>45'
    
    print '------------------------------'
    for n,c in wwflow.iteritems():
        print '%-30s: %s'% (n,c)

    topflow = CutFlow()

    topflow['base']  = wwflow.string()
    topflow['bctrl'] = 'jettche1>2.1 && bveto_mu'
    topflow['btag']  = 'jettche2>2.1'

    print '-'*80

    for n,s in samples.iteritems():
        print '%-10s'%n,s
    analysers = hwwlatino.makeanalysers(samples,orchard,topflow,opt.lumi)

    for n,a in analysers.iteritems():
        old = a.worker.entries()
        a.bufferentries()
        print '  ',n,':',old,'>>',a.selectedentries(),'...done'
    

    print '-'*80
    imgext = ['pdf','png']

    # ---
    # variables
    jetptbins = range(10,30,2)+range(30, 100, 5)+range(100, 201, 20)

    # block A: control plots for the standard estimate
    vars = {
#         'mll_1j'        : ('mll'           , 'njet >= 1', (100,  0, 600) , 'm_{ll} [GeV]' ),
#         'tche2_1j'      : ('jettche2'      , 'njet >= 1', ( 40,-20,  20) , 'TCHE_{j2}' ),
        'jetpt2_1j'     : ('jetpt2'        , 'njet >= 1', ( jetptbins, ) , 'pt_{j2}' ) ,
        'jeteta2_1j'    : ('fabs(jeteta2)' , 'njet == 1', ( 12,  0,  3 ) , 'eta_{j2}' ),

    }

    prefix=''
    if opt.out:
        hwwtools.ensuredir(opt.out)
        prefix = opt.out+'/'
   
    # prepare 
    bplots = AlienDict()

    for v,(expr,cut,bins,xaxis) in vars.iteritems():
        for n,a in analysers.iteritems():
            print v,':',n,

            pf = a.plotsflow(n+'_'+v,expr,extra=cut,bins=bins)
            print '...done'

            for c,h in pf.iteritems():
                bplots[v][c][n] = h


    
    for v,(expr,cut,bins,xaxis) in vars.iteritems():
        print 'Printing bplots:',v
        hwwlatino.printplots(bplots[v]['bctrl'],prefix+'bplot_bctrl_%s' % v,xaxis=xaxis, label='b_{CTRL}, %s' % cut, lumi=opt.lumi, exts=imgext)
        hwwlatino.printplots(bplots[v]['btag'] ,prefix+'bplot_btag_%s'  % v,xaxis=xaxis, label='b_{TAG} , %s' % cut, lumi=opt.lumi, exts=imgext)

    # block B: control shapes
    vars = {
#         'njet'            : ('njet'                                               , 'jettche1 > 2.1', (10,  0,10), 'n^{jets}'),
#         'tche1'           : ('jettche1'                                           , ''              , (40,-20,20), 'TCHE_{j1}' ),
#         'softjetpt_2j'    : ('jetpt3'                                             , 'jettche1 > 2.1 && njet == 2' , (10,10,30) , 'pt_{j3}' )                     ,
#         'softjetpt_3j'    : ('jetpt4'                                             , 'jettche1 > 2.1 && njet == 3' , (10,10,30) , 'pt_{j4}' )                     ,

        'softjet_1j_btotal_eta2' : ('fabs(jeteta2)' , 'bveto_mu && njet == 1 && jettche1 > 2.1'                      , ( 12,  0,  3 ) , 'eta_{j2}' ),
        'softjet_1j_btag_eta2'   : ('fabs(jeteta2)' , 'bveto_mu && njet == 1 && jettche1 > 2.1 && jettche2 >  2.1'   , ( 12,  0,  3 ) , 'eta_{j2}' ),
        'softjet_0j_btotal_eta1' : ('fabs(jeteta1)' , 'bveto_mu && njet == 0'                                        , ( 12,  0,  3 ) , 'eta_{j2}' ),
        'softjet_0j_btag_eta1'   : ('fabs(jeteta1)' , 'bveto_mu && njet == 0 && jettche1 >  2.1'                     , ( 12,  0,  3 ) , 'eta_{j2}' ),

        'softjet_1j_btotal_pt'   : ('jetpt2'        , 'bveto_mu && njet == 1 && jettche1 > 2.1'                      , ( 10, 10, 30 ) , 'pt_{j2}' ) ,
        'softjet_1j_btag_pt'     : ('jetpt2'        , 'bveto_mu && njet == 1 && jettche1 > 2.1 && jettche2 >  2.1'   , ( 10, 10, 30 ) , 'pt_{j2}' ) ,
#         'softjet_1j_bveto_pt'    : ('jetpt2'        , 'bveto_mu && njet == 1 && jettche1 > 2.1 && jettche2 <= 2.1'   , ( 10, 10, 30 ) , 'pt_{j2}' ) ,
        'softjet_0j_btotal_pt'   : ('jetpt1'        , 'bveto_mu && njet == 0'                                        , ( 10, 10, 30 ) , 'pt_{j1}' ) ,
        'softjet_0j_btag_pt'     : ('jetpt1'        , 'bveto_mu && njet == 0 && jettche1 >  2.1'                     , ( 10, 10, 30 ) , 'pt_{j1}' ) ,
#         'softjet_0j_bveto_pt'    : ('jetpt1'        , 'bveto_mu && njet == 0 && jettche1 <= 2.1'                     , ( 10, 10, 30 ) , 'pt_{j1}' ) ,
    }

    xplots = AlienDict()
    for v,(expr,cut,bins,xaxis) in vars.iteritems():
        for n,a in analysers.iteritems():
            print v,':',n,
            xplots[v][n] = a.views['base'].plot('xplots_%s_%s' % (n,v), expr,extra=cut,bins=bins)
            print '...done'
        hwwlatino.printplots(xplots[v],prefix+'xplots_base_%s' % v, xaxis=xaxis, label='base, %s' % cut, lumi=opt.lumi, exts=imgext)

        
    if 'Top' not in analysers:
        print 'No top, no party'
        return

    from ginger.plotter import H1RatioPlotter
    others = ['WW','ggWW','WJet','DYLL','DYTT','VV','Vg','VgS']

    colors  = [ROOT.kRed+1      , ROOT.kAzure-5   ]
    markers = [ROOT.kFullCircle , ROOT.kFullCircle]

    # ---
    # by ptj2
    if 'jetpt2_1j' in bplots:
        ptj2_bctrl = bplots['jetpt2_1j']['bctrl']
        ptj2_btag  = bplots['jetpt2_1j']['btag']

        ptj2_bctrl_ds = ptj2_bctrl['Data'].Clone('ptj2_bctrl_ds')
        for p in others:  ptj2_bctrl_ds -= ptj2_bctrl[p]

        ptj2_btag_ds  = ptj2_btag['Data'].Clone('ptj2_btag_ds')
        for p in others:  ptj2_btag_ds  -= ptj2_btag[p]

        heff_ds_ptj2 = ptj2_btag_ds.Clone('heff_ds_ptj2')
        heff_ds_ptj2.Divide(ptj2_btag_ds,ptj2_bctrl_ds,1,1,'b')
        heff_ds_ptj2.SetTitle('data - mc_{others}')

        heff_mc_ptj2 = ptj2_btag['Top'].Clone('heff_mc_ptj2')
        heff_mc_ptj2.Divide(ptj2_btag['Top'],ptj2_bctrl['Top'],1,1,'b')
        heff_mc_ptj2.SetTitle('mc (tW/t#bar{t});pt_{j2};tag/ctrl')


        hratio_pt2j = H1RatioPlotter(colors=colors,markers=markers)
        hratio_pt2j.scalemax = 1.1
        hratio_pt2j.legalign = ('r','b')
        hratio_pt2j.ltitle = 'top tag efficiency'
        hratio_pt2j.rtitle = 'L=%.2f fb^{-1}' % opt.lumi
        hratio_pt2j.ytitle2 = 'data/mc'
        hratio_pt2j.set(heff_mc_ptj2, heff_ds_ptj2)
        c = hratio_pt2j.plot()
        
        for ext in imgext:
            c.Print(prefix+'heff_ptj2.'+ext)
        
        del c
        hratio_pt2j.userrange = (10.,30-0.01)
        hratio_pt2j.scalemax  = 0.6
        c = hratio_pt2j.plot()

        for ext in imgext:
            c.Print(prefix+'heff_ptj2_zoom.'+ext)
         

        trans_ds_pt2j  = heff_ds_ptj2.Clone('htrans_ds_ptj2')
        trans_mc_pt2j  = heff_mc_ptj2.Clone('htrans_mc_ptj2')

        for i in xrange(0,trans_mc_pt2j.GetNbinsX()+2):
            e  = trans_ds_pt2j.GetBinContent(i)
            ee = trans_ds_pt2j.GetBinError(i)
            trans_ds_pt2j.SetBinContent(i,(1-e)/e)
            trans_ds_pt2j.SetBinError(i,ee/e**2)
            e = trans_mc_pt2j.GetBinContent(i)
            ee = trans_ds_pt2j.GetBinError(i)
            trans_mc_pt2j.SetBinContent(i,(1-e)/e)
            trans_mc_pt2j.SetBinError(i,ee/e**2)

        htrans_pt2j = H1RatioPlotter(colors=colors,markers=markers)
        htrans_pt2j.scalemax = 2
        htrans_pt2j.legalign = ('r','t')
        htrans_pt2j.ltitle = 'transfer factor'
        hratio_pt2j.rtitle = 'L=%.2f fb^{-1}' % opt.lumi
        htrans_pt2j.ytitle2 = 'data/mc'
        htrans_pt2j.set(trans_mc_pt2j, trans_ds_pt2j)
        htrans_pt2j.userrange = (10.,30-0.01)

        c = htrans_pt2j.plot()
        for ext in imgext:
            c.Print(prefix+'htrans_ptj2_zoom.'+ext)
    # ---
    # by etaj2
    if 'jeteta2_1j' in bplots:
        etaj2_bctrl = bplots['jeteta2_1j']['bctrl']
        etaj2_btag  = bplots['jeteta2_1j']['btag']

        etaj2_bctrl_ds = etaj2_bctrl['Data'].Clone('etaj2_bctrl_ds')
        for p in others:  etaj2_bctrl_ds -= etaj2_bctrl[p]

        etaj2_btag_ds  = etaj2_btag['Data'].Clone('etaj2_btag_ds')
        for p in others:  etaj2_btag_ds  -= etaj2_btag[p]

        heff_ds_etaj2 = etaj2_btag_ds.Clone('heff_ds_etaj2')
        heff_ds_etaj2.Divide(etaj2_btag_ds, etaj2_bctrl_ds,1,1,'b')
        heff_ds_etaj2.SetTitle('data - mc_{others}')

        heff_mc_etaj2 = etaj2_btag['Top'].Clone('heff_mc_etaj2')
        heff_mc_etaj2.Divide(etaj2_btag['Top'],etaj2_bctrl['Top'],1,1,'b')
        heff_mc_etaj2.SetTitle('mc (tW/t#bar{t});eta_{j2};tag/ctrl')


        hratio_eta2j = H1RatioPlotter(colors=colors,markers=markers)
        hratio_eta2j.scalemax = 1.2
        hratio_eta2j.ltitle = 'top tag efficiency'
        hratio_eta2j.rtitle = 'L=%.2f fb^{-1}' % opt.lumi
        hratio_eta2j.ytitle2 = 'data/mc'
        hratio_eta2j.legalign = ('r','t')
        hratio_eta2j.set(heff_mc_etaj2, heff_ds_etaj2)
        c = hratio_eta2j.plot()
        
        for ext in imgext:
            c.Print(prefix+'heff_etaj2.'+ext)

    if 'softjet_0j_btag_pt' in xplots and 'softjet_1j_btag_pt' in xplots: 

        ptj2_btag  = xplots['softjet_1j_btag_pt']
        ptj1_btag  = xplots['softjet_0j_btag_pt']

        ptj2_btag_mc  = ptj2_btag['Top']
        ptj2_btag_ds  = ptj2_btag['Data'].Clone('ptj2_btag_ds')
        for p in others:  ptj2_btag_ds  -= ptj2_btag[p]

        ptj1_btag_mc  = ptj1_btag['Top']
        ptj1_btag_ds  = ptj1_btag['Data'].Clone('ptj1_btag_ds')
        for p in others:  ptj1_btag_ds  -= ptj1_btag[p]

        # btag region
        ptj2_btag_ds.Scale(1./ptj2_btag_ds.Integral())
        ptj2_btag_mc.Scale(1./ptj2_btag_mc.Integral())
        ptj1_btag_ds.Scale(1./ptj1_btag_ds.Integral())
        ptj1_btag_mc.Scale(1./ptj1_btag_mc.Integral())

        ptj2_btag_mc.SetTitle('mc_{top} n_{jet}=1;pt_{j-soft}')
        ptj1_btag_mc.SetTitle('mc_{top} n_{jet}=0;pt_{j-soft}')
        ptj2_btag_ds.SetTitle('data-mc_{others} n_{jet}=1;pt_{j-soft}')
        ptj1_btag_ds.SetTitle('data-mc_{others} n_{jet}=0;pt_{j-soft}')

        hratio_0j1j_btag = H1RatioPlotter(colors=colors,markers=markers)
#         hratio_0j1j_btag.set(ptj2_btag_mc,ptj1_btag_mc,ptj2_btag_ds,ptj1_btag_ds)
        hratio_0j1j_btag.set(ptj2_btag_ds,ptj1_btag_ds)
        hratio_0j1j_btag.ltitle  = '"hardest" pt_{j-soft}'
        hratio_0j1j_btag.rtitle  = '0j vs 1j'
        hratio_0j1j_btag.ytitle2 = 'x/mc_{top} n_{jet}=1'
        hratio_0j1j_btag.legtextsize = 25
        hratio_0j1j_btag.legboxsize  = 30
        hratio_0j1j_btag.scalemax = 1.2
        hratio_0j1j_btag.markersize = 16

        c = hratio_0j1j_btag.plot()
        for ext in imgext:
            c.Print(prefix+'ratio_softjet_pt_0j1j_btag.'+ext)

    if 'softjet_0j_btag_eta1' in xplots and 'softjet_1j_btag_eta2' in xplots: 

        etaj2_btag  = xplots['softjet_1j_btag_eta2']
        etaj1_btag  = xplots['softjet_0j_btag_eta1']

        etaj2_btag_mc  = etaj2_btag['Top']
        etaj2_btag_ds  = etaj2_btag['Data'].Clone('etaj2_btag_ds')
        for p in others:  etaj2_btag_ds  -= etaj2_btag[p]

        etaj1_btag_mc  = etaj1_btag['Top']
        etaj1_btag_ds  = etaj1_btag['Data'].Clone('etaj1_btag_ds')
        for p in others:  etaj1_btag_ds  -= etaj1_btag[p]

        # btag region
        etaj2_btag_ds.Scale(1./etaj2_btag_ds.Integral())
        etaj2_btag_mc.Scale(1./etaj2_btag_mc.Integral())
        etaj1_btag_ds.Scale(1./etaj1_btag_ds.Integral())
        etaj1_btag_mc.Scale(1./etaj1_btag_mc.Integral())

        etaj2_btag_mc.SetTitle('mc_{top} n_{jet}=1;#eta_{j-soft}')
        etaj1_btag_mc.SetTitle('mc_{top} n_{jet}=0;#eta_{j-soft}')
        etaj2_btag_ds.SetTitle('data-mc_{others} n_{jet}=1;#eta_{j-soft}')
        etaj1_btag_ds.SetTitle('data-mc_{others} n_{jet}=0;#eta_{j-soft}')

        hratio_0j1j_btag = H1RatioPlotter(colors=colors,markers=markers)
        hratio_0j1j_btag.set(etaj2_btag_ds,etaj1_btag_ds)
        hratio_0j1j_btag.ltitle  = '"hardest" #eta_{j-soft}'
        hratio_0j1j_btag.rtitle  = '0j vs 1j'
        hratio_0j1j_btag.ytitle2 = 'x/mc_{top} n_{jet}=1'
        hratio_0j1j_btag.legtextsize = 25
        hratio_0j1j_btag.legboxsize  = 30
        hratio_0j1j_btag.scalemax    = 1.2
        hratio_0j1j_btag.markersize  = 16
        hratio_0j1j_btag.legalign    = ('r','t')

        c = hratio_0j1j_btag.plot()
        for ext in imgext:
            c.Print(prefix+'ratio_softjet_eta_0j1j_btag.'+ext)


# ---
if __name__ == '__main__':
    import optparse
    import hwwtools
    import sys
    import bdb

    
    parser = optparse.OptionParser()
    parser.add_option('-d', '--debug'    , dest='debug'       , help='Debug level'            , default=0 , action='count' )
    parser.add_option('-l', '--lumi'     , dest='lumi'        , help='Luminosity'             , default=19.468 )
    parser.add_option('-o', '--out'      , dest='out'         , help='Output'                 , default=None )
    

    (opt, args) = parser.parse_args()

    hwwtools.setDebugLevel(opt)

    import os.path
    import ROOT
    shape_path = os.path.join(os.getenv('CMSSW_BASE'),'src/HWWAnalysis/ShapeAnalysis')
    print 'Shape directory is',shape_path
    ROOT.gInterpreter.ExecuteMacro(shape_path+'/macros/LatinoStyle2.C')

    try:
        main( opt )
    except SystemExit:
        pass
    except bdb.BdbQuit:
        pass
    except:
        import traceback
        exc_type, exc_value, exc_traceback = sys.exc_info()

        print '-'*80
        print '--> Exception:'
        formatted_lines = traceback.format_exc().splitlines()
        print '   ',formatted_lines[-1]
        print '-'*80
        traceback.print_exc()
        print '-'*80
    finally:
        print 'over and out!'

    try:
        __IPYTHON__
    except NameError:
        print 'Cleaning up'
  
