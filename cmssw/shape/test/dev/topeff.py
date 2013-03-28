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

    wwflow = CutFlow(wwcuts.wwcommon)

    del wwflow['bveto_mu']
    del wwflow['bveto_ip']
    wwflow['ptll'] = 'ptll>45'
    
    print '------------------------------'
    for n,c in wwflow.iteritems():
        print '%-30s: %s'% (n,c)

    topflow = CutFlow()

    topflow['base']  = wwflow.string()
    topflow['bctrl'] = 'jettche1>2.1 && njet>=1'
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

    # ---
    # variables
    vars = {
#         'mll'  :  ('mll',      (100 , 0   , 600), 'm_{ll} [GeV]' ),
#         'tche2':  ('jettche2', (40  , -20 , 20) , 'TCHE_{j2}' ),
#         'jetpt2': ('jetpt2',   (40,10,200), 'pt_{j2}' ),
        'jetpt2': ('jetpt2',  ( range(10,100,5)+range(100,201,20),) , 'pt_{j2}' ),
        'jeteta2': ('fabs(jeteta2)',   (12,0,3), 'eta_{j2}' ),
    }
   
    # prepare 
    plots = AlienDict()

    for v,(expr,bins,xaxis) in vars.iteritems():
        for n,a in analysers.iteritems():
            print v,':',n,

            pf = a.plotsflow(n+'_'+v,expr,bins=bins)
            print '...done'

            for c,h in pf.iteritems():
                plots[v][c][n] = h
                
    prefix=''
    if opt.out:
        hwwtools.ensuredir(opt.out)
        prefix = opt.out+'/'

    
#     for v,(expr,bins,xaxis) in vars.iteritems():
#         hwwlatino.printplots(plots[v]['bctrl'],prefix+'plot_bctrl_%s' % v,xaxis=xaxis, label='8TeV', lumi=opt.lumi, exts=['pdf','png'])
#         hwwlatino.printplots(plots[v]['btag'] ,prefix+'plot_btag_%s'  % v,xaxis=xaxis, label='8TeV', lumi=opt.lumi, exts=['pdf','png'])


    from ginger.plotter import H1RatioPlotter
    others = ['WW','ggWW','WJet','DYLL','DYTT','VV','Vg','VgS']

    # ---
    # by ptj2
    ptj2_bctrl = plots['jetpt2']['bctrl']
    ptj2_btag  = plots['jetpt2']['btag']

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


    hratio_pt2j = H1RatioPlotter()
    hratio_pt2j.scalemax = 1.3
    hratio_pt2j.ltitle = 'top tag efficiency'
    hratio_pt2j.rtitle = ''
    hratio_pt2j.ytitle2 = 'data/mc'
    hratio_pt2j.set(heff_mc_ptj2, heff_ds_ptj2)
    c = hratio_pt2j.draw()
    
    c.Print(prefix+'heff_ptj2.pdf')
    c.Print(prefix+'heff_ptj2.png')

    
    # ---
    # by etaj2
    etaj2_bctrl = plots['jeteta2']['bctrl']
    etaj2_btag  = plots['jeteta2']['btag']

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


    hratio_eta2j = H1RatioPlotter()
    hratio_eta2j.scalemax = 1.3
    hratio_eta2j.ltitle = 'top tag efficiency'
    hratio_eta2j.rtitle = ''
    hratio_eta2j.ytitle2 = 'data/mc'
    hratio_eta2j.legalign = 'r'
    hratio_eta2j.set(heff_mc_etaj2, heff_ds_etaj2)
    c = hratio_eta2j.draw()
    
    c.Print(prefix+'heff_etaj2.pdf')
    c.Print(prefix+'heff_etaj2.png')



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
        print 'GOODBYE!'

    try:
        __IPYTHON__
    except NameError:
        print 'Cleaning up'
  
