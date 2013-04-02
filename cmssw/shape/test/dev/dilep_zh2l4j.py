#!/usr/bin/env python
from ginger.analysis import TreeAnalyser,Cut,CutFlow,Latino,Sample
import HWWAnalysis.Misc.odict as odict

orchard = '/shome/mtakahashi/HWW/Tree/ShapeAna/53x_195fb/tree_skim_wwmin/'
orchard = '/shome/thea/HWW/work/shapeMoriond/trees/dileptons'


order=[
    'ggH',  
    'vbfH', 
    'wzttH',
    'VH',   
    'wH',   
    'zH',   


    'WW',   
    'ggWW', 
    'VV',   
    'WJet', 
    'Top',  
    'Vg',   
    'VgS',  
    'DYTT', 
    'DYLL', 
    'DYee', 
    'DYmm', 
]
# ---
def main( opt ):
    print opt
    import hwwlatino
    samples = hwwlatino.samples(125,'8TeV','Data2012','SM','zh4j_mm')

    samples.pop('VgS-template')
    samples.pop('Vg-template')

    cutflow = CutFlow()
    cutflow['os_sf']  = '(ch1*ch2)<0. && sameflav & trigger == 1' # don't use 0
    cutflow['pt12']   = 'pt1>20 && pt2 > 10'
    cutflow['mllmin'] = 'mll > 12'
    cutflow['met']    = 'pfmet < 20'
    cutflow['nextra'] = 'nextra == 0'
    cutflow['njet4']  = 'jetpt3>20'
    
    for n,s in samples.iteritems():
        print '%-10s'%n,s
    del samples['VVV']
    analysers = hwwlatino.makeanalysers(samples,orchard,cutflow,opt.lumi)

    print '-'*80

    for n,a in analysers.iteritems():
        old = a.worker.entries()
        a.bufferentries()
        print '  ',n,':',old,'>>',a.selectedentries(),'...done'
    

    print '-'*80

    plots_mm_mll = {}
    plots_ee_mll = {}
    plots_mm = odict.OrderedDict(zip(cutflow.keys(),[{} for _ in cutflow]))
    plots_ee = odict.OrderedDict(zip(cutflow.keys(),[{} for _ in cutflow]))

    for n,a in analysers.iteritems():
        print n,
        pf_mm = a.plotsflow(n,'mll','channel==0',bins=(100,0,500))
        pf_ee = a.plotsflow(n,'mll','channel==1',bins=(100,0,500))
        print '...done'
            
        # re-order
        for c in cutflow:
            plots_mm[c][n] = pf_mm[c]
            plots_ee[c][n] = pf_ee[c]
        

    prefix=''
    if opt.out:
        hwwtools.ensuredir(opt.out)
        prefix = opt.out+'/'

    for i,c in enumerate(cutflow):
        hwwlatino.printplots(plots_mm[c],prefix+'plot_%d_mm_mll_%s' % (i,c),xaxis='m_{ll} [GeV]', lumi=opt.lumi, logy=True, exts=['pdf','png'], order=order)
        hwwlatino.printplots(plots_ee[c],prefix+'plot_%d_ee_mll_%s' % (i,c),xaxis='m_{ll} [GeV]', lumi=opt.lumi, logy=True, exts=['pdf','png'], order=order)

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
  
