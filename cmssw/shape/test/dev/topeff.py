#!/usr/bin/env python
from ginger.analysis import TreeAnalyser,Cut,CutFlow,Latino,Sample
import HWWAnalysis.Misc.odict as odict
from hwwinfo2g import wwnamedcuts as wwcuts
from ginger.painter import Canvas,Pad,Legend
import sys
import hwwlatino
from ginger.plotter import H1RatioPlotter

# orchard = '/shome/mtakahashi/HWW/Tree/ShapeAna/53x_195fb/tree_skim_wwmin/'
# orchard = '/shome/thea/HWW/work/shapeMoriond/trees/dileptons'
orchard = '/shome/thea/HWW/work/dds/trees/top'

# useful lists
others = ['WW','ggWW','WJet','DYLL','DYTT','VV','Vg','VgS']
tops   = ['ttbar','tW']

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
def sumplots( plots, newname, processes ):

#     print plots
#     print newname
#     print processes
    if not processes: raise ValueError('Processes is emtty!')

    hnew = plots[processes[0]].Clone(newname)
    hnew.Reset()
    for p in processes: 
        hnew += plots[p]

    return hnew

class histcollector(odict.OrderedDict):

    def __add__(self,obj):
        if not isinstance(onj,ROOT.TH1): raise TypeError('I only eat TH1')

        self[obj.GetName()] = obj


#     ______________
#    / ____/ __/ __/
#   / __/ / /_/ /_  
#  / /___/ __/ __/  
# /_____/_/ /_/     
                  

# ---
def makeefficiency( name, bplots, x, xlabel, scalemax, legalign, lumi, imgext, prefix, save=None ):

    colors  = [ROOT.kRed+1      , ROOT.kAzure-5   ]
    markers = [ROOT.kFullCircle , ROOT.kFullCircle]

    x_bctrl = bplots[x]['bctrl']
    x_btag  = bplots[x]['btag']

    x_bctrl_ds  = x_bctrl['Data'].Clone('bctrl_%s_ds' % x)
    x_bctrl_ds -= sumplots(x_bctrl, 'bctrl_others', others)
#     for p in others:  x_bctrl_ds -= x_bctrl[p]

    x_btag_ds  = x_btag['Data'].Clone('btag_%s_ds' % x)
    x_btag_ds -= sumplots(x_btag, 'bctrl_others', others)
#     for p in others:  x_btag_ds  -= x_btag[p]

    heff_ds_x = x_btag_ds.Clone('heff_ds_%s' % x)
    heff_ds_x.Divide(x_btag_ds,x_bctrl_ds,1,1,'b')
    heff_ds_x.SetTitle('data - mc_{others}')

    x_bctrl_mc = sumplots(x_bctrl, 'bctrl_top_%s' % x, tops)
    x_btag_mc  = sumplots(x_btag,  'btag_top_%s'  % x, tops)

    heff_mc_x = x_bctrl_mc.Clone('heff_mc_%s' % x)
    heff_mc_x.Reset()
    heff_mc_x.Divide(x_btag_mc,x_bctrl_mc,1,1,'b')
    heff_mc_x.SetTitle('mc (tW/t#bar{t});%s;tag/ctrl' % xlabel)


    hratio = H1RatioPlotter(colors=colors,markers=markers)
    hratio.scalemax   = scalemax
    hratio.legalign   = legalign
    hratio.ltitle     = 'top tag efficiency' 
    hratio.rtitle     = 'L=%.2f fb^{-1}' % lumi
    hratio.ytitle2    = 'data/mc'
    hratio.markersize = 16
    hratio.set(heff_mc_x, heff_ds_x)
    c = hratio.plot()
    
    for ext in imgext:
        c.Print(prefix+name+'.'+ext )

    if save is not None:
        save[heff_ds_x.GetName()] = heff_ds_x
        save[heff_mc_x.GetName()] = heff_mc_x

    return hratio

        # to check
#         trans_ds_pt2j  = heff_ds_ptj2.Clone('htrans_ds_ptj2')
#         trans_mc_pt2j  = heff_mc_ptj2.Clone('htrans_mc_ptj2')

#         for i in xrange(0,trans_mc_pt2j.GetNbinsX()+2):
#             e  = trans_ds_pt2j.GetBinContent(i)
#             ee = trans_ds_pt2j.GetBinError(i)
#             trans_ds_pt2j.SetBinContent(i,(1-e)/e)
#             trans_ds_pt2j.SetBinError(i,ee/e**2)
#             e = trans_mc_pt2j.GetBinContent(i)
#             ee = trans_ds_pt2j.GetBinError(i)
#             trans_mc_pt2j.SetBinContent(i,(1-e)/e)
#             trans_mc_pt2j.SetBinError(i,ee/e**2)

#         htrans_pt2j = H1RatioPlotter(colors=colors,markers=markers)
#         htrans_pt2j.scalemax = 2
#         htrans_pt2j.legalign = ('r','t')
#         htrans_pt2j.ltitle = 'transfer factor'
#         hratio_pt2j.rtitle = 'L=%.2f fb^{-1}' % opt.lumi
#         htrans_pt2j.ytitle2 = 'data/mc'
#         htrans_pt2j.set(trans_mc_pt2j, trans_ds_pt2j)
#         htrans_pt2j.userrange = (10.,30-0.01)

#         c = htrans_pt2j.plot()
#         for ext in imgext:
#             c.Print(prefix+'htrans_ptj2_zoom.'+ext)


# ---
def doefficiencies( analysers, imgext, prefix, opt ):
    # ---
    # variables
    softjet_ptbins = range(10,  30, 2) + range( 30, 70, 5) + range( 70, 150, 20) + range(150,201, 50)
    jet_ptbins     =                     range( 30, 70, 5) + range( 70, 150, 20) + range(150,201, 50)
    

    # block A: control plots for the standard estimate
    vars = {
#         'mll_1j'        : ('mll'           , 'njet >= 1', (100,  0, 600) , 'm_{ll} [GeV]' ),
#         'tche2_1j'      : ('jettche2'      , 'njet >= 1', ( 40,-20,  20) , 'TCHE_{j2}' ),

        'nj_jetpt2'            : ('jetpt2'        , 'njet >= 1'                                                 , ( softjet_ptbins, ) , 'pt_{j2}' ) ,
        'nj_jetpt2-jeteta2-b0' : ('jetpt2'        , 'njet >= 1 && fabs(jeteta2) < 0.75'                         , ( softjet_ptbins, ) , 'pt_{j2}' ) ,
        'nj_jetpt2-jeteta2-b1' : ('jetpt2'        , 'njet >= 1 && fabs(jeteta2) >=0.75 && fabs(jeteta2)< 1.5  ' , ( softjet_ptbins, ) , 'pt_{j2}' ) ,
        'nj_jetpt2-jeteta2-b2' : ('jetpt2'        , 'njet >= 1 && fabs(jeteta2) >=1.5  && fabs(jeteta2)< 2.8'   , ( softjet_ptbins, ) , 'pt_{j2}' ) ,
        'nj_jetpt2-jeteta2-b3' : ('jetpt2'        , 'njet >= 1 && fabs(jeteta2) >=2.8  && fabs(jeteta2)< 5. '   , ( softjet_ptbins, ) , 'pt_{j2}' ) ,

        '1j_jeteta2'           : ('fabs(jeteta2)' , 'njet == 1' ,  ( 12,  0,  3 )      , '#eta_{j2}' ),
        '1j_jeteta2-bins'      : ('fabs(jeteta2)' , 'njet == 1' ,  ( [0.,0.75,1.5,2.8,5], ) , '#eta_{j2}' ),
        '1j_jetpt2-fine'       : ('jetpt2'        , 'njet == 1'                                                  , ( 10, 10, 30 ) , 'pt_{j2}' ) ,

        '1j_jetpt2'            : ('jetpt2'        , 'njet == 1'                                                            , (  5, 10, 30 ) , 'pt_{j2}' ) ,
        '1j_jetpt2-jeteta2-b0' : ('jetpt2'        , 'njet == 1 && fabs(jeteta2) < 0.75'                                    , (  5, 10, 30 ) , 'pt_{j2}' ) ,
        '1j_jetpt2-jeteta2-b1' : ('jetpt2'        , 'njet == 1 && fabs(jeteta2) >=0.75 && fabs(jeteta2) < 1.5  '           , (  5, 10, 30 ) , 'pt_{j2}' ) ,
        '1j_jetpt2-jeteta2-b2' : ('jetpt2'        , 'njet == 1 && fabs(jeteta2) >=1.5  && fabs(jeteta2) < 2.8'             , (  5, 10, 30 ) , 'pt_{j2}' ) ,
        '1j_jetpt2-jeteta2-b3' : ('jetpt2'        , 'njet == 1 && fabs(jeteta2) >=2.8  && fabs(jeteta2) < 5. '             , (  5, 10, 30 ) , 'pt_{j2}' ) ,

        '2j_jeteta2'           : ('fabs(jeteta2)' , 'bveto_ip && njet == 2'                                                , ( 12,  0,  3 )      , '#eta_{j2}' ),
        '2j_jeteta2-bins'      : ('fabs(jeteta2)' , 'bveto_ip && njet == 2'                                                , ( [0.,1.,2.,2.8, 5.], ) , '#eta_{j2}' ),

        '2j_jetpt2'            : ('jetpt2'        , 'bveto_ip && njet == 2'                                                , ( jet_ptbins, ) , 'pt_{j2}' ) ,
        '2j_jetpt2-jeteta2-b0' : ('jetpt2'        , 'bveto_ip && njet == 2 && fabs(jeteta2) < 0.75'                        , ( jet_ptbins, ) , 'pt_{j2}' ) ,
        '2j_jetpt2-jeteta2-b1' : ('jetpt2'        , 'bveto_ip && njet == 2 && fabs(jeteta2) >=0.75 && fabs(jeteta2) < 1.5' , ( jet_ptbins, ) , 'pt_{j2}' ) ,
        '2j_jetpt2-jeteta2-b2' : ('jetpt2'        , 'bveto_ip && njet == 2 && fabs(jeteta2) >=1.5  && fabs(jeteta2) < 2.8' , ( jet_ptbins, ) , 'pt_{j2}' ) ,
        '2j_jetpt2-jeteta2-b3' : ('jetpt2'        , 'bveto_ip && njet == 2 && fabs(jeteta2) >=2.8  && fabs(jeteta2) < 5. ' , ( jet_ptbins, ) , 'pt_{j2}' ) ,
    }

    # prepare 
    bplots = AlienDict()

    for v,(expr,cut,bins,xaxis) in vars.iteritems():
        print '%-30s: [' % v,
        for n,a in analysers.iteritems():
            print '%s,' % n,

            pf = a.plotsflow(n+'_'+v,expr,extra=cut,bins=bins)
            sys.stdout.flush()

            for c,h in pf.iteritems():
                bplots[v][c][n] = h
        print ']'

    bplots.lock()

    
    if opt.datamc:
        for v,(expr,cut,bins,xaxis) in vars.iteritems():
            print 'Printing bplots:',v
            hwwlatino.printplots(bplots[v]['bctrl'],prefix+'bplot_bctrl_%s' % v,xaxis=xaxis, label='b_{CTRL}, %s' % cut, lumi=opt.lumi, exts=imgext)
            hwwlatino.printplots(bplots[v]['btag'] ,prefix+'bplot_btag_%s'  % v,xaxis=xaxis, label='b_{TAG} , %s' % cut, lumi=opt.lumi, exts=imgext)

    heffs = odict.OrderedDict()

    # common paramete
    commons = {
        'lumi'   : opt.lumi,
        'imgext' : imgext,
        'prefix' : prefix,
        'save'   : heffs,
    }

#     print '\n'.join(bplots)
    makeefficiency( 'heff_inc_pt-full'        , bplots , 'nj_jetpt2'            ,  'pt_{j2}'   , 1.1 , ('r','b') , **commons )
    makeefficiency( 'heff_inc_pt-full-eta-b0' , bplots , 'nj_jetpt2-jeteta2-b0' ,  'pt_{j2}'   , 1.1 , ('r','b') , **commons )
    makeefficiency( 'heff_inc_pt-full-eta-b1' , bplots , 'nj_jetpt2-jeteta2-b1' ,  'pt_{j2}'   , 1.1 , ('r','b') , **commons )
    makeefficiency( 'heff_inc_pt-full-eta-b2' , bplots , 'nj_jetpt2-jeteta2-b2' ,  'pt_{j2}'   , 1.3 , ('r','b') , **commons )
    makeefficiency( 'heff_inc_pt-full-eta-b3' , bplots , 'nj_jetpt2-jeteta2-b3' ,  'pt_{j2}'   , 1.3 , ('r','b') , **commons )
                                                
    makeefficiency( 'heff_0j_pt-fine'         , bplots , '1j_jetpt2-fine'       ,  'pt_{j2}'   , 1.1 , ('l','t') , **commons )
    makeefficiency( 'heff_0j_eta'             , bplots , '1j_jeteta2'           ,  '#eta_{j2}' , 1.3 , ('r','t') , **commons )
    makeefficiency( 'heff_0j_eta-bin'         , bplots , '1j_jeteta2-bins'      ,  '#eta_{j2}' , 1.3 , ('r','t') , **commons )
                                                
    makeefficiency( 'heff_0j_pt'              , bplots , '1j_jetpt2'            ,  'pt_{j2}'   , 1.1 , ('l','t') , **commons )
    makeefficiency( 'heff_0j_pt-eta-b0'       , bplots , '1j_jetpt2-jeteta2-b0' ,  'pt_{j2}'   , 1.1 , ('l','t') , **commons )
    makeefficiency( 'heff_0j_pt-eta-b1'       , bplots , '1j_jetpt2-jeteta2-b1' ,  'pt_{j2}'   , 1.1 , ('l','t') , **commons )
    makeefficiency( 'heff_0j_pt-eta-b2'       , bplots , '1j_jetpt2-jeteta2-b2' ,  'pt_{j2}'   , 1.3 , ('l','t') , **commons )
    makeefficiency( 'heff_0j_pt-eta-b3'       , bplots , '1j_jetpt2-jeteta2-b3' ,  'pt_{j2}'   , 1.3 , ('l','t') , **commons )
                                                
    makeefficiency( 'heff_1j_eta'             , bplots , '2j_jeteta2'           ,  '#eta_{j2}' , 1.3 , ('r','t') , **commons )
    makeefficiency( 'heff_1j_eta-bin'         , bplots , '2j_jeteta2-bins'      ,  '#eta_{j2}' , 1.3 , ('r','t') , **commons )
                                                
    makeefficiency( 'heff_1j_pt'              , bplots , '2j_jetpt2'            ,  'pt_{j2}'   , 1.1 , ('r','b') , **commons )
    makeefficiency( 'heff_1j_pt-eta-b0'       , bplots , '2j_jetpt2-jeteta2-b0' ,  'pt_{j2}'   , 1.1 , ('r','b') , **commons )
    makeefficiency( 'heff_1j_pt-eta-b1'       , bplots , '2j_jetpt2-jeteta2-b1' ,  'pt_{j2}'   , 1.1 , ('r','b') , **commons )
    makeefficiency( 'heff_1j_pt-eta-b2'       , bplots , '2j_jetpt2-jeteta2-b2' ,  'pt_{j2}'   , 1.3 , ('r','b') , **commons )
    makeefficiency( 'heff_1j_pt-eta-b3'       , bplots , '2j_jetpt2-jeteta2-b3' ,  'pt_{j2}'   , 1.3 , ('r','b') , **commons )

    markers = [ROOT.kFullCircle , ROOT.kFullCircle , ROOT.kFullCircle , ROOT.kFullCircle]

    heff_fj = [ heffs['heff_ds_%s' % n] for n in [ 'nj_jetpt2','nj_jetpt2-jeteta2-b0','nj_jetpt2-jeteta2-b1','nj_jetpt2-jeteta2-b2']] 
    heff_fj[0].SetTitle('#eta-inclusive;p_{T} [GeV]; efficiency')
    heff_fj[1].SetTitle('0.   < |#eta| < 0.75')
    heff_fj[2].SetTitle('0.75 < |#eta| < 1.5')
    heff_fj[3].SetTitle('1.5  < |#eta| < 2.8')
    hratio = H1RatioPlotter(markers=markers)
    hratio.set(*heff_fj)
    hratio.markersize = 12
    hratio.scalemax = 1.2
    hratio.legtextsize = 25
    hratio.legboxsize  = 30
    hratio.legalign    = ('r','b')

    c = hratio.plot()
    for ext in imgext:
        c.Print(prefix+'heff_ratio_fj.%s' % ext )
    

    # --0jet---
    heff_0j = [ heffs['heff_ds_%s' % n] for n in [ '1j_jetpt2','1j_jetpt2-jeteta2-b0','1j_jetpt2-jeteta2-b1','1j_jetpt2-jeteta2-b2']] 
    heff_0j[0].SetTitle('#eta-inclusive;p_{T} [GeV]; efficiency')
    heff_0j[1].SetTitle('0.   < |#eta| < 0.75')
    heff_0j[2].SetTitle('0.75 < |#eta| < 1.5')
    heff_0j[3].SetTitle('1.5  < |#eta| < 2.8')
    hratio = H1RatioPlotter(markers=markers)
    hratio.set(*heff_0j)
    hratio.markersize = 16
    hratio.scalemax = 1.2
    hratio.legtextsize = 25
    hratio.legboxsize  = 30

    c = hratio.plot()
    for ext in imgext:
        c.Print(prefix+'heff_ratio_0j.%s' % ext )

    # --1jet---
    heff_1j = [ heffs['heff_ds_%s' % n] for n in [ '2j_jetpt2','2j_jetpt2-jeteta2-b0','2j_jetpt2-jeteta2-b1','2j_jetpt2-jeteta2-b2']] 
    heff_1j[0].SetTitle('#eta-inclusive;p_{T} [GeV]; efficiency')
    heff_1j[1].SetTitle('0.   < |#eta| < 0.75')
    heff_1j[2].SetTitle('0.75 < |#eta| < 1.5')
    heff_1j[3].SetTitle('1.5  < |#eta| < 2.8')
    hratio = H1RatioPlotter(markers=markers)
    hratio.set(*heff_1j)
    hratio.markersize = 16
    hratio.scalemax = 1.
    hratio.legtextsize = 25
    hratio.legboxsize  = 30
    hratio.legalign    = ('r','b')

    c = hratio.plot()
    for ext in imgext:
        c.Print(prefix+'heff_ratio_1j.%s' % ext )

def doraemon(analysers, imgext, prefix, opt ):
    
    n ='ttbar'
    a = analysers[n]
    old = a.worker.entries()
#     a.bufferentries()
    a.worker.setalias('bjet1','(jettche1 > 2.1)')
    a.worker.setalias('bjet2','(jettche2 > 2.1)')
    a.worker.setalias('bjet3','(jettche3 > 2.1)')
    a.worker.setalias('bjet4','(jettche4 > 2.1)')
    a.worker.setalias('myjets','(njet >= 1)')
    a.worker.setalias('nbtags','bjet1 + bjet2 + bjet3 + bjet4')

    print '  ',n,':',old,'>>',a.selectedentries(),'...done'
    print a.worker.aliases()

    print 'bctrl:' ,a.views['bctrl'].entries('njet <= 2')
    print 'btag:'  ,a.views['btag'].entries('njet <= 2')

    print 'btag X00:',a.views['bctrl'].entries('njet<=2 && !bjet3 && !bjet4')
    print 'btag 0X0:',a.views['bctrl'].entries('njet<=2 && !bjet2 && !bjet4')
    print 'btag 00X:',a.views['bctrl'].entries('njet<=2 && !bjet2 && !bjet3')

    print 'btag 100:',a.views['bctrl'].entries('njet<=2 &&  bjet2 && !bjet3 && !bjet4')
    print 'btag 010:',a.views['bctrl'].entries('njet<=2 && !bjet2 &&  bjet3 && !bjet4')
    print 'btag 001:',a.views['bctrl'].entries('njet<=2 && !bjet2 && !bjet3 &&  bjet4')

    print 'btag 011:',a.views['bctrl'].entries('njet<=2 && !bjet2 &&  bjet3 &&  bjet4')
    print 'btag 101:',a.views['bctrl'].entries('njet<=2 &&  bjet2 && !bjet3 &&  bjet4')
    print 'btag 110:',a.views['bctrl'].entries('njet<=2 &&  bjet2 &&  bjet3 && !bjet4')

    print 'btag 111:',a.views['bctrl'].entries('njet<=2 &&  bjet2 &&  bjet3 &&  bjet4')
    print 'btag 000:',a.views['bctrl'].entries('njet<=2 && !bjet2 && !bjet3 && !bjet4')

    print '- bctrl:' ,a.views['bctrl'].entries('myjets')
    print '- btag:'  ,a.views['btag'].entries('myjets')

    print '- btag X00:',a.views['bctrl'].entries('myjets && jetpt1 > 50 && !bjet3 && !bjet4')
    print '- btag 0X0:',a.views['bctrl'].entries('myjets && jetpt1 > 50 && !bjet2 && !bjet4')
    print '- btag 00X:',a.views['bctrl'].entries('myjets && jetpt1 > 50 && !bjet2 && !bjet3')

    print '- btag 100:',a.views['bctrl'].entries('myjets && jetpt1 > 50 &&  bjet2 && !bjet3 && !bjet4')
    print '- btag 010:',a.views['bctrl'].entries('myjets && jetpt1 > 50 && !bjet2 &&  bjet3 && !bjet4')
    print '- btag 001:',a.views['bctrl'].entries('myjets && jetpt1 > 50 && !bjet2 && !bjet3 &&  bjet4')

    print '- btag 011:',a.views['bctrl'].entries('myjets && jetpt1 > 50 && !bjet2 &&  bjet3 &&  bjet4')
    print '- btag 101:',a.views['bctrl'].entries('myjets && jetpt1 > 50 &&  bjet2 && !bjet3 &&  bjet4')
    print '- btag 110:',a.views['bctrl'].entries('myjets && jetpt1 > 50 &&  bjet2 &&  bjet3 && !bjet4')

    print '- btag 111:',a.views['bctrl'].entries('myjets && jetpt1 > 50 &&  bjet2 &&  bjet3 &&  bjet4')
    print '- btag 000:',a.views['bctrl'].entries('myjets && jetpt1 > 50 && !bjet2 && !bjet3 && !bjet4')

    h = a.views['bctrl'].plot('ziogatto','nbtags','njet >= 1', bins=(5,0,5))

    for i in xrange(h.GetNbinsX()):
        print i,h.GetBinContent(i+1)

#     print 'btag X00:',a.views['bctrl'].entries('jettche3 < 2.1 && jettche4 < 2.1')
#     print 'btag 0X0:',a.views['bctrl'].entries('jettche2 < 2.1 && jettche4 < 2.1')
#     print 'btag 00X:',a.views['bctrl'].entries('jettche2 < 2.1 && jettche3 < 2.1')

#     print 'btag 100:',a.views['bctrl'].entries('jettche2 > 2.1 && jettche3 < 2.1 && jettche4 < 2.1')
#     print 'btag 010:',a.views['bctrl'].entries('jettche2 < 2.1 && jettche3 > 2.1 && jettche4 < 2.1')
#     print 'btag 001:',a.views['bctrl'].entries('jettche2 < 2.1 && jettche3 < 2.1 && jettche4 > 2.1')

#     print 'btag 011:',a.views['bctrl'].entries('jettche2 < 2.1 && jettche3 > 2.1 && jettche4 > 2.1')
#     print 'btag 101:',a.views['bctrl'].entries('jettche2 > 2.1 && jettche3 < 2.1 && jettche4 > 2.1')
#     print 'btag 110:',a.views['bctrl'].entries('jettche2 > 2.1 && jettche3 > 2.1 && jettche4 < 2.1')

#     print 'btag 111:',a.views['bctrl'].entries('jettche2 > 2.1 && jettche3 > 2.1 && jettche4 > 2.1')
#     print 'btag 000:',a.views['bctrl'].entries('jettche2 < 2.1 && jettche3 < 2.1 && jettche4 < 2.1')



#     ____        __ 
#    / __ \____ _/ /_
#   / /_/ / __ `/ __/
#  / _, _/ /_/ / /_  
# /_/ |_|\__,_/\__/  
#                    
#---

commons = {
    'x1tag':'n_{jet}=1',
    'x2tag':'n_{jet}=0',
    'ltitle':'n_{jet}=0',
    'rtitle':'n_{jet}=0',

} 


def plotbtags(name, x1_plots,x2_plots, **kwargs): #, ltitle, rtitle, legtextsize, legboxsize, scalemax, markersize ):

    colors  = [ROOT.kRed+1      , ROOT.kAzure-5   ]
    markers = [ROOT.kFullCircle , ROOT.kFullCircle]
    x1_mc  = sumplots(x1_plots, 'x1_top', tops)
    x1_ds  = x1_plots['Data'].Clone('x1_ds')
    x1_ds -= sumplots(x1_plots, 'x1_others', others)

    x2_mc  = sumplots(x2_plots, 'x2_top', tops)
    x2_ds  = x2_plots['Data'].Clone('x2_ds')
    x2_ds -= sumplots(x2_plots, 'x2_others', others)

    # btag region
    x1_ds.Scale(1./x1_ds.Integral())
    x1_mc.Scale(1./x1_mc.Integral())
    x2_ds.Scale(1./x2_ds.Integral())
    x2_mc.Scale(1./x2_mc.Integral())

    x1_tag = kwargs.get('x1tag','')
    x2_tag = kwargs.get('x2tag','')
    xtitle = kwargs.get('xtitle','')
    ytitle = kwargs.get('ytitle','')

    x1_mc.SetTitle('mc_{top} %s;%s;%s' % (x1_tag,xtitle,ytitle))
    x2_mc.SetTitle('mc_{top} %s;%s;%s' % (x2_tag,xtitle,ytitle))
    x1_ds.SetTitle('data-mc_{others} %s;%s;%s' % (x1_tag,xtitle,ytitle))
    x2_ds.SetTitle('data-mc_{others} %s;%s;%s' % (x2_tag,xtitle,ytitle))

#         x1_mc.SetTitle('mc_{top} n_{jet}=1;'+xlabel)
#         x2_mc.SetTitle('mc_{top} n_{jet}=0;'+xlabel)
#         x1_ds.SetTitle('data-mc_{others} n_{jet}=1;'+xlabel)
#         x2_ds.SetTitle('data-mc_{others} n_{jet}=0;'+xlabel)

    hratio = H1RatioPlotter(colors=colors,markers=markers)
    hratio.set(x1_ds,x2_ds)
    hratio.ltitle      = kwargs.get('ltitle','')
    hratio.rtitle      = kwargs.get('rtitle','')
    hratio.ytitle2     = kwargs.get('ytitle2','ratio')
    hratio.legtextsize = kwargs.get('legtextsize',25)
    hratio.legboxsize  = kwargs.get('legboxsize',30)
    hratio.scalemax    = kwargs.get('scalemax',1.)
    hratio.markersize  = kwargs.get('markersize',16)
    hratio.legalign    = kwargs.get('legalign',('l','t'))

    c = hratio.plot()
    imgext = kwargs.get('imgext',[])
    prefix = kwargs.get('prefix',[])
    for ext in imgext:
        c.Print(prefix+name+'.'+ext)

def doratios(analysers, imgext, prefix, opt ):

    # block B: control shapes
    vars = {
        'njet'                      : ('njet'          , 'base'     ,  'jettche1 > 2.1'                                  , (10,  0,10), 'n^{jets}'),
        'tche1'                     : ('jettche1'      , 'base'     ,  ''                                                , (40,-20,20), 'TCHE_{j1}' ),

        'softjet_1j_btotal_jeteta2' : ('fabs(jeteta2)' , 'bveto-mu' ,  'njet == 1 && jettche1 > 2.1'                      , ( 12,  0,  3 ) , 'eta_{j2}' ),
        'softjet_1j_btag_jeteta2'   : ('fabs(jeteta2)' , 'bveto-mu' ,  'njet == 1 && jettche1 > 2.1 && jettche2 >  2.1'   , ( 12,  0,  3 ) , 'eta_{j2}' ),
        'softjet_1j_bveto_jeteta2'  : ('fabs(jeteta2)' , 'bveto-mu' ,  'njet == 1 && jettche1 > 2.1 && jettche2 <= 2.1'   , ( 12,  0,  3 ) , 'eta_{j2}' ),
        'softjet_0j_btotal_jeteta1' : ('fabs(jeteta1)' , 'bveto-mu' ,  'njet == 0'                                        , ( 12,  0,  3 ) , 'eta_{j1}' ),
        'softjet_0j_btag_jeteta1'   : ('fabs(jeteta1)' , 'bveto-mu' ,  'njet == 0 && jettche1 >  2.1'                     , ( 12,  0,  3 ) , 'eta_{j1}' ),
        'softjet_0j_bveto_jeteta1'  : ('fabs(jeteta1)' , 'bveto-mu' ,  'njet == 0 && jettche1 <= 2.1'                     , ( 12,  0,  3 ) , 'eta_{j1}' ),

        'softjet_1j_btotal_jetpt2'  : ('jetpt2'        , 'bveto-mu' ,  'njet == 1 && jettche1 > 2.1'                      , ( 10, 10, 30 ) , 'pt_{j2}' ) ,
        'softjet_1j_btag_jetpt2'    : ('jetpt2'        , 'bveto-mu' ,  'njet == 1 && jettche1 > 2.1 && jettche2 >  2.1'   , ( 10, 10, 30 ) , 'pt_{j2}' ) ,
        'softjet_1j_bveto_jetpt2'   : ('jetpt2'        , 'bveto-mu' ,  'njet == 1 && jettche1 > 2.1 && jettche2 <= 2.1'   , ( 10, 10, 30 ) , 'pt_{j2}' ) ,
        'softjet_0j_btotal_jetpt1'  : ('jetpt1'        , 'bveto-mu' ,  'njet == 0'                                        , ( 10, 10, 30 ) , 'pt_{j1}' ) ,
        'softjet_0j_btag_jetpt1'    : ('jetpt1'        , 'bveto-mu' ,  'njet == 0 && jettche1 >  2.1'                     , ( 10, 10, 30 ) , 'pt_{j1}' ) ,
        'softjet_0j_bveto_jetpt1'   : ('jetpt1'        , 'bveto-mu' ,  'njet == 0 && jettche1 <= 2.1'                     , ( 10, 10, 30 ) , 'pt_{j1}' ) ,


        'jet_2j_btotal_jeteta2'     : ('fabs(jeteta2)' , 'bveto-mu' ,  'bveto_ip && njet == 2 && jettche1 > 2.1'                    , ( 12,  0,  3 ) , 'eta_{j2}' ),
        'jet_2j_btag_jeteta2'       : ('fabs(jeteta2)' , 'bveto-mu' ,  'bveto_ip && njet == 2 && jettche1 > 2.1 && jettche2 >  2.1' , ( 12,  0,  3 ) , 'eta_{j2}' ),
        'jet_2j_bveto_jeteta2'      : ('fabs(jeteta2)' , 'bveto-mu' ,  'bveto_ip && njet == 2 && jettche1 > 2.1 && jettche2 <= 2.1' , ( 12,  0,  3 ) , 'eta_{j2}' ),
        'jet_1j_btotal_jeteta1'     : ('fabs(jeteta1)' , 'bveto-mu' ,  'bveto_ip && njet == 1'                                      , ( 12,  0,  3 ) , 'eta_{j1}' ),
        'jet_1j_btag_jeteta1'       : ('fabs(jeteta1)' , 'bveto-mu' ,  'bveto_ip && njet == 1 && jettche1 >  2.1'                   , ( 12,  0,  3 ) , 'eta_{j1}' ),
        'jet_1j_bveto_jeteta1'      : ('fabs(jeteta1)' , 'bveto-mu' ,  'bveto_ip && njet == 1 && jettche1 <= 2.1'                   , ( 12,  0,  3 ) , 'eta_{j1}' ),

        'jet_2j_btotal_jetpt2'      : ('jetpt2'        , 'bveto-mu' ,  'bveto_ip&& njet == 2 && jettche1 > 2.1'                     , ( 34, 30, 200 ) , 'pt_{j2}' ) ,
        'jet_2j_btag_jetpt2'        : ('jetpt2'        , 'bveto-mu' ,  'bveto_ip&& njet == 2 && jettche1 > 2.1 && jettche2 >  2.1'  , ( 34, 30, 200 ) , 'pt_{j2}' ) ,
        'jet_2j_bveto_jetpt2'       : ('jetpt2'        , 'bveto-mu' ,  'bveto_ip&& njet == 2 && jettche1 > 2.1 && jettche2 <= 2.1'  , ( 34, 30, 200 ) , 'pt_{j2}' ) ,
        'jet_1j_btotal_jetpt1'      : ('jetpt1'        , 'bveto-mu' ,  'bveto_ip&& njet == 1'                                       , ( 34, 30, 200 ) , 'pt_{j1}' ) ,
        'jet_1j_btag_jetpt1'        : ('jetpt1'        , 'bveto-mu' ,  'bveto_ip&& njet == 1 && jettche1 >  2.1'                    , ( 34, 30, 200 ) , 'pt_{j1}' ) ,
        'jet_1j_bveto_jetpt1'       : ('jetpt1'        , 'bveto-mu' ,  'bveto_ip&& njet == 1 && jettche1 <= 2.1'                    , ( 34, 30, 200 ) , 'pt_{j1}' ) ,
    }

    xplots = AlienDict()
    for v,(expr,lvl,cut,bins,xaxis) in vars.iteritems():
        print '%-30s: [' % v,
        for n,a in analysers.iteritems():
            print '%s,' % n,
            xplots[v][n] = a.views[lvl].plot('xplots_%s_%s' % (n,v), expr,extra=cut,bins=bins)
            sys.stdout.flush()
        print ']'
        if opt.datamc:
            hwwlatino.printplots(xplots[v],prefix+'xplots_%s_%s' % (lvl,v), xaxis=xaxis, label='base, %s' % cut, lumi=opt.lumi, exts=imgext)

    xplots.lock()
        
#     if 'Top' not in analysers:
#         print 'No top, no party'
#         return



    colors  = [ROOT.kRed+1      , ROOT.kAzure-5   ]
    markers = [ROOT.kFullCircle , ROOT.kFullCircle]

    # 
    #  xplots
    # 

    commons = {
        'imgext'     : imgext,
        'prefix'     : prefix,
        'ytitle'     : 'normalised',
        'legboxsize' : 40,
    } 

    #    ____        _      __ 
    #   / __ \      (_)__  / /_
    #  / / / /_____/ / _ \/ __/
    # / /_/ /_____/ /  __/ /_  
    # \____/   __/ /\___/\__/  
    #         /___/            
    decorations = {
        'x1tag'      : 'n_{jet}=1',
        'x2tag'      : 'n_{jet}=0',
        'ltitle'     : '"hardest" jet pt_{j}^{soft}',
        'rtitle'     : '0j vs 1j',
        'xtitle'     : 'pt_{j}^{soft}',
        'ytitle2'    : '0j / 1j',
        'scalemax'   : 1.2
    } 

    kwargs = dict(commons,**decorations)
    plotbtags('ratio_softjet_pt_0j1j_btag_ziogano' , xplots['softjet_1j_btag_jetpt2'], xplots['softjet_0j_btag_jetpt1'], **kwargs)

    decorations = {
        'x1tag'      : 'n_{jet}=1',
        'x2tag'      : 'n_{jet}=0',
        'ltitle'     : '"hardest" jet #eta_{j}^{soft}',
        'rtitle'     : '0j vs 1j',
        'xtitle'     : '#eta_{j}^{soft}',
        'ytitle2'    : '0j / 1j',
        'scalemax'   : 1.2,
        'legalign'   : ('r','t'),
    } 
    kwargs = dict(commons,**decorations)
    plotbtags('ratio_softjet_eta_0j1j_btag_ziogano', xplots['softjet_1j_btag_jeteta2'], xplots['softjet_0j_btag_jeteta1'], **kwargs)

    #    ___      _      __ 
    #   <  /     (_)__  / /_
    #   / /_____/ / _ \/ __/
    #  / /_____/ /  __/ /_  
    # /_/   __/ /\___/\__/  
    #      /___/            
    decorations = {
        'x1tag'      : 'n_{jet}=2',
        'x2tag'      : 'n_{jet}=1',
        'ltitle'     : 'second jet pt_{j}',
        'rtitle'     : '1j vs 2j',
        'xtitle'     : 'pt_{j}',
        'ytitle2'    : '1j / 2j',
        'scalemax'   : 1.2,
        'legalign'   : ('r','t'),
    } 

    kwargs = dict(commons,**decorations)
    plotbtags('ratio_jet_pt_1j2j_btag_ziogano' , xplots['jet_2j_btag_jetpt2'], xplots['jet_1j_btag_jetpt1'], **kwargs)


    decorations = {
        'x1tag'      : 'n_{jet}=2',
        'x2tag'      : 'n_{jet}=1',
        'ltitle'     : 'second jet #eta_{j}',
        'rtitle'     : '1j vs 2j',
        'xtitle'     : '#eta_{j}',
        'ytitle2'    : '1j / 2j',
        'scalemax'   : 1.2,
        'legalign'   : ('r','t'),
    } 
    kwargs = dict(commons,**decorations)
    plotbtags('ratio_jet_eta_1j2j_btag_ziogano', xplots['jet_2j_btag_jeteta2'], xplots['jet_1j_btag_jeteta1'], **kwargs)
    #    ____        _      __ 
    #   / __ \      (_)__  / /_
    #  / / / /_____/ / _ \/ __/
    # / /_/ /_____/ /  __/ /_  
    # \____/   __/ /\___/\__/  
    #         /___/            

#     colors  = [ROOT.kRed+1      , ROOT.kAzure-5   , ROOT.kRed+1      , ROOT.kAzure-5   ]
#     markers = [ROOT.kOpenCircle , ROOT.kOpenCircle, ROOT.kFullCircle , ROOT.kFullCircle]
    if 'softjet_0j_btag_jetpt1' in xplots and 'softjet_1j_btag_jetpt2' in xplots: 

        ptj2_btag  = xplots['softjet_1j_btag_jetpt2']
        ptj1_btag  = xplots['softjet_0j_btag_jetpt1']

        ptj2_btag_mc  = sumplots(ptj2_btag, 'ptj2_btag_top', tops)
        ptj2_btag_ds  = ptj2_btag['Data'].Clone('ptj2_btag_ds')
        ptj2_btag_ds -= sumplots(ptj2_btag, 'ptj2_btag_others', others)

        ptj1_btag_mc  = sumplots(ptj1_btag, 'ptj1_btag_top', tops)
        ptj1_btag_ds  = ptj1_btag['Data'].Clone('ptj1_btag_ds')
        ptj1_btag_ds -= sumplots(ptj1_btag, 'ptj1_btag_others', others)

#         ptj2_btag_mc  = ptj2_btag['Top']
#         ptj2_btag_ds  = ptj2_btag['Data'].Clone('ptj2_btag_ds')
#         for p in others:  ptj2_btag_ds  -= ptj2_btag[p

#         ptj1_btag_mc  = ptj1_btag['Top']
#         ptj1_btag_ds  = ptj1_btag['Data'].Clone('ptj1_btag_ds')
#         for p in others:  ptj1_btag_ds  -= ptj1_btag[p]

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

    if 'softjet_0j_btag_jeteta1' in xplots and 'softjet_1j_btag_jeteta2' in xplots: 

        etaj2_btag  = xplots['softjet_1j_btag_jeteta2']
        etaj1_btag  = xplots['softjet_0j_btag_jeteta1']

        etaj2_btag_mc  = sumplots(etaj2_btag, 'etaj2_btag_top', tops)
        etaj2_btag_ds  = etaj2_btag['Data'].Clone('etaj2_btag_ds')
        etaj2_btag_ds -= sumplots(etaj2_btag, 'etaj2_btag_others', others)

        etaj1_btag_mc  = sumplots(etaj1_btag, 'etaj1_btag_top', tops)
        etaj1_btag_ds  = etaj1_btag['Data'].Clone('etaj1_btag_ds')
        etaj1_btag_ds -= sumplots(etaj1_btag, 'etaj1_btag_others', others)
#         etaj2_btag_mc  = etaj2_btag['Top']
#         etaj2_btag_ds  = etaj2_btag['Data'].Clone('etaj2_btag_ds')
#         for p in others:  etaj2_btag_ds  -= etaj2_btag[p]

#         etaj1_btag_mc  = etaj1_btag['Top']
#         etaj1_btag_ds  = etaj1_btag['Data'].Clone('etaj1_btag_ds')
#         for p in others:  etaj1_btag_ds  -= etaj1_btag[p]

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
#         hratio_0j1j_btag.set(etaj2_btag_mc,etaj1_btag_mc,etaj2_btag_ds,etaj1_btag_ds)
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

    #    ___      _      __ 
    #   <  /     (_)__  / /_
    #   / /_____/ / _ \/ __/
    #  / /_____/ /  __/ /_  
    # /_/   __/ /\___/\__/  
    #      /___/            

    # 1-jet bin estimate
    if 'jet_1j_btag_jetpt1' in xplots and 'jet_2j_btag_jetpt2' in xplots: 

        ptj2_btag  = xplots['jet_2j_btag_jetpt2']
        ptj1_btag  = xplots['jet_1j_btag_jetpt1']

        ptj2_btag_mc  = sumplots(ptj2_btag, 'ptj2_btag_top', tops)
        ptj2_btag_ds  = ptj2_btag['Data'].Clone('ptj2_btag_ds')
        ptj2_btag_ds -= sumplots(ptj2_btag, 'ptj2_btag_others', others)

        ptj1_btag_mc  = sumplots(ptj1_btag, 'ptj1_btag_top', tops)
        ptj1_btag_ds  = ptj1_btag['Data'].Clone('ptj1_btag_ds')
        ptj1_btag_ds -= sumplots(ptj1_btag, 'ptj1_btag_others', others)

#         ptj2_btag_mc  = ptj2_btag['Top']
#         ptj2_btag_ds  = ptj2_btag['Data'].Clone('ptj2_btag_ds')
#         for p in others:  ptj2_btag_ds  -= ptj2_btag[p]

#         ptj1_btag_mc  = ptj1_btag['Top']
#         ptj1_btag_ds  = ptj1_btag['Data'].Clone('ptj1_btag_ds')
#         for p in others:  ptj1_btag_ds  -= ptj1_btag[p]

        # btag region
        ptj2_btag_ds.Scale(1./ptj2_btag_ds.Integral())
        ptj2_btag_mc.Scale(1./ptj2_btag_mc.Integral())
        ptj1_btag_ds.Scale(1./ptj1_btag_ds.Integral())
        ptj1_btag_mc.Scale(1./ptj1_btag_mc.Integral())

        ptj2_btag_mc.SetTitle('mc_{top} n_{jet}=2;pt_{j}')
        ptj1_btag_mc.SetTitle('mc_{top} n_{jet}=1;pt_{j}')
        ptj2_btag_ds.SetTitle('data-mc_{others} n_{jet}=2;pt_{j}')
        ptj1_btag_ds.SetTitle('data-mc_{others} n_{jet}=1;pt_{j}')

        hratio_1j2j_btag = H1RatioPlotter(colors=colors,markers=markers)
#         hratio_1j2j_btag.set(ptj2_btag_mc, ptj1_btag_mc, ptj2_btag_ds, ptj1_btag_ds)
        hratio_1j2j_btag.set(ptj2_btag_ds,ptj1_btag_ds)
        hratio_1j2j_btag.ltitle      = 'pt_{j-soft}'
        hratio_1j2j_btag.rtitle      = '1j vs 2j'
        hratio_1j2j_btag.ytitle2     = 'x/mc_{top} n_{jet} = 2'
        hratio_1j2j_btag.legtextsize = 25
        hratio_1j2j_btag.legboxsize  = 30
        hratio_1j2j_btag.scalemax    = 1.2
        hratio_1j2j_btag.markersize  = 16
        hratio_1j2j_btag.legalign    = ('r','t')

        c = hratio_1j2j_btag.plot()
        for ext in imgext:
            c.Print(prefix+'ratio_jet_pt_1j2j_btag.'+ext)

    if 'jet_1j_btag_jeteta1' in xplots and 'jet_2j_btag_jeteta2' in xplots: 

        etaj2_btag  = xplots['jet_2j_btag_jeteta2']
        etaj1_btag  = xplots['jet_1j_btag_jeteta1']

        etaj2_btag_mc  = sumplots(etaj2_btag, 'etaj2_btag_top', tops)
        etaj2_btag_ds  = etaj2_btag['Data'].Clone('etaj2_btag_ds')
        etaj2_btag_ds -= sumplots(etaj2_btag, 'etaj2_btag_others', others)

        etaj1_btag_mc  = sumplots(etaj1_btag, 'etaj1_btag_top', tops)
        etaj1_btag_ds  = etaj1_btag['Data'].Clone('etaj1_btag_ds')
        etaj1_btag_ds -= sumplots(etaj1_btag, 'etaj1_btag_others', others)

#         etaj2_btag_mc  = etaj2_btag['Top']
#         etaj2_btag_ds  = etaj2_btag['Data'].Clone('etaj2_btag_ds')
#         for p in others:  etaj2_btag_ds  -= etaj2_btag[p]

#         etaj1_btag_mc  = etaj1_btag['Top']
#         etaj1_btag_ds  = etaj1_btag['Data'].Clone('etaj1_btag_ds')
#         for p in others:  etaj1_btag_ds  -= etaj1_btag[p]

        # btag region
        etaj2_btag_ds.Scale(1./etaj2_btag_ds.Integral())
        etaj2_btag_mc.Scale(1./etaj2_btag_mc.Integral())
        etaj1_btag_ds.Scale(1./etaj1_btag_ds.Integral())
        etaj1_btag_mc.Scale(1./etaj1_btag_mc.Integral())

        etaj2_btag_mc.SetTitle('mc_{top} n_{jet}=1;#eta_{j}')
        etaj1_btag_mc.SetTitle('mc_{top} n_{jet}=0;#eta_{j}')
        etaj2_btag_ds.SetTitle('data-mc_{others} n_{jet}=1;#eta_{j}')
        etaj1_btag_ds.SetTitle('data-mc_{others} n_{jet}=0;#eta_{j}')

        hratio_1j2j_btag = H1RatioPlotter(colors=colors,markers=markers)
#         hratio_1j2j_btag.set(etaj2_btag_mc, etaj1_btag_mc, etaj2_btag_ds, etaj1_btag_ds)
        hratio_1j2j_btag.set(etaj2_btag_ds,etaj1_btag_ds)
        hratio_1j2j_btag.ltitle      = '#eta_{j}'
        hratio_1j2j_btag.rtitle      = '1j vs 2j'
        hratio_1j2j_btag.ytitle2     = 'x/mc_{top} n_{jet} = 1'
        hratio_1j2j_btag.legtextsize = 25
        hratio_1j2j_btag.legboxsize  = 30
        hratio_1j2j_btag.scalemax    = 1.2
        hratio_1j2j_btag.markersize  = 16
        hratio_1j2j_btag.legalign    = ('r','t')

        c = hratio_1j2j_btag.plot()
        for ext in imgext:
            c.Print(prefix+'ratio_jet_eta_1j2j_btag.'+ext)
# ---
def main( opt ):
    samples = hwwlatino.samples(125,'8TeV','Data2012','SM','topestimate')
    samples = hwwlatino.samples(125,'8TeV','Data2012','SM','topplots')

    wwflow = CutFlow(wwcuts.wwcommon)

    del wwflow['bveto_mu']
    del wwflow['bveto_ip']
    wwflow['ptll'] = 'ptll>45'
    
    print '-'*80
    for n,c in wwflow.iteritems():
        print '%-30s: %s'% (n,c)

    print '-'*80
    print wwflow.string()
    print '-'*80

    topflow = CutFlow()

    topflow['base']     = wwflow.string()
    topflow['bveto-mu'] = 'bveto_mu'
    topflow['bctrl']    = 'jettche1>2.1'
    topflow['btag']     = 'jettche2>2.1'

    for n,s in samples.iteritems():
        print '%-10s'%n,s
    analysers = hwwlatino.makeanalysers(samples,orchard,topflow,opt.lumi)

    if not opt.nobuf:
        for n,a in analysers.iteritems():
            old = a.worker.entries()
            a.bufferentries()
            print '  ',n,':',old,'>>',a.selectedentries(),'...done'
    

    print '-'*80
    imgext = ['pdf','png']
    prefix=''
    if opt.out:
        hwwtools.ensuredir(opt.out)
        prefix = opt.out+'/'
   
    if opt.eff: doefficiencies( analysers, imgext, prefix, opt )
    if opt.rat: doratios(analysers, imgext, prefix, opt )
    if opt.dora: doraemon(analysers, imgext, prefix, opt )


# ---
if __name__ == '__main__':
    import optparse
    import hwwtools
    import sys
    import bdb

    
    parser = optparse.OptionParser()
    parser.add_option('-d', '--debug'    , dest = 'debug'       , help='Debug level'            , default=0 , action='count' )
    parser.add_option('-l', '--lumi'     , dest = 'lumi'        , help='Luminosity'             , default=19.468 )
    parser.add_option('-o', '--out'      , dest = 'out'         , help='Output'                 , default=None )
    parser.add_option('--no-buff'        , dest = 'nobuf'       , help='Don\'t pre-buffer'      , action='store_true' ,default=False )
    parser.add_option('--datamc-plots'   , dest = 'datamc'      , help='print datamc-plots'     , type='int',default=False )
    parser.add_option('--eff'            , dest = 'eff'         , help='do efficiency'          , action='store_true',default=False )
    parser.add_option('--rat'            , dest = 'rat'         , help='do ratios'              , action='store_true',default=False )
    parser.add_option('--dora'           , dest = 'dora'        , help='doraemon'               , action='store_true',default=False )

    

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
  
