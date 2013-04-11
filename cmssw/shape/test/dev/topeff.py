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

others = ['WW','ggWW','WJet','DYLL','DYTT','VV','Vg','VgS','tW']
tops   = ['ttbar']

from base import AlienDict

# ---
def thsum( plots, newname, processes ):

#     print plots
#     print newname
#     print processes
    if not processes: raise ValueError('Processes is emtty!')

    hnew = plots[processes[0]].Clone(newname)
    hnew.Reset()
    for p in processes: 
        hnew += plots[p]

    return hnew

#     ______________
#    / ____/ __/ __/
#   / __/ / /_/ /_  
#  / /___/ __/ __/  
# /_____/_/ /_/     

# _______________________________________________________________________________
def vars2plots( analysers, vars, opt ):
    # prepare 
    plots = AlienDict()

    morevars = {}

    for v,(expr,cut,bins,div,xlabel,units) in vars.iteritems():
        print '%-30s: [' % v,
        for n,a in analysers.iteritems():
            print '%s,' % n,

            pf = a.plotsflow(n+'_'+v,expr,extra=cut,bins=bins)
            sys.stdout.flush()
            

            for c,h in pf.iteritems():
                if isinstance(h,ROOT.TH2):
                    #make the slices
                    for i in xrange(1,h.GetNbinsY()+1):
                        hp = h.ProjectionX('%s_b%d' % (h.GetName(),i),i,i,'e')
                        nv = '%s_by%d' % (v,i)
                        plots[nv][c][n] = hp
                        morevars[nv] = (expr,cut,bins,div,xlabel,units)

                    # and the inclusive
                    hp = h.ProjectionX('%s_inc' % h.GetName(),1,h.GetNbinsY(),'e')
                    nv = '%s_inc' % v
                    plots[nv][c][n] = hp
                    morevars[nv] = (expr,cut,bins,div,xlabel,units)

                plots[v][c][n] = h
        print ']'


    plots.lock()
    vars.update(morevars)

    if opt.datamc:
        for v,(expr,cut,bins,div,xlabel,units) in vars.iteritems():
            print 'Printing plots:',v
            # but only if they are 1D
            if isinstance(plots[v]['base']['Data'],ROOT.TH2): continue
            print v
            hwwlatino.printplots(plots[v]['bctrl'],prefix+'bplot_bctrl_%s' % v,xlabel=xlabel, units=units, label='b_{CTRL}, %s' % cut, div=div, lumi=opt.lumi, exts=imgext)
            hwwlatino.printplots(plots[v]['btag'] ,prefix+'bplot_btag_%s'  % v,xlabel=xlabel, units=units, label='b_{TAG} , %s' % cut, div=div, lumi=opt.lumi, exts=imgext)
    
    
    return plots

# _______________________________________________________________________________
def makeefficiency( name, bplots, x, xlabel, scalemax, legalign, lumi, imgext, prefix,  logx=False,save=None ):

    colors  = [ROOT.kRed+1      , ROOT.kAzure-5   ]
    markers = [ROOT.kFullCircle , ROOT.kFullCircle]

    try:
        x_bctrl = bplots[x]['bctrl']
        x_btag  = bplots[x]['btag']
    except KeyError:
        return

    x_bctrl_ds  = x_bctrl['Data'].Clone('bctrl_%s_ds' % x)
    x_bctrl_ds -= thsum(x_bctrl, 'bctrl_others', others)
#     for p in others:  x_bctrl_ds -= x_bctrl[p]

    x_btag_ds  = x_btag['Data'].Clone('btag_%s_ds' % x)
    x_btag_ds -= thsum(x_btag, 'bctrl_others', others)
#     for p in others:  x_btag_ds  -= x_btag[p]

    heff_ds_x = x_btag_ds.Clone('heff_ds_%s' % x)
    heff_ds_x.Divide(x_btag_ds,x_bctrl_ds,1,1,'b')
    heff_ds_x.SetTitle('data - mc_{others}')

    x_bctrl_mc = thsum(x_bctrl, 'bctrl_top_%s' % x, tops)
    x_btag_mc  = thsum(x_btag,  'btag_top_%s'  % x, tops)

    heff_mc_x = x_bctrl_mc.Clone('heff_mc_%s' % x)
    heff_mc_x.Reset()
    heff_mc_x.Divide(x_btag_mc,x_bctrl_mc,1,1,'b')
#     heff_mc_x.SetTitle('mc (tW/t#bar{t});%s;tag/ctrl' % xlabel)
    heff_mc_x.SetTitle('mc (t#bar{t});%s;tag/ctrl' % xlabel)


    hratio = H1RatioPlotter(colors=colors,markers=markers)
    hratio.scalemax   = scalemax
    hratio.legalign   = legalign
    hratio.ltitle     = 'top tag efficiency' 
    hratio.rtitle     = 'L=%.2f fb^{-1}' % lumi
    hratio.ytitle2    = 'data/mc'
    hratio.markersize = 16
    hratio.logx       = logx
    hratio.morelogx   = logx
    hratio.yrange     = (0.,0.9)
    hratio.set(heff_mc_x, heff_ds_x)
    c = hratio.plot()
    
    for ext in imgext:
        c.Print(prefix+name+'.'+ext )

    if save is not None:
        save['%s_ds' % name] = heff_ds_x
        save['%s_mc' % name] = heff_mc_x

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
#     jet_ptbins = range(10,  30, 4) + range( 30, 70, 5) + range( 70, 150, 20) + range(150,201, 50)

    softjets_ptbins = range(5, 31, 5)
    jet_ptbins  = range(5, 70, 5) + range( 70, 150, 20) + range(150,251, 50)
    hardjet_ptbins =                   range( 30, 70, 5) + range( 70, 150, 20) + range(150,251, 50)

    etabins = [0.,0.75,1.5,2.8,5]
    
    eta_b1 = 'fabs(jeteta2) < 0.75'                        
    eta_b2 = 'fabs(jeteta2) >=0.75 && fabs(jeteta2)< 1.5  '
    eta_b3 = 'fabs(jeteta2) >=1.5  && fabs(jeteta2)< 2.8'  
    eta_b4 = 'fabs(jeteta2) >=2.8  && fabs(jeteta2)< 5. '  
    
    hp50 = 'jetpt1 > 50 && ( njet <= 1 || (bveto_ip && njet >= 2) ) && jettche3<=2.1 && jettche4<=2.1'
    hp   =                '( njet <= 1 || (bveto_ip && njet >= 2) ) && jettche3<=2.1 && jettche4<=2.1'
    
    divbins = False


    # block A: distributions used to calculate the efficiencies
    vars = {
#         'mll_1j'        : ('mll'           , 'njet >= 1', (100,  0, 600) , 'm_{ll} [GeV]' ),
#         'tche2_1j'      : ('jettche2'      , 'njet >= 1', ( 40,-20,  20) , 'TCHE_{j2}' ),
        'kj_jet_eta2-pt2-hp50' : ('fabs(jeteta2):jetpt2', 'njet >= 1 && %s' % (hp50,) , ( jet_ptbins, etabins )     , divbins , 'p_{T}^{j2}', 'GeV' ) ,

        'kj_jet_eta2-pt2-hp'   : ('fabs(jeteta2):jetpt2', 'njet >= 1 && %s' % (hp,)   , ( jet_ptbins, etabins )     , divbins , 'p_{T}^{j2}', 'GeV' ) ,

        'kj_jet_eta2-pt2'      : ('fabs(jeteta2):jetpt2', 'njet >= 1'                 , ( jet_ptbins, etabins )     , divbins , 'p_{T}^{j2}', 'GeV' ) ,

        '1j_jeteta2'           : ('fabs(jeteta2)'       , 'njet == 1'                 , ( 12,  0,  3 )              , False   , '#eta_{j2}' , '' ),
        '1j_jeteta2-bins'      : ('fabs(jeteta2)'       , 'njet == 1'                 , ( etabins, )                , False   , '#eta_{j2}' , '' ),
        '1j_jetpt2-fine'       : ('jetpt2'              , 'njet == 1'                 , ( 10, 10, 30 )              , False   , 'p_{T}^{j2}', 'GeV' ) ,

        '1j_jet_eta2-pt2'      : ('fabs(jeteta2):jetpt2', 'njet == 1'                 , ( softjets_ptbins, etabins ) , divbins , 'p_{T}^{j2}', 'GeV' ) ,

        '2j_jeteta2'          : ('fabs(jeteta2)'        , 'bveto_ip && njet == 2'     , ( 12,  0,  3 )              , False   , '#eta_{j2}', '' ),
        '2j_jeteta2-bins'     : ('fabs(jeteta2)'        , 'bveto_ip && njet == 2'     , ( etabins, )                , False   , '#eta_{j2}', '' ),
                                                                                                                    
        '2j_jet_eta2-pt2'     : ('fabs(jeteta2):jetpt2' , 'bveto_ip && njet == 2'     , ( hardjet_ptbins, etabins ) , divbins , 'p_{T}^{j2}', 'GeV' ) ,
    }

    # use this for tests
    if opt.testmode:
        vars = {
            'kj_jet_eta2-pt2-hp50' : ('fabs(jeteta2):jetpt2', 'njet >= 1 && %s' % (hp50,) , ( jet_ptbins, etabins ) , divbins , 'p_{T}^{j2}', 'GeV' ) ,
            'kj_jet_eta2-pt2-hp'   : ('fabs(jeteta2):jetpt2', 'njet >= 1 && %s' % (hp,)   , ( jet_ptbins, etabins ) , divbins , 'p_{T}^{j2}', 'GeV' ) ,
        }


    bplots = vars2plots( analysers, vars, opt)

    heffs = odict.OrderedDict()

    # common paramete
    commons = {
        'lumi'   : opt.lumi,
        'imgext' : imgext,
        'prefix' : prefix,
        'save'   : heffs,
    }


    # calls to plot the efficiencies
    for x in ['inc','by1','by2','by3','by4']:

        makeefficiency( 'heff_kj_pt-hp50_%s' %x , bplots , 'kj_jet_eta2-pt2-hp50_%s' %x ,  'pt_{j2}' , 1.1 , ('r','b') , logx=True , **commons )
        makeefficiency( 'heff_kj_pt-hp_%s'   %x , bplots , 'kj_jet_eta2-pt2-hp_%s'   %x ,  'pt_{j2}' , 1.1 , ('r','b') , logx=True , **commons )
        makeefficiency( 'heff_kj_pt_%s' %x      , bplots , 'kj_jet_eta2-pt2_%s'      %x ,  'pt_{j2}' , 1.1 , ('r','b') , logx=True , **commons )
        makeefficiency( 'heff_0j_pt_%s' %x      , bplots , '1j_jet_eta2-pt2_%s'      %x ,  'pt_{j2}' , 1.1 , ('l','t') , logx=False, **commons )
        makeefficiency( 'heff_1j_pt_%s' %x      , bplots , '2j_jet_eta2-pt2_%s'      %x ,  'pt_{j2}' , 1.1 , ('r','b') , logx=False, **commons )
   
    makeefficiency( 'heff_0j_pt-fine'       , bplots , '1j_jetpt2-fine'          ,  'pt_{j2}'   , 1.1 , ('l','t') , **commons )
    makeefficiency( 'heff_0j_eta'           , bplots , '1j_jeteta2'              ,  '#eta_{j2}' , 1.3 , ('r','t') , **commons )
    makeefficiency( 'heff_0j_eta-bin'       , bplots , '1j_jeteta2-bins'         ,  '#eta_{j2}' , 1.3 , ('r','t') , **commons )
                                                                                 
    makeefficiency( 'heff_1j_eta'           , bplots , '2j_jeteta2'              ,  '#eta_{j2}' , 1.3 , ('r','t') , **commons )
    makeefficiency( 'heff_1j_eta-bin'       , bplots , '2j_jeteta2-bins'         ,  '#eta_{j2}' , 1.3 , ('r','t') , **commons )

    markers = [ROOT.kFullCircle , ROOT.kFullCircle , ROOT.kFullCircle , ROOT.kFullCircle]

    try:
#         heff_fj = [ heffs['heff_ds_%s' % n] for n in [ 'kj_jetpt2-hp50','kj_jetpt2-hp50-jeteta2-b0','kj_jetpt2-hp50-jeteta2-b1','kj_jetpt2-hp50-jeteta2-b2']] 
        heff = [heffs['heff_kj_pt-hp50_%s_ds' %x] for x in ['inc','by1','by2','by3','by4']]
        heff[0].SetTitle('#eta-inclusive;p_{T} [GeV]; efficiency')
        heff[1].SetTitle('0.   < |#eta| < 0.75')
        heff[2].SetTitle('0.75 < |#eta| < 1.5')
        heff[3].SetTitle('1.5  < |#eta| < 2.8')
        hratio = H1RatioPlotter(markers=markers)
        hratio.set(*heff)
        hratio.markersize  = 12
        hratio.scalemax    = 1.2
        hratio.ltitle      = 'eff in eta bins'
        hratio.rtitle      = 'high purity (p_{T}^{jet} > 50 GeV)'
        hratio.legtextsize = 25
        hratio.legboxsize  = 30
        hratio.legalign    = ('r','b')
        hratio.logx        = True
        hratio.yrange      = (0.,0.9)

        c = hratio.plot()
        for ext in imgext:
            c.Print(prefix+'heff_ratio_kj-hp50.%s' % ext )
    except KeyError as ke:
        print 'Not found:',ke
           

    try:
#         heff_fj = [ heffs['heff_ds_%s' % n] for n in [ 'kj_jetpt2-hp','kj_jetpt2-hp-jeteta2-b0','kj_jetpt2-hp-jeteta2-b1','kj_jetpt2-hp-jeteta2-b2']] 
        heff = [heffs['heff_kj_pt-hp_%s_ds' %x] for x in ['inc','by1','by2','by3','by4']]
        heff[0].SetTitle('#eta-inclusive;p_{T} [GeV]; efficiency')
        heff[1].SetTitle('0.   < |#eta| < 0.75')
        heff[2].SetTitle('0.75 < |#eta| < 1.5')
        heff[3].SetTitle('1.5  < |#eta| < 2.8')
        hratio = H1RatioPlotter(markers=markers)
        hratio.set(*heff)
        hratio.markersize  = 12
        hratio.scalemax    = 1.2
        hratio.ltitle      = 'eff in eta bins'
        hratio.rtitle      = 'high purity'
        hratio.legtextsize = 25
        hratio.legboxsize  = 30
        hratio.legalign    = ('r','b')
        hratio.logx        = True
        hratio.yrange      = (0.,0.9)

        c = hratio.plot()
        for ext in imgext:
            c.Print(prefix+'heff_ratio_kj-hp.%s' % ext )
    except KeyError as ke:
        print 'Not found:',ke

    try:
#         heff_fj = [ heffs['heff_ds_%s' % n] for n in [ 'kj_jetpt2','kj_jetpt2-jeteta2-b0','kj_jetpt2-jeteta2-b1','kj_jetpt2-jeteta2-b2']] 
        heff = [heffs['heff_kj_pt_%s_ds' %x] for x in ['inc','by1','by2','by3','by4']]
        heff[0].SetTitle('#eta-inclusive;p_{T} [GeV]; efficiency')
        heff[1].SetTitle('0.   < |#eta| < 0.75')
        heff[2].SetTitle('0.75 < |#eta| < 1.5')
        heff[3].SetTitle('1.5  < |#eta| < 2.8')
        hratio = H1RatioPlotter(markers=markers)
        hratio.set(*heff)
        hratio.markersize  = 12
        hratio.scalemax    = 1.2
        hratio.ltitle      = 'eff in eta bins'
        hratio.legtextsize = 25
        hratio.legboxsize  = 30
        hratio.legalign    = ('r','b')
        hratio.logx        = True
        hratio.yrange      = (0.,0.9)

        c = hratio.plot()
        for ext in imgext:
            c.Print(prefix+'heff_ratio_kj.%s' % ext )
    except KeyError as ke:
        print 'Not found:',ke

    try:
        # --0jet---
#         heff = [ heffs['heff_ds_%s' % n] for n in [ '1j_jetpt2','1j_jetpt2-jeteta2-b0','1j_jetpt2-jeteta2-b1','1j_jetpt2-jeteta2-b2']] 
        heff = [heffs['heff_0j_pt_%s_ds' %x] for x in ['inc','by1','by2','by3','by4']]
        heff[0].SetTitle('#eta-inclusive;p_{T} [GeV]; efficiency')
        heff[1].SetTitle('0.   < |#eta| < 0.75')
        heff[2].SetTitle('0.75 < |#eta| < 1.5')
        heff[3].SetTitle('1.5  < |#eta| < 2.8')
        hratio = H1RatioPlotter(markers=markers)
        hratio.set(*heff)
        hratio.markersize = 16
        hratio.scalemax = 1.2
        hratio.legtextsize = 25
        hratio.legboxsize  = 30

        c = hratio.plot()
        for ext in imgext:
            c.Print(prefix+'heff_ratio_0j.%s' % ext )
    except KeyError as ke:
        print 'Not found:',ke

    try:
        # --1jet---
#         heff = [ heffs['heff_ds_%s' % n] for n in [ '2j_jetpt2','2j_jetpt2-jeteta2-b0','2j_jetpt2-jeteta2-b1','2j_jetpt2-jeteta2-b2']] 
        heff = [heffs['heff_1j_pt_%s_ds' %x] for x in ['inc','by1','by2','by3','by4']]
        heff[0].SetTitle('#eta-inclusive;p_{T} [GeV]; efficiency')
        heff[1].SetTitle('0.   < |#eta| < 0.75')
        heff[2].SetTitle('0.75 < |#eta| < 1.5')
        heff[3].SetTitle('1.5  < |#eta| < 2.8')
        hratio = H1RatioPlotter(markers=markers)
        hratio.set(*heff)
        hratio.markersize = 16
        hratio.scalemax = 1.
        hratio.legtextsize = 25
        hratio.legboxsize  = 30
        hratio.legalign    = ('r','b')

        c = hratio.plot()
        for ext in imgext:
            c.Print(prefix+'heff_ratio_1j.%s' % ext )
    except KeyError as ke:
        print 'Not found:',ke

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


def makeratio(name, x1_plots,x2_plots, **kwargs): #, ltitle, rtitle, legtextsize, legboxsize, scalemax, markersize ):

    colors  = [ROOT.kRed+1      , ROOT.kAzure-5   ]
    markers = [ROOT.kFullCircle , ROOT.kFullCircle]
    x1_mc  = thsum(x1_plots, 'x1_top', tops)
    x1_ds  = x1_plots['Data'].Clone('x1_ds')
    x1_ds -= thsum(x1_plots, 'x1_others', others)

    x2_mc  = thsum(x2_plots, 'x2_top', tops)
    x2_ds  = x2_plots['Data'].Clone('x2_ds')
    x2_ds -= thsum(x2_plots, 'x2_others', others)

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

#---
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

        # pt ratio
        'jet_2j_btotal_jetpt1'      : ('jetpt1'        , 'bveto-mu' ,  'bveto_ip && njet == 2 && jettche2 > 2.1'                     , ( 34, 30, 200 ) , 'pt_{j2}' ) ,
        'jet_2j_btag_jetpt1'        : ('jetpt1'        , 'bveto-mu' ,  'bveto_ip && njet == 2 && jettche2 > 2.1 && jettche1 >  2.1'  , ( 34, 30, 200 ) , 'pt_{j2}' ) ,
        'jet_2j_bveto_jetpt1'       : ('jetpt1'        , 'bveto-mu' ,  'bveto_ip && njet == 2 && jettche2 > 2.1 && jettche1 <= 2.1'  , ( 34, 30, 200 ) , 'pt_{j2}' ) ,

        'jet_2j_btotal_jetpt2'      : ('jetpt2'        , 'bveto-mu' ,  'bveto_ip && njet == 2 && jettche1 > 2.1'                     , ( 34, 30, 200 ) , 'pt_{j2}' ) ,
        'jet_2j_btag_jetpt2'        : ('jetpt2'        , 'bveto-mu' ,  'bveto_ip && njet == 2 && jettche1 > 2.1 && jettche2 >  2.1'  , ( 34, 30, 200 ) , 'pt_{j2}' ) ,
        'jet_2j_bveto_jetpt2'       : ('jetpt2'        , 'bveto-mu' ,  'bveto_ip && njet == 2 && jettche1 > 2.1 && jettche2 <= 2.1'  , ( 34, 30, 200 ) , 'pt_{j2}' ) ,

        'jet_1j_btotal_jetpt1'      : ('jetpt1'        , 'bveto-mu' ,  'bveto_ip && njet == 1'                                       , ( 34, 30, 200 ) , 'pt_{j1}' ) ,
        'jet_1j_btag_jetpt1'        : ('jetpt1'        , 'bveto-mu' ,  'bveto_ip && njet == 1 && jettche1 >  2.1'                    , ( 34, 30, 200 ) , 'pt_{j1}' ) ,
        'jet_1j_bveto_jetpt1'       : ('jetpt1'        , 'bveto-mu' ,  'bveto_ip && njet == 1 && jettche1 <= 2.1'                    , ( 34, 30, 200 ) , 'pt_{j1}' ) ,

        # secondary jet, measurement area
        'jet_kj_btotal_jetpt2'      : ('jetpt2'        , 'bveto-mu' ,  'njet >= 1 && jettche1 > 2.1'                     , ( 39, 5, 200 ) , 'pt_{j2}' ) ,
        'jet_kj_btag_jetpt2'        : ('jetpt2'        , 'bveto-mu' ,  'njet >= 1 && jettche1 > 2.1 && jettche2 >  2.1'  , ( 39, 5, 200 ) , 'pt_{j2}' ) ,
        'jet_kj_bveto_jetpt2'       : ('jetpt2'        , 'bveto-mu' ,  'njet >= 1 && jettche1 > 2.1 && jettche2 <= 2.1'  , ( 39, 5, 200 ) , 'pt_{j2}' ) ,
        # leading jet, application area
        'jet_kj_btotal_jetpt1'      : ('jetpt1'        , 'bveto-mu' ,  'njet >= 0'                                       , ( 39, 5, 200 ) , 'pt_{j1}' ) ,
        'jet_kj_btag_jetpt1'        : ('jetpt1'        , 'bveto-mu' ,  'njet >= 0 && jettche1 >  2.1'                    , ( 39, 5, 200 ) , 'pt_{j1}' ) ,
        'jet_kj_bveto_jetpt1'       : ('jetpt1'        , 'bveto-mu' ,  'njet >= 0 && jettche1 <= 2.1'                    , ( 39, 5, 200 ) , 'pt_{j1}' ) ,


        # secondary jet, measurement area
        'jet_kj_btag_jetpt2-b0'     : ('jetpt2'        , 'bveto-mu' ,  'fabs(jeteta2) < 0.75                        && njet >= 1 && jettche1 > 2.1 && jettche2 >  2.1' , ( 39, 5, 200 ) , 'pt_{j2}' ) ,
        'jet_kj_btag_jetpt2-b1'     : ('jetpt2'        , 'bveto-mu' ,  'fabs(jeteta2) >=0.75 && fabs(jeteta2) < 1.5 && njet >= 1 && jettche1 > 2.1 && jettche2 >  2.1' , ( 39, 5, 200 ) , 'pt_{j2}' ) ,
        'jet_kj_btag_jetpt2-b2'     : ('jetpt2'        , 'bveto-mu' ,  'fabs(jeteta2) >=1.5  && fabs(jeteta2) < 2.8 && njet >= 1 && jettche1 > 2.1 && jettche2 >  2.1' , ( 39, 5, 200 ) , 'pt_{j2}' ) ,
        'jet_kj_btag_jetpt2-b3'     : ('jetpt2'        , 'bveto-mu' ,  'fabs(jeteta2) >=2.8  && fabs(jeteta2) < 5.  && njet >= 1 && jettche1 > 2.1 && jettche2 >  2.1' , ( 39, 5, 200 ) , 'pt_{j2}' ) ,
        # leading jet, application area
        'jet_kj_btag_jetpt1-b0'     : ('jetpt1'        , 'bveto-mu' ,  'fabs(jeteta1) < 0.75                        && njet >= 0 && jettche1 > 2.1'                    , ( 39, 5, 200 ) , 'pt_{j1}' ) ,
        'jet_kj_btag_jetpt1-b1'     : ('jetpt1'        , 'bveto-mu' ,  'fabs(jeteta1) >=0.75 && fabs(jeteta1) < 1.5 && njet >= 0 && jettche1 > 2.1'                    , ( 39, 5, 200 ) , 'pt_{j1}' ) ,
        'jet_kj_btag_jetpt1-b2'     : ('jetpt1'        , 'bveto-mu' ,  'fabs(jeteta1) >=1.5  && fabs(jeteta1) < 2.8 && njet >= 0 && jettche1 > 2.1'                    , ( 40, 5, 200 ) , 'pt_{j1}' ) ,
        'jet_kj_btag_jetpt1-b3'     : ('jetpt1'        , 'bveto-mu' ,  'fabs(jeteta1) >=2.8  && fabs(jeteta1) < 5.  && njet >= 0 && jettche1 > 2.1'                    , ( 39, 5, 200 ) , 'pt_{j1}' ) ,
    }

    # do some filtering
    btags = [
        'jet_2j_btag_jetpt1', 'jet_2j_btag_jetpt2','jet_1j_btag_jetpt1',
        'jet_kj_btag_jetpt2', 'jet_kj_btag_jetpt1',

        'jet_kj_btag_jetpt2-b0','jet_kj_btag_jetpt2-b1', 'jet_kj_btag_jetpt2-b2', 'jet_kj_btag_jetpt2-b3',
        'jet_kj_btag_jetpt1-b0','jet_kj_btag_jetpt1-b1', 'jet_kj_btag_jetpt1-b2', 'jet_kj_btag_jetpt1-b3',
    ]

#     btags = [ n for n in vars if 'btag' in n]

    vars = vars.fromkeys(btags)

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

    #     ____           __
    #    /  _/___  _____/ /
    #    / // __ \/ ___/ / 
    #  _/ // / / / /__/ /  
    # /___/_/ /_/\___/_/   
    #                      

    decorations = {
        'x1tag'      : 'trailing',
        'x2tag'      : 'leading',
#         'ltitle'     : '"hardest" jet pt_{j}^{soft}',
        'rtitle'     : 'lead vs trail',
        'xtitle'     : 'p_{T}^{jet}',
#         'ytitle2'    : '0j / 1j',
        'scalemax'   : 1.2,
        'markersize' : 12,
        'legalign'   : ('r','t'),
    } 

    kwargs = dict(commons,**decorations)
    makeratio('ratio_inc_pt_1st2nd_btag'    , xplots['jet_kj_btag_jetpt2']    , xplots['jet_kj_btag_jetpt1']   , **kwargs)
    makeratio('ratio_inc_pt_1st2nd_btag-b0' , xplots['jet_kj_btag_jetpt2-b0'] , xplots['jet_kj_btag_jetpt1-b0'], **kwargs)
    makeratio('ratio_inc_pt_1st2nd_btag-b1' , xplots['jet_kj_btag_jetpt2-b1'] , xplots['jet_kj_btag_jetpt1-b1'], **kwargs)
    makeratio('ratio_inc_pt_1st2nd_btag-b2' , xplots['jet_kj_btag_jetpt2-b2'] , xplots['jet_kj_btag_jetpt1-b2'], **kwargs)
#     makeratio('ratio_inc_pt_1st2nd_btag-b3' , xplots['jet_kj_btag_jetpt2-b3'] , xplots['jet_kj_btag_jetpt1-b3'], **kwargs)

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
    makeratio('ratio_softjet_pt_0j1j_btag' , xplots['softjet_1j_btag_jetpt2'], xplots['softjet_0j_btag_jetpt1'], **kwargs)

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
    makeratio('ratio_softjet_eta_0j1j_btag', xplots['softjet_1j_btag_jeteta2'], xplots['softjet_0j_btag_jeteta1'], **kwargs)

    #    ___      _      __ 
    #   <  /     (_)__  / /_
    #   / /_____/ / _ \/ __/
    #  / /_____/ /  __/ /_  
    # /_/   __/ /\___/\__/  
    #      /___/            
    decorations = {
        'x1tag'      : 'n_{jet}=2',
        'x2tag'      : 'n_{jet}=1',
        'ltitle'     : 'leading jet pt_{j}',
        'rtitle'     : '1j vs 2j',
        'xtitle'     : 'pt_{j}',
        'ytitle2'    : '1j / 2j',
        'scalemax'   : 1.2,
        'legalign'   : ('r','t'),
    } 

    kwargs = dict(commons,**decorations)
    makeratio('ratio_leadjet_pt_1j2j_btag' , xplots['jet_2j_btag_jetpt1'], xplots['jet_1j_btag_jetpt1'], **kwargs)

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
    makeratio('ratio_jet_pt_1j2j_btag' , xplots['jet_2j_btag_jetpt2'], xplots['jet_1j_btag_jetpt1'], **kwargs)


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
    makeratio('ratio_jet_eta_1j2j_btag', xplots['jet_2j_btag_jeteta2'], xplots['jet_1j_btag_jeteta1'], **kwargs)

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
    print '-'*80

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
    parser.add_option('--eff'            , dest = 'eff'         , help='do efficiency'          , action='store_true' , default=False )
    parser.add_option('--rat'            , dest = 'rat'         , help='do ratios'              , action='store_true' , default=False )
    parser.add_option('--test'           , dest = 'testmode'    , help='run fast to test'       , action='store_true' , default=False )
    parser.add_option('--dora'           , dest = 'dora'        , help='doraemon'               , action='store_true' , default=False )

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
  
