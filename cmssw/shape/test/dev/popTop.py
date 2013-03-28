#!/usr/bin/env python

import HWWAnalysis.Misc.odict as odict
from ginger.analysis import TreeAnalyser,Cut,CutFlow,Latino,Sample
from ginger.tree import Yield
from hwwinfo2g import wwnamedcuts as wwcuts
import hwwlatino
import ROOT
import sys
import logging
import math
import copy
import pdb
import bdb

#orchard = '/shome/mtakahashi/HWW/Tree//ShapeAna/tree_skim_wwmin'
orchard = '/shome/mtakahashi/HWW/Tree/all_2012_53x_195fb_skim/'
orchard = '/shome/mtakahashi/HWW/Tree/ShapeAna/53x_195fb/tree_skim_wwmin/'
orchard = '/shome/thea/HWW/work/dds/trees'
orchard = '/shome/thea/HWW/work/dds/trees/top'

#_______________________________________________________________________
def printbins(*hists):
    if len(hists) == 0:
        return

    names = [h.GetName() for h in hists]
    l = max(map(len,names))

    for i in xrange(hists[0].GetNbinsX()):
        yb = ['%.3f +/- %.3f'% (h.GetBinContent(i+1),h.GetBinError(i+1)) for h in hists]
        print '|'.join([' %-10s: %-30s' %(n,s) for n,s in zip(names,yb)])

#_______________________________________________________________________
#   ______               __  ___      _     
#  /_  __/___  ____     /  |/  /___ _(_)___ 
#   / / / __ \/ __ \   / /|_/ / __ `/ / __ \
#  / / / /_/ / /_/ /  / /  / / /_/ / / / / /
# /_/  \____/ .___/  /_/  /_/\__,_/_/_/ /_/ 
#          /_/                              

#_______________________________________________________________________
def topestimate( opt ):
    print ' -Pop the top 2'

    samples = hwwlatino.samples(125,'8TeV','Data2012','SM','topestimate')
    samples['Data'].weight = '1'
    samples['WJet'].weight = 'baseW*fakeW'
    samples['DYLL'].weight = 'baseW*puW*effW*triggW'
    samples['WW'].weight   = 'baseW*puW*effW*triggW*(1+(mjj>500)*(detajj>3.5))'

    wwflow = CutFlow(wwcuts.wwcommon)

    del wwflow['bveto_mu']
    del wwflow['bveto_ip']
    wwflow['ptll'] = 'ptll>45'

    print '------------------------------'
    for n,c in wwflow.iteritems():
        print '%-30s: %s'% (n,c)

    topflow = CutFlow()

    topflow['base']   = wwflow.string()

    print '-'*80
    print 'Baseline cut'
    print topflow['base']
    print '-'*80
#     topflow['base']=''
    
    analysers = hwwlatino.makeanalysers(samples,orchard,topflow,opt.lumi)
    
    if opt.buffer:
        print '--Buffering-------'
        for n,a in analysers.iteritems():
            a.bufferentries()
            print '  ',n,':',a.entries(),'>>',a.selectedentries(),'...done'


#     eff0j = efftop0j( copy.deepcopy(analysers) )
#     eff0j = (0.489415, 0.048752)
#     estimation0j( copy.deepcopy(analysers), eff0j )

#     eff1j = efftop1j( copy.deepcopy(analysers) )
#     eff1j = (0.647059, 0.004662)
#     estimation1j( copy.deepcopy(analysers), eff1j )


    eff2jdata, eff2jmc = efftopvbf( copy.deepcopy(analysers) )
#     eff2j.Print()
#     eff2j = None
    estimationvbf( copy.deepcopy(analysers), eff2jdata, eff2jmc)

#_______________________________________________________________________
def extrtop(eff,yldtag,yldbkg):
    eff_val, eff_err = eff

    yldsub = yldtag-yldbkg
    val = yldsub.value * (1-eff_val)/eff_val;
    err = math.sqrt(
        ( yldsub.value * eff_err / eff_val**2)**2 +
        ( (1-eff_val)/eff_val * yldbkg.error)**2
    );

    return Yield(val,err)

#_______________________________________________________________________
#    ____        _      __      
#   / __ \      (_)__  / /______
#  / / / /_____/ / _ \/ __/ ___/
# / /_/ /_____/ /  __/ /_(__  ) 
# \____/   __/ /\___/\__/____/  
#         /___/                 

#_______________________________________________________________________
def efftop0j( analysers ):
    #---
    # - measure the btag efficiency for pt<30 GeV

    print '- Top efficiency 0j'

    print '--Updating buffers-------'
    for n,a in analysers.iteritems():
        old = a.selectedentries()
        a.append('ctrl', 'jettche1>2.1')
        a.append('tag',  'softtche>2.1 || !bveto_munj30')
        a.bufferentries()
        print '  ',n,':',old,'>>',a.selectedentries(),'...done'
        
    # ---

    ttdataset = 'dataset==19'
    njet1 = 'njet==1'
    njet0 = 'njet==0'

    print '-'*80

    y_data = analysers['Data'].yieldsflow(njet1)
    effb = y_data['tag'].value/y_data['ctrl'].value

    print 'N_ctrl',y_data['ctrl'],'  N_tag',y_data['tag']
    print 'softtag efficiency:',effb

#     print 'ctrl data: ',analysers['Data'].cutstring('ctrl',extra=njet1)
#     print 'ctrl singlet: ',analysers['Top'].cutstring('ctrl',extra='!('+ttdataset+') && '+njet1)


    y_t = analysers['Top'].yieldsflow(extra='!('+ttdataset+') && '+njet1)
    print 'singlet ctrl = %30s tag = %30s' % (y_t['ctrl'],y_t['tag'])

    y_ww = analysers['WW'].yieldsflow(njet1)
    print 'WW      ctrl = %30s tag = %30s' % (y_ww['ctrl'],y_ww['tag'])

    y_ggww = analysers['ggWW'].yieldsflow(njet1)
    print 'ggWW    ctrl=%30s tag = %30s' % (y_ggww['ctrl'],y_ggww['tag'])

    y_dyll = analysers['DYLL'].yieldsflow(njet1)
    print 'DYLL    ctrl = %30s tag = %30s' % ( y_dyll['ctrl'],y_dyll['tag'])

    y_dytt = analysers['DYTT'].yieldsflow(njet1)
    print 'DYTT    ctrl = %30s tag = %30s' % ( y_dytt['ctrl'],y_dytt['tag'])

    print '-'*80
    n_ctrl = y_data['ctrl']-(y_t['ctrl']+y_ww['ctrl']+y_ggww['ctrl']+y_dyll['ctrl'])
    n_tag  = y_data['tag'] -(y_t['tag'] +y_ww['tag'] +y_ggww['tag'] +y_dyll['tag'] )
    print 'Nctrl_sub = %s  Ntag_sub = %s ' % (n_ctrl,n_tag)

    eff_softtoptag = n_tag.value/n_ctrl.value
    eff_softtoptag_err = math.sqrt(eff_softtoptag*(1-eff_softtoptag)/n_ctrl.value)
    print 'Btag efficiency for 10 < jpt < 30 (bkg sub): %f +/- %f' % (eff_softtoptag,eff_softtoptag_err) 

    y_top = analysers['Top'].yieldsflow(njet0)
    print 'top0j   ctrl = %30s tag = %30s' % (y_top['tag'],y_top['ctrl'])

    y_tt  = analysers['Top'].yieldsflow(ttdataset+' && '+njet0)
    print 'tt0j    ctrl = %30s tag = %30s' % (y_tt['tag'],y_tt['ctrl'])

    fttbar = y_tt['base'].value/y_top['base'].value
    fsinglet = 1.0 - fttbar;

    fttbar_over_singlet_err = 0.17; # generator uncertainty, this is a relative error

#     fttbar_err = 1.0/fttbar * fttbar_over_singlet_err
#     fsinglet_err = 1.0/fsinglet * fttbar_over_singlet_err
# ---
#     fttbar_err = (1-fttbar)**2 * fttbar_over_singlet_err
#     fsinglet_err = fsinglet**2 * fttbar_over_singlet_err
# ---
#     fttbar_err   = (1-fttbar) * fttbar * fttbar_over_singlet_err
#     fsinglet_err = fsinglet * (1-fsinglet) * fttbar_over_singlet_err
# ---
#     # assume that the x-section error is absolute and not relative
#     fttbar_err = math.sqrt(fttbar**2 * (1-fttbar)**2/(fttbar**4-fttbar**2+1))*fttbar_over_singlet_err*(fttbar/fsinglet)
#     fsinglet_err = fttbar_err
# ---
    # assume that the x-section error is absolute and not relative
    fttbar_err   = math.fabs( fttbar_over_singlet_err * (1-fttbar)**2 * fttbar/fsinglet )
    fsinglet_err = fttbar_err
# ---

    print '-'*80
    print 'ftt =   %10.3f +/- %10.3f ' % ( fttbar,fttbar_err )
    print 'fst =   %10.3f +/- %10.3f ' % ( fsinglet,fsinglet_err )

    eff2b_tt = 1 - (1 - eff_softtoptag)**2;
    eff2b_tt_err = 2 * (1-eff2b_tt) * eff_softtoptag_err;

    print 'effsofttoptag %10.3s +/- %10.3s ' % ( eff_softtoptag,eff_softtoptag )
    print 'eff2b_tt      %10.3s +/- %10.3s ' % ( eff2b_tt,eff2b_tt_err )

    print '-'*80

    x=0.37
    eff_top_0j = (fttbar+(fsinglet)*x) * eff2b_tt + fsinglet * (1-x) * eff_softtoptag;
    eff_top_0j_err = math.sqrt( 
        fttbar_err**2 * eff2b_tt**2 + 
        eff2b_tt_err**2 * fttbar**2 + 
        (1-x)**2 * (fsinglet_err**2 * eff_softtoptag**2 +
        eff_softtoptag_err**2 * fsinglet**2)
    )

    print ' error propagation'
    print '  fttbar_err**2 * eff2b_tt**2 =', math.sqrt(fttbar_err**2 * eff2b_tt**2)
    print '  eff2b_tt_err**2 * fttbar**2 =', math.sqrt(eff2b_tt_err**2 * fttbar**2)
    print '  (1-x)**2 * (fsinglet_err**2 * eff_softtoptag**2) =',math.sqrt((1-x)**2 * (fsinglet_err**2 * eff_softtoptag**2))
    print '  (1-x)**2 * (eff_softtoptag_err**2 * fsinglet**2) =',math.sqrt((1-x)**2*(eff_softtoptag_err**2 * fsinglet**2))

    print '-'*80
    print 'eff_top %f +/- %f ' % ( eff_top_0j,eff_top_0j_err )
    print '-'*80

    return (eff_top_0j, eff_top_0j_err)


#_______________________________________________________________________
def estimation0j( analysers, eff ):
    '''TODO: DY data-driven data/mc factor'''

    print '- Top estimation 0j'

    atop = analysers['Top'].clone()
    atop.append('bveto','njet == 0 && bveto')

    for n,a in analysers.iteritems():
        old = a.selectedentries()
        a.append('ctrltop','njet==0 && !bveto')
        a.bufferentries()
        print '  ',n,':',old,'>>',a.selectedentries(),'...done'
    print '- buffering complete'

    nctrl_data  = analysers['Data'].yields()
    nctrl_top   = analysers['Top'].yields()
#     nctrl_dy    = analysers['DYLL'].yields() + analysers['DYTT'].yields()
    nctrl_dy    =  6.91159*analysers['DYLL'].yields()
    nctrl_other = analysers['VV'].yields()   + analysers['Vg'].yields()  + analysers['VgS'].yields()
    nctrl_ww    = analysers['WW'].yields() + analysers['ggWW'].yields()
    nctrl_wjet  = analysers['WJet'].yields()

    nctrl_bkg = nctrl_dy+nctrl_ww+nctrl_wjet+nctrl_other

    nctrl_mm    = analysers['Data'].yields(extra='channel == 0')
    nctrl_ee    = analysers['Data'].yields(extra='channel == 1')
    nctrl_em    = analysers['Data'].yields(extra='channel == 2')
    nctrl_me    = analysers['Data'].yields(extra='channel == 3')
    
    print 'tagged events in data'
    print 'nctrl_mm   ',nctrl_mm
    print 'nctrl_ee   ',nctrl_ee
    print 'nctrl_em   ',nctrl_em
    print 'nctrl_me   ',nctrl_me
    print '-'*80
    print 'nctrl      ',nctrl_data
    print 'nctrl_top  ',nctrl_top
    print 'nctrl_dy   ',nctrl_dy
    print 'nctrl_ww   ',nctrl_ww
    print 'nctrl_wjet ',nctrl_wjet
    print 'nctrl_other',nctrl_other
    print 'nctrl_bkg  ',nctrl_bkg
    
    print 'nctrl_sub  ',nctrl_data-nctrl_bkg
    
    print '-'*80
    print 'top events from data:', extrtop(eff,nctrl_data,nctrl_bkg)
    print 'top events from mc:  ', atop.yields()

    print '-'*80

#_______________________________________________________________________
#    ___      _      __ 
#   <  /     (_)__  / /_
#   / /_____/ / _ \/ __/
#  / /_____/ /  __/ /_  
# /_/   __/ /\___/\__/  
#      /___/            

#_______________________________________________________________________
def efftop1j( analysers ):
    '''
    Measure the top-tag efficiency for the 1j bin case (pt > 30 GeV)
    '''

    print '- Top efficiency 1j'

    print '--Updating buffers-------'
    for n,a in analysers.iteritems():
        old = a.selectedentries()
        a.append('ctrl','njet==2 && bveto_munj30 && softtche<=2.1 && jettche2>2.1')
        a.append('tag','jettche1>2.1')
        a.bufferentries()
        print '  ',n,':',old,'>>',a.selectedentries(),'...done'


    print '-'*80
    y_data = analysers['Data'].yieldsflow()

    print 'N_ctrl',y_data['ctrl'],'  N_tag',y_data['tag']

    eff_top_1j = y_data['tag'].value/y_data['ctrl'].value
    eff_top_1j_err = math.sqrt(eff_top_1j*(1-eff_top_1j)/y_data['ctrl'].value)

    print '-'*80
    print 'eff_top %f +/- %f ' % ( eff_top_1j,eff_top_1j_err )
    print '-'*80

    return (eff_top_1j, eff_top_1j_err)

#_______________________________________________________________________
def estimation1j( analysers, eff ):

    print '- Top estimation 1j'

    atop = analysers['Top'].clone()
    atop.append('bveto','njet == 1 && bveto_mu && bveto_ip && nbjettche==0')

    for n,a in analysers.iteritems():
        old = a.selectedentries()
        a.append('ctrltop','njet==1 && bveto_munj30 && jettche1>2.1 && softtche<=2.1')
        a.bufferentries()
        print '  ',n,':',old,'>>',a.selectedentries(),'...done'
    print '- buffering complete'

    nctrl_data  = analysers['Data'].yields()
    nctrl_top   = analysers['Top'].yields()
#     nctrl_dy    = analysers['DYLL'].yields() + analysers['DYTT'].yields()
    nctrl_dy    =  3.73366*analysers['DYLL'].yields()
    nctrl_other = analysers['VV'].yields()   + analysers['Vg'].yields()  + analysers['VgS'].yields()
    nctrl_ww    = analysers['WW'].yields() + analysers['ggWW'].yields()
    nctrl_wjet  = analysers['WJet'].yields()

    nctrl_bkg = nctrl_dy+nctrl_ww+nctrl_wjet+nctrl_other

    nctrl_mm    = analysers['Data'].yields(extra='channel == 0')
    nctrl_ee    = analysers['Data'].yields(extra='channel == 1')
    nctrl_em    = analysers['Data'].yields(extra='channel == 2')
    nctrl_me    = analysers['Data'].yields(extra='channel == 3')
    
    print 'tagged events in data'
    print 'nctrl_mm   ',nctrl_mm
    print 'nctrl_ee   ',nctrl_ee
    print 'nctrl_em   ',nctrl_em
    print 'nctrl_me   ',nctrl_me
    print '-'*80
    print 'nctrl      ',nctrl_data
    print 'nctrl_top  ',nctrl_top
    print 'nctrl_dy   ',nctrl_dy
    print 'nctrl_ww   ',nctrl_ww
    print 'nctrl_wjet ',nctrl_wjet
    print 'nctrl_other',nctrl_other
    print 'nctrl_bkg  ',nctrl_bkg
    
    print 'nctrl_sub  ',nctrl_data-nctrl_bkg
    
    print '-'*80
    print 'top events from data:', extrtop(eff,nctrl_data,nctrl_bkg)
    print 'top events from mc:  ', atop.yields()

    print '-'*80

#_______________________________________________________________________
#  _    ______  ______            ___     _ 
# | |  / / __ )/ ____/           |__ \   (_)
# | | / / __  / /_     ______    __/ /  / / 
# | |/ / /_/ / __/    /_____/   / __/  / /  
# |___/_____/_/                /____/_/ /   
#                                  /___/    

#_______________________________________________________________________
def estimationvbf( analysers, heff_ds, heff_mc ):
    print '- Top estimation vbf'
    
    data = analysers['Data']
    top  = analysers['Top']

    print data.cuts.keys()[-1]


    # b-veto+b-tag
    bctrl  = 'nbjettche==0   || (nbjettche==1 && ((abs(jeteta1)<abs(jeteta2)  && jettche1>2.10) || (abs(jeteta1)>=abs(jeteta2) && jettche2>2.10)))'

    # b-tag
    btag   = 'nbjettche==1 && ( ( abs(jeteta1)<abs(jeteta2)  && jettche1>2.10 ) || (abs(jeteta1)>=abs(jeteta2) && jettche2>2.10) )'

    # b-veto
    bveto  = 'nbjettche==0 && (!(( abs(jeteta1)<abs(jeteta2)  && jettche1>2.10 ) || (abs(jeteta1)>=abs(jeteta2) && jettche2>2.10)) )'

    vbfflow = CutFlow( wwcuts.vbfcb( 160 ) )

    vbfbctrl = CutFlow()
    vbfbctrl['softbveto'] = 'bveto_mu && bveto_ip'
    vbfbctrl['vbf-base']  = wwcuts.vbf
    vbfbctrl['vbf-level'] = vbfflow.string()
    vbfbctrl['bctrl']     = bctrl
    
    for n,a in analysers.iteritems():
        old = a.selectedentries()
        a.update(vbfbctrl)
        a.bufferentries()
        print '  ',n,':',old,'>>',a.selectedentries(),'...done'
    print '- buffering complete'



#     import itertools
#     for i,((kd,yd),(kt,yt)) in enumerate(itertools.izip(data.yieldsflow().iteritems(),top.yieldsflow().iteritems())):
#         print '%2d %10s| %-30s | %-30s ' % (i,kd,yd,yt)
        
    # clever binning
    etabins = ([0,0.5,1.,1.5,2.5,5.0],)

    
    regions = ['btag','bveto']
    plots = odict.OrderedDict(zip(regions,[{} for _ in xrange(len(regions))]))

    cjetaexpr = 'abs(jeteta1)*( abs(jeteta1)<abs(jeteta2) ) + abs(jeteta2)*( abs(jeteta1)>=abs(jeteta2) )'
    for p,a in analysers.iteritems():
        print '->',p
        cjeta_btag  = a.views['bctrl'].plot(p+'_btag' ,cjetaexpr,extra=btag, bins=etabins)
        cjeta_bveto = a.views['bctrl'].plot(p+'_bveto',cjetaexpr,extra=bveto,bins=etabins)

        plots['btag'][p]  = cjeta_btag
        plots['bveto'][p] = cjeta_bveto

    hothers_tag = plots['btag']['Data'].Clone('other_tag')
    hothers_tag.Reset()
    hothers_tag.Add(plots['btag']['WW'])
    hothers_tag.Add(plots['btag']['ggWW'])
    hothers_tag.Add(plots['btag']['WJet'])
    hothers_tag.Add(plots['btag']['DYLL'])
    hothers_tag.Add(plots['btag']['DYTT'])
    hothers_tag.Add(plots['btag']['VV'])
    hothers_tag.Add(plots['btag']['Vg'])
    hothers_tag.Add(plots['btag']['VgS']) 

    hdsub_tag = plots['btag']['Data'].Clone('datasub_tag')
    hdsub_tag.Add(hothers_tag,-1)

    print'-'*80
    printbins(plots['btag']['Data'],  hothers_tag, hdsub_tag, plots['btag']['Top'], heff_ds)
    print'-'*80
    printbins(plots['bveto']['Data'], plots['bveto']['Top'], heff_ds)
    print'-'*80


    # da fare

#     hdsub_tag*(1-heff)/(heff)


#_______________________________________________________________________
def efftopvbf( analysers ):
    '''
    Measure the top-tag efficiency for the 2j bin case (pt > 30 GeV)
    '''


    print '- Top efficiency vbf'

    # b-veto+b-tag
    bctrl  = 'nbjettche==0   || (nbjettche==1 && ((abs(jeteta1)<abs(jeteta2)  && jettche1>2.10) || (abs(jeteta1)>=abs(jeteta2) && jettche2>2.10)))'

    # b-tag
    btag   = 'nbjettche==1 && ( ( abs(jeteta1)<abs(jeteta2)  && jettche1>2.10 ) || (abs(jeteta1)>=abs(jeteta2) && jettche2>2.10) )'

    # b-veto
    bveto  = 'nbjettche==0 && (!(( abs(jeteta1)<abs(jeteta2)  && jettche1>2.10 ) || (abs(jeteta1)>=abs(jeteta2) && jettche2>2.10)) )'


    vbfbeff = CutFlow()
    vbfbeff['df']        = '!sameflav'
    vbfbeff['softbveto'] = 'bveto_mu && bveto_ip'
    vbfbeff['vbf']       = wwcuts.vbf
    vbfbeff['bctrl']     = bctrl
    vbfbeff['btag']      = btag

    print '--Updating buffers-------'
    for n,a in analysers.iteritems():
        old = a.selectedentries()
        a.update(vbfbeff)
        a.bufferentries()
        print '  ',n,':',old,'>>',a.selectedentries(),'...done'

    print '-'*80
    y_data = analysers['Data'].yieldsflow()
    print '          N_ctrl',y_data['bctrl'],'  N_tag',y_data['btag']


    print analysers['Data'].cuts.string()

    etabins = ([0,0.5,1.,1.5,2.5,5.0],)
#     etabins = ([0,0.1,0.2,0.3,0.4,0.5,1.,1.5,2.5],)
#     etabins = (6,0.,3.)

    cuts = analysers['Data'].cuts
    plots = odict.OrderedDict(zip(cuts.keys(),[{} for _ in xrange(len(cuts))]))

    cjetaexpr = 'abs(jeteta1)*( abs(jeteta1)<abs(jeteta2) ) + abs(jeteta2)*( abs(jeteta1)>=abs(jeteta2) )'
    for p,a in analysers.iteritems():
        print '->',p
        cjeta = a.plotsflow(p,cjetaexpr,bins=etabins)

        for c,h in cjeta.iteritems():
            plots[c][p] = h

    for i,(k,p) in enumerate(plots.iteritems()):
        hwwlatino.printplots(p,'plot_%d_%s' % (i,k), xlabel='#eta_{CJ}',lumi=19.468)

    heff_data = plots['btag']['Data'].Clone('heff_data')
    heff_data.Divide(plots['bctrl']['Data'])

    hctrl_ds  =  plots['bctrl']['Data'].Clone('hctrl_ds')
    hctrl_ds.Add(plots['bctrl']['WW'],-1)
    hctrl_ds.Add(plots['bctrl']['ggWW'],-1)
    hctrl_ds.Add(plots['bctrl']['WJet'],-1)
    hctrl_ds.Add(plots['bctrl']['DYLL'],-1)
    hctrl_ds.Add(plots['bctrl']['DYTT'],-1)
    hctrl_ds.Add(plots['bctrl']['VV'],-1)
    hctrl_ds.Add(plots['bctrl']['Vg'],-1)
    hctrl_ds.Add(plots['bctrl']['VgS'],-1)

    htag_ds  =  plots['btag']['Data'].Clone('htag_ds')
    htag_ds.Add(plots['btag']['WW'],-1)
    htag_ds.Add(plots['btag']['ggWW'],-1)
    htag_ds.Add(plots['btag']['WJet'],-1)
    htag_ds.Add(plots['btag']['DYLL'],-1)
    htag_ds.Add(plots['btag']['DYTT'],-1)
    htag_ds.Add(plots['btag']['VV'],-1)
    htag_ds.Add(plots['btag']['Vg'],-1)
    htag_ds.Add(plots['btag']['VgS'],-1)

    heff_ds = htag_ds.Clone('heff_ds')
    heff_ds.Divide(hctrl_ds)

    print '-'*80
    print 'bctrl ds',hctrl_ds.Integral()
    print 'bctrl mc',plots['bctrl']['Top'].Integral()

    heff_mc = plots['btag']['Top'].Clone('heff_mc')
    heff_mc.Divide(plots['bctrl']['Top'])


    hothers_ctrl = plots['bctrl']['Data'].Clone('other_ctrl')
    hothers_ctrl.Reset()
    hothers_ctrl.Add(plots['bctrl']['WW'])
    hothers_ctrl.Add(plots['bctrl']['ggWW'])
    hothers_ctrl.Add(plots['bctrl']['WJet'])
    hothers_ctrl.Add(plots['bctrl']['DYLL'])
    hothers_ctrl.Add(plots['bctrl']['DYTT'])
    hothers_ctrl.Add(plots['bctrl']['VV'])
    hothers_ctrl.Add(plots['bctrl']['Vg'])
    hothers_ctrl.Add(plots['bctrl']['VgS'])

    hothers_tag = plots['btag']['Data'].Clone('other_tag')
    hothers_tag.Reset()
    hothers_tag.Add(plots['btag']['WW'])
    hothers_tag.Add(plots['btag']['ggWW'])
    hothers_tag.Add(plots['btag']['WJet'])
    hothers_tag.Add(plots['btag']['DYLL'])
    hothers_tag.Add(plots['btag']['DYTT'])
    hothers_tag.Add(plots['btag']['VV'])
    hothers_tag.Add(plots['btag']['Vg'])
    hothers_tag.Add(plots['bctrl']['VgS'])

    hdsub_ctrl = plots['bctrl']['Data'].Clone('datasub_ctrl')
    hdsub_ctrl.Add(hothers_ctrl,-1)

    hdsub_tag = plots['btag']['Data'].Clone('datasub_tag')
    hdsub_tag.Add(hothers_tag,-1)

    print'-'*80
    printbins(plots['bctrl']['Data'],hothers_ctrl,hdsub_ctrl, plots['bctrl']['Top'])

    print'-'*80
    printbins(plots['btag']['Data'], hothers_tag,hdsub_tag, plots['btag']['Top'])

    print'-'*80
    printbins(heff_data, heff_ds, heff_mc)

    return heff_ds,heff_mc
    

#_______________________________________________________________________
def plotctrlregion():

    import os
    cmsbase = os.getenv('CMSSW_BASE')
    mypath = os.path.join(cmsbase,'src/HWWAnalysis/ShapeAnalysis/macros')

    ROOT.gInterpreter.ExecuteMacro(mypath+'/LatinoStyle2.C')
    loaded = -1
    try:
        loaded = ROOT.gROOT.LoadMacro(mypath+'/MWLPlot.C+g')
    except RuntimeError:
        loaded = ROOT.gROOT.LoadMacro(mypath+'/MWLPlot.C++g')



    wwflow = CutFlow(wwcuts.wwcommon)
    wwflow.insert(0,'of','!sameflav')

    del wwflow['bveto_mu']
    del wwflow['bveto_ip']
    wwflow.collapse('base')
    wwflow['2jet'] = 'njet == 2'
    # wwflow['top bctrl'] ='(njet==1 && jettche1>2.1)'
    #---
#     wwflow['top bctrl'] ='(njet==1 && jettche1>2.1)'
#     wwflow.collapse('top den')
#     wwflow['top num']  ='(softtche>2.1 || !bveto_munj30)'
    for n,c in wwflow.iteritems():
        print '%-30s: %s'% (n,c)
    
    latinos = hwwsamples2g.samples(125,'8TeV','Data2012','SM','topestimate')

    analysers = makeanalysers(latinos,orchard,wwflow)

    name = 'histo'
    varexp  = 'jettche1'
    binning = (100,-20,30)


    vars = [
        ( 'jettche1', (100,-20,30) ),
        ( 'jettche2', (100,-20,30) ),
    ]

    import time

    for varexp,binning in vars:
        start = time.time()

        plots = odict.OrderedDict() 

        for c in wwflow:
            plots[c] = {}

        lumi=19.468;

        for n,a in analysers.iteritems():
            a.lumi = lumi if n != 'Data' else 1
            print n
            hists = a.plotsflow(name+'_'+n,varexp, bins=binning)
            for c,h in hists.iteritems():
                plots[c][n] = h
                print '%s %-20s >> %10d : %-10f' % (c,n,h.GetEntries(), h.Integral())


        for c,l in plots.iteritems():
            printplots(l,'plot_%s_%s'% (varexp,c), lumi=lumi, xlabel=varexp)
        
        elapsed = (time.time()-start)

        print 'Elapsed',elapsed


if __name__ == '__main__':
    import optparse

    
    parser = optparse.OptionParser()
    parser.add_option('-d', '--debug'    , dest='debug'       , help='Debug level'            , default=0 , action='count' )
    parser.add_option('-n', '--nobuffer' , dest='buffer'      , help='buffer common entries'  , default=0 , action='store_false' )
    parser.add_option('-l', '--lumi'     , dest='lumi'        , help='Luminosity'             , default=19.468 )
    (opt, args) = parser.parse_args()

    if not opt.debug:
        pass
    elif opt.debug >= 2:
        print 'Logging level set to DEBUG (%d)' % opt.debug
        logging.basicConfig(level=logging.DEBUG)
    elif opt.debug == 1:
        print 'Logging level set to INFO (%d)' % opt.debug
        logging.basicConfig(level=logging.INFO)

    import os.path
    shape_path = os.path.join(os.getenv('CMSSW_BASE'),'src/HWWAnalysis/ShapeAnalysis')
    print 'Shape directory is',shape_path
    ROOT.gInterpreter.ExecuteMacro(shape_path+'/macros/LatinoStyle2.C')

    try:
        topestimate( opt )
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


