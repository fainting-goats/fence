#!/usr/bin/env python

import HWWAnalysis.Misc.odict as odict
#from ginger.analysis import TreeAnalyser,Cut,CutFlow,Latino,Sample
from ginger.analysis import CutFlow
from ginger.core import ValErr,Yield
from hwwinfo2g import wwnamedcuts as wwcuts
import hwwlatino
import ROOT
import sys
import logging
import math
import copy
import pdb
import bdb
import hwwtools
import ctypes
import os

from HWWAnalysis.Misc.ROOTAndUtils import TStyleSentry,PadPrinter,Tee

from ginger.plotter import H1RatioPlotter
from ginger.painter import Canvas,Pad,Legend,EmbedPad
from ginger.filters import UnderOverTucker

# move somewhere safe
from commons import AlienDict

#orchard = '/shome/mtakahashi/HWW/Tree//ShapeAna/tree_skim_wwmin'
orchard = '/shome/mtakahashi/HWW/Tree/all_2012_53x_195fb_skim/'
orchard = '/shome/mtakahashi/HWW/Tree/ShapeAna/53x_195fb/tree_skim_wwmin/'
orchard = '/shome/thea/HWW/work/dds/trees'
orchard = '/shome/thea/HWW/work/dds/trees/top'

'''
--- some definitions
The top estimate is based on the ABCD schema. Here we define the regions according to the following sketch
            _________
            |   |   |
measurement | A | B | (top enriched)
            |___|___|
            |   |   |
estimation  | C | D | WW level, pre-btag
            |___|___|

The definition above applies to all the estimate algorihms in this script.
'''

#_______________________________________________________________________
#     __  __     __
#    / / / /__  / /___  ___  __________
#   / /_/ / _ \/ / __ \/ _ \/ ___/ ___/
#  / __  /  __/ / /_/ /  __/ /  (__  )
# /_/ /_/\___/_/ .___/\___/_/  /____/
#             /_/
#_______________________________________________________________________
def printbins(*hists):
    if len(hists) == 0:
        return

    names = [h.GetName() for h in hists]

    for i in xrange(hists[0].GetNbinsX()):
        yb = ['%.3f +/- %.3f'% (h.GetBinContent(i+1),h.GetBinError(i+1)) for h in hists]
        print '|'.join([' %-10s: %-30s' %(n,s) for n,s in zip(names,yb)])

#_______________________________________________________________________
def thsameminmax( *hists, **kwargs):
    nh = len(hists)
    # appli max min or calculate the global limits
    themax = kwargs.get('max',max( h.GetMaximum() for h in hists ) )
    themin = kwargs.get('min',max( h.GetMinimum() for h in hists ) )
    map(ROOT.TH1.SetMaximum, hists, [themax]*nh)
    map(ROOT.TH1.SetMinimum, hists, [themin]*nh)

    #print 'minmax',' '.join(['%s(%f,%f)' % (h.GetName(),h.GetMinimum(),h.GetMaximum()) for h in hists]),themin,themax
    return themin,themax

#_______________________________________________________________________
def thsum( plots, newname, processes ):

#     print plots
#     print newname
#     print processes
    if not processes: raise ValueError('Processes is empty!')

    hnew = plots[processes[0]].Clone(newname)
    hnew.Reset()
    for p in processes:
        hnew += plots[p]

    return hnew

#_______________________________________________________________________
def th2yield( h, binx=None, biny=None ):
    err = ctypes.c_double(0.)
    binx1,binx2 = binx if binx else (1,h.GetNbinsX())
    biny1,biny2 = biny if biny else (1,h.GetNbinsY())
    return Yield(h.IntegralAndError(binx1,binx2,biny1,biny2,err),err.value)

#_______________________________________________________________________
def extrtop(eff,yldtag,yldbkg):
    #eff_val, eff_err = eff

    yldsub = yldtag-yldbkg
    val = yldsub.value * (1-eff.value)/eff.value;
    err = math.sqrt(
        ( yldsub.value * eff.error / eff.value**2)**2 +
        ( (1-eff.value)/eff.value * yldbkg.error)**2
    )
#     err = math.sqrt(
#         ( yldsub.value * eff_err / eff_val**2)**2 +
#         ( (1-eff_val)/eff_val )**2 * (yldbkg.error**2+yldtag**2)
#     )

    return Yield(val,err)

#_______________________________________________________________________
def eff2alpha( eff, name='alpha', title=None ):

    if title is None: title = name
    alpha = eff.Clone()
    alpha.Reset()
    alpha.SetNameTitle(name,title)

    for i in xrange(eff.GetSize()):

        bc = eff.GetBinContent(i)
        be = eff.GetBinError(i)

        if bc == 0: continue

        a = (1-bc)/bc
        e = math.fabs(be/bc**2)

        alpha.SetBinContent(i,a)
        alpha.SetBinError(i,e)


    return alpha

#_______________________________________________________________________
def applyalpha( alpha, btag ):

    bveto = btag.Clone('bveto')
    bveto.SetTitle('bveto estimation')
    bveto.Reset()

    ptax = btag.GetXaxis()
    etax = btag.GetYaxis()

    bin_pt  = ctypes.c_int()
    bin_eta = ctypes.c_int()
    bin_z   = ctypes.c_int()

    pt  = ctypes.c_double()
    eta = ctypes.c_double()


    #print 'bveto empty:',th2yield(bveto)

    for i in xrange(btag.GetSize()):
        # bin to bin_pt bin_eta
        btag.GetBinXYZ(i,bin_pt,bin_eta,bin_z)

        # bins to pt,eta
        pt  = ptax.GetBinCenter(bin_pt.value)
        eta = etax.GetBinCenter(bin_eta.value)

        # pt,eta to alpha, alphaerr
        #ac = alpha.GetBinContent(bin_pt.value, bin_eta.value)
        #ae = alpha.GetBinError(bin_pt.value, bin_eta.value)
        bin_a = alpha.FindBin(pt,eta)
        ac = alpha.GetBinContent(bin_a)
        ae = alpha.GetBinError(bin_a)

        # btag content
        tc = btag.GetBinContent(i)
        te = btag.GetBinError(i)

        # bveto estimation
        vc = ac*tc
        ve = math.sqrt((ae/ac)**2 + (te/tc)**2 )*vc if (ac!=0 and tc!=0) else 0.
        #print ','.join( str(x).ljust(20) for x in [vc,ve,ac,ae,tc,te])

        # bveto = alpha*btag
        # bveto_err = ((e_alpha/alpha)**2 + (e_btag/btag)**2 )*bveto

        bveto.SetBinContent(i,vc)
        bveto.SetBinError(i,ve)

    #print 'bveto estimate:',th2yield(bveto)

    return bveto

#_______________________________________________________________________
def th2slices( th2, drawopt='',**opts ):

    #print th2.GetName(),th2.GetMinimumStored(),th2.GetMaximumStored()

    slices = []

    yax = th2.GetYaxis()
    for i in xrange(1,yax.GetNbins()+1):
        h = th2.ProjectionX('%s_%d' % (th2.GetName(),i),i,i,'e')
        h.SetTitle( '%.1f < #eta < %.1f' % (yax.GetBinLowEdge(i),yax.GetBinUpEdge(i)) )
        h.SetMinimum(th2.GetMinimumStored())
        h.SetMaximum(th2.GetMaximumStored())
        slices.append(h)

    # if not defined as kwarg, maintain the minmax from the th2
    #if 'yrange' not in opts:
        #opts['yrange'] = (th2.GetMinimum(),th2.GetMaximum())
        #print 'yrange',opts['yrange']

    hratio = H1RatioPlotter(**opts)
    hratio.set(*slices)
    c = hratio.plot(drawopt)
    c.parent = hratio
    return c


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

#     samples = hwwlatino.samples(125,'8TeV','Data2012','SM','topestimate')
    samples = hwwlatino.samples(125,'8TeV','Data2012','SM','topplots')
    # there are for 2j
#     samples['Data'].weight = '1'
#     samples['WJet'].weight = 'baseW*fakeW'
#     samples['DYLL'].weight = 'baseW*puW*effW*triggW'
#     samples['WW'].weight   = 'baseW*puW*effW*triggW*(1+(mjj>500)*(detajj>3.5))'

    if opt.closure:
        montecarlos = [ s for s in samples if s != 'Data' ]
        pseudo = []
        for s in montecarlos:
            pseudo.extend(samples[s])
        samples['PseudoData'] = pseudo

    if opt.onetop:
        tops = ['ttbar','tW']
        mytop = []
        for s in tops:
            mytop.extend(samples[s])
            samples.pop(s)

        samples['Top'] = mytop

    wwflow = CutFlow(wwcuts.wwcommon)

    del wwflow['bveto_mu']
    del wwflow['bveto_ip']
    wwflow['ptll'] = 'ptll>45'

    print '-'*80
    for n,c in wwflow.iteritems():
        print '%-30s: %s'% (n,c)

    topflow = CutFlow()

    topflow['base']   = wwflow.string()

    print '-'*80
    print 'Baseline cut'
    print topflow['base']
    print '-'*80

    analysers = hwwlatino.makeanalysers(samples,orchard,topflow,opt.lumi)

    # replace after building to ensure the lumi scaling to be applied
    if opt.closure: analysers['Data'] = analysers.pop('PseudoData')

    if opt.buffer:
        print '--Buffering-------'
        for n,a in analysers.iteritems():
            a.bufferentries()
            print '  ',n,':',a.entries(),'>>',a.selectedentries(),'...done'

    # are there further options to implement:
    # cuts as ABCD(mu): a dictionary of cuts for the standard regions
    # it's another story for the definitions of tops and others, because they
    # change between estimation and
    abcd1j = {}
    abcd1j['A']   = 'njet==2 && bveto_munj30 && softtche<=2.1 && jettche2>2.1'
    abcd1j['B']   = 'jettche1>2.1'
    abcd1j['C']   = 'njet==1 && bveto_munj30 && softtche<=2.1'
    abcd1j['D']   = 'jettche1>2.1 '
    abcd1j['C-D'] = 'nbjettche==0'
    abcd1j['mu']  = 'bveto_mu'

    abcd1j_mu = {}
    abcd1j_mu['A']   = 'njet==2 && bveto_mu && softtche<=2.1 && jettche2>2.1'
    abcd1j_mu['B']   = 'jettche1>2.1'
    abcd1j_mu['C']   = 'njet==1 && bveto_mu && softtche<=2.1'
    abcd1j_mu['D']   = 'jettche1>2.1 '
    abcd1j_mu['C-D'] = 'nbjettche==0'
    abcd1j_mu['mu']  = 'bveto_mu'

    # select what ABCD definition to use
    regions = abcd1j_mu

    if opt.chan in ['0j','all']:
        eff0j = efftop0j( copy.deepcopy(analysers) )
    #     eff0j = (0.489415, 0.048752)
        if opt.estimate:
            estimation0j( copy.deepcopy(analysers), eff0j )

    if opt.chan in ['1j','all']:
        eff1j = efftop1j( copy.deepcopy(analysers), regions )
    #     eff1j = (0.647059, 0.004662)
        if opt.estimate:
            estimation1j( copy.deepcopy(analysers), eff1j, regions )

    if opt.chan in ['1j2g_lead','all']:
        eff1j = efftop1j2g_lead( copy.deepcopy(analysers), regions, opt )
    #     eff1j = (0.647059, 0.004662)
        if opt.estimate:
            estimation1j2g( copy.deepcopy(analysers), eff1j, regions, opt )

    if opt.chan in ['1j2g_sublead','all']:
        eff1j = efftop1j2g_sublead( copy.deepcopy(analysers), opt )
    #     eff1j = (0.647059, 0.004662)
        if opt.estimate:
            estimation1j2g( copy.deepcopy(analysers), eff1j, opt )

    if opt.chan in ['vbf','all']:
        eff2jdata, eff2jmc = efftopvbf( copy.deepcopy(analysers) )
    #     eff2j.Print()
    #     eff2j = None
        if opt.estimate:
            estimationvbf( copy.deepcopy(analysers), eff2jdata, eff2jmc)

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


    y_t = analysers['Top'].yieldsflow('!('+ttdataset+') && '+njet1)
    print 'tW   ctrl = %15s tag = %15s' % (y_t['ctrl'],y_t['tag'])

    y_ww = analysers['WW'].yieldsflow(njet1)
    print 'WW   ctrl = %15s tag = %15s' % (y_ww['ctrl'],y_ww['tag'])

    y_ggww = analysers['ggWW'].yieldsflow(njet1)
    print 'ggWW ctrl = %15s tag = %15s' % (y_ggww['ctrl'],y_ggww['tag'])

    y_dyll = analysers['DYLL'].yieldsflow(njet1)
    print 'DYLL ctrl = %15s tag = %15s' % ( y_dyll['ctrl'],y_dyll['tag'])

    y_dytt = analysers['DYTT'].yieldsflow(njet1)
    print 'DYTT ctrl = %15s tag = %15s' % ( y_dytt['ctrl'],y_dytt['tag'])

    print '-'*80
    n_ctrl = y_data['ctrl']-(y_t['ctrl']+y_ww['ctrl']+y_ggww['ctrl']+y_dyll['ctrl'])
    n_tag  = y_data['tag'] -(y_t['tag'] +y_ww['tag'] +y_ggww['tag'] +y_dyll['tag'] )
    print 'Nctrl_sub = %s  Ntag_sub = %s ' % (n_ctrl,n_tag)

    eff_softtoptag = n_tag.value/n_ctrl.value
    eff_softtoptag_err = math.sqrt(eff_softtoptag*(1-eff_softtoptag)/n_ctrl.value)
    print 'Btag eff for 10 < jpt < 30 (bkg sub): %s' % Yield(eff_softtoptag,eff_softtoptag_err)

    y_top = analysers['Top'].yieldsflow(njet0)
    print 'top0j   ctrl = %15s tag = %15s' % (y_top['tag'],y_top['ctrl'])

    y_tt  = analysers['Top'].yieldsflow(ttdataset+' && '+njet0)
    print 'tt0j    ctrl = %15s tag = %15s' % (y_tt['tag'],y_tt['ctrl'])

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
    print 'ftt           = %s ' % Yield( fttbar,fttbar_err )
    print 'fst           = %s ' % Yield( fsinglet,fsinglet_err )

    eff2b_tt = 1 - (1 - eff_softtoptag)**2;
    eff2b_tt_err = 2 * (1-eff2b_tt) * eff_softtoptag_err;

    print 'effsofttoptag = %s ' % Yield( eff_softtoptag,eff_softtoptag )
    print 'eff2b_tt      = %s ' % Yield( eff2b_tt,eff2b_tt_err )

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

    return ValErr(eff_top_0j, eff_top_0j_err)


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

    nctrl_mm    = analysers['Data'].yields('channel == 0')
    nctrl_ee    = analysers['Data'].yields('channel == 1')
    nctrl_em    = analysers['Data'].yields('channel == 2')
    nctrl_me    = analysers['Data'].yields('channel == 3')

    print '-'*80
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
def efftop1j2g_sublead( analysers, opt ):
    '''
    Measure the top-tag efficiency for the 1j bin case (pt > 30 GeV)
    '''

    print '- Top efficiency 1j (on subleading) in pt/eta bins'

    print '--Updating buffers-------'
    for n,a in analysers.iteritems():
        old = a.selectedentries()
#         a.append('ctrl','njet==2 && bveto_munj30 && softtche<=2.1 && jettche1>2.1 && ( njet <= 1 || (bveto_ip && njet >= 2) ) && jettche3<=2.1 && jettche4<=2.1')
#         a.append('ctrl','njet==2 && bveto_munj30 && softtche<=2.1 && jettche1>2.1 && jetpt1>50 && jetpt3 < 10')
        #a.append('ctrl','njet==2 && bveto_munj30 && softtche<=2.1 && jettche1>2.1 && jetpt1>50 && jetpt3 < 10')
        #a.append('ctrl','njet==2 && bveto_munj30 && softtche<=2.1 && jettche1>2.1 && jetpt3 < 10')
        #a.append('ctrl','njet==2 && bveto_munj30 && softtche<=2.1 && jettche1>2.1')
        a.append('ctrl','njet==2 && bveto_mu && softtche<=2.1 && jettche1>2.1 && jetpt3 < 10')
        a.append('tag','jettche2>2.1')
        a.bufferentries()
        print '  ',n,':',old,'>>',a.selectedentries(),'...done'

    print '-'*80

    ptbins  = range(30, 70, 10) + range( 70, 150, 20) + range(150,251, 50)
    etabins = [0.,0.75,1.5,2.8,5]
    #etabins = [0.,2.8,5]

    plots = AlienDict()
    for n,a in analysers.iteritems():
        # plot the 2d histograms for the efficiency
        pf = a.plotsflow(n+'_etapt','fabs(jeteta2):jetpt2',bins=(ptbins,etabins))

        for c,h in pf.iteritems():
            # rearrange the plots as cut,process
            plots[c][n] = h

    plots.lock()
    if 'Top' in analysers:
        others = ['DYLL', 'DYTT', 'VV', 'Vg', 'VgS', 'WJet', 'WW', 'ggH', 'ggWW', 'qqH', 'wzttH']
        tops   = ['Top']
    else:
        others = ['DYLL', 'DYTT', 'VV', 'Vg', 'VgS', 'WJet', 'WW', 'ggH', 'ggWW', 'tW', 'qqH', 'wzttH']
        tops   = ['ttbar']

    # assure a proper subtraction
    assert( set(analysers.iterkeys()) == set(others+tops+['Data']) )


    # calculate the efficiency
    plt_A_ds = plots['ctrl']['Data'].Clone('plt_A_ds')
    plt_A_ds.Add( thsum(plots['ctrl'],'ctrl_others',others) ,-1)

    plt_B_ds = plots['tag']['Data'].Clone('plt_B_ds')
    plt_B_ds.Add( thsum(plots['tag'],'tag_others',others) ,-1)

    plt_A_mc = thsum(plots['ctrl'], 'plt_A_mc', tops)
    plt_B_mc  = thsum(plots['tag'] , 'plt_B_mc' , tops)

    eff_ds = plt_B_ds.Clone('eff_1j2g')
    eff_ds.Reset()

    eff_ds.Divide(plt_B_ds,plt_A_ds,1,1,'b')
    eff_ds.SetTitle('b-tag efficiency [data-mc_{other}]')

    eff_mc = plt_B_mc.Clone('eff_1j2g')
    eff_mc.Reset()

    eff_mc.Divide(plt_B_mc,plt_A_mc,1,1,'b')
    eff_mc.SetTitle('b-tag efficiency [mc_{top}]')

    if opt.prefix:
        hwwtools.ensuredir(opt.prefix)

        thsameminmax(plt_A_ds, plt_A_mc)
        thsameminmax(plt_B_ds, plt_B_mc)
        thsameminmax(eff_ds, eff_mc, max=1)


        with TStyleSentry( '2Dplots' ) as sentry:
            ROOT.gStyle.SetPaintTextFormat('2.2f')
            ROOT.gStyle.SetOptLogx()

            padstoprint = {}

            c = ROOT.TCanvas('stoca','stoca',500,750)
            c.Divide(2,3)

            c.cd(1)
            plt_A_ds.Draw('text45 colz')

            c.cd(2)
            plt_A_mc.Draw('text45 colz')

            c.cd(3)
            plt_B_ds.Draw('text45 colz')
            c.cd(4)
            plt_B_mc.Draw('text45 colz')

            padstoprint['eff_ds'] = c.cd(5)
            eff_ds.Draw('e text45 colz')

            padstoprint['eff_mc'] = c.cd(6)
            eff_mc.Draw('e text45 colz')

            c.Print(opt.prefix+'/stocazzo_sub.png')
            c.Print(opt.prefix+'/stocazzo_sub.pdf')


            padstoprint = { 'yield_ds_A':1, 'yield_ds_B':3, 'eff_ds':5, 'eff_mc':6, }

            printer = PadPrinter(opt.prefix)
            printer.savefromcanvas(c,**padstoprint)

    return eff_ds, eff_mc

#_______________________________________________________________________
def efftop1j2g_lead( analysers, regions, opt ):
    '''
    Measure the top-tag efficiency for the 1j bin case (pt > 30 GeV)
    '''
    print '- Top efficiency 1j (on leading) in pt/eta bins'

    if 'Top' in analysers:
        others = ['DYLL', 'DYTT', 'VV', 'Vg', 'VgS', 'WJet', 'WW', 'ggH', 'ggWW', 'qqH', 'wzttH']
        tops   = ['Top']
    else:
        if   opt.mcsub == 0:
            # classic set
            tops   = ['ttbar','tW','DYLL', 'DYTT', 'VV', 'Vg', 'VgS', 'WJet', 'WW', 'ggH', 'ggWW', 'qqH', 'wzttH']
            others = []
        elif opt.mcsub == 1:
            # classic + mc subtraction
            tops   = ['ttbar','tW']
            others = ['DYLL', 'DYTT', 'VV', 'Vg', 'VgS', 'WJet', 'WW', 'ggH', 'ggWW', 'qqH', 'wzttH']
        elif opt.mcsub == 2:
            # mc & tW subtraction
            tops   = ['ttbar']
            others = ['DYLL', 'DYTT', 'VV', 'Vg', 'VgS', 'WJet', 'WW', 'ggH', 'ggWW', 'tW', 'qqH', 'wzttH']
        else:
            raise ValueError('What\'s this shit? '+opt.mcsub)

    # drop the analyzers which are not needed
    analysers = odict.OrderedDict([ (n,a) for n,a in analysers.iteritems() if n in (tops+others+['Data']) ])

    # specific cuts
    cuts = CutFlow()
    cuts['A'] = regions['A']
    cuts['B'] = regions['B']

    # print the cuts because it's nice
    print '--1j cuts--'.ljust(80,'-')
    print '\n'.join(['%-30s: %s'% i for i in cuts.iteritems()])

    # and we want to tuck the under and overflows in
    uo = UnderOverTucker()

    print '--Updating buffers--'.ljust(80,'-')
    for n,a in analysers.iteritems():
        old = a.selectedentries()
        a.update(cuts)
        a.bufferentries()
        # and tuck in underflow/overflows by default
        a.filters.append(uo)
        print '  ',n,':',old,'>>',a.selectedentries(),'...done'

    print '-'*80

    ptbins  = range(30, 70, 10) + range( 70, 150, 20) + range(150,250, 50) + [250,500,1000]
    etabins = [0.,0.75,1.5,2.8,5]
    etabins = [0.,1.4,2.8,5]
    etabins = [0.,2.8,5]
    #ptbins  = range(30,300,10)#+range(300,1001,100)
    ptbins  = range(30,300,10)#+range(300,1001,100)
    #ptbins  = [30,1000]
    #etabins = [0.,2.8,5]

    print 'plots: [' ,
    plots = AlienDict()
    for n,a in analysers.iteritems():
        print '%s,' % n,
        # plot the 2d histograms for the efficiency
        pf = a.plotsflow(n+'_etapt','fabs(jeteta1):jetpt1',bins=(ptbins,etabins))

        for c,h in pf.iteritems():
            # rearrange the plots as cut,process
            plots[c][n] = h
        sys.stdout.flush()

    plots.lock()
    print ']'

    # assure a proper subtraction
    #assert( set(analysers.iterkeys()) == set(others+tops+['Data']) )

    # calculate the efficiency
    # regions for data subtracted
    plt_A_ds = plots['A']['Data'].Clone('plt_A_ds')
    if others: plt_A_ds.Add( thsum(plots['A'],'ctrl_others',others) ,-1)
    plt_A_ds.SetTitle('Data etapt [A]')

    plt_B_ds = plots['B']['Data'].Clone('plt_B_ds')
    if others: plt_B_ds.Add( thsum(plots['B'],'tag_others',others) ,-1)
    plt_B_ds.SetTitle('Data etapt [B]')

    # regions for mc
    plt_A_mc = thsum(plots['A'], 'plt_A_mc', tops)
    plt_B_mc = thsum(plots['B'] , 'plt_B_mc' , tops)

    # eff data-sub
    eff_ds = plt_B_ds.Clone('eff_1j2g')
    eff_ds.Reset()

    eff_ds.Divide(plt_B_ds,plt_A_ds,1,1,'b')
    eff_ds.SetTitle('b-tag efficiency [data-mc_{other}]')

    # and montecarlo
    eff_mc = plt_B_mc.Clone('eff_1j2g')
    eff_mc.Reset()

    eff_mc.Divide(plt_B_mc,plt_A_mc,1,1,'b')
    eff_mc.SetTitle('b-tag efficiency [mc_{top}]')

    # for each component
    eff_tops = {}
    for t in tops:
        eff_t = plots['A'][t].Clone('plt_BA_'+t)
        eff_t.Reset()
        eff_t.Divide(plots['B'][t],plots['A'][t],1,1,'b')
        eff_tops[t] = eff_t

    # x-check on the efficiency, diff -> inc to compare with pure mc
    # calculation
    print '-- diff->inc --'.ljust(80,'-')
    test_ds = plt_A_ds.Clone('test_ds')
    test_ds.Scale(1/test_ds.Integral())
    test_ds.Multiply(eff_ds)
    print 'eff_ds :',th2yield(test_ds)
    del test_ds

    test_mc = plt_A_mc.Clone('test_mc')
    test_mc.Scale(1/test_mc.Integral())
    test_mc.Multiply(eff_mc)
    print 'eff_mc :',th2yield(test_mc)
    del test_mc

    for t in tops:
        test_t = plots['A'][t].Clone('test_'+t)
        if test_t.Integral() == 0:
            eff_t = Yield()
        else:
            test_t.Scale(1/test_t.Integral())
            test_t.Multiply(eff_tops[t])
            eff_t = th2yield(test_t)
        print 'eff_%s :' %t, eff_t
        del test_t


    nA_ds = th2yield(plt_A_ds)
    nB_ds = th2yield(plt_B_ds)
    nA_mc = th2yield(plt_A_mc)
    nB_mc = th2yield(plt_B_mc)

    print '-'*80
    print '-- inclusive efficiency --'.ljust(80,'-')
    print 'N_A_data  :', th2yield( plots['A']['Data'] ),plots['A']['Data'].GetEntries()
    print 'N_B_data  :', th2yield( plots['B']['Data'] ),plots['B']['Data'].GetEntries()

    print 'N_A_others:', th2yield(thsum(plots['A'],'ctrl_others',others) ) if others else Yield()
    print 'N_B_others:', th2yield(thsum(plots['B'],'ctrl_others',others) ) if others else Yield()
    print '---'
    print 'N_A_ds = ',nA_ds
    print 'N_B_ds = ',nB_ds
    print 'N_A_mc = ',nA_mc
    print 'N_B_mc = ',nB_mc

    eff_inc_val = (nB_ds/nA_ds).value
    eff_inc_err = math.sqrt(eff_inc_val*(1-eff_inc_val)/nA_ds.value)
    eff_inc= Yield(eff_inc_val, eff_inc_err)

    print 'eff (inc) =',eff_inc

    print '-'*80

    if opt.prefix:
        with TStyleSentry( '2Dplots' ) as sentry:
            ROOT.gStyle.SetPaintTextFormat('2.2f')
            hwwtools.ensuredir(opt.prefix)

            # save ranges for later
            thsameminmax(plt_A_ds, plt_A_mc)
            thsameminmax(plt_B_ds, plt_B_mc)
            thsameminmax(eff_ds, eff_mc, *(eff_tops.itervalues()),max=1)


            # main plotting style
            style  = {
                    'plotratio'  : False,
                    #'logx'       : True,
                    'morelogx'   : True,
                    'legalign'   : ('r','t'),
                    'rtitle'     : '%.2f fb^{-1}' % opt.lumi,
                    #'legtextsize': 20,
                    'legboxsize' : 30,
                    }

            # plotting style for efficiencies
            styleeff = style.copy()
            styleeff.update({
                    'drawopt'    : 'text',
                    'markersize' : 10,     # need it for text
                    'markers'    : [1]*10, # with a small marker
            })

            canvases = {}
            canvases[0,0] = th2slices(plt_A_ds, scalemax=1.5, ltitle='A region (ds)', **style )
            canvases[1,0] = th2slices(plt_A_mc, scalemax=1.5, ltitle='A region (mc)', **style )
            canvases[0,1] = th2slices(plt_B_ds, scalemax=1.5, ltitle='B region (ds)', **style )
            canvases[1,1] = th2slices(plt_B_mc, scalemax=1.5, ltitle='B region (mc)', **style )
            canvases[0,2] = th2slices(eff_ds,   scalemax=1.2, ltitle='efficiency (ds)', **styleeff )
            canvases[1,2] = th2slices(eff_mc,   scalemax=1.2, ltitle='efficiency (mc)', **styleeff )

            for i,t in enumerate(tops):
                canvases[i+2,0] = th2slices(plots['A'][t], scalemax=1.5, ltitle='A region (%s)' % t, **style)
                canvases[i+2,1] = th2slices(plots['B'][t], scalemax=1.5, ltitle='B region (%s)' % t, **style)
                canvases[i+2,2] = th2slices(eff_tops[t],   scalemax=1.2, ltitle='efficiency (%s)' % t, **styleeff)

            padstoprint = {
                    'yield_ds_A':(0,0),
                    'yield_ds_B':(0,1),
                    'eff_ds':(0,2),
                    'eff_mc':(1,2)
                    }
            for name,p in padstoprint.iteritems():
                c = canvases[p]
                c.Print('%s/%s.pdf' % (opt.prefix,name))
                c.Print('%s/%s.png' % (opt.prefix,name))

            call = Canvas()

            for p,c in canvases.iteritems():
                call[p] = EmbedPad(c)

            call.makecanvas()
            call.applystyle()

            call.SetCanvasSize(call.GetWw()/2,call.GetWh()/2)
            print 'Scaling PS for',max(call.gridsize())
            ROOT.gStyle.SetLineScalePS(1/max(call.gridsize()))
            call.Print(opt.prefix+'/eff_summary.pdf')
            call.Print(opt.prefix+'/eff_summary.png')


    return eff_ds, eff_mc

#_______________________________________________________________________
def estimation1j2g( analysers, eff, regions,  opt ):
    # sanity check
    # note: the ctrl region contains both ttbar and tW, therefore tW goes in tops and not in ttbar
    others = ['DYLL', 'DYTT', 'VV', 'Vg', 'VgS', 'WJet', 'WW', 'ggH', 'ggWW', 'qqH', 'wzttH']
    tops   = ['tW', 'ttbar']
    #others = ['DYLL', 'DYTT', 'VV', 'Vg', 'VgS', 'WJet', 'WW', 'ggH', 'ggWW', 'qqH', 'wzttH','tW']
    #tops   = ['ttbar']

    # insure a proper subtraction
    assert( set(analysers.iterkeys()) == set(others+tops+['Data']) )

    eff_ds, eff_mc = eff

    # specific cuts
    cuts = CutFlow()
    cuts['C'] = regions['C']
    cuts['D'] = regions['D']

    print '--1j cuts--'.ljust(80,'-')
    print '\n'.join(['%-30s: %s'% i for i in cuts.iteritems()])

    # tuck under and overflows in
    uo = UnderOverTucker()

    print '---Updating buffers-'.ljust(80,'-')
    for n,a in analysers.iteritems():
        old = a.selectedentries()
        a.update(cuts)
        a.bufferentries()
        a.filters.append(uo)
        print '  ',n,':',old,'>>',a.selectedentries(),'...done'
    print '- buffering complete'

    xbins = eff_ds.GetXaxis().GetXbins()
    ybins = eff_ds.GetYaxis().GetXbins()

    ptbins  = [ xbins.At(i) for i in xrange(xbins.GetSize()) ]
    etabins = [ ybins.At(i) for i in xrange(ybins.GetSize()) ]

    #npt  = eff_ds.GetNbinsX()
    #neta = eff_ds.GetNbinsY()

    #ptbins  = range(30, 70, 10) + range( 70, 150, 20) + range(150,251, 50)
    #etabins = [0.,0.75,1.5,2.8,5]
    #ptbins  = [30,1000]
    #etabins = [0.,2.8,5]

    npt  = len(ptbins)-1
    neta = len(etabins)-1


    print 'plots: [' ,
    plots = AlienDict()
    for n,a in analysers.iteritems():
        print '%s,' % n,
        # plot the 2d histograms for the efficiency
        pf = a.plotsflow(n+'_etapt','fabs(jeteta1):jetpt1',bins=(ptbins,etabins))

        for c,h in pf.iteritems():
            # rearrange the plots as cut,process
            plots[c][n] = h

        # add the bveto plots, for later
        p = a.views['C'].plot(n+'_etapt','fabs(jeteta1):jetpt1',cut=regions['C-D'],bins=(ptbins,etabins))
        # tuck under and overflows in
        uo(p)
        plots['C-D'][n] = p
        sys.stdout.flush()

    plots.lock()
    print ']'

    print plots.keys()

    plt_C_ds = plots['C']['Data'].Clone('plt_C_ds')
    if others: plt_C_ds.Add( thsum(plots['C'],'pretag_others',others), -1)
    plt_C_ds.SetTitle('WW-lvl - b_{veto} [data-mc_{other}]')

    plt_C_mc = thsum(plots['C'],'plt_C_mc',tops)
    plt_C_mc.SetTitle('WW-lvl - b_{veto} [mc_{top}]')

    plt_D_ds = plots['D']['Data'].Clone('plt_D_ds')
    if others: plt_D_ds.Add( thsum(plots['D'],'btag_others',others), -1)
    plt_D_ds.SetTitle('top control region [btagged,data-mc_{other}]')

    plt_D_mc = thsum(plots['D'],'plt_D_mc',tops)
    plt_D_mc.SetTitle('top control region [btagged,mc_{top}]')

    # same for the btag-ctrl region

    # clean negative bins
    for i in xrange(plt_D_ds.GetSize()):
        bc = plt_D_ds.GetBinContent(i)
        if bc > 0: continue

        plt_D_ds.SetAt(0,i)
        plt_D_ds.SetBinError(i,0)

    # plots of alpha factors
    alpha_ds = eff2alpha( eff_ds, 'alpha_ds' )
    alpha_mc = eff2alpha( eff_mc, 'alpha_mc' )

    # applied to the b-tagged region
    plt_CD_ds = applyalpha( alpha_ds, plt_D_ds)
    plt_CD_ds.SetTitle('WW-lvl t#bar{t}+tW estimate [data-mc_{other}]')

    plt_CD_mc = applyalpha( alpha_mc, plt_D_mc)
    plt_CD_mc.SetTitle('WW-lvl t#bar{t}+tW estimate [mc_{top}]')

    # while extracting the yield in the CD region directly
    plt_CD_mc_topwwlvl = thsum(plots['C-D'],'plt_CD_mc',tops)
    plt_CD_mc_topwwlvl.SetTitle('top mc, WW-level yield [C-D,mc_{top}]')

    # so this is the transfer factor calculated on the mc from eta 0-2.8 -> 0.5
    beta= th2yield(plt_C_mc, biny=(neta,neta))/th2yield(plt_C_mc, biny=(1,neta-1))
    print '-'*80
    print 'beta',beta
    print '-'*80

    y_plt_D_ds  = th2yield(plt_D_ds)
    y_plt_CD_ds = th2yield(plt_CD_ds)

    print '1 --D yields [0,5]--'.ljust(80,'-')
    print '[raw]   '.ljust(15), th2yield(plots['D']['Data'])
    print '[others]'.ljust(15), th2yield(thsum(plots['D'],'pretag_others',others))
    print '[ds]    '.ljust(15), y_plt_D_ds
    for s in tops:
        print ('['+s+']').ljust(15), th2yield(plots['D'][s])
    print '[mc]'.ljust(15), th2yield(plt_D_mc)
    print '2 --C-D yields [0,2.8]--'.ljust(80,'-')
    print '[ds]'.ljust(15), y_plt_CD_ds
    print '[mc]'.ljust(15), th2yield(plt_CD_mc_topwwlvl,biny=(1,neta-1))
    print '3 --C yields [0,5]--'.ljust(80,'-')
    print '[est,ds]'.ljust(15), (y_plt_D_ds+y_plt_CD_ds)*(1+beta)
    print '[puremc]'.ljust(15), th2yield(plt_C_mc)
    print '4 --C-D yields [0,5]--'.ljust(80,'-')
    print '[est,ds]'.ljust(15), beta*y_plt_D_ds+(1+beta)*y_plt_CD_ds
    print '[puremc]'.ljust(15), th2yield(plt_CD_mc_topwwlvl)
    print '-'*80


    if opt.prefix:
        with TStyleSentry( '2Dplots' ) as sentry:
            ROOT.gStyle.SetPaintTextFormat('2.2f')
            hwwtools.ensuredir(opt.prefix)

            thsameminmax( plt_C_ds, plt_C_mc )
            thsameminmax( plt_D_ds, plt_D_mc )
            thsameminmax( alpha_ds, alpha_mc, max=2. )
            thsameminmax( plt_CD_ds,plt_CD_mc,plt_CD_mc_topwwlvl )

            plt_DC_mc = plt_D_mc.Clone('plt_DC_mc')
            plt_DC_mc.SetTitle('C/D regions vs p_{T}')
            plt_DC_mc.Reset()
            plt_DC_mc.Divide(plt_D_mc,plt_C_mc,1,1,'b')
            plt_DC_mc.SetMaximum(1)

            plt_DC_tops = {}
            for t in tops:
                plt_t = plots['C'][t].Clone('plt_DC_'+t)
                plt_t.Reset()
                plt_t.Divide(plots['D'][t],plots['C'][t],1,1,'b')
                plt_DC_tops[t] = plt_t

            #save minmax for later
            thsameminmax(plt_DC_mc, *(plt_DC_tops.itervalues()))

            heff_rat = plt_DC_mc.Clone('heff_rat')
            heff_rat.SetTitle('efficiency / (btag/all)')
            heff_rat.Reset()
            heff_rat.Divide(eff_mc,plt_DC_mc, 1,1)
            heff_rat.SetMaximum(1.5)
            heff_rat.SetMinimum(0.5)

            plt_CD_mc_rat = plt_CD_mc_topwwlvl.Clone('rat')
            plt_CD_mc_rat.SetTitle('Closure: mc_{top}^{est}/mc_{top}^{ww-lvl}')
            plt_CD_mc_rat.Reset()
            plt_CD_mc_rat.Divide(plt_CD_mc,plt_CD_mc_topwwlvl)
            plt_CD_mc_rat.Draw('e text45 colz')

            # main plotting style
            style  = {
                    'plotratio'  : False,
                    #'logx'       : True,
                    'morelogx'   : True,
                    'legalign'   : ('r','t'),
                    'rtitle'     : '%.2f fb^{-1}' % opt.lumi,
                    'legboxsize' : 30,
                    }
            # plotting style for efficiencies
            styleeff = style.copy()
            styleeff.update({
                    'drawopt'    : 'text',
                    'markersize' : 10,     # need it for text
                    'markers'    : [1]*10, # with a small marker
            })

            canvases = {}
            # Cs
            # first row
            canvases[0,0] = th2slices(plt_C_ds,  scalemax=1.5, ltitle='C region(ds)',   **style )
            canvases[1,0] = th2slices(plt_C_mc,  scalemax=1.5, ltitle='C region(mc)',   **style )

            canvases[2,0] = th2slices(plt_DC_mc, scalemax=1.5, ltitle='D/C(mc)', **styleeff )
            canvases[3,0] = th2slices(heff_rat,  yrange=(0,2), ltitle='#frac{B/A}{D/C}', **styleeff )

            # second row
            canvases[0,1] = th2slices(plt_D_ds,  scalemax=1.5, ltitle='D region(ds)',   **style )
            canvases[1,1] = th2slices(plt_D_mc,  scalemax=1.5, ltitle='D region(mc)',   **style )
            canvases[2,1] = th2slices(alpha_ds, yrange=(0,2), ltitle='#alpha (ds)',    **styleeff )
            canvases[3,1] = th2slices(alpha_mc, yrange=(0,2), ltitle='#alpha (mc)',    **styleeff )

            # third row
            canvases[0,2] = th2slices(plt_CD_ds, scalemax=1.5, ltitle='C-D region(ds)', **style )
            canvases[1,2] = th2slices(plt_CD_mc, scalemax=1.5, ltitle='C-D region(mc)', **style )
            canvases[2,2] = th2slices(plt_CD_mc_topwwlvl, scalemax=1.5, ltitle='C-D t#bar{t}+tW', **style )
            canvases[3,2] = th2slices(plt_CD_mc_rat, yrange=(0,2), ltitle='Closure: mc_{top}^{est}/mc_{top}^{ww-lvl}', **styleeff )

            for i,t in enumerate(tops):
                # tops addition to first row
                canvases[4+i,0] = th2slices(plt_DC_tops[t],scalemax=1.5, ltitle='C/D (%s)' % t, **styleeff)
                # tops addition to second row
                canvases[4+i,1] = th2slices(plots['D'][t], scalemax=1.5, ltitle='D (%s)' % t, **style)

            padstoprint = { 'yield_ds_D': (0,1), 'alpha_ds': (2,1), 'yield_ds_C-D':(0,2), 'yield_top_C-D':(2,2), 'effBA_vs_effDC':(3,0),'closure':(3,2) }
            for name,p in padstoprint.iteritems():
                c = canvases[p]
                c.Print('%s/%s.pdf' % (opt.prefix,name))
                c.Print('%s/%s.png' % (opt.prefix,name))

            call = Canvas()

            for p,c in canvases.iteritems():
                call[p] = EmbedPad(c)

            call.makecanvas()
            call.applystyle()

            call.SetCanvasSize(call.GetWw()/2,call.GetWh()/2)
            ROOT.gStyle.SetLineScalePS(1/max(call.gridsize()))
            #call.Print(opt.prefix+'/est_summary.eps')
            call.Print(opt.prefix+'/est_summary.png')
            call.Print(opt.prefix+'/est_summary.pdf')
            #os.system('epstopdf %s/est_summary.eps' % opt.prefix)
            #gStyle.SetLineScalePS(1.)


    print 'Ok Gringo!'

#_______________________________________________________________________
def efftop1j( analysers, regions={} ):
    '''
    Measure the top-tag efficiency for the 1j bin case (pt > 30 GeV)
    '''

    print '- Top efficiency 1j'
    if 'Top' in analysers:
        others = ['DYLL', 'DYTT', 'VV', 'Vg', 'VgS', 'WJet', 'WW', 'ggH', 'ggWW', 'qqH', 'wzttH']
        tops   = ['Top']
    else:
        if   opt.mcsub == 0:
            # classic set, with no subtraction
            tops   = ['ttbar','tW','DYLL', 'DYTT', 'VV', 'Vg', 'VgS', 'WJet', 'WW', 'ggH', 'ggWW', 'qqH', 'wzttH']
            others = []
        elif opt.mcsub == 1:
            # classic + mc subtraction
            tops   = ['ttbar','tW']
            others = ['DYLL', 'DYTT', 'VV', 'Vg', 'VgS', 'WJet', 'WW', 'ggH', 'ggWW', 'qqH', 'wzttH']
        elif opt.mcsub == 2:
            # mc & tW subtraction
            tops   = ['ttbar']
            others = ['DYLL', 'DYTT', 'VV', 'Vg', 'VgS', 'WJet', 'WW', 'ggH', 'ggWW', 'tW', 'qqH', 'wzttH']
        else:
            raise ValueError('What\'s this shit? '+opt.mcsub)
        #others = ['DYLL', 'DYTT', 'VV', 'Vg', 'VgS', 'WJet', 'WW', 'ggH', 'ggWW', 'tW', 'qqH', 'wzttH']
        #tops   = ['ttbar']

        ## switch this on to run in classic mode
        ## with no others
        #others = []
        ## and 2 samples
        ##others = ['DYLL', 'DYTT', 'VV', 'Vg', 'VgS', 'WJet', 'WW', 'ggH', 'ggWW', 'qqH', 'wzttH']
        #tops   = ['ttbar','tW']

    if  set(analysers.iterkeys()) != set(others+tops+['Data']):
        print 'Warning: not all the samples are used'

    # drop the analyzers which are not needed
    analysers = odict.OrderedDict([ (n,a) for n,a in analysers.iteritems() if n in (tops+others+['Data']) ])

    cuts = CutFlow()
    cuts['A'] = regions['A']
    cuts['B'] = regions['B']
    print '---1j cuts'.ljust(80,'-')
    print '\n'.join(['%-30s: %s'% i for i in cuts.iteritems()])

    print '---Updating buffers-'.ljust(80,'-')
    for n,a in analysers.iteritems():
        old = a.selectedentries()
        a.update(cuts)
        a.bufferentries()
        print '  ',n,':',old,'>>',a.selectedentries(),'...done'


    print '-'*80
    #y_data = analysers['Data'].yieldsflow()
    ys = {n:a.yieldsflow() for n,a in analysers.iteritems() }
    #es = {n:a.entriesflow() for n,a in analysers.iteritems() }
    y_data = ys['Data']
    y_A_o = sum(ys[o]['A'] for o in others) if others else Yield()
    y_B_o = sum(ys[o]['B'] for o in others) if others else Yield()

    print 'N_A_data  :', y_data['A']
    print 'N_B_data  :', y_data['B']

    print 'N_A_others:', y_A_o
    print 'N_B_others:', y_B_o

    eff_ds_val = ((y_data['B']-y_B_o)/(y_data['A']-y_A_o)).value   # other subtraction
    eff_ds_err = math.sqrt(eff_ds_val*(1-eff_ds_val)/y_data['A'].value) # other error subtraction?
    eff_ds = ValErr(eff_ds_val, eff_ds_err)

    y_A_t = sum(ys[t]['A'] for t in tops) if tops else Yield()
    y_B_t = sum(ys[t]['B'] for t in tops) if tops else Yield()

    eff_mc_val = (y_B_t/y_A_t).value   # other subtraction
    eff_mc_err = math.sqrt(eff_mc_val*(1-eff_mc_val)/y_A_t.value) # other error subtraction?
    eff_mc = ValErr(eff_mc_val, eff_mc_err)

    print '-'*80
    print 'eff_ds:',eff_ds
    print 'eff_mc:',eff_mc
    print '-'*80

    return eff_ds,eff_mc

#_______________________________________________________________________
def estimation1j( analysers, eff, regions={} ):

    others = ['DYLL', 'DYTT', 'VV', 'Vg', 'VgS', 'WJet', 'WW', 'ggH', 'ggWW', 'qqH', 'wzttH']
    if 'Top' in analysers:
        tops   = ['Top']
    else:
        tops   = ['tW','ttbar']

    assert( set(analysers.iterkeys()) == set(others+tops+['Data']) )

    eff_ds,eff_mc = eff

    print '- Top estimation 1j'

    # cuts for the mc closure
    cuts = CutFlow()
    cuts['C']   = regions['C']
    cuts['C-D'] = regions['C-D']
    cuts['mu']  = regions['mu']

    print '--1j cuts for closure--'.ljust(80,'-')
    print '\n'.join(['%-30s: %s'% i for i in cuts.iteritems()])
    atops = odict.OrderedDict( copy.deepcopy( [ (n,analysers[n]) for n in tops ] ) )
    for a in atops.itervalues():
        a.update(cuts)

    # specific cuts
    cuts = CutFlow()
    cuts['C'] = regions['C']
    cuts['D'] = regions['D']
    print '--1j cuts--'.ljust(80,'-')
    print '\n'.join(['%-30s: %s'% i for i in cuts.iteritems()])
    print '-'*80

    for n,a in analysers.iteritems():
        old = a.selectedentries()
        a.update(cuts)
        a.bufferentries()
        print '  ',n,':',old,'>>',a.selectedentries(),'...done'
    print '- buffering complete'

    nctrl_data  = analysers['Data'].yields()
    nctrl_top   = sum( a.yields() for n,a in analysers.iteritems() if n in tops)
#     nctrl_dy    = analysers['DYLL'].yields() + analysers['DYTT'].yields()
    nctrl_dy    =  3.73366*analysers['DYLL'].yields()
    nctrl_other = analysers['VV'].yields()   + analysers['Vg'].yields()  + analysers['VgS'].yields()
    nctrl_ww    = analysers['WW'].yields() + analysers['ggWW'].yields()
    nctrl_wjet  = analysers['WJet'].yields()

    nctrl_bkg = nctrl_dy+nctrl_ww+nctrl_wjet+nctrl_other

    nctrls = {n:analysers[n].yields() for n in others}
    print '-'*80
    print '\n'.join(['%s = %s' % n for n in nctrls.iteritems() ])
    print
    print 'nctrl_bkg = ',sum(nctrls.itervalues())
    print '-'*80


    #nctrl_mm    = analysers['Data'].yields('channel == 0')
    #nctrl_ee    = analysers['Data'].yields('channel == 1')
    #nctrl_em    = analysers['Data'].yields('channel == 2')
    #nctrl_me    = analysers['Data'].yields('channel == 3')

    #print '-'*80
    #print 'tagged events in data'
    #print 'nctrl_mm   ',nctrl_mm
    #print 'nctrl_ee   ',nctrl_ee
    #print 'nctrl_em   ',nctrl_em
    #print 'nctrl_me   ',nctrl_me
    print '-- D region --'.ljust(80,'-')
    print 'data ',nctrl_data
    print 'top  ',nctrl_top
    print 'dy   ',nctrl_dy
    print 'ww   ',nctrl_ww
    print 'wjet ',nctrl_wjet
    print 'other',nctrl_other
    print 'bkg  ',nctrl_bkg

    print 'nctrl_sub  ',nctrl_data-nctrl_bkg

    y_CD_tops = sum(a.views['C-D'].yields() for a in atops.itervalues() )
    y_mu_tops = sum(a.views['mu'].yields() for a in atops.itervalues() )

    #y_CD_data = extrtop(eff,nctrl_data,nctrl_bkg)
    y_CD_ds   = extrtop(eff_ds,nctrl_data,sum(nctrls.itervalues()))

    print '-- top events in C-D--'.ljust(80,'-')
    #print 'data :', y_CD_data
    print 'ds :', y_CD_ds
    print 'mc :', y_CD_tops
    print '-'*80

    eff2mu_val = y_mu_tops.value/y_CD_tops.value
    eff2mu_err = math.sqrt(eff2mu_val*(1-eff2mu_val)/y_CD_tops.value)
    eff2mu = ValErr(eff2mu_val,eff2mu_err)

    print 'eff2mu', eff2mu
    print '-- up to bveto_mu --'.ljust(80,'-')

    #print 'data :', eff2mu*y_CD_data
    print 'ds :', eff2mu*y_CD_ds
    print 'mc :', y_mu_tops


    print '-'*80
    if opt.prefix:
        with TStyleSentry( 'LatinosStyle' ) as sentry:
            plots_mll = {}
            plots_dphi = {}
            for n,a in analysers.iteritems():
                plots_mll[n]=a.plot('mll','mll',bins=(27,30,300))
                plots_dphi[n]=a.plot('dphill','dphill*180/pi',bins=(36,0,180))

            hwwlatino.printplots(plots_mll,opt.prefix+'/stack_D_mll',xlabel='mll', units='GeV', lumi=opt.lumi, exts=['pdf','png'])
            hwwlatino.printplots(plots_dphi,opt.prefix+'/stack_D_dphi',xlabel='dphi', units='deg', lumi=opt.lumi, exts=['pdf','png'])



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

    # take only the vbf specific part, from lepcnt1 onwards
    vbfflow = CutFlow( wwcuts.vbfcbfull( 160 ) )['lepcnt1':]

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

    print '-'*80
    print 'vbfbctrl cut',vbfbctrl.string()


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
        cjeta_btag  = a.views['bctrl'].plot(p+'_btag' ,cjetaexpr,cut=btag, bins=etabins)
        cjeta_bveto = a.views['bctrl'].plot(p+'_bveto',cjetaexpr,cut=bveto,bins=etabins)

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

    print '-'*80
    print 'vbfeff cut',vbfbeff.string()

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

    print '-'*80
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
    parser.add_option('-d', '--debug'   , dest='debug'   , help='Debug level'                     , default=0 , action='count' )
    parser.add_option('-o', '--out'     , dest='prefix'  , help='out dir'                         , default=None , )

    parser.add_option('-b', '--buf'     , dest='buffer'  , help='buffer common entries'           , default=False, action='store_true' )
    parser.add_option('-l', '--lumi'    , dest='lumi'    , help='Luminosity'                      , default=19.468 )
    parser.add_option('-c', '--chan'    , dest='chan'    , help='Channel [0j,1j,vbf,all]'         , default='all' )
    parser.add_option('-C', '--closure' , dest='closure' , help='Run the closure test'            , default=False, action='store_true' )
    parser.add_option('-T', '--onetop'  , dest='onetop'  , help='Top = ttbar+tW'                  , default=False, action='store_true' )
    parser.add_option(      '--mcsub'   , dest='mcsub'   , help='MC subtraction level (%default)' , default=2, type=int)
    parser.add_option(      '--logx2d'  , dest='logx2d'  , help='LogX in 2d'                      , default=False, action='store_true' )
    parser.add_option(      '--noest'   , dest='estimate', help='Dont\' perform estimate'         , default=True, action='store_false' )
    (opt, args) = parser.parse_args()

    if opt.prefix:
        # ensure the target dir
        hwwtools.ensuredir(opt.prefix)
        # add the index.php file
        os.system('addindex.sh '+opt.prefix)
        # attache the logfile
        logname = os.path.splitext(os.path.basename(sys.argv[0]))[0]
        logpath = '%s/%s.log' % (opt.prefix,logname)
        print 'Logging to', logpath
        tee = Tee(logpath)

    print '='*80
    print '  Command line:'
    print '  '+' '.join(sys.argv)
    print '='*80

    if not opt.debug:
        pass
    elif opt.debug >= 2:
        print 'Logging level set to DEBUG (%d)' % opt.debug
        logging.basicConfig(level=logging.DEBUG)
    elif opt.debug == 1:
        print 'Logging level set to INFO (%d)' % opt.debug
        logging.basicConfig(level=logging.INFO)

    # this is for ROOT
    sys.argv.append('-b')

    ROOT.gROOT.SetBatch()
    teststyle = ROOT.gROOT.GetStyle("Modern").Clone('2Dplots')
    teststyle.SetOptStat(0)
    teststyle.SetPadRightMargin(0.2)
    ROOT.gROOT.GetListOfStyles().Add(teststyle)

    shape_path = os.path.join(os.getenv('CMSSW_BASE'),'src/HWWAnalysis/ShapeAnalysis')
    print 'Shape directory is',shape_path
    ROOT.gInterpreter.ExecuteMacro(shape_path+'/macros/LatinoStyle2.C')

    print opt
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
        print '\n'.join(formatted_lines)
        print '-'*80
    finally:
        print 'GOODBYE!'

    try:
        __IPYTHON__
    except NameError:
        print 'Cleaning up'


