#!/usr/bin/env python

import HWWAnalysis.Misc.odict as odict
#from ginger.analysis import TreeAnalyser,Cut,CutFlow,Latino,Sample
from ginger.analysis import CutFlow
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
import hwwtools
import ctypes

from HWWAnalysis.Misc.ROOTAndUtils import TStyleSentry,PadPrinter,Tee

from ginger.plotter import H1RatioPlotter
from ginger.painter import Canvas,Pad,Legend,EmbedPad

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

#_______________________________________________________________________
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

#_______________________________________________________________________
def th2yield( h, binx=None, biny=None ):
    err = ctypes.c_double(0.)
    binx1,binx2 = binx if binx else (1,h.GetNbinsX())
    biny1,biny2 = biny if biny else (1,h.GetNbinsY())
    #print binx1,binx2,biny1,biny2
    return Yield(h.IntegralAndError(binx1,binx2,biny1,biny2,err),err.value)

#_______________________________________________________________________
def extrtop(eff,yldtag,yldbkg):
    eff_val, eff_err = eff

    yldsub = yldtag-yldbkg
    val = yldsub.value * (1-eff_val)/eff_val;
    err = math.sqrt(
        ( yldsub.value * eff_err / eff_val**2)**2 +
        ( (1-eff_val)/eff_val * yldbkg.error)**2
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
def th2slices( th2, **opts ):

    slices = []
    yax = th2.GetYaxis()
    for i in xrange(1,yax.GetNbins()+1):
        h = th2.ProjectionX('%s_%d' % (th2.GetName(),i),i,i,'e')
        h.SetTitle( '%.1f < #eta < %.1f' % (yax.GetBinLowEdge(i),yax.GetBinUpEdge(i)) )
        slices.append(h)

    hratio = H1RatioPlotter(**opts)
    hratio.set(*slices)
    c = hratio.plot()
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

    print '------------------------------'
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


    if opt.chan in ['0j','all']:
        eff0j = efftop0j( copy.deepcopy(analysers) )
    #     eff0j = (0.489415, 0.048752)
        estimation0j( copy.deepcopy(analysers), eff0j )

    if opt.chan in ['1j','all']:
        eff1j = efftop1j( copy.deepcopy(analysers) )
    #     eff1j = (0.647059, 0.004662)
        estimation1j( copy.deepcopy(analysers), eff1j )

    if opt.chan in ['1j2g_lead','all']:
        eff1j = efftop1j2g_lead( copy.deepcopy(analysers), opt )
    #     eff1j = (0.647059, 0.004662)
        estimation1j2g( copy.deepcopy(analysers), eff1j, opt )

    if opt.chan in ['1j2g_sublead','all']:
        eff1j = efftop1j2g_sublead( copy.deepcopy(analysers), opt )
    #     eff1j = (0.647059, 0.004662)
        estimation1j2g( copy.deepcopy(analysers), eff1j, opt )

    if opt.chan in ['vbf','all']:
        eff2jdata, eff2jmc = efftopvbf( copy.deepcopy(analysers) )
    #     eff2j.Print()
    #     eff2j = None
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
def efftop1j2g_lead( analysers, opt ):
    '''
    Measure the top-tag efficiency for the 1j bin case (pt > 30 GeV)
    '''

    print '- Top efficiency 1j (on leading) in pt/eta bins'

    print '--Updating buffers-------'
    for n,a in analysers.iteritems():
        old = a.selectedentries()
        #a.append('A','njet==2 && bveto_munj30 && softtche<=2.1 && jettche2>2.1 && ( njet <= 1 || (bveto_ip && njet >= 2) ) && jettche3<=2.1 && jettche4<=2.1')
        #a.append('A','njet==2 && bveto_munj30 && softtche<=2.1 && jettche2>2.1')
        #a.append('A','njet==2 && bveto_munj30 && softtche<=2.1 && jettche2>2.1 && jetpt3<10')
        #a.append('A','jetpt2>20 && bveto_munj30 && jettche2>2.1 && jetpt3<10')
        #a.append('A','njet==2 && bveto_mu  && softtche<=2.1 && jettche2>2.1')
        #a.append('A','jetpt2>20 && bveto_mu && jettche2>2.1 && jetpt3<10')
        #a.append('A','njet==2 && bveto_munj30 && jettche2>2.1 && jetpt3<10')
        a.append('A','njet==2 && bveto_munj30 && softtche<=2.1 && jettche2>2.1') # classic
        a.append('B','jettche1>2.1')
        a.bufferentries()
        print '  ',n,':',old,'>>',a.selectedentries(),'...done'

    print '-'*80

    ptbins  = range(30, 70, 10) + range( 70, 150, 20) + range(150,250, 50) + [250,500,1000]
    etabins = [0.,0.75,1.5,2.8,5]
    etabins = [0.,2.8,5]
    ptbins  = range(30,200,10)
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

    if 'Top' in analysers:
        others = ['DYLL', 'DYTT', 'VV', 'Vg', 'VgS', 'WJet', 'WW', 'ggH', 'ggWW', 'qqH', 'wzttH']
        tops   = ['Top']
    else:
        others = ['DYLL', 'DYTT', 'VV', 'Vg', 'VgS', 'WJet', 'WW', 'ggH', 'ggWW', 'tW', 'qqH', 'wzttH']
        tops   = ['ttbar']
        others = ['DYLL', 'DYTT', 'VV', 'Vg', 'VgS', 'WJet', 'WW', 'ggH', 'ggWW', 'qqH', 'wzttH']
        tops   = ['ttbar','tW']

    # assure a proper subtraction
    #assert( set(analysers.iterkeys()) == set(others+tops+['Data']) )

    # calculate the efficiency
    plt_A_ds = plots['A']['Data'].Clone('plt_A_ds')
    plt_A_ds.Add( thsum(plots['A'],'ctrl_others',others) ,-1)
    plt_A_ds.SetTitle('Data etapt [A]')

    plt_B_ds = plots['B']['Data'].Clone('plt_B_ds')
    plt_B_ds.Add( thsum(plots['B'],'tag_others',others) ,-1)
    plt_B_ds.SetTitle('Data etapt [B]')

    plt_A_mc = thsum(plots['A'], 'plt_A_mc', tops)
    plt_B_mc = thsum(plots['B'] , 'plt_B_mc' , tops)

    eff_ds = plt_B_ds.Clone('eff_1j2g')
    eff_ds.Reset()

    eff_ds.Divide(plt_B_ds,plt_A_ds,1,1,'b')
    eff_ds.SetTitle('b-tag efficiency [data-mc_{other}]')

    eff_mc = plt_B_mc.Clone('eff_1j2g')
    eff_mc.Reset()

    eff_mc.Divide(plt_B_mc,plt_A_mc,1,1,'b')
    eff_mc.SetTitle('b-tag efficiency [mc_{top}]')



    test_ds = plt_A_ds.Clone('test_ds')
    test_ds.Reset()
    test_ds.Multiply(eff_ds,plt_A_ds)
    test_ds.Scale(1/ plt_A_ds.Integral())

    test_mc = plt_A_mc.Clone('test_mc')
    test_mc.Reset()
    test_mc.Multiply(eff_mc,plt_A_mc)
    test_mc.Scale(1/ plt_A_mc.Integral())

    print 'eff_ds (diff->inc):',test_ds.Integral()
    print 'eff_mc (diff->inc):',test_mc.Integral()

    nA_ds = th2yield(plt_A_ds)
    nB_ds = th2yield(plt_B_ds)
    nA_mc = th2yield(plt_A_mc)
    nB_mc = th2yield(plt_B_mc)

    print '-'*80
    print '---inclusive efficiency'

    print 'N_A_data  :', th2yield( plots['A']['Data'] )
    print 'N_B_data  :', th2yield( plots['B']['Data'] )
    print 'N_A_others:', th2yield(thsum(plots['A'],'ctrl_others',others) )
    print 'N_B_others:', th2yield(thsum(plots['B'],'ctrl_others',others) )
    print '---'


    print 'N^ds_A = ',nA_ds
    print 'N^ds_B = ',nB_ds
    print 'N^mc_A = ',nA_mc
    print 'N^mc_B = ',nB_mc

    eff_inc_val = (nB_ds/nA_ds).value
    eff_inc_err = math.sqrt(eff_inc_val*(1-eff_inc_val)/nA_ds.value)
    eff_inc= Yield(eff_inc_val, eff_inc_err)

    print 'eff (inc) =',eff_inc



    print '-'*80

    if opt.prefix:
        hwwtools.ensuredir(opt.prefix)

        thsameminmax(plt_A_ds, plt_A_mc)
        thsameminmax(plt_B_ds, plt_B_mc)
        thsameminmax(eff_ds, eff_mc, max=1)

        #with TStyleSentry( '2Dplots' ) as sentry:
            #ROOT.gStyle.SetPaintTextFormat('2.2f')
            #if opt.logx2d:
                #ROOT.gStyle.SetOptLogx()
            #else:
                #ROOT.gStyle.SetOptLogx(False)
    ##         ROOT.gStyle.SetMarkerColor(ROOT.kWhite)
    ##         ROOT.gStyle.SetMarkerSize(1.)

    ##         map(ROOT.TH2.UseCurrentStyle,[plt_A_ds,plt_B_ds,eff_ds,plt_A_mc,plt_B_mc,eff_mc])

            #c = ROOT.TCanvas('stoca','stoca',500,750)
            #c.Divide(2,3)

            #c.cd(1)
            #plt_A_ds.Draw('text45 colz')
            #c.cd(2)
            #plt_A_mc.Draw('text45 colz')

            #c.cd(3)
            #plt_B_ds.Draw('text45 colz')
            #c.cd(4)
            #plt_B_mc.Draw('text45 colz')

            #c.cd(5)
            #eff_ds.Draw('e text45 colz')
            #c.cd(6)
            #eff_mc.Draw('e text45 colz')

            #c.Print(opt.prefix+'/stocazzo.png')
            #c.Print(opt.prefix+'/stocazzo.pdf')

            #padstoprint = { 'yield_ds_A':1, 'yield_ds_B':3, 'eff_ds':5, 'eff_mc':6, }

            #printer = PadPrinter(opt.prefix)
            #printer.savefromcanvas(c,**padstoprint)



        canvases = {}
        options  = {
                'plotratio'  : False,
                #'logx'       : True,
                'morelogx'   : True,
                'legalign'   : ('r','t'),
                'rtitle'     : '19 fb^{-1}',
                #'legtextsize': 20,
                'legboxsize' : 30,
                }
        print 'here'
        canvases[0,0] = th2slices(plt_A_ds, scalemax=1.5, ltitle='A region (ds)', **options )
        canvases[1,0] = th2slices(plt_A_mc, scalemax=1.5, ltitle='A region (mc)', **options )
        canvases[0,1] = th2slices(plt_B_ds, scalemax=1.5, ltitle='B region (ds)', **options )
        canvases[1,1] = th2slices(plt_B_mc, scalemax=1.5, ltitle='B region (mc)', **options )
        canvases[0,2] = th2slices(eff_ds,   scalemax=1.8, ltitle='efficiency (ds)', **options )
        canvases[1,2] = th2slices(eff_mc,   scalemax=1.8, ltitle='efficiency (mc)', **options )
        print 'there'

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
        call.Print(opt.prefix+'/eff_summary.pdf')
        call.Print(opt.prefix+'/eff_summary.png')

    return eff_ds, eff_mc

#_______________________________________________________________________
def estimation1j2g( analysers, eff, opt ):
    # sanity check
    # note: the ctrl region contains both ttbar and tW, therefore tW goes in tops and not in ttbar
    others = ['DYLL', 'DYTT', 'VV', 'Vg', 'VgS', 'WJet', 'WW', 'ggH', 'ggWW', 'qqH', 'wzttH']
    tops   = ['tW', 'ttbar']
    others = ['DYLL', 'DYTT', 'VV', 'Vg', 'VgS', 'WJet', 'WW', 'ggH', 'ggWW', 'qqH', 'wzttH','tW']
    tops   = ['ttbar']
    # insure a proper subtraction
    assert( set(analysers.iterkeys()) == set(others+tops+['Data']) )

    eff_ds, eff_mc = eff

    for n,a in analysers.iteritems():
        old = a.selectedentries()
        a.append('pre-tag','njet==1 && bveto_munj30 && softtche<=2.1') # classic
        #a.append('pre-tag','njet==1 && bveto_mu && softtche<=2.1')
        a.append('ctrltop','jettche1>2.1')
        a.bufferentries()
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
#         plots['bveto'][n] = a.views['pre-tag'].plot(n+'_etapt','fabs(jeteta1):jetpt1',cut='bveto_ip && bveto_mu && nbjettche==0',bins=(ptbins,etabins))
        plots['bveto'][n] = a.views['pre-tag'].plot(n+'_etapt','fabs(jeteta1):jetpt1',cut='nbjettche==0',bins=(ptbins,etabins))
        sys.stdout.flush()

    plots.lock()
    print ']'

    print plots.keys()

    pretag_ds = plots['pre-tag']['Data'].Clone('pretag_ds')
    pretag_ds.Add( thsum(plots['pre-tag'],'pretag_others',others), -1)
    pretag_ds.SetTitle('WW-lvl - b_{veto} [data-mc_{other}]')

    pretag_mc = thsum(plots['pre-tag'],'pretag_mc',tops)
    pretag_mc.SetTitle('WW-lvl - b_{veto} [mc_{top}]')

    btag_ds = plots['ctrltop']['Data'].Clone('btag_ds')
    btag_ds.Add( thsum(plots['ctrltop'],'btag_others',others), -1)
    btag_ds.SetTitle('top control region [btagged,data-mc_{other}]')

    btag_mc = thsum(plots['ctrltop'],'btag_mc',tops)
    btag_mc.SetTitle('top control region [btagged,mc_{top}]')

    # same for the btag-ctrl region

    # clean negative bins
    for i in xrange(btag_ds.GetSize()):
        bc = btag_ds.GetBinContent(i)
        if bc > 0: continue

        btag_ds.SetAt(0,i)
        btag_ds.SetBinError(i,0)

    alpha_ds = eff2alpha( eff_ds, 'alpha_ds' )
    alpha_mc = eff2alpha( eff_mc, 'alpha_mc' )


    #bveto_ds = btag_ds.Clone('bveto_ds')
    #bveto_ds.Reset()
    #bveto_ds.SetTitle('WW-lvl t#bar{t}+tW estimate [data-mc_{other}]')
    #bveto_ds.Multiply(alpha_ds,btag_ds)


    #bveto_mc = btag_mc.Clone('bveto_mc')
    #bveto_mc.Reset()
    #bveto_mc.SetTitle('WW-lvl t#bar{t}+tW estimate [mc_{top}]')
    #bveto_mc.Multiply(alpha_mc,btag_mc)

    bveto_ds = applyalpha( alpha_ds, btag_ds)
    bveto_ds.SetTitle('WW-lvl t#bar{t}+tW estimate [data-mc_{other}]')

    bveto_mc = applyalpha( alpha_mc, btag_mc)
    bveto_mc.SetTitle('WW-lvl t#bar{t}+tW estimate [mc_{top}]')

    bveto_mc_topwwlvl = thsum(plots['bveto'],'bveto_mc',tops)
    bveto_mc_topwwlvl.SetTitle('top mc, WW-level yield [bveto,mc_{top}]')


    print '-'*80

    err = ctypes.c_double(0.)

    print '-'*80
    beta= th2yield(pretag_ds, biny=(neta,neta))/th2yield(pretag_ds, biny=(1,neta-1))

    print 'beta',beta

    y_btag_ds  = th2yield(btag_ds)
    y_bveto_ds = th2yield(bveto_ds)

    print '-'*80
    print 'ww btag  yield [rw,0.0-5.0]', th2yield(plots['ctrltop']['Data'])
    print 'ww btag  yield [ot,0.0-5.0]', th2yield(thsum(plots['ctrltop'],'pretag_others',others))
    print 'ww btag  yield [ds,0.0-5.0]', y_btag_ds
    for s in tops:
        print s,'yield [xx,0.0,5.0]:', th2yield(plots['ctrltop'][s])

    print '='*4
    print 'ww btag  yield [mc,0.0-5.0]', th2yield(btag_mc)
    print 'ww bveto yield [ds,0.0-5.0]', y_bveto_ds

    print '-'*80

    print 'ww-pre bveto [est,ds,0-5]', (y_btag_ds+y_bveto_ds)*(1+beta)
    print 'ww-pre bveto [puremc,0-5]', th2yield(pretag_mc)

    print '-'*80
    print 'ww bveto [est,ds,0-5]', beta*y_btag_ds+(1+beta)*y_bveto_ds
    print 'ww bveto [puremc,0-5]', th2yield(bveto_mc_topwwlvl)
    print '-'*80


    if opt.prefix:
        hwwtools.ensuredir(opt.prefix)

        thsameminmax( pretag_ds, pretag_mc )
        thsameminmax( btag_ds, btag_mc )
        thsameminmax( alpha_ds, alpha_mc, max=2. )
        thsameminmax( bveto_ds,bveto_mc,bveto_mc_topwwlvl )

        #with TStyleSentry( '2Dplots' ) as sentry:
            #ROOT.gStyle.SetPaintTextFormat('2.2f')
            #if opt.logx2d:
                #ROOT.gStyle.SetOptLogx()
            #else:
                #ROOT.gStyle.SetOptLogx(False)

            #padstoprint = {}

            #c = ROOT.TCanvas('stami','stami',1000,750)

            #c.SetWindowSize(2*c.GetWw(),2*c.GetWh())
            #c.Modified()
            #c.Update()

            #c.Divide(4,3)

            ## pre-tag (WW-pre-bveto)
            #c.cd(1)
            #pretag_ds.Draw('text45 colz')
            #c.cd(2)
            #pretag_mc.Draw('text45 colz')

            #c.cd(3)
    ##         ROOT.gPad.SetLogz()
            #heff_mc = btag_mc.Clone('heff_mc')
            #heff_mc.SetTitle('btag/preveto on 1j')
            #heff_mc.Reset()
            #heff_mc.Divide(btag_mc,pretag_mc,1,1,'b')
            #heff_mc.SetMaximum(1)
            #heff_mc.Draw('e text45 colz')

            #c.cd(4)
            #heff_rat = heff_mc.Clone('heff_rat')
            #heff_rat.SetTitle('efficiency / (btag/all)')
            #heff_rat.Reset()
            #heff_rat.Divide(eff_mc,heff_mc, 1,1)
            #heff_rat.SetMaximum(1.5)
            #heff_rat.SetMinimum(0.5)
            #heff_rat.Draw('e text45 colz')


    ##         ROOT.gPad.SetLogz()
            #c.cd(5)
            #btag_ds.Draw('e text45 colz')
            #c.cd(6)
            #btag_mc.Draw('e text45 colz')


            ## alphas
            #c.cd(7)
            #alpha_ds.Draw('e text45 colz')
            #c.cd(8)
            #alpha_mc.Draw('e text45 colz')

            ## bvetoed
            #c.cd(9)
            #bveto_ds.Draw('e text45 colz')
            #c.cd(10)
            #bveto_mc.Draw('e text45 colz')
            #c.cd(11)
            #bveto_mc_topwwlvl.Draw('e text45 colz')


            #padstoprint['closure_mc'] = c.cd(12)
            #bveto_mc_rat = bveto_mc_topwwlvl.Clone('rat')
            #bveto_mc_rat.SetTitle('Closure: mc_{top}^{est}/mc_{top}^{ww-lvl}')
            #bveto_mc_rat.Reset()
            #bveto_mc_rat.Divide(bveto_mc,bveto_mc_topwwlvl)
            #bveto_mc_rat.Draw('e text45 colz')
            #bveto_mc_rat.SetMaximum(2)

            #c.Print(opt.prefix+'/stameng.png')
            #c.Print(opt.prefix+'/stameng.pdf')

            #padstoprint = { 'yield_ds_D':5, 'alpha_ds':7, 'yield_ds_C-D':9, 'yield_top_C-D':11, 'closure':12 }

            #printer = PadPrinter(opt.prefix)
            #printer.savefromcanvas(c,**padstoprint)

        heff_mc = btag_mc.Clone('heff_mc')
        heff_mc.SetTitle('btag/preveto on 1j')
        heff_mc.Reset()
        heff_mc.Divide(btag_mc,pretag_mc,1,1,'b')
        heff_mc.SetMaximum(1)

        heff_rat = heff_mc.Clone('heff_rat')
        heff_rat.SetTitle('efficiency / (btag/all)')
        heff_rat.Reset()
        heff_rat.Divide(eff_mc,heff_mc, 1,1)
        heff_rat.SetMaximum(1.5)
        heff_rat.SetMinimum(0.5)

        bveto_mc_rat = bveto_mc_topwwlvl.Clone('rat')
        bveto_mc_rat.SetTitle('Closure: mc_{top}^{est}/mc_{top}^{ww-lvl}')
        bveto_mc_rat.Reset()
        bveto_mc_rat.Divide(bveto_mc,bveto_mc_topwwlvl)
        bveto_mc_rat.Draw('e text45 colz')

        canvases = {}
        options  = {
                'plotratio'  : False,
                #'logx'       : True,
                'morelogx'   : True,
                'legalign'   : ('r','t'),
                'rtitle'     : '19 fb^{-1}',
                'legboxsize' : 30,
                }
        print 'here'
        # Cs
        canvases[0,0] = th2slices(pretag_ds,scalemax=1.5, ltitle='C region(ds)',   **options )
        canvases[1,0] = th2slices(pretag_mc,scalemax=1.5, ltitle='C region(mc)',   **options )

        canvases[2,0] = th2slices(heff_mc,  scalemax=1.5, ltitle='btag/preveto on 1j (C/D)', **options )
        canvases[3,0] = th2slices(heff_rat, scalemax=1.5, ltitle='efficiency / (C/D)', **options )

        canvases[0,1] = th2slices(btag_ds,  scalemax=1.5, ltitle='D region(ds)',   **options )
        canvases[1,1] = th2slices(btag_mc,  scalemax=1.5, ltitle='D region(mc)',   **options )
        canvases[2,1] = th2slices(alpha_ds, yrange=(0,2), ltitle='#alpha (ds)',    **options )
        canvases[3,1] = th2slices(alpha_mc, yrange=(0,2), ltitle='#alpha (mc)',    **options )

        canvases[0,2] = th2slices(bveto_ds, scalemax=1.5, ltitle='C-D region(ds)', **options )
        canvases[1,2] = th2slices(bveto_mc, scalemax=1.5, ltitle='C-D region(mc)', **options )
        canvases[2,2] = th2slices(bveto_mc_topwwlvl, scalemax=1.5, ltitle='#alpha (ds)', **options )
        canvases[3,2] = th2slices(bveto_mc_rat, yrange=(0,2), ltitle='Closure: mc_{top}^{est}/mc_{top}^{ww-lvl}', **options )
        print 'there'

        padstoprint = { 'yield_ds_D': (0,1), 'alpha_ds': (2,1), 'yield_ds_C-D':(0,2), 'yield_top_C-D':(2,2), 'closure':(3,2) }
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
        call.Print(opt.prefix+'/est_summary.pdf')
        call.Print(opt.prefix+'/est_summary.png')
    print 'Ok Gringo!'

#_______________________________________________________________________
def efftop1j( analysers ):
    '''
    Measure the top-tag efficiency for the 1j bin case (pt > 30 GeV)
    '''

    print '- Top efficiency 1j'
    if 'Top' in analysers:
        others = ['DYLL', 'DYTT', 'VV', 'Vg', 'VgS', 'WJet', 'WW', 'ggH', 'ggWW', 'qqH', 'wzttH']
        tops   = ['Top']
    else:
        others = ['DYLL', 'DYTT', 'VV', 'Vg', 'VgS', 'WJet', 'WW', 'ggH', 'ggWW', 'tW', 'qqH', 'wzttH']
        tops   = ['ttbar']
    assert( set(analysers.iterkeys()) == set(others+tops+['Data']) )



    print '--Updating buffers-------'
    for n,a in analysers.iteritems():
        old = a.selectedentries()
        a.append('ctrl','njet==2 && bveto_munj30 && softtche<=2.1 && jettche2>2.1')  # classic
        #a.append('ctrl','njet==2 && bveto_munj30 && jettche2>2.1')
        a.append('tag','jettche1>2.1')
        a.bufferentries()
        print '  ',n,':',old,'>>',a.selectedentries(),'...done'


    print '-'*80
    #y_data = analysers['Data'].yieldsflow()
    ys = {n:a.yieldsflow() for n,a in analysers.iteritems() }
    y_data = ys['Data']
    y_A = sum(ys[o]['ctrl'] for o in others)
    y_B = sum(ys[o]['tag'] for o in others)

    print 'N_A_data  :', y_data['ctrl']
    print 'N_B_data  :', y_data['tag']
    print 'N_A_others:', y_A
    print 'N_B_others:', y_B
    print 'eff_ds    :', (y_data['tag']-y_B)/(y_data['ctrl']-y_A)

    eff_top_1j = ((y_data['tag']-y_B)/(y_data['ctrl']-y_A)).value   # other subtraction
    #eff_top_1j = y_data['tag'].value/y_data['ctrl'].value          # classic
    eff_top_1j_err = math.sqrt(eff_top_1j*(1-eff_top_1j)/y_data['ctrl'].value)

    print '-'*80
    print 'eff_top %f +/- %f ' % ( eff_top_1j,eff_top_1j_err )
    print '-'*80

    return (eff_top_1j, eff_top_1j_err)

#_______________________________________________________________________
def estimation1j( analysers, eff ):

    others = ['DYLL', 'DYTT', 'VV', 'Vg', 'VgS', 'WJet', 'WW', 'ggH', 'ggWW', 'qqH', 'wzttH']
    if 'Top' in analysers:
        tops   = ['Top']
    else:
        tops   = ['tW','ttbar']

    assert( set(analysers.iterkeys()) == set(others+tops+['Data']) )

    print '- Top estimation 1j'

#     atop = analysers['Top'].clone()
#     atop.append('bveto','njet == 1 && bveto_munj30 && bveto_ip && nbjettche==0')
    atops = odict.OrderedDict( copy.deepcopy( [ (n,analysers[n]) for n in tops ] ) )
    for a in atops.itervalues():
        a.append('bveto','njet == 1 && bveto_munj30 && bveto_ip && nbjettche==0')

    for n,a in analysers.iteritems():
        old = a.selectedentries()
        a.append('ctrltop','njet==1 && bveto_munj30 && jettche1>2.1 && softtche<=2.1')
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
    print 'top events from data :', extrtop(eff,nctrl_data,nctrl_bkg)
    print 'top events from dataX:', extrtop(eff,nctrl_data,sum(nctrls.itervalues()))
#     print 'top events from mc:  ', atop.yields()
    print 'top events from mcX  : ', sum(a.yields() for a in atops.itervalues() )

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
    parser.add_option('-d', '--debug'  , dest='debug'   , help='Debug level'            , default=0 , action='count' )
    parser.add_option('-o', '--out'    , dest='prefix'  , help='out dir'                , default=None , )

    parser.add_option('-b', '--buf'    , dest='buffer'  , help='buffer common entries'  , default=False, action='store_true' )
    parser.add_option('-l', '--lumi'   , dest='lumi'    , help='Luminosity'             , default=19.468 )
    parser.add_option('-c', '--chan'   , dest='chan'    , help='Channel [0j,1j,vbf,all]', default='all' )
    parser.add_option('-C', '--closure', dest='closure' , help='Run the closure test'   , default=False, action='store_true' )
    parser.add_option('-T', '--onetop' , dest='onetop'  , help='Top = ttbar+tW'         , default=False, action='store_true' )
    parser.add_option(      '--logx2d' , dest='logx2d'  , help='LogX in 2d'             , default=False, action='store_true' )
    (opt, args) = parser.parse_args()

    if opt.prefix:
        hwwtools.ensuredir(opt.prefix)
        import os
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

    teststyle = ROOT.gROOT.GetStyle("Modern").Clone('2Dplots')
    teststyle.SetOptStat(0)
    teststyle.SetPadRightMargin(0.2)
    ROOT.gROOT.GetListOfStyles().Add(teststyle)

    import os.path
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


