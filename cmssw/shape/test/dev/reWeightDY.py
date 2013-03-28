#!/usr/bin/env python


import hwwinfo
import hwwsamples
import hwwtools

import os.path

path_main='/shome/thea/HWW/work/shape2012/trees/'
path_bdt='/shome/thea/HWW/work/shape2012/trees/bdt_skim/mva_MH{mass}_{category}'
treename = 'latino'

# _____________________________________________________________________________
def _buildchain(treeName,files):
    import ROOT
    tree = ROOT.TChain(treeName)
    for path in files:
        if not os.path.exists(path):
            raise RuntimeError('File '+path+' doesn\'t exists')
        tree.Add(path) 

    return tree


# _____________________________________________________________________________
def _connect( samples, tree, path, mask=None, friends=[] ):
    print mask
    print tree, path
    print friends

    inputs = {}
    for process,filenames in samples.iteritems():
        if mask and process not in mask:
            continue
        chain = _buildchain(tree,[ os.path.join(path,f) for f in filenames])
        inputs[process] = chain
        chain.GetEntries()
        for ftree,fpath in friends:
            fchain =  _buildchain(ftree,[ os.path.join(fpath,f) for f in filenames])
            if chain.GetEntriesFast() != fchain.GetEntries():
                raise RuntimeError('Mismatching number of entries: '
                                   +tree.GetName()+'('+str(tree.GetEntriesFast())+'), '
                                   +ftree.GetName()+'('+str(ftree.GetEntriesFast())+')')
            chain.AddFriend(fchain)


        print process, '>>', chain.GetEntriesFast()
    
    return inputs
    
# _____________________________________________________________________________
def _disconnect(inputs):
    import ROOT
    for n in inputs.keys():
        friends = inputs[n].GetListOfFriends()
        if friends.__nonzero__():
            for fe in friends:
                friend = fe.GetTree()
                inputs[n].RemoveFriend(friend)
                ROOT.SetOwnership(friend,True)
                del friend
        del inputs[n]


def getminmax(tree,var,binsize):
    import math
    xmin,xmax = dyll.GetMinimum(var),dyll.GetMaximum(var)
    
    xmin,xmax = math.floor(xmin/binsize)*binsize,math.ceil(xmax/binsize)*binsize

    return xmin,xmax

if __name__ == '__main__':
    print 'start!'

    import optparse
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)

    parser.add_option('-o',dest='out', help='output path', default='')
    (opt, args) = parser.parse_args()

    if opt.out and not os.path.exists(opt.out):
        os.system('mkdir -p '+opt.out)

    import ROOT
    mypath = os.path.dirname(os.path.abspath(__file__))
    ROOT.gInterpreter.ExecuteMacro(mypath+'/LatinoStyle2.C')

    inputs = _connect(hwwsamples.backgrounds, 'latino',path_main, mask=['DYLL'] ) 

    print 'inputs =',inputs

    dyll = inputs['DYLL']

    ROOT.TH1.SetDefaultSumw2()

    import numpy
    edges = range(0,60,5)+range(60,100,10)+range(100,200,20)+range(200,320,40)+range(320,800,80)
    bins = numpy.array(edges,dtype='d')
    h = ROOT.TH1F('ptll','ptll;ptll',len(bins)-1,bins)
    h.SetMarkerStyle(24)

    h.SetLineWidth(2)
    h.SetBit(ROOT.TH1.kNoStats)

    hlo = h.Clone('ptll_lo')
    hhi = h.Clone('ptll_hi')
    
    hlo.SetLineColor(ROOT.kRed)
    hlo.SetMarkerColor(ROOT.kRed)
    hhi.SetLineColor(ROOT.kBlue)
    hhi.SetMarkerColor(ROOT.kBlue)


    formats=['pdf','png']

#     c.SetLogy()
    basic = 'trigger==1 && nextra==0 && !zveto && njet==0'
    cutlo = basic+'&& ( 20 < mpmet  && mpmet < 45 )'
    cuthi = basic+' && mpmet > 45 '

    print 'Using basic selection:',basic

    hlo.SetTitle(cutlo)
    hhi.SetTitle(cuthi)
    print 'lomet:',dyll.Draw('ptll >> ptll_lo',cutlo,'goff')

    print 'himet:',dyll.Draw('ptll >> ptll_hi',cuthi,'goff')


    hs = ROOT.THStack('ptll_hilomet','')
    hs.Add(hlo)
    hs.Add(hhi)

    c = ROOT.TCanvas()
    hs.Draw('nostack')
    c.BuildLegend()

    c.SaveAs(os.path.join(opt.out,'check_ptll.pdf') )
    c.SaveAs(os.path.join(opt.out,'check_ptll.png') )

    ratio = hhi.Clone('ptll_hilomet_ratio')
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetMarkerColor(ROOT.kBlack)
    ratio.Divide(hlo)
    ratio.SetTitle('')

#     fitfun = ROOT.TF1('ptll_hilomet_fun','pol1',ratio.GetXaxis().GetXmin(), ratio.GetXaxis().GetXmax())
    fitfun = ROOT.TF1('ptll_hilomet_fun','pol1',45., ratio.GetXaxis().GetXmax())
    fitfun.SetLineWidth(1)
    fitfun.SetLineColor(ROOT.kOrange+7)
    c.Clear()
    ratio.Draw('hist')
    # fit from 45 up only
    r = ratio.Fit(fitfun,'','',45,ratio.GetXaxis().GetXmax())
    c.SaveAs(os.path.join(opt.out,'ratio.pdf') )
    c.SaveAs(os.path.join(opt.out,'ratio.png') )

    ratio.GetXaxis().SetRangeUser(0,200)
    c.SetLogy()
    c.SaveAs(os.path.join(opt.out,'ratio_zoom.pdf') )
    c.SaveAs(os.path.join(opt.out,'ratio_zoom.png') )


    of = ROOT.TFile.Open('hilomet_ratio.root','recreate')
    ratio.Write()
    fitfun.Write()
    of.Write()
    of.Close()
    del of

    _disconnect(inputs)



