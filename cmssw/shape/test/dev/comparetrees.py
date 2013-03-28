#!/usr/bin/env python


import hwwinfo
import hwwsamples
import hwwtools

import os.path
import HWWAnalysis.Misc.ROOTAndUtils as utils
import HWWAnalysis.Misc.odict as odict
import ROOTplotter


import ROOTtree


        
# _____________________________________________________________________________
if __name__ == '__main__':
    print 'start!'


    import optparse
    import ROOT

    ROOT.TH1.SetDefaultSumw2()
    mypath = os.path.dirname(os.path.abspath(__file__))
    ROOT.gInterpreter.ExecuteMacro(mypath+'/LatinoStyle2.C')

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)

    parser.add_option('-o',dest='out', help='output path', default='')
    parser.add_option('-c',dest='cfg', help='cfg path', default='cfg.py')
    (opt, args) = parser.parse_args()

    filename=opt.cfg
    if os.path.exists(filename):
        handle = open(filename,'r')
        pars = {}
        exec(handle,pars)
        handle.close()
    
    print pars.keys()

    if opt.out and not os.path.exists(opt.out):
        os.system('mkdir -p '+opt.out)

    sA = pars['sampleA']
    sB = pars['sampleB']

    
    tA = ROOTtree.TreeWorker(sA.tree,sA.files)
    tA.setselection( sA.selection )

    tB = ROOTtree.TreeWorker(sB.tree,sB.files)
    tB.setselection( sB.selection )
    tB.setweight(sB.weight)


#     if True:
#         print '-'*80
#         import sys

#         testH = ROOT.TH1F('xx','yy',18,0,360)
#         testH.SetDirectory(0x0)

#         tA.fill(testH,'dphill*180/pi')



#         sys.exit(0)

    if False:
        print '-'*80
        y =  tA.yieldsflow( hwwinfo.wwcutsB.wwhi )
        for name,value in y.iteritems():
            print name, '%.2f' % value
        print


        print '-'*80
        y =  tB.yieldsflow( hwwinfo.wwcutsB.wwhi )
        for name,value in y.iteritems():
            print name, '%.2f' % value
        print

#     cut = ' && '.join( '(%s)' % s for s in hwwinfo.wwcutsB.wwhi.values()[:-2])
#     cut = '!zveto && njet==1 && channel==0' 
    cut = pars['cut']

    print 'Active cut',cut
    canvas = ROOT.TCanvas('xxx','xxx') 
    for p in pars['plots'].itervalues():
        canvas.Clear()

        hA = tA.plot('h_'+p.name,p.formula,cut,bins=p.bins, options=p.options)
        hA.SetTitle(sA.title+';'+p.title)

        hB = tB.plot('h_'+p.name,p.formula,cut,bins=p.bins, options=p.options)
        hB.SetTitle(sB.title+';'+p.title)

        diff = ROOTplotter.H1DiffPlotter( logx = p.logx, logy = p.logy)
        diff.set(hA,hB)
        diff._ltitle = p.title
        diff._rtitle = sA.title+' vs. '+sB.title
        diff._rtitle = pars['ltitle'] # sA.title+' vs. '+sB.title
        diff.draw()

        map( lambda fmt: canvas.SaveAs(os.path.join(opt.out,p.name+'.'+fmt)), pars['formats'] )


    by2 = {}
    by2['split']    = ('_split',lambda h,t: h)
    by2['profileY'] = ('_py',ROOT.TH2.ProfileY)
    by2['profileX'] = ('_px',ROOT.TH2.ProfileX)

    for name,(ext,func) in by2.iteritems():
        if name not in pars: continue
        for p in pars[name].itervalues():
            canvas.Clear()
            canvas.Divide(1,2)

            canvas.cd(1)

            hA = tA.plot('h_'+p.name,p.formula,cut,bins=p.bins, options=p.options)
            hA.SetTitle(sA.title+';'+p.title)
            hA.SetMarkerColor(ROOT.kBlue)
            hA.SetLineColor(ROOT.kBlue)
            func(hA,sA.title+ext).Draw()
            

            canvas.cd(2)
            hB = tB.plot('h_'+p.name,p.formula,cut,bins=p.bins, options=p.options)
            hB.SetTitle(sB.title+';'+p.title)
            hB.SetMarkerColor(ROOT.kRed)
            hB.SetLineColor(ROOT.kRed)
            func(hB,sB.title+ext).Draw()

            map( lambda fmt: canvas.SaveAs(os.path.join(opt.out,p.name+ext+'.'+fmt)), pars['formats'] )



    

