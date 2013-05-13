#!/usr/bin/env python


import ROOT


def maketrees():

    datat1='''
    2 1 2 3
    2 2 3 4
    2 3 4 5
    '''

    datat2='''
    1 6 7 8
    1 8 9 10
    1 9 10 11
    '''

    # file 1
    with file('datat1.dat','w') as filet1: 
        print 'printing to datat1.dat'
        filet1.write(datat1)

    tfilet1 = ROOT.TFile.Open('datat1.root','recreate')
    t1 = ROOT.TTree('latino','t1')
    n = t1.ReadFile('datat1.dat','baseW:mll:mth:dphill')
    print 'Read',n,'lines'
    t1.Scan()
    tfilet1.Write()
    tfilet1.Close()

    # file 2
    with file('datat2.dat','w') as filet2: 
        print 'printing to datat2.dat'
        filet2.write(datat2)

    tfilet2 = ROOT.TFile.Open('datat2.root','recreate')
    t2 = ROOT.TTree('latino','t2')
    n = t2.ReadFile('datat2.dat','baseW:mll:mth:dphill')
    print 'Read',n,'lines'
    t2.Scan()
    tfilet2.Write()
    tfilet2.Close()


def testchain():
    from ginger.analysis import TreeAnalyser,Cut,CutFlow,Latino
    from ginger.tree import TreeWorker,ChainWorker,Sample

    orchard = '/shome/thea/HWW/work/dds/trees/top'

    s1 = Sample('latino',[orchard+'/nominals/latino_1250_ggToH250toWWTo2LAndTau2Nu.root'],  weight='baseW')
    s2 = Sample('latino',[orchard+'/nominals/latino_2250_vbfToH250toWWTo2LAndTau2Nu.root'], weight='baseW')

    s3A = Sample('latino',[orchard+'/nominals/latino_3250_wzttH250ToWW.root'],weight='baseW',preselection='(mctruth == 26)')
    s3B = Sample('latino',[orchard+'/nominals/latino_3250_wzttH250ToWW.root'],weight='baseW',preselection='(mctruth == 24)')
    s3C = Sample('latino',[orchard+'/nominals/latino_3250_wzttH250ToWW.root'],weight='baseW',preselection='(mctruth == 121)')

    print s1
    print s2

    t1  = TreeWorker.fromsample(s1)
    t2  = TreeWorker.fromsample(s2)
    t3A = TreeWorker.fromsample(s3A)
    t3B = TreeWorker.fromsample(s3B)
    t3C = TreeWorker.fromsample(s3C)

    testcut = 'mth < 100'
    testcut = '!sameflav'

    varexp='detajj'


    t1.scale = 2.
    t2.scale = 2.
    t1.scale = 1.
    t2.scale = 1.

    print 'y t1:',t1.yields(testcut)
    print 'y t2:',t2.yields(testcut)
    print 'y t3A:',t3A.yields(testcut)
    print 'y t3B:',t3B.yields(testcut)
    print 'y t3C:',t3C.yields(testcut)

    # chains
    c1 = ChainWorker.fromsamples( s1, s2, s3A, s3B, s3C )
    c1.scale = 2.
    c1.scale = 1.

    c2 = ChainWorker.fromsamples( s1, s2 )
    c3 = ChainWorker.fromsamples( s3A, s3B, s3C )

    c4 = ChainWorker( t1, t2, t3A, t3B, t3C )
    c5 = ChainWorker( c2, c3)

    print 'y c1:',c1.yields(testcut)
    print 'y c2:',c2.yields(testcut)
    print 'y c3:',c3.yields(testcut)
    print 'y c4:',c4.yields(testcut)
    print 'y c5:',c5.yields(testcut)

    bins = (15,0,15)
    bins = None
    bins = ([0,1,2,4,8],)
    bins = (20,0,200)
    bins = (20,0,10)
    
    # single tree histograms
    ht1  = t1.plot('ht1' ,varexp,testcut,bins=bins)
    ht2  = t2.plot('ht2' ,varexp,testcut,bins=bins)
    ht3A = t3A.plot('ht1',varexp,testcut,bins=bins)
    ht3B = t3B.plot('ht2',varexp,testcut,bins=bins)
    ht3C = t3C.plot('ht2',varexp,testcut,bins=bins)


    ht1.SetFillColor(ROOT.kRed)
    ht2.SetFillColor(ROOT.kBlue)
    ht3A.SetFillColor(ROOT.kGreen)
    ht3B.SetFillColor(ROOT.kGreen+1)
    ht3C.SetFillColor(ROOT.kGreen+2)

    hs1 = ROOT.THStack('chainme1','chainme1')
    hs1.Add(ht3C.Clone(),'hist')
    hs1.Add(ht3B.Clone(),'hist')
    hs1.Add(ht3A.Clone(),'hist')
    hs1.Add(ht2.Clone(),'hist')
    hs1.Add(ht1.Clone(),'hist')

    # overall histogram
    hc1 = c1.plot('hc1',varexp,testcut,bins=bins)

    # subchains
    hc2 = c2.plot('hc2',varexp,testcut,bins=bins)
    hc3 = c3.plot('hc3',varexp,testcut,bins=bins)

    hc2.SetFillColor(ROOT.kOrange)
    hc3.SetFillColor(ROOT.kOrange+1)

    hc1bis = hc1.Clone('hc1bis')
    hc1bis.SetMarkerStyle(20)
    hc1bis.SetMarkerSize(1)

    hs2 = ROOT.THStack('chainme2','chainme2')
    hs2.Add(hc3.Clone(),'hist')
    hs2.Add(hc2.Clone(),'hist')


    canv = ROOT.TCanvas('x','x',1500,2000)
    canv.Divide(3,4)

    canv.cd(1)
    ROOT.gPad.SetLogy()
    ht3A.Draw('hist')
    canv.cd(2)
    ROOT.gPad.SetLogy()
    ht3B.Draw('hist')
    canv.cd(3)
    ROOT.gPad.SetLogy()
    ht3C.Draw('hist')
    canv.cd(4)
    ROOT.gPad.SetLogy()
    ht1.Draw('hist')
    canv.cd(5)
    ROOT.gPad.SetLogy()
    ht2.Draw('hist')
    canv.cd(6)
    ROOT.gPad.SetLogy()
    hc1.Draw('hist')
    canv.cd(7)
    ROOT.gPad.SetLogy()
    hc2.Draw('hist')
    canv.cd(8)
    ROOT.gPad.SetLogy()
    hc3.Draw('hist')

    canv.cd(9)
    ROOT.gPad.SetLogy()
    hs1.Draw()
    hc1bis.Draw('same')

    canv.cd(10)
    ROOT.gPad.SetLogy()
    hs2.Draw()
    hc1bis.Draw('same')
    canv.cd(11)
    
    hdiff = hc1.Clone('hdiff')
    hdiff.Add(hs1.GetStack().Last(),-1)

    hdiff.SetMarkerStyle(25)
    hdiff.SetMarkerSize(2)

    hdiff.Draw('hist E0 L')

    hdiff2 = hc1.Clone('hdiff')
    hdiff2.Add(hs2.GetStack().Last(),-1)
    hdiff2.SetMarkerStyle(24)
    hdiff2.SetMarkerSize(2)

    hdiff.Draw('hist E0 L same')
    
#     canv.Print('www/test/chained.png')
    canv.Print('www/test/chained.pdf')

def testviews():
    from ginger.analysis import TreeAnalyser,Cut,CutFlow,Latino
    from ginger.tree import TreeWorker,ChainWorker,Sample,TreeView,ChainView

    orchard = '/shome/thea/HWW/work/dds/trees/top'

    s1 = Sample('latino',[orchard+'/nominals/latino_1250_ggToH250toWWTo2LAndTau2Nu.root'],  weight='baseW')
    s2 = Sample('latino',[orchard+'/nominals/latino_2250_vbfToH250toWWTo2LAndTau2Nu.root'], weight='baseW')

    s3A = Sample('latino',[orchard+'/nominals/latino_3250_wzttH250ToWW.root'],weight='baseW',preselection='(mctruth == 26)')
    s3B = Sample('latino',[orchard+'/nominals/latino_3250_wzttH250ToWW.root'],weight='baseW',preselection='(mctruth == 24)')
    s3C = Sample('latino',[orchard+'/nominals/latino_3250_wzttH250ToWW.root'],weight='baseW',preselection='(mctruth == 121)')

    testcut = 'mth < 100'
    testcut = '!sameflav'

    varexp='detajj'

    # workers
    tw1  = TreeWorker.fromsample(s1)
    tw2  = TreeWorker.fromsample(s2)
    tw3A = TreeWorker.fromsample(s3A)
    tw3B = TreeWorker.fromsample(s3B)
    tw3C = TreeWorker.fromsample(s3C)

    t1  = TreeView( tw1 , testcut ) 
    t2  = TreeView( tw2 , testcut )
    t3A = TreeView( tw3A, testcut )
    t3B = TreeView( tw3B, testcut )
    t3C = TreeView( tw3C, testcut )

    print 'y t1:',t1.yields()
    print 'y t2:',t2.yields()
    print 'y t3A:',t3A.yields()
    print 'y t3B:',t3B.yields()
    print 'y t3C:',t3C.yields()

    # chains
    cw1 = ChainWorker.fromsamples( s1, s2, s3A, s3B, s3C )
    cw2 = ChainWorker.fromsamples( s1, s2 )
    cw3 = ChainWorker.fromsamples( s3A, s3B, s3C )

    cw4 = ChainWorker( tw1, tw2, tw3A, tw3B, tw3C )
    cw5 = ChainWorker( cw2, cw3)

    c1 = ChainView( cw1, testcut ) 
    c2 = ChainView( cw2, testcut )
    c3 = ChainView( cw3, testcut )
    c4 = ChainView( cw4, testcut )
    c5 = ChainView( cw5, testcut )


    print 'y c1:',c1.yields()
    print 'y c2:',c2.yields()
    print 'y c3:',c3.yields()
    print 'y c4:',c4.yields()
    print 'y c5:',c5.yields()

    bins = (15,0,15)
    bins = None
    bins = ([0,1,2,4,8],)
    bins = (20,0,200)
    bins = (20,0,10)
    
    # single tree histograms
    ht1  = t1.plot('ht1' ,varexp,bins=bins)
    ht2  = t2.plot('ht2' ,varexp,bins=bins)
    ht3A = t3A.plot('ht1',varexp,bins=bins)
    ht3B = t3B.plot('ht2',varexp,bins=bins)
    ht3C = t3C.plot('ht2',varexp,bins=bins)


    ht1.SetFillColor(ROOT.kRed)
    ht2.SetFillColor(ROOT.kBlue)
    ht3A.SetFillColor(ROOT.kGreen)
    ht3B.SetFillColor(ROOT.kGreen+1)
    ht3C.SetFillColor(ROOT.kGreen+2)

    hs1 = ROOT.THStack('chainme1','chainme1')
    hs1.Add(ht3C.Clone(),'hist')
    hs1.Add(ht3B.Clone(),'hist')
    hs1.Add(ht3A.Clone(),'hist')
    hs1.Add(ht2.Clone() ,'hist')
    hs1.Add(ht1.Clone() ,'hist')

    # overall histogram
    hc1 = c1.plot('hc1',varexp,bins=bins)

    # subchains
    hc2 = c2.plot('hc2',varexp,bins=bins)
    hc3 = c3.plot('hc3',varexp,bins=bins)

    hc2.SetFillColor(ROOT.kOrange)
    hc3.SetFillColor(ROOT.kOrange+1)

    hc1bis = hc1.Clone('hc1bis')
    hc1bis.SetMarkerStyle(20)
    hc1bis.SetMarkerSize(1)

    hs2 = ROOT.THStack('chainme2','chainme2')
    hs2.Add(hc3.Clone(),'hist')
    hs2.Add(hc2.Clone(),'hist')


    canv = ROOT.TCanvas('x','x',1500,2000)
    canv.Divide(3,4)

    canv.cd(1)
    ROOT.gPad.SetLogy()
    ht3A.Draw('hist')
    canv.cd(2)
    ROOT.gPad.SetLogy()
    ht3B.Draw('hist')
    canv.cd(3)
    ROOT.gPad.SetLogy()
    ht3C.Draw('hist')
    canv.cd(4)
    ROOT.gPad.SetLogy()
    ht1.Draw('hist')
    canv.cd(5)
    ROOT.gPad.SetLogy()
    ht2.Draw('hist')
    canv.cd(6)
    ROOT.gPad.SetLogy()
    hc1.Draw('hist')
    canv.cd(7)
    ROOT.gPad.SetLogy()
    hc2.Draw('hist')
    canv.cd(8)
    ROOT.gPad.SetLogy()
    hc3.Draw('hist')

    canv.cd(9)
    ROOT.gPad.SetLogy()
    hs1.Draw()
    hc1bis.Draw('same')

    canv.cd(10)
    ROOT.gPad.SetLogy()
    hs2.Draw()
    hc1bis.Draw('same')
    canv.cd(11)
    
    hdiff = hc1.Clone('hdiff')
    hdiff.Add(hs1.GetStack().Last(),-1)

    hdiff.SetMarkerStyle(25)
    hdiff.SetMarkerSize(2)

    hdiff.Draw('hist E0 L')

    hdiff2 = hc1.Clone('hdiff')
    hdiff2.Add(hs2.GetStack().Last(),-1)
    hdiff2.SetMarkerStyle(24)
    hdiff2.SetMarkerSize(2)

    hdiff.Draw('hist E0 L same')
    
#     canv.Print('www/test/chained.png')
    canv.Print('www/test/chained_views.pdf')

def testviews2D():
    from ginger.analysis import TreeAnalyser,Cut,CutFlow,Latino
    from ginger.tree import TreeWorker,ChainWorker,Sample,TreeView,ChainView

    orchard = '/shome/thea/HWW/work/dds/trees/top'

    s1 = Sample('latino',[orchard+'/nominals/latino_1250_ggToH250toWWTo2LAndTau2Nu.root'],  weight='baseW')
    s2 = Sample('latino',[orchard+'/nominals/latino_2250_vbfToH250toWWTo2LAndTau2Nu.root'], weight='baseW')

    s3A = Sample('latino',[orchard+'/nominals/latino_3250_wzttH250ToWW.root'],weight='baseW',preselection='(mctruth == 26)')
    s3B = Sample('latino',[orchard+'/nominals/latino_3250_wzttH250ToWW.root'],weight='baseW',preselection='(mctruth == 24)')
    s3C = Sample('latino',[orchard+'/nominals/latino_3250_wzttH250ToWW.root'],weight='baseW',preselection='(mctruth == 121)')

    testcut = 'mth < 100'
    testcut = '!sameflav'

    varexp='mll:detajj'

    # workers
    tw1  = TreeWorker.fromsample(s1)
    tw2  = TreeWorker.fromsample(s2)
    tw3A = TreeWorker.fromsample(s3A)
    tw3B = TreeWorker.fromsample(s3B)
    tw3C = TreeWorker.fromsample(s3C)

    t1  = TreeView( tw1 , testcut ) 
    t2  = TreeView( tw2 , testcut )
    t3A = TreeView( tw3A, testcut )
    t3B = TreeView( tw3B, testcut )
    t3C = TreeView( tw3C, testcut )

    print 'y t1:',t1.yields()
    print 'y t2:',t2.yields()
    print 'y t3A:',t3A.yields()
    print 'y t3B:',t3B.yields()
    print 'y t3C:',t3C.yields()

    # chains
    cw1 = ChainWorker.fromsamples( s1, s2, s3A, s3B, s3C )
    cw2 = ChainWorker.fromsamples( s1, s2 )
    cw3 = ChainWorker.fromsamples( s3A, s3B, s3C )

    cw4 = ChainWorker( tw1, tw2, tw3A, tw3B, tw3C )
    cw5 = ChainWorker( cw2, cw3)

    c1 = ChainView( cw1, testcut ) 
    c2 = ChainView( cw2, testcut )
    c3 = ChainView( cw3, testcut )
    c4 = ChainView( cw4, testcut )
    c5 = ChainView( cw5, testcut )


    print 'y c1:',c1.yields()
    print 'y c2:',c2.yields()
    print 'y c3:',c3.yields()
    print 'y c4:',c4.yields()
    print 'y c5:',c5.yields()

    bins = (15,0,15)
    bins = None
    bins = ([0,1,2,4,8],)
    bins = (20,0,200)
    bins = (20,0,10,20,0,300)
    
    # single tree histograms
    ht1  = t1.plot('ht1' ,varexp,bins=bins)
    ht2  = t2.plot('ht2' ,varexp,bins=bins)
    ht3A = t3A.plot('ht1',varexp,bins=bins)
    ht3B = t3B.plot('ht2',varexp,bins=bins)
    ht3C = t3C.plot('ht2',varexp,bins=bins)

    ht1.SetFillColor(ROOT.kRed)
    ht2.SetFillColor(ROOT.kBlue)
    ht3A.SetFillColor(ROOT.kGreen)
    ht3B.SetFillColor(ROOT.kGreen+1)
    ht3C.SetFillColor(ROOT.kGreen+2)

    hs1 = ROOT.THStack('chainme1','chainme1')
    hs1.Add(ht3C.Clone(),'hist')
    hs1.Add(ht3B.Clone(),'hist')
    hs1.Add(ht3A.Clone(),'hist')
    hs1.Add(ht2.Clone() ,'hist')
    hs1.Add(ht1.Clone() ,'hist')

    # overall histogram
    hc1 = c1.plot('hc1',varexp,bins=bins)

    # subchains
    hc2 = c2.plot('hc2',varexp,bins=bins)
    hc3 = c3.plot('hc3',varexp,bins=bins)

    hc2.SetFillColor(ROOT.kOrange)
    hc3.SetFillColor(ROOT.kOrange+1)

    hc1bis = hc1.Clone('hc1bis')
    hc1bis.SetMarkerStyle(20)
    hc1bis.SetMarkerSize(1)

    hs2 = ROOT.THStack('chainme2','chainme2')
    hs2.Add(hc3.Clone(),'hist')
    hs2.Add(hc2.Clone(),'hist')


    canv = ROOT.TCanvas('x','x',1500,2000)
    canv.Divide(3,4)

    canv.cd(1)
    ROOT.gPad.SetLogz()
    ht3A.Draw('colz')
    canv.cd(2)
    ROOT.gPad.SetLogz()
    ht3B.Draw('colz')
    canv.cd(3)
    ROOT.gPad.SetLogz()
    ht3C.Draw('colz')
    canv.cd(4)
    ROOT.gPad.SetLogz()
    ht1.Draw('colz')
    canv.cd(5)
    ROOT.gPad.SetLogz()
    ht2.Draw('colz')
    canv.cd(6)
    ROOT.gPad.SetLogz()
    hc1.Draw('colz')
    canv.cd(7)
    ROOT.gPad.SetLogz()
    hc2.Draw('colz')
    canv.cd(8)
    ROOT.gPad.SetLogz()
    hc3.Draw('colz')

#     canv.cd(9)
#     ROOT.gPad.SetLogy()
#     hs1.Draw()
#     hc1bis.Draw('same')

#     canv.cd(10)
#     ROOT.gPad.SetLogy()
#     hs2.Draw()
#     hc1bis.Draw('same')
    canv.cd(11)
    
    hdiff = hc1.Clone('hdiff')
    hdiff.Add(hs1.GetStack().Last(),-1)

    hdiff.SetMarkerStyle(25)
    hdiff.SetMarkerSize(2)

    hdiff.Draw('text 0')

    canv.cd(12)
    hdiff2 = hc1.Clone('hdiff')
    hdiff2.Add(hs2.GetStack().Last(),-1)
    hdiff2.SetMarkerStyle(24)
    hdiff2.SetMarkerSize(2)

    hdiff.Draw('text')
    
#     canv.Print('www/test/chained.png')
    canv.Print('www/test/chained_views_2D.pdf')


def testviews_old():
    from ginger.analysis import TreeAnalyser,Cut,CutFlow,Latino
    from ginger.tree import TreeWorker,ChainWorker,Sample,TreeView,ChainView

#     s1 = Sample('latino',['datat1.root'], weight='baseW')
#     s2 = Sample('latino',['datat2.root'], weight='baseW')
    
    orchard = '/shome/thea/HWW/work/dds/trees/top'
    s1 = Sample('latino',[orchard+'/nominals/latino_1250_ggToH250toWWTo2LAndTau2Nu.root'],  weight='baseW')
    s2 = Sample('latino',[orchard+'/nominals/latino_2250_vbfToH250toWWTo2LAndTau2Nu.root'], weight='baseW')

    print s1
    print s2

    t1 = TreeWorker.fromsample(s1)
    t2 = TreeWorker.fromsample(s2)
    c1 = ChainWorker.fromsamples( s1, s2 )

    viewcut = 'mth < 200'
    varexp  = 'mth'
    bins = (30,10,310)

    v1A = TreeView(t1)
    v2A = TreeView(t2)
    w1A = ChainView(c1)

    print '#----'
    print 'Y vt1A:',v1A.yields()
    print 'Y vt2A:',v2A.yields()
    print 'Y vc1A:',w1A.yields()

    print '#----'
    v1B = v1A.spawn(viewcut)
    v2B = v2A.spawn(viewcut)
    w1B = w1A.spawn(viewcut)

    print 'Y vt1B:',v1B.yields()
    print 'Y vt2B:',v2B.yields()
    print 'Y vc1B:',w1B.yields()

    hv1 = v1B.plot('hv1B',varexp,bins=bins)
    hv2 = v2B.plot('hv2B',varexp,bins=bins)

    hv1.SetFillColor(ROOT.kRed)
    hv2.SetFillColor(ROOT.kBlue)

    hw1 = w1B.plot('hw1B',varexp,bins=bins)
    hw1bis = hw1.Clone('hc1bis')
    hw1bis.SetMarkerStyle(20)
    hw1bis.SetMarkerSize(1)

    hs = ROOT.THStack('chainme','chainme')
    hs.Add(hv1.Clone(),'hist')
    hs.Add(hv2.Clone(),'hist')

    canv = ROOT.TCanvas('x','x',1500,1000)
    canv.Divide(3,2)

    canv.cd(1)
    hv1.Draw('hist')
    canv.cd(2)
    hv2.Draw('hist')
    canv.cd(3)
    hw1.Draw('hist')
    canv.cd(4)
    hs.Draw()
    hw1bis.Draw('same')
    canv.cd(5)

    hdiff = hw1.Clone('hdiff')
    hdiff -= hs.GetStack().Last()

    hdiff.SetMarkerStyle(20)
    hdiff.SetMarkerSize(1)

    hdiff.Draw('hist E0')

#     canv.Print('www/test/chained_views.png')
    canv.Print('www/test/chained_views.pdf')

if __name__ == '__main__':
    import optparse
    import hwwtools

    parser = optparse.OptionParser()
    parser.add_option('--maketrees', dest = 'maketrees', help='Create test trees' , default=False , action='store_true' )
    parser.add_option('--testchain', dest = 'testchain', help='Test chain objs'   , default=False , action='store_true' )
    parser.add_option('--testviews', dest = 'testviews', help='Test views objs'   , default=False , action='store_true' )
    parser.add_option('--test2D'   , dest = 'test2D'   , help='Test views 2D'     , default=False , action='store_true' )
    parser.add_option('-d', '--debug'    , dest = 'debug'    , help='Debug level'       , default=0     , action='count' )

    (opt, args) = parser.parse_args()

    hwwtools.setDebugLevel(opt)


    if opt.maketrees:
        maketrees()

    if opt.testchain:
        testchain()

    if opt.testviews:
        testviews()
        
    if opt.test2D:
        testviews2D()

