#!/bin/env python

from ginger.painter import Pad,Canvas,Legend
import sys
import ROOT

if __name__ == '__main__':

    sys.argv.append('-b')
    ROOT.gROOT.SetBatch()

    def probe():
        print '-'*80
        pad = ROOT.gPad.func()
        pw = pad.GetWw()*pad.GetWNDC()
        ph = pad.GetWh()*pad.GetHNDC()

        uxmin = pad.GetUxmin()
        uxmax = pad.GetUxmax()
        uymin = pad.GetUymin()
        uymax = pad.GetUymax()

        dx = (uxmax-uxmin)/(1-pad.GetLeftMargin()-pad.GetRightMargin())
        dy = (uymax-uymin)/(1-pad.GetTopMargin()-pad.GetBottomMargin())

        print pw,ph
        dummy = ROOT.TLatex(0,0,'stocazzo')
        dummy.SetTextFont(44)
        dummy.SetTextSize(20)
        dummy.Print()
        print dummy.GetTitle(),dummy.GetXsize()*pw/(ph*dx),dummy.GetYsize()/dy
        dummy.Draw()
        ROOT.SetOwnership(dummy,False)
        print '-'*80

    c = Canvas(minsize = (0,0))

    axsty = {
    'labelfamily'       : 4,
    'labelsize'         : 20,
    'labeloffset'       : 5,
    'titlefamily'       : 4,
    'titlesize'         : 20,
    'titleoffset'       : 50,
    'ticklength'        : 100
    }
#     axsty = {'label-family':42, 'label-size':0.04, }
    nolabs = axsty.copy()
    nolabs['labelsize'] = 0
    nolabs['titlesize'] = 0

    p0 = Pad('p0',500,200, margins=(60,60,60,60), xaxis = axsty,  yaxis = axsty, align=('l','m')  )
    p1 = Pad('p1',500,500, margins=(80,20,20,20), xaxis = nolabs, yaxis = axsty )
    p2 = Pad('p2',500,200, margins=(80,20,20,20), xaxis = nolabs, yaxis = axsty )
    p3 = Pad('p3',500,200, margins=(80,20,20,80), xaxis = axsty,  yaxis = axsty )
    p4 = Pad('p4',200,200, margins=(20,80,20,80), xaxis = axsty,  yaxis = axsty )
    p5 = Pad('p5',200,200, margins=(20,80,80,20), xaxis = nolabs, yaxis = axsty, align=('l','m') )
    p1._xaxis['label-size'] = 0.0
    p2._xaxis['label-size'] = 0.0
#     p3._xaxis['title-offset'] = 2.5
#     p4._xaxis['title-offset'] = 2.5

    c[1,0] = p0
    c[0,0] = p1
    c[0,1] = p2
    c[0,2] = p3
    c[1,2] = p4
    c[1,1] = p5

    tc = c.makecanvas()
    tc.SetName('aaa')

    p5.SetFillColor(ROOT.kRed)

    bins = 10
    hdummy = ROOT.TH1F('dummy','',bins,0,bins)
    hs = ROOT.THStack('stack','stocazz')
    hcols = []
    for i in xrange(bins):
        h = hdummy.Clone('col%d' % i)
        h.SetTitle(h.GetName())
        h.SetFillColor(i+ROOT.kOrange)
        h.SetLineColor(i+ROOT.kOrange)
        h.SetFillStyle(3001)
        h.SetLineWidth(2)
        h.Fill(i,i)
        ROOT.SetOwnership(h,False)
        hs.Add(h)
        hcols.append(h)

#     hcols[0].SetTitle('stocazz')

    p1.cd()
    p1.SetTicks()
    hs1 = hs.Clone('xxx')
    hs1.Draw()
    p1.Modified()
    p1.Update()
    hs1.GetYaxis().SetTitle('y-axis')
    hs1.GetXaxis().SetTitle('x-axis')
#     hs.Draw()
    leg = Legend(4,4, 30, anchor=(90,30),align=('l','t') )
    sequence = leg.sequence
    sequence.remove( (1,3) )
    sequence.remove( (2,3) )
    sequence.remove( (2,2) )
    leg.sequence = sequence
    leg.addentry(hcols[0],'f')
    leg.addentry(hcols[1],'f')
    leg.addentry(hcols[2],'f')
    leg.addentry(hcols[3],'f')
    leg.addentry(hcols[4],'f')
    leg.addentry(hcols[5],'f')
    leg.addentry(hcols[6],'f')
    leg.addentry(hcols[7],'f')
    leg.addentry(hcols[8],'f')
    leg.addentry(hcols[9],'f')
    leg.draw()
    probe()

    p2.cd()
    hs2 = hs.Clone('yyy')
    hs2.Draw()
    hs2.GetYaxis().SetTitle('y-axis')
#     hs.Draw()

    pad = c.get(0,2)
    pad.cd()
    d3 = hdummy.Clone('d3')
    d3.GetYaxis().SetTitle('y-axis')
    d3.GetXaxis().SetTitle('x-axis')
    d3.Draw()

    pad = c.get(1,2)
    pad.cd()
    d4 = hdummy.Clone('d4')
    d4.GetXaxis().SetTitle('x-axis')
    d4.Draw('Y+')

    pad = c.get(1,1)
    pad.cd()
    d5 = hdummy.Clone('d5')
    d5.GetXaxis().SetTitle('x-axis')
    d5.Draw('Y+')

#     p0.cd()
#     d0 = hdummy.Clone('d0')
#     d0.GetXaxis().SetTitle('x-axis')
#     d0.Draw('Y+')

    pad = c.get(1,0)
    pad.cd()

    hdummy = ROOT.TH1D("pippo","pippo;topolino;paperino",10,0,1)
    hdummy.Draw()





    c.applystyle()
    c.SetFillColor(ROOT.kOrange)
    #p0.SetFillStyle(0)
    #p1.SetFillStyle(0)
    #p2.SetFillStyle(0)
    #p3.SetFillStyle(0)
    #p4.SetFillStyle(0)
    #p5.SetFillStyle(0)

#     tc.ls()
    ROOT.gSystem.ProcessEvents()


#     tc.Print('des.png')
    tc.Print('des.pdf')
    tc.Print('des.png')

    ROOT.gStyle.SetPaperSize(200,240)
    tc.Print('des_paper.pdf')
    tc.Print('des_paper.png')

    tc.SetCanvasSize(tc.GetWw()/2,tc.GetWh()/2)
    tc.Modified()
    tc.Update()
    tc.Print('des_small.pdf')
    tc.Print('des_small.png')




