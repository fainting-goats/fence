#!/bin/env python

import ROOT
ROOT.TH1.SetDefaultSumw2()

def resize(x,ratio):
    x.SetLabelSize(x.GetLabelSize()*ratio/(1-ratio))
    x.SetTitleSize(x.GetTitleSize()*ratio/(1-ratio))


ratio = 0.7
outer = 0.1
inner = 0.02
marg  = 0.1
ltitle = 'la rava'
rtitle = 'e la fava'
ytitle2 = '(h1-h2)/h1'
legSize = (0.2,0.2)

h1 = ROOT.TH1F('aa','aaaaaaaaaaaaa;ax;ayo\'',100,0,100)
h2 = ROOT.TH1F('bb','bbbbbbbbbbbbb;bx;by',100,0,100)

hFill = ROOT.TH1F.Fill
gaus = ROOT.gRandom.Gaus

entries = 100000

for i in xrange(entries):
    hFill(h1,gaus(50,10))
    hFill(h2,gaus(51,10))

c = ROOT.TCanvas('xx','xx')

#- pad1 ---
pad1 = ROOT.TPad('pad1','pad1',0.,(1-ratio),1.,1.)
pad1.SetLeftMargin(marg)
pad1.SetRightMargin(marg)
pad1.SetTopMargin(outer/ratio)
pad1.SetBottomMargin(inner/ratio)
pad1.SetTicks()
pad1.Draw()

pad1.cd()


hists = [h1,h2]
map(lambda h: ROOT.TH1.SetLineWidth(h,2),hists)
h1.SetLineColor(ROOT.kRed)
h2.SetLineColor(ROOT.kBlue)

ymax = max([ h.GetMaximum() for h in hists ])
ymin = min([ h.GetMinimum() for h in hists ])
xmax = max([ h.GetXaxis().GetXmax() for h in hists ])
xmin = min([ h.GetXaxis().GetXmin() for h in hists ])

ymin = ymin if ymin == 0 else ymin*1.1
ymax = ymax if ymax == 0 else ymax*1.1

print xmin,xmax
print ymin,ymax

hframe = pad1.DrawFrame(xmin,ymin,xmax,ymax)
ROOT.SetOwnership(hframe,False)
hframe.GetXaxis().SetLabelSize(0.00)
hframe.GetYaxis().SetTitle(h1.GetYaxis().GetTitle())

h1.Draw('same hist')
h2.Draw('same hist')

anchor = (1-marg,1-outer/ratio)
leg = ROOT.TLegend(anchor[0]-legSize[0],anchor[1]-legSize[1],anchor[0],anchor[1],'','NDC')
leg.SetFillColor(ROOT.kWhite)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
# leg.SetNColumns(2)
leg.AddEntry(h1,'','l')
leg.AddEntry(h2,'','l')

leg.Draw()


l = ROOT.TLatex()
l.SetNDC()
l.SetTextAlign(12)
l.DrawText(marg,1-(0.5*outer/ratio),ltitle)
l.SetTextAlign(32)
l.DrawText(1-marg,1-(0.5*outer/ratio),rtitle)


#- pad2 ---

c.cd()
pad2 = ROOT.TPad('pad1','pad1',0.,0.0,1.,(1-ratio))
pad2.SetTopMargin(inner/(1-ratio))
pad2.SetBottomMargin(outer/(1-ratio))
pad2.SetLeftMargin(marg)
pad2.SetRightMargin(marg)
pad2.SetTicks()
pad2.SetGridy()
pad2.Draw()

pad2.cd()

hdiff = h1.Clone('diff')
hdiff.Add(h2,-1)
hdiff.Divide(h1)

hdiff.SetMarkerStyle(20)
hdiff.SetLineWidth(1)
hdiff.SetLineColor(ROOT.kBlack)
hframe = pad2.DrawFrame(xmin,-0.15,xmax,0.15)
ROOT.SetOwnership(hframe,False)


ax = hframe.GetXaxis()
ay = hframe.GetYaxis()
ax.SetTitle(hdiff.GetXaxis().GetTitle())
ay.SetTitle(ytitle2)
ay.SetTitleOffset(ay.GetTitleOffset()/ratio*(1-ratio) )
# ax.SetLabelSize(ax.GetLabelSize()*ratio/(1-ratio))
# ay.SetLabelSize(ay.GetLabelSize()*ratio/(1-ratio))
resize(ax,ratio)
resize(ay,ratio)

hdiff.Draw('same')

# l10p,l05p,l05m,l10m = [ ROOT.TGraph(2) for i in range(4)]
# l10p.SetPoint(0,xmin,0.1)
# l10p.SetPoint(1,xmax,0.1)
# l10p.Draw()







