#!/bin/env python

import ROOT
import sys
from ginger.painter import Pad,Canvas,Legend,Latex
from HWWAnalysis.Misc.ROOTAndUtils import TH1AddDirSentry
import numpy
import array

class Directory:

    @property
    def name(self):
        return self._name

    @property
    def objs(self):
        return dict([ (k,a)  for k,a in self.__dict__.iteritems() if (not isinstance(a,Directory) and not k[0] == '_') ] )

    @property
    def dirs(self):
        return dict([ (k,a)  for k,a in self.__dict__.iteritems() if (isinstance(a,Directory) and not k[0] == '_') ] )

    #---
    def get(self, name):
        if name[0] == '_': return None
        return self.__dict__[name]

    #---
    def __getitem__(self, name):
        return self.get(name)

    #---
    def show(self, indent=0, level=-1):
        print ' '*indent+'+',self._name
        for k,a in self.__dict__.iteritems():
            if k[0] == '_': continue
            if isinstance(a,Directory):
                a.show(indent+1)
            else:
                print ' '*(indent+1)+'-',k,'=',a

    #---
    def __init__(self,tdir):
        if not isinstance(tdir,ROOT.TDirectory):
            raise TypeError('Can\t build a Directory from a %s' % tdir.__class__.__name__)

#         dir = Directory(tdir.GetName())
        self._name = tdir.GetName()

        for k in tdir.GetListOfKeys():
            o = None
            if k.GetClassName().startswith('TDirectory'):
                o = Directory(k.ReadObj())
            else:
                o = k.ReadObj()

            setattr(self, k.GetName(), o)


# bin 2 axis definition
def bins2axis( axis ):

    if isinstance(axis,tuple):
        if   len(axis) == 1:
#             n = axis[0]
#             axdef = (n,0,n)
            nx = axis[0]
            xbins = range(nx+1)
        elif len(axis) == 3:
#             n = axis[0]
#             axdef = axis
            nx,xmin,xmax = axis
            xbins = numpy.arange(xmin,xmax,(xmax-xmin)/nx)
        else:
            raise ValueError('What\'s this?')
    elif isinstance(axis,list):
#         n = len(axis)-1
#         axdef = (n, array.array('d',axis))
        nx = len(axis)-1
        xbins = axis

    return nx,array.array('d',xbins)


#---
def fold(h,xaxis,yaxis,folding='xy'):

    nx,xdef = bins2axis(xaxis)
    ny,ydef = bins2axis(yaxis)

    if h.GetNbinsX() != (nx*ny):
        raise ValueError('Can\'t fit %d bins in (%dx%d)' % (h.GetNbinsX(),nx,ny))

    hclass = getattr(ROOT,h.__class__.__name__.replace('TH1','TH2'))
    hfolded = hclass(h.GetName()+'_folded',h.GetTitle()+'_folded',nx,xdef,ny,ydef)

    if   folding == 'xy':
        for k in xrange(h.GetNbinsX()):
            i,j = k / ny, k % ny
            val = h.GetBinContent(k+1)
            hfolded.SetBinContent(i+1,j+1,val)
    elif folding == 'yx':
        for k in xrange(h.GetNbinsX()):
            i,j = k % nx, k / nx
            val = h.GetBinContent(k+1)
            hfolded.SetBinContent(i,j,val)
    else:
        raise ValueError('%s is not a valid options')

    return hfolded

class RollProjector:

    def __init__(self, axis,proj,nfold=None):
        n,axdef = bins2axis(axis)
        self._nbins = n
        self._axdef = axdef

        # this is the number of folds contained in the original object
        # for instance in the mth x mll case (13x9 bins)
        #  mth=1    |mth=2  |mth=3    
        #  123456789123456789123456789...
        # 
        self._nfold = nfold if nfold else n
        self._proj = proj
        # projection axis
        if proj == 'x':
            self._pax = 1
        elif proj == 'y':
            self._pax = 0

    def project(self,o):
        
        if isinstance(o,ROOT.TH1):
            return self._projTH1(o)
        elif isinstance(o,ROOT.TGraphAsymmErrors):
            return self._projTGraph(o)

    def _projTH1(self,h):
        sentry = TH1AddDirSentry()

        typestr = h.__class__.__name__[3:]
        if   typestr == 'D':
            htype = numpy.double
        elif typestr == 'F':
            htype = numpy.float
        elif typestr == 'I':
            htype = numpy.int32
        else:
            raise ValueError('Boh?!??!')
        
        # add a 2 because the array contains under and overflow
        harray = numpy.ndarray( (h.GetNbinsX()+2,),dtype=htype, buffer=h.GetArray() )
#         print len(harray[1:-1])
#         print self._nfold,h.GetNbinsX()
        hsplit =  numpy.hsplit(harray[1:-1],self._nfold)
        newarray = numpy.sum(hsplit,axis=self._pax) 

        print 'Entries',numpy.sum(hsplit)

#         print self._axdef
        hproj = (h.__class__)('%s_proj_%s' % (h.GetName(),self._proj), '%s proj %s' % (h.GetTitle(),self._proj),self._nbins, self._axdef )

#         print newarray
        
        for i,x in enumerate(newarray):
            hproj.SetAt(x,i+1)
            
        ROOT.TAttLine.Copy(h,hproj)
        ROOT.TAttFill.Copy(h,hproj)
        ROOT.TAttMarker.Copy(h,hproj)

        return hproj

    # ---
    def _projTGraph(self,g):

        sentry = TH1AddDirSentry()
        
        y      = numpy.ndarray( (g.GetN(),),dtype=numpy.double, buffer=g.GetY() )
        ysplit = numpy.hsplit(y,self._nfold)
        p_y    = numpy.sum(ysplit,axis=self._pax) 

        eyh = numpy.ndarray( (g.GetN(),),dtype=numpy.double, buffer=g.GetEYhigh() )
        eyh2_split = numpy.hsplit( (eyh**2) ,self._nfold)
        p_eyh      = numpy.sqrt( numpy.sum(eyh2_split,axis=self._pax) )

        eyl = numpy.ndarray( (g.GetN(),),dtype=numpy.double, buffer=g.GetEYlow() )
        eyl2_split = numpy.hsplit( (eyl**2) ,self._nfold)
        p_eyl      = numpy.sqrt( numpy.sum(eyl2_split,axis=self._pax) )

        x = array.array('d',[0]*self._nbins)
        exh = array.array('d',[0]*self._nbins)
        exl = array.array('d',[0]*self._nbins)
        for i in xrange(self._nbins):
            x[i]            = (self._axdef[i+1]+self._axdef[i])/2.
            exh[i] = exl[i] = (self._axdef[i+1]-self._axdef[i])/2.

        p_g = ROOT.TGraphAsymmErrors(self._nbins, x, p_y, exl, exh, p_eyl, p_eyh)
        p_g.SetNameTitle('%s_proj_%s' % (g.GetName(),self._proj),'%s proj %s' % (g.GetTitle(),self._proj))

        ROOT.TAttLine.Copy(g,p_g)
        ROOT.TAttFill.Copy(g,p_g)
        ROOT.TAttMarker.Copy(g,p_g)
        return p_g


# ---
def tester():
    
    dumpfile = ROOT.TFile.Open('dumplings.root')
    
    d = Directory(dumpfile)

    d.info.show()

    processes = [str(p) for p in d.info.processes]
    signals   = [str(s) for s in d.info.signals]
    nuisances = [str(n) for n in d.info.nuisances]

#     print nuisances

    nuismap = dict([(str(k),[str(s) for s in d.info.map_binnuisances.GetValue(k)]) for k in d.info.map_binnuisances])
    bins = nuismap.keys()


    b = bins[0]
    p = 'Top'
    n = nuismap[b][0]
    n = 'all'

    print b,n,p

    sentry = TH1AddDirSentry()
    
    def setstyle(obj, **opts):
        methods = [
            ('linewidth'  , obj.SetLineWidth), 
            ('linecolor'  , obj.SetLineColor), 
            ('linestyle'  , obj.SetLineStyle), 
            ('fillstyle'  , obj.SetFillStyle), 
            ('fillcolor'  , obj.SetFillColor), 
            ('markercolor', obj.SetMarkerColor), 
            ('markerstyle', obj.SetMarkerStyle), 
            ('markersize' , obj.SetMarkerSize), 
        ]

        for l,m in methods:
            x = opts.get(l,None)
            if not x is None: m(x)#; print x,m.__name__ 

    h_data = d.init[b]['histo_Data']
    h_init = d.init[b]['histo_%s' % p]
    e_init = d.init[b][p]['%s_errs_%s' % (p,n)]
    h_sig  = d.sig[b]['histo_%s' % p]
    e_sig  = d.sig[b][p]['%s_errs_%s' % (p,n)]
    h_bkg  = d.bkg[b]['histo_%s' % p]
    e_bkg  = d.bkg[b][p]['%s_errs_%s' % (p,n)]

    setstyle(h_data,markerstyle=20,markersize=1)
    setstyle(h_init,linewidth=2,linecolor=ROOT.kRed)
    setstyle(e_init,linewidth=2,linecolor=ROOT.kRed,fillstyle=3001,fillcolor=ROOT.kRed)
    setstyle(h_sig ,linewidth=2,linecolor=ROOT.kGreen)
    setstyle(e_sig ,linewidth=2,linecolor=ROOT.kGreen, fillstyle=3005,fillcolor=ROOT.kGreen)
    setstyle(h_bkg ,linewidth=2,linecolor=ROOT.kBlue)
    setstyle(e_bkg ,linewidth=2,linecolor=ROOT.kBlue, fillstyle=3004,fillcolor=ROOT.kBlue)

    mth = [60,70,80,90,100,110,120,140,160,180,200,220,240,260,280]
    mll = [12,30,45,60,75,100,125,150,175,200]

    print 'xxx'
    hf = fold(h_init,mth,mll)
    rp = RollProjector(mth,'x',nfold=(len(mth)-1))

    hp_init = rp.project(h_init)
    ep_init = rp.project(e_init)
#     hp = project(h_init,mth)

#     ROOT.gStyle.SetPalette(1)
#     hp.Draw()
#     ep.Draw('2')
#     ROOT.gPad.Print('hproj.pdf')

#     return
    noaxis = { 'labelsize':0, }
    xaxis = { 'labelsize':15, 'titlesize':20}
    yaxis = { 'titleoffset':60, 'titlesize':20, 'labeloffset':5, 'labelsize':15}

    p0 = Pad('main'  ,500,300,margins=(80,20,50, 0), xaxis=noaxis, yaxis=yaxis)
    p1 = Pad('diffs' ,500,100,margins=(80,20,0 , 0), xaxis=noaxis, yaxis=yaxis)
    p2 = Pad('priors',500,100,margins=(80,20,0 , 0), xaxis=noaxis, yaxis=yaxis)
    p3 = Pad('post'  ,500,150,margins=(80,20,0 ,60), xaxis= xaxis, yaxis=yaxis)

    c = Canvas(1,4)
    c[0,0] = p0
    c[0,1] = p1
    c[0,2] = p2
    c[0,3] = p3

    canv = c.makecanvas(n)

    hp_data = rp.project(h_data) 
    hp_init = rp.project(h_init) 
    hp_sig  = rp.project(h_sig)  
    hp_bkg  = rp.project(h_bkg)  
    ep_init = rp.project(e_init) 
    ep_sig  = rp.project(e_sig)  
    ep_bkg  = rp.project(e_bkg)  

    print 'Kolmogorov init-init',hp_init.KolmogorovTest(hp_init)
    print 'Kolmogorov init-sig',hp_init.KolmogorovTest(hp_sig)
    print 'Kolmogorov init-bkg',hp_init.KolmogorovTest(hp_bkg)
    print 'Kolmogorov sig-bkg',hp_sig.KolmogorovTest(hp_bkg)


    p0.cd()
    hs = ROOT.THStack('postfit','postfit;;entries')
    hframe  = hp_init.Clone('hp_frame')
    hframe.SetLineColor(ROOT.kBlack)
    hframe.SetXTitle('m^{ll,#slash{E}}_{T}')
    hframe.Reset()

    hs.Add( hp_init.Clone() )
    hs.Add( hp_sig.Clone()  )
    hs.Add( hp_bkg.Clone()  )
#     hs.Add( hp_data.Clone(),'e' )
    hs.Draw('nostack')
#     hs.SetYTitle('init errs')

    txt = Latex('Sample %s (syst \'%s\')' % (p,n) ,anchor=(80,40), align=('l','b'), style={'textsize':20})
    txt.draw()
#     ROOT.gPad.ls()

    leg = Legend(1,4,100,20,anchor=(360,60) )
    leg.addentry(hp_init,'l',title='init')
    leg.addentry(hp_sig,'l',title='s+b')
    leg.addentry(hp_bkg,'l',title='b')
    leg.draw()

    p1.cd()
    hp_sig_diff = hp_sig.Clone()
    hp_sig_diff.Add(hp_init,-1) 
    hp_sig_diff.Divide(hp_init)

    hp_bkg_diff = hp_bkg.Clone()
    hp_bkg_diff.Add(hp_init,-1) 
    hp_bkg_diff.Divide(hp_init)

    hsB = ROOT.THStack('diffs','diffs;;s-b diffs')
    hsB.Add( hp_sig_diff )
    hsB.Add( hp_bkg_diff )
    hsB.Draw('nostack')

    

    p2.cd()
    import ctypes
    er_init = ep_init.Clone()
    for i in xrange(er_init.GetN()):
        x,y = ctypes.c_double(),ctypes.c_double()
        er_init.GetPoint(i,x,y)
        er_init.SetPoint(i,x.value,0.)
        er_init.SetPointEYhigh(i,er_init.GetErrorYhigh(i)/y.value)
        er_init.SetPointEYlow (i,er_init.GetErrorYlow(i)/y.value)
    hf2 = hframe.Clone('frame2')
    hf2.SetMinimum(er_init.GetHistogram().GetMinimum())
    hf2.SetMaximum(er_init.GetHistogram().GetMaximum())
    hf2.SetYTitle('init errs')
    hf2.Draw()
    er_init.Draw('2')

    p3.cd()
    mg = ROOT.TMultiGraph()
    er_sig  = ep_sig.Clone()
    er_bkg  = ep_bkg.Clone()

    for i in xrange(er_init.GetN()):
        x,y = ctypes.c_double(),ctypes.c_double()

        er_sig.GetPoint(i,x,y)
        er_sig.SetPoint(i,x.value,0.)
        er_sig.SetPointEYhigh(i,er_sig.GetErrorYhigh(i)/y.value)
        er_sig.SetPointEYlow (i,er_sig.GetErrorYlow(i)/y.value)

        er_bkg.GetPoint(i,x,y)
        er_bkg.SetPoint(i,x.value,0.)
        er_bkg.SetPointEYhigh(i,er_bkg.GetErrorYhigh(i)/y.value)
        er_bkg.SetPointEYlow (i,er_bkg.GetErrorYlow(i)/y.value)

    mg.Add(er_sig,'2')
    mg.Add(er_bkg,'2')

    mg.Paint('a')
    p3.Clear()

    hf3 = hframe.Clone()
    hf3.SetMinimum(mg.GetHistogram().GetMinimum())
    hf3.SetMaximum(mg.GetHistogram().GetMaximum())
    hf3.SetYTitle('post errs')
    hf3.Draw()
    mg.Draw()

    c.applystyle()
    c.Modified()
    c.Update()
    canv.Print('zzz.pdf')
    

    # for each bin and process need to get the 3 histograms

    # for each nuisance, the errorgraph

#     for b in bins:
#         binfold = d.init.get(b)
#         
#         bkg = []
#         sig = []
#         for p in processes:
#             o = binfold.get('histo_'+p)

#             if p in signals:
#                 sig.append(o)
#             else:
#                 bkg.append(o)

#         break
    

#     print d.info.objs
#     print d.init.dirs



#     d.show()

#     dump = container()
#     info = container()
#     
#     info.nuisances = dumpfile.Get('info/nuisances')
#     info.nuisances = dumpfile.Get('info/processes')
#     info.nuisances = dumpfile.Get('info/signals')
#     info.nuisances = dumpfile.Get('info/nuisbybin')

#     dump.info = info

# trigger the tester
if __name__ == '__main__':

    sys.argv.append('-b')
    ROOT.gROOT.SetBatch()
    tester()
