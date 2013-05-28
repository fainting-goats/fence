#!/usr/bin/env python

#----
class UnderOverTucker(object):
    def __init__(self):
        pass
    @staticmethod
    def _tuckbin(p,frombin,tobin):
        '''
        move the content of bin 'frombin' to bin 'tobin'.
        '''
        print frombin,tobin,':',p.GetAt(frombin),p.GetAt(tobin),'=',
        sumw2 = p.GetSumw2()
        # bincontent
        p.SetAt (p.GetAt (frombin)+p.GetAt (tobin),tobin)
        p.SetAt (0.,frombin)

        # sum of weights
        if p.GetSumw2N() != 0:
            sumw2.SetAt(sumw2.GetAt(frombin)+sumw2.GetAt(tobin),tobin)
            sumw2.SetAt(0.,frombin)

        print p.GetAt(tobin)

    def __call__(self,p):
        # TH1 histograms
        if p.GetDimension() == 1:
            nBins = h.GetNbinsX()
            entries = h.GetEntries()
            underFlow = h.GetBinContent(0)
            overFlow  = h.GetBinContent(nBins+1)
            bin1      = h.GetBinContent(1)
            binN      = h.GetBinContent(nBins)

            h.SetAt(0.,0)
            h.SetAt(underFlow+bin1,1)
            h.SetAt(overFlow+binN, nBins)
            h.SetAt(0.,nBins+1,)
        elif p.GetDimension() == 2:
            nx    = p.GetNbinsX()
            ny    = p.GetNbinsY()
            entries = p.GetEntries()

            print 'tucking',p.GetName(),nx,ny
            print 'xsscan'
            # x-scan to tuck in the over/underflows
            for i in xrange(nx+2):
                bu = p.GetBin(i,0)
                b1 = p.GetBin(i,1)
                self._tuckbin(p,bu,b1)
            for i in xrange(nx+2):
                bn = p.GetBin(i,ny)
                bo = p.GetBin(i,ny+1)
                self._tuckbin(p,bo,bn)

            print 'yscan'
            # y-scan to tuck in the over/underflows
            for j in xrange(ny+2):
                bu = p.GetBin(0   ,j)
                b1 = p.GetBin(1   ,j)
                self._tuckbin(p,bu,b1)

            for j in xrange(ny+2):
                bn = p.GetBin(ny  ,j)
                bo = p.GetBin(ny+1,j)

                self._tuckbin(p,bo,bn)

        elif p.GetDimension() == 3:
            pass
        else:
            raise ValueError('Don\' know about histograms with dimension '+str(p.GetDimension()))


uo = UnderOverTucker()

import ROOT

h2 = ROOT.TH2D('aa','aa',3,0,3,3,0,3)

h2.Fill(0.,0.)
h2.Fill(-1,-1)

print '-'*80
h2.Print('all')
print '-'*80
uo(h2)
h2.Print('all')
print '-'*80
