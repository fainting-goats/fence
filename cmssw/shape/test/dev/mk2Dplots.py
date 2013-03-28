#!/usr/bin/env python

import array
import re
import ROOT
import hwwtools
import optparse
import HWWAnalysis.Misc.ROOTAndUtils as utils
import sys
import os
'''

===================================
accumPalette.resize(nAccum);

double rgbAccum[3][2] = {{0.70, 0.00}, {0.90, 0.10}, {0.90, 0.90}};

for(int i(0); i < nAccum; i++){
 double r((rgbAccum[0][1] - rgbAccum[0][0]) / (nAccum - 1) * i + rgbAccum[0][0]);
 double g((rgbAccum[1][1] - rgbAccum[1][0]) / (nAccum - 1) * i + rgbAccum[1][0]);
 double b((rgbAccum[2][1] - rgbAccum[2][0]) / (nAccum - 1) * i + rgbAccum[2][0]);
 new TColor(iCol + i, r, g, b);
 accumPalette[i] = iCol + i;
}
===================================

accumPalette is a vector<int>. nAccum is 50.

You would then do
gStyle->SetPalette(accumPalette.size(), &(accumPalette[0]));
'''

def monochrome():
    # pretty monochrome palette
    Number = 2;
    Red    = array.array( 'd', [ 0.70, 0.00 ] )
    Green  = array.array( 'd', [ 0.90, 0.10 ] )
    Blue   = array.array( 'd', [ 0.90, 0.90 ] )
    Length = array.array( 'd', [ 0.00, 1.00 ] )
    nb=50;
    ROOT.TColor.CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);


def prettypal():
    NRGBs = 5;
    NCont = 255;

    stops   = array.array('d',[ 0.00, 0.34, 0.61, 0.84, 1.00 ] )
    red     = array.array('d',[ 0.00, 0.00, 0.87, 1.00, 0.51 ] )
    green   = array.array('d',[ 0.00, 0.81, 1.00, 0.20, 0.00 ] )
    blue    = array.array('d',[ 0.51, 1.00, 0.12, 0.00, 0.00 ] )
    TColor.CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont)
    gStyle.SetNumberContours(NCont)


class ROOTFile(ROOT.TFile):
    def __init__(self, **kwargs):

        # run the initialiser
        ROOT.TFile.__init__(self,**kwargs)

    def Get(self,name):
        o = ROOT.TFile.Get(self,name)
        if not o.__nonzero__():
            raise RuntimeError('Object %s not found' % name)
        return o



def getFiles(dir):
    list = os.listdir(dir)
    list.sort()
    ## remove non root files
    cleaned = [l for l in list if '.root' in l]
    return cleaned

def getNominals(file):
    finder = utils.ObjFinder('TH1D')
    names = finder.findRecursive(file)
    nomRegex  = re.compile('^histo_([^_]+)_2d$')

    nominals = {}
    for name in names:
        if not nomRegex.match(name):
            continue
        h = getHist(file,name)
        nominals[name.replace('histo_','')] = h
                 
    return nominals

if __name__ == '__main__':
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i' , '--input'       , dest='inputdir'    , help='Input dir')
    parser.add_option('-o' , '--output'      , dest='outputdir'   , help='Output dir'                          , default='.' )

    hwwtools.addOptions(parser)
    hwwtools.loadOptDefaults(parser)

    (opt, args) = parser.parse_args()
    sys.argv.append('-b')

    if opt.inputdir is None:
        parser.error('No input file defined')

    indir = opt.inputdir
    outdir = opt.outputdir
    mass = opt.mass

    filenames = getFiles(indir)
    for file in filenames:
        path = indir+'/'+file
        if str(mass) not in file and mass > 0:
            continue
        print file


        ROOTFile('pippo')
        tf = ROOTFile(file)
        

    print 'Used options'
    print ', '.join([ '{0} = {1}'.format(a,b) for a,b in opt.__dict__.iteritems()])
