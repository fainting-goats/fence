#!/usr/bin/env python

import ROOT
import hwwinfo
import optparse
import os.path

usage = 'usage: %prog mass cut'
parser = optparse.OptionParser(usage)
parser.add_option('-e','--extra',dest='extra',default=None)
parser.add_option('-s','--scan',dest='scan',default=None)
parser.add_option('-o','--opts',dest='opts',default='')
parser.add_option('--bdt',dest='bdtnjet',default=None)
parser.add_option('-l','--log',dest='log',default=None)

(opt, args) = parser.parse_args()

if not args or len(args) < 2:
    parser.error('Il cut, deficiente!')

mass = int(args[0])
cut = args[1]


processes = hwwinfo.samples(mass)

treepath = '/shome/thea/HWW/ShapeAnalysis/trees/latino_skim/all/'
bdtpathtmpl = '/shome/thea/HWW/ShapeAnalysis/trees/bdt_skim/ntupleMVA_MH{mass}_njet{njet}/all/'

selections = hwwinfo.massSelections(mass)
sel = selections[cut]
if opt.extra:
    sel = '('+sel+')'+opt.extra

print '-'*80
print selections[cut]
print '-'*80
chains = {}
for process in sorted(processes):
    if process != 'Data': continue
    latchain = ROOT.TChain('latino')
    for f in processes[process]:
        print 'latino::adding',f
        latchain.Add(treepath+f)
    if opt.bdt:
        bdtchain = ROOT.TChain('latinobdt')
        bdtpath = bdtpathtmpl.format(mass=mass,njet=opt.bdtnjet)
        for f in processes[process]:
            print 'bdt::adding',os.path.join(bdtpath,f)
            bdtchain.Add( bdtpath+f)
        latchain.AddFriend(bdtchain)
   

    if opt.scan:
        if opt.log:
            print 'Logging to',opt.log
            latchain.GetPlayer().SetScanRedirect(True)
            latchain.GetPlayer().SetScanFileName(opt.log)
        entries = latchain.Scan(opt.scan,sel,opt.opts)
    else:
        entries = latchain.Draw('trigger>>htmp',sel,'goff'+opt.opts)
        htmp = ROOT.gDirectory.Get('htmp')
        print '{0:<15} {1:<10} {2:<10}'.format(process,entries,htmp.Integral())
print '-'*80

