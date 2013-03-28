#!/bin/env python

from ginger.tree import TreeWorker,TreeView
from ginger.analysis import CutFlow,TreeAnalyser

import HWWAnalysis.Misc.odict as odict
import hwwsamples2g
import uuid
import logging
import itertools
import pdb





def test():
    latinos= hwwsamples2g.samples(125,'8TeV','Data2012','SM','topestimate')
    s =  latinos['Top'].makeSample('/shome/thea/HWW/work/dds/trees/top')
    s.preselection = 'njet > 2'

    w = TreeWorker.fromsample(s)

    all =  w.entries()
    print 'all:',all
#     print '7',w.entries('njet <= 7')
#     print '6',w.entries('njet <= 6')
#     print '5',w.entries('njet <= 5')
#     print '4',w.entries('njet <= 4')
#     print '3',w.entries('njet <= 3')

#     r1 = w.entries('njet <= 4')
#     r2 = w.entries('njet > 4')

#     print 'R1',r1
#     print 'R2',r2
#     print r1+r2,all


#     m7 = TreeView(w,'njet <= 7')
#     print 'm7',m7.entries(),m7.cut()
#     m6 = m7.spawn('njet <= 6')
#     print 'm6',m6.entries(),m6.cut()
#     m5 = m6.spawn('njet <= 5')
#     print 'm5',m5.entries(),m5.cut()
#     m4 = m5.spawn('njet <= 4')
#     print 'm4',m4.entries(),m4.cut()

#     k1 = TreeView(w,'njet <= 4')
#     k2 = TreeView(w,'njet  > 4')

#     print 'k1',k1.entries()
#     print 'k2',k2.entries()

#     x1 = k1.entries('ptll > 60')
#     x2 = k2.entries('ptll > 60')

#     print x1,'+',x2

#     print x1+x2,w.entries('ptll > 60')

#     x1 = k1.yields('ptll > 60')
#     x2 = k2.yields('ptll > 60')

#     print x1,'+',x2

#     print x1+x2,w.yields('ptll > 60')



#     print misteries
#     print odict.OrderedDict([ (n,m.yields()) for n,m in misteries.iteritems() ])
#     xyz =  odict.OrderedDict([ (n,m.plot('ziogano%i' % i,'mll')) for i,(n,m) in enumerate(misteries.iteritems()) ])
#     for w in xyz.itervalues(): w.Print()

    print 'nosel:',w.GetEntries()
    print 'elist',w._elist.GetN()
    last = TreeView(w)
    print 'mistery:',last.entries()
    print 'xcheck:',w.entries(last.cut)

    flow = CutFlow()
    flow['n4'] = 'ptll > 45'
    flow['n5'] = 'ptll > 50'
    flow['n6'] = 'ptll > 55'
    flow['n7'] = 'ptll > 60'
    
    #---
    def grow( cutflow, views ):
        # explect cutflow to be longer than 
    
        nv = len(views)
        nc = len(cutflow)
        if nv == nc : return
        elif nv > nc : raise ValueError('WTF!')

        last = views.values()[-1] if nv > 0 else TreeView(w) # <<-- w is to be replaced)

        newcuts = cutflow[nv:]
        print newcuts
        for i,(n,c) in enumerate(newcuts.iteritems()):
            m = last.spawn(c,'elist%d' % (i+nv) )
            print m.cut
            last = m  
            views[n] = m


    misteries = odict.OrderedDict()
    grow(flow,misteries)

    for n,m in misteries.iteritems():
        print n,m.entries(),m.yields()

    del flow['n4']
    flow['n8'] = 'ptll > 70'

    
    def purge( cutflow, views ):
        # check the keys
        nv = len(views)

        print 'cutflow',cutflow.keys()
        print 'view',views.keys()

        if len(cutflow) == 0:
            # the new list is empty! what do we do?
            numok = 0
#             mlast = TreeView(w)
        else:
            k = 0
            matches = False
#             lastmatch = None
            # find the first mismatching cut
            for k, ( n,m ) in enumerate(itertools.izip(cutflow.iterkeys(),views.iterkeys())):
                matches =  ( n == m ) and ( str(cutflow[n]) == str(views[m].cut) )
                if not matches: break
#                 lastmatch = n
            if k==0 and not matches:
                # disagreement at the first match
                numok = 0
#                 mlast = TreeView(w)
            elif not matches:
                numok = k
#                 mlast = views[lastmatch]
            else:
                numok = k+1
#                 mlast = views[n] 

#         print numok,mlast
        print views.keys()[numok:nv]

        # purge the rest of elists
        if numok < nv:
            for n in views.keys()[numok:nv]:
                logging.debug('Deleting %s',n)
                del views[n]

        for i,(n,l) in enumerate(views.iteritems()):
            logging.debug('- %d %s,%d', i,n,l.entries())
    
        return numok

    n = purge(flow, misteries)
    
    print n, len(flow)


    grow(flow,misteries)

    for n,m in misteries.iteritems():
        print n,m.entries(),m.yields()






    # setting the last common
#     self._chain.SetEntryList(elast)

#     newcuts = cuts[numok:]
#     # create the missing entrllists 
#     newlist = self._makeentrylists(newcuts,numok)
#     elists.update(newlist)

#     # restore the preselection
#     self._chain.SetEntryList(self._elist if self._elist else 0x0)


def testanalyser():

    lumi = 19.601
    latinos= hwwsamples2g.samples(125,'8TeV','Data2012','SM','topestimate')
    s =  latinos['Top'].makeSample('/shome/thea/HWW/work/dds/trees/top')
    s.preselection = 'njet > 2'

    flow = CutFlow()
    flow['n4'] = 'ptll > 45'
    flow['n5'] = 'ptll > 50'
#     flow['n6'] = 'ptll > 55'
#     flow['n7'] = 'ptll > 60'

    a = TreeAnalyser(s,flow)
    a.lumi = lumi
    a.bufferentries()

    print '# before','-'*50
    for i,(n,l) in enumerate(a._entrylists.iteritems()):
        print i,n,l.entries(),l.yields()

#     a.remove('n4')
#     a.append('n8','ptll > 70')

#     a.bufferentries()
    
#     print '# after ','-'*50
#     for i,(n,l) in enumerate(a._entrylists.iteritems()):
#         print i,yields(extra)n,l.entries(),l.yields()
    print a._entrylists

    print a.yieldsflow()
    print a.plotsflow('ziofatto','mll',bins=(100,0,500))
    

if __name__ == '__main__':
    import optparse
    import logging
    import bdb,sys
    parser = optparse.OptionParser()
    parser.add_option('-d', '--debug'    , dest='debug'       , help='Debug level'                           , default=0      , action='count' )
    (opt, args) = parser.parse_args()

    if not opt.debug:
        pass
    elif opt.debug == 2:
        print 'Logging level set to DEBUG (%d)' % opt.debug
        logging.basicConfig(level=logging.DEBUG)
    elif opt.debug == 1:
        print 'Logging level set to INFO (%d)' % opt.debug
        logging.basicConfig(level=logging.INFO)

    try:
        testanalyser()
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
    #         print formatted_lines[0]
        print '   ',formatted_lines[-1]
        print '-'*80
        traceback.print_exc()
    #         print "*** format_exception:"
    #         print repr(traceback.format_exception(exc_type, exc_value,
    #                                               exc_traceback))
        print '-'*80
    finally:
        print 'GOODBYE!'


