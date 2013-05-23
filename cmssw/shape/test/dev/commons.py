class AlienDict(dict):
    """Implementation of perl's autovivification feature."""
    def __init__(self,*args, **kwargs):
        # init the dict
        super(self.__class__,self).__init__(self, *args, **kwargs)
        self._lock = False

    def lock(self):
        self._lock = True
        for a in self.itervalues():
            if type(a) == type(self):
                a.lock()

    def unlock(self):
        self._lock = False
        for a in self.itervalues():
            if type(a) == type(self):
                a.unlock()

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            if self._lock:
                raise
            else:
                value = self[item] = type(self)()
                return value

from ginger.plotter import H1RatioPlotter
import ROOT

def diffembed(th2, options='',markers=None, colors=None):

    thepad = ROOT.gPad.func()

    from HWWAnalysis.Misc.ROOTAndUtils import TH1AddDirSentry, TH1Sumw2Sentry
    sentrys2 = TH1Sumw2Sentry()
    sentry = TH1AddDirSentry()

    dummy1 = ROOT.TH1D('dummyA','dummy',10,0,10)
    dummy2 = ROOT.TH1D('dummyB','dummy',10,0,10)
    map(dummy1.Fill,[1,1,2,7,4,4,5])
    map(dummy2.Fill,[1,1,2,7,4,6,5])

    dummy1.Print("V")

    h1 = H1RatioPlotter(plotratio=False)
    h1.set(dummy1,dummy2)
    h1._style.scale(0.5)
    #for attr in h1._style.__dict__:
        #x = getattr(h1,attr)
        #if isinstance(x,int):
            #setattr(h1,attr,x/2)

    c = h1.plot()


    thepad.Clear()
    thepad.cd()
    ROOT.gROOT.SetSelectedPad(thepad)
    c.DrawClonePad()
    thepad.ls()
    thepad.cd()

