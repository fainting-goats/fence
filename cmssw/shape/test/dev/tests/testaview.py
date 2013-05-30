
def run( path, globals=None, locals=None ):
    import os
    if path.startswith('./'):
        d = os.path.dirname(__file__)
        path = d+path[1:]
    elif path.startswith('../'):
        d = os.path.dirname(__file__)
        path = d+'/'+path[:]
    execfile(path, globals, locals)

#run /shome/thea/HWW/git/fainting-goats/cmssw/shape/test/dev/loadsamples.py

run('../loadsamples.py',globals(),locals())

from ginger.filters import UnderOverTucker
atW = analysers['tW']
uo = UnderOverTucker()
#atW.filters.append(uo)
print atW.filters
print atW._aview._filters
av2 = atW._aview.spawn('mll > 30')
av2.entries()
atW._aview.entries()
av2
h1 = av2.plot('xx','mll',bins=(1,30,40))
h2 = av2.plot('xx','mll',bins=(1,30,40))
