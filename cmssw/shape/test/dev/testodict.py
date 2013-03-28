#!/usr/bin/env python

from collections import OrderedDict
from ginger.analysis import Cut





#______________________________________________________________________________
class CutFlow(OrderedDict):

    #---
    def __init__(self, init_val=()):
        OrderedDict.__init__(self,()) 
        self._import(init_val)
    
    #---
    def _import(self,l):
        if isinstance(l,list):
            for item in l:
                if ( isinstance(item,tuple) and isinstance(item,list) ) and len(item)==2:
                    self[item[0]] = item[1]
                elif isinstance(item,str):
                    self[item] = item
                else:
                    raise ValueError('CutFlow: list entry not supported: '+item.__class__.__name__)
        elif isinstance(l,OrderedDict):
            for key,val in l.iteritems():
                self[key] = val
        elif isinstance(l,tuple):
            pass
        else:
            raise ValueError('CutFlow: init value not supported: '+l.__class__.__name__)


    #---
    def __setitem__(self, key, val):
        if isinstance(val,str):
            val = Cut(val)    
        elif not isinstance(val,Cut):
            raise TypeError('CutFlow need CutSteps')
        
        val.name = key

        OrderedDict.__setitem__(self,key,val)

    #---
    def __getitem__(self, key):
        """
        Allows slicing. Returns an OrderedDict if you slice.
        >>> b = OrderedDict([(7, 0), (6, 1), (5, 2), (4, 3), (3, 4), (2, 5), (1, 6)])
        >>> b[::-1]
        OrderedDict([(1, 6), (2, 5), (3, 4), (4, 3), (5, 2), (6, 1), (7, 0)])
        >>> b[2:5]
        OrderedDict([(5, 2), (4, 3), (3, 4)])
        >>> type(b[2:4])
        <class '__main__.OrderedDict'>
        """
        item = OrderedDict.__getitem__(self,key)
        if isinstance(item, OrderedDict):
            return CutFlow(item)
        else:
            return item

    #---
    def __repr__(self):
        return '%s([%s])' % (self.__class__.__name__, ', '.join( 
            ['(%r, %r)' % (key, self[key].cut) for key in self.iterkeys()]))
    
    __str__ = __repr__


    #---
    def at(self,index):
        return self.__getitem__(self._sequence[index])

    #---
    def collapse(self,name):
        cut = self.string()
        self.clear()
        self.__setitem__(name,cut)

    #---
    def insert(self, index, key, value):

        if isinstance(value,str):
            value = Cut(value)

        OrderedDict.insert(self, index, key, value)

    #---
    def rename(self, old_key, new_key):

        OrderedDict.rename(self,old_key,new_key)

        self.__getitem__(new_key).name = new_key

    #---
    def string(self):
        return ' && '.join( [ '(%s)' % step.cut for step in self.itervalues() ] ) 

    #---
    def list(self):
        return [ (step.name,step.cut) for step in self.itervalues() ]

    #---
    def rawlist(self):
        return [ step.cut for step in self.itervalues() ]


c = CutFlow()

c['a'] = '1'
c['b'] = '2'

