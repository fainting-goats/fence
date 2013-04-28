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
