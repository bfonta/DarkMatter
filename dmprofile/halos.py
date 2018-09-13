import os, warnings
import pynbody as pn

class Halos():
    """
    Halos management
    """
    def __init__(self, filename, N):
        self.filename = filename
        s = pn.load(os.path.join(self.filename))
        s['eps'] = 200.*pn.units.pc
        s.physical_units()

        if N>len(s): 
            warnings.warn("The specified halo number is larger than the number of halos in the simulation box.")
        self.N = N
        self.halos = s.halos()[:N]

    def get_halo(self, n):
        if n<=self.N: 
            return self.halos[n]
        else:
            raise IndexError
            return self.halos[-1]

    def get_halos(self):
        return self.halos
