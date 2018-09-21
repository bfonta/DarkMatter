import os, warnings
import pynbody as pn
import pynbody.units as u
from dmprofile.fit import nfw_fit
from dmprofile.move import centering_com

class Halos():
    """
    Halos management
    """
    def __init__(self, filename, N):
        self.filename = filename
        s = pn.load(os.path.join(self.filename))
        s['eps'] = 200.*pn.units.pc
        s.physical_units() #concentration will be wrong if this line is removed

        if N>len(s): 
            warnings.warn("The specified halo number is larger than the number of halos in the simulation box.")
        self._N = N
        self._halos = s.halos()[:self._N]

    def _check_subhalo_exists(self, halo_idx):
        if self._halos[halo_idx].properties['NsubPerHalo'] == 0:
            return False
        else:
            return True

    def get_halo(self, halo_idx):
        if halo_idx<=self._N: 
            return self._halos[halo_idx]
        else:
            raise IndexError("The halo you requested was not stored.")
            return self._halos[-1]

    def get_subhalo(self, halo_idx, sub_idx=0):
        """
        By default returns the largest subhalo in the specified halo.
        """
        if not self._check_subhalo_exists(halo_idx):
            raise RuntimeError('This halo has no subhalos!')
        return self._halos[halo_idx].sub[sub_idx]

    def get_subhalos(self, sub_idx=0):
        """
        By default returns all the largest subhalo of each halo.
        """
        return [self.get_subhalo(i, sub_idx) for i in range(self._N)]

    def get_halos(self, n=0):
        """Returns the first 'n' halos."""
        if n==0: n=self._N
        if n>self._N:
            raise("You are asking for too many halos.")
        return self._halos[:n]

    def concentration_200(self, idx, sub='False', sub_idx=0):
        """
        Calculates the concentration at r200 
        """
        _current_halo = self.get_halo(idx) if sub=='False' else self.get_subhalo(idx, sub_idx)
        with centering_com(_current_halo):
            #we always use the main halo for the r200 distance, even for subhalos
            _r200 = self.get_halo(idx).properties['Halo_R_Crit200'].in_units('kpc')
            _profile = pn.analysis.profile.Profile(_current_halo.dm, 
                                                   ndim=3, nbins=250, min=0.5, max=35)
            _rs = nfw_fit(_profile)[1]
            print(_r200, _rs)
        return _r200/_rs
      

    def concentrations_200(self, sub='False', sub_idx=0, n=None):
        """
        Calculates the first 'n' concentrations at r200. Returns a list.
        """
        if n==None: 
            n = self._N
        return [self.concentration_200(i, sub, sub_idx) for i in range(n)]
