import os, warnings
import numpy as np
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

    def get_halos(self):
        return self._halos

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

    def get_profile(self, idx, component, bins, bin_type='linear', normalize=False, centering=True):
        """
        Obtain the profile of a single halo.
        Arguments:
        'idx': which halo to get the profile from
        'component': 'dm', 'stars', 'gas' or 'all'
        'bins': tuple with (start,stop,nbins)
        'bin_type': either 'linear' or 'log'
        'normalize': 'None' for no normalization and 'r200' for normalization
        'centering': whether to centre the halo before obtaining the profile or not
        """
        _h = self.get_halo(idx)
        if centering:
            centering_com(_h)
        
        components = {'dm':_h.dm, 'stars':_h.s, 'gas':_h.g, 'all':_h}
        if component not in components:
            raise ValueError('The specified component does not exist.')
        if type(bins) is not tuple or len(bins) != 3:
            raise TypeError('Please insert bins with the following pattern: (start, stop, nbins).')
        
        if normalize:
            r200 = _h.properties['Halo_R_Crit200'].in_units('kpc')
            calc_x = lambda x: x['r']/r200
        else:
            calc_x = lambda x: x['r']
            
        if bin_type=='linear':
            return pn.analysis.profile.Profile(components[component],
                                               calc_x = calc_x,
                                               ndim=3, 
                                               bins=np.linspace(bins[0],
                                                                bins[1],
                                                                bins[2]))
        elif bin_type=='log':
            return pn.analysis.profile.Profile(components[component], 
                                               calc_x = calc_x,
                                               ndim=3, 
                                               bins=np.logspace(np.log10(bins[0]),
                                                                np.log10(bins[1]),
                                                                bins[2]))
        else:
            raise ValueError('The specified bin type does not exist.')

    def relaxation(self, idx, profile=None, component='dm', 
                   bins=(4,50,30), bin_type='linear', centering=True):
        """
        Creates a t_relax(r) / t_circ(r200) array for the provided halo.
        If you do not provide a profile, the method creates one for you. 
        The profile is needed to obtain the density as a function of r.
        """
        _h = self.get_halo(idx)
        _rho_crit_z = pn.array.SimArray(pn.analysis.cosmology.rho_crit(_h, unit="Msol kpc**-3"), 
                                       "Msol kpc**-3")
        if profile==None:
            _prof = self.get_profile(idx, component=component, bins=bins, bin_type=bin_type,
                                     normalize=False, centering=centering)
            _density = _prof['density']
            _rbins = _prof['rbins']
            print(_density)
            print(_rbins)
        else:
            _density = profile['density']
            _rbins= profile['rbins']

        _nparts = pn.array.SimArray([len(_h[_h['r']<_rbins[i]]) for i in range(len(_rbins))],'')
        if 0. in _density:
            raise ValueError('At least one of the values of the density is zero. The relaxation time cannot be calculated when this happens.')
        return (np.sqrt(200)*_nparts)/(8*np.log(_nparts))*np.sqrt(_rho_crit_z/_density)
