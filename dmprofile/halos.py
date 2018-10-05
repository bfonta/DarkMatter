#MISSING: Add methods to retrieve the halos and subhalos stored in backup

import os, warnings
import numpy as np
import pynbody as pn
import pynbody.units as u
from pynbody.analysis.halo import center_of_mass as com
from dmprofile.fit import nfw_fit
from dmprofile.utilities import intersect
from dmprofile.move import centering_com

class Halos():
    """
    Halos management
    """
    def __init__(self, filename, N=1):
        self.filename = filename
        _s = pn.load(os.path.join(self.filename))
        _s['eps'] = 200.*pn.units.pc
        _s.physical_units() #concentration will be wrong if this line is removed
        if N>len(_s): 
            warnings.warn("The specified halo number is larger than the number of halos in the simulation box.")
        self._N = N
        self._halos = _s.halos()[:self._N]
        self._halos_bak = self._halos #halos backup

        #backup for subhalos (id, subhalo): 
        #the id is helpful because each halo can have more than one subhalo
        #'0' is the id for the main subhalo, '1' is for the second largest, etc
        self._subhalos = [[] for _ in range(self._N)] #one row per halo in descending size order
        for _i in range(self._N):
            length = int(self._halos[_i].properties['NsubPerHalo'])#number subhalos in a halo
            #columns represent the subhalos in descending size order
            self._subhalos[_i].extend((_isub, self._halos[_i].sub[_isub]) for _isub in range(length))
        self._subhalos_bak = self._subhalos
    
    def _get_r200(self, idx):
        return self._halos[idx].properties['Halo_R_Crit200'].in_units('kpc')

    def _check_backup(self, idx, sub_idx=-1):
        if sub_idx < 0:
            return self._halos[idx]==self._halos_bak[idx]
        else:
            return self._subhalos[idx][sub_idx]==self._subhalos_bak[idx][sub_idx]

    def restore(self, idx, sub_idx=-1):
        if idx >= self._N: raise ValueError('The halo index is too large.')
        if type(idx) != int or type(sub_idx) != int: raise TypeError('The indexes must be integers.')
        if sub_idx<0:
            self._halos[idx], self._halos_bak[idx] = self._halos_bak[idx], self._halos[idx]
        else:
            self._subhalos[idx][sub_idx], self._subhalos_bak[idx][sub_idx] = self._subhalos_bak[idx][sub_idx], self._subhalos[idx][sub_idx]

    def filter_(self, halo_idxs, sub_idx=-1, filter_str='Sphere_1.2', option=None):
        """
        halo_idxs: indexes of the halos to filter; if more than one, it must be a list
        option: None, 'half1', 'half2' or 'all'
        filters:
        Sphere_radius, where radius is a number relative to r200. Example: Sphere_1
        BandPass_prop_min_max. Example: BandPass_y_1 kpc_2 kpc. It must have units.
        The inputs for each filter are separated using the underscore '_' carachter
        More filters can be added to 'filter_dict' in the same fashion with variable number of arguments
        """
        if type(halo_idxs) != int and type(halo_idxs) != list:
            raise TypeError('The introduced indexes have the wrong format')
        if type(halo_idxs) is int:
            halo_idxs = [halo_idxs]

        filter_dict = {'Sphere': lambda _r,_factors,_halo: 
                       pn.filt.Sphere(float(_factors[0])*_r, com(_halo)),
                       'BandPass': lambda _r,_factors,_halo: 
                       pn.filt.BandPass(_factors[0], _factors[1], _factors[2])}
        for item in filter_dict.keys():
            if item in filter_str:
                break
            else:
                if item==[*filter_dict.keys()][-1]:
                    raise ValueError('The specified filtering option is not supported.')
        _split = filter_str.split('_')
        _split_str = _split[0]
        _split_values = {i: s for i,s in enumerate(_split[1:])}
        print("---FILTER--- Type:", _split_str, "-- Parameters:", *_split_values.values())
        _filter_func = filter_dict[_split_str]

        if option=='all': halo_idxs=range(self._N) 
        elif option=='half1': halo_idxs=range(int(self._N/2))
        elif option=='half2': halo_idxs=range(int(self._N/2),self._N)
        elif option==None: pass
        else: raise ValueError('That option is not supported')

        for i in halo_idxs:
            _r200 = self._get_r200(i)
            
            if sub_idx < 0: #filter the halo
                if not self._check_backup(i):
                    print(self._check_backup(i))
                    raise warnings.warn('This halo backup is already being used!')
                self._halos_bak[i] = self._halos[i]
                with centering_com(self._halos[i]):
                    self._halos[i] = self._halos[i][_filter_func(_r200, _split_values, 
                                                                 self._halos[i])]
            else: #filter the halo but centered in the subhalo
                if not self._check_backup(i):
                    raise warnings.warn('This subhalo backup is already being used!')
                if self._subhalos[i][sub_idx][0] != sub_idx:
                    raise RuntimeError('Different subhalos are being mixed!')
                self._subhalos_bak[i][sub_idx] = self._subhalos[i][sub_idx]
                with centering_com(self._subhalos[i][sub_idx][1]):
                    _subhalo = self._halos[i][_filter_func(_r200, _split_values, 
                                                           self._subhalos[i][sub_idx][1])]
                    self._subhalos[i][sub_idx] = (sub_idx,_subhalo)


    def is_relaxed(self, halo_idx, sub_idx=-1):
        if sub_idx < 0:
            if len(self._subhalos[halo_idx])==0:
                raise RuntimeError('This halo has no subhalos!')
            _h = self._halos[halo_idx]
        else:
            _h = self._halos[halo_idx].sub[sub_idx]
        with centering_com(self._halos[halo_idx]): #centering the main halo
            _com = pn.analysis.halo.center_of_mass(_h)
            _r200 = self._get_r200(halo_idx)
            _dist = np.sqrt( np.power(_h[0]['pos']-_com,2).sum() )
            #The minimum potential point is obtained with:
        #halo[halo_index].subhalo[subhalo_index][0]['pos']
        return bool(_dist < 0.07*_r200)

    def is_resolved(self, halo_idx, intersect_value=1.):
        with centering_com(self._halos[halo_idx]): #centering the main halo
            _r200 = self._get_r200(halo_idx)
            _prof = self.get_profile(halo_idx, 'dm', bins=(2.4,25,30), bin_type='log', normalize=False)
            _relax = self.relaxation(halo_idx, _prof)
            _inter = intersect(_prof['rbins'], _relax, intersect_value=intersect_value)
            if _inter[1]==True:            
                print(_inter[0], " ---- ", _r200*np.power(10,-1.25))
                return _inter[0] < _r200*np.power(10,-1.25)
            else:
                raise RuntimeError('No intersection was found for this halo. The resolution criteria cannot thus be checked.')
                return False

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
        return self._subhalos[halo_idx][sub_idx][1]

    def get_subhalos(self, sub_idx=-1):
        """
        Either returns all subhalos or just the subhalos with a specific index for all halos.
        The subhalos will be wrapped in a tuple with their respective indexes.
        """
        if sub_idx<0:
            return self._subhalos
        else:
            return [self._subhalos[i][sub_idx] for i in range(self._N)]

    def get_halos(self):
        return self._halos

    def concentration_200(self, idx, sub_idx=-1):
        """
        Calculates the concentration at r200 
        """
        if sub_idx<0:
            _current_halo = self._halos[idx] 
        else:
            _current_halo = self._subhalos[idx][sub_idx][1]

        with centering_com(self._halos[idx]):
            #we always use the main halo for the r200 distance, even for subhalos
            _r200 = _get_r200(idx)
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

    def get_profile(self, idx, sub_idx=-1, component='dm', 
                    bins=(1,45,30), bin_type='linear', normalize=False):
        """
        Obtain the profile of a single halo.
        Arguments:
        'idx': which halo to get the profile from
        'sub_idx': which subhalo to get the profile from
        'component': 'dm', 'stars', 'gas' or 'all'
        'bins': tuple with (start,stop,nbins)
        'bin_type': either 'linear' or 'log'
        'normalize': 'None' for no normalization and 'r200' for normalization
        'centering': whether to centre the halo before obtaining the profile or not
        """
        _h = self._halos[idx] if sub_idx<0 else self._subhalos[idx][sub_idx][1]
        with centering_com(_h):
            components = {'dm':_h.dm, 'stars':_h.s, 'gas':_h.g, 'all':_h}
            if component not in components:
                raise ValueError('The specified component does not exist.')
            if type(bins) is not tuple or len(bins) != 3:
                raise TypeError('Please insert bins with the following pattern: (start, stop, nbins).')
        
            if normalize:
                r200 = self._get_r200(idx)
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

    def relaxation(self, idx, profile=None, component='dm', bins=(4,50,30), bin_type='linear'):
        """
        Creates a t_relax(r) / t_circ(r200) array for the provided halo.
        If you do not provide a profile, the method creates one for you. 
        The profile is needed to obtain the density as a function of r.
        """
        _h = self._halos[idx]
        _rho_crit_z = pn.array.SimArray(pn.analysis.cosmology.rho_crit(_h, unit="Msol kpc**-3"), 
                                        "Msol kpc**-3")
        with centering_com(_h):
            if profile==None:
                _prof = self.get_profile(idx, component=component, bins=bins, 
                                         bin_type=bin_type, normalize=False)
                _density = _prof['density']
                _rbins = _prof['rbins']
                print(_density)
                print(_rbins)
            else:
                _density = profile['density']
                _rbins= profile['rbins']

            _nparts = pn.array.SimArray([len(_h[_h['r']<_rbins[i]]) for i in range(len(_rbins))],'')

        if 0. in _density:
            print(_density)
            raise ValueError('At least one of the values of the density is zero. The relaxation time cannot be calculated when this happens.')
        return (np.sqrt(200)*_nparts)/(8*np.log(_nparts))*np.sqrt(_rho_crit_z/_density)
