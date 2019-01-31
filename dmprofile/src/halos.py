import os, warnings
import numpy as np
import copy
import pynbody as pn
import pynbody.units as u
from pynbody.analysis.halo import center_of_mass as com

from src.fit import nfw_fit
from src.utilities import intersect, memory_usage_resource as mem
from src.move import centering_com

class Halos():
    """
    Halos management
    """
    __slots__ = ['_sim', '_N', '_halos', '_halos_bak', 
                 '_subhalos', '_subhalos_bak', '_vars', '_mbp']

    def __init__(self, halos_filename, sim_filename='', N=None, min_size=300):
        """
        Tasks:
        1. Loads the halos and full simulation. Make sure they correspond to the same particles.
        2. Creates backup copies of all halos and subhalos

        Note: The full simulaton includes the particles that were not assigned to any halo
              by the Friends of Friends algorithm.
        """
        if type(min_size)!=int:
            raise ValueError("The 'min_size' parameter has to be an integer")
        if min_size<1 or (type(N)==int and N<1):
            raise ValueError('The number of selected halos has to be larger than zero.')

        _s_halos = pn.load(halos_filename)
        _s_halos['eps'] = 200.*pn.units.pc
        _s_halos.physical_units()
        
        self._sim = None
        if sim_filename!='':
            self._sim = pn.load(sim_filename)
            self._sim['eps'] = 200.*pn.units.pc
            self._sim.physical_units()

        if N!=None:
            if N>len(_s_halos): 
                warnings.warn("The specified halo number is larger than the number of halos in the simulation box.")
        _halos_tot = _s_halos.halos()
        self._halos = [item for item in _halos_tot if len(item)>=min_size]
        
        if N!=None:
            self._halos = self._halos[:N]
        self._N = len(self._halos)
        self._halos_bak = self._halos.copy() #halos backup

        #backup for subhalos (id, subhalo): 
        #the id is helpful because each halo can have more than one subhalo
        #'0' is the id for the main subhalo, '1' is for the second largest, etc
        self._subhalos = [[] for _ in range(self._N)] #one row per halo in descending size order
        for _i in range(self._N):
            length = int(self._halos[_i].properties['NsubPerHalo'])#number subhalos in a halo
            #columns represent the subhalos in descending size order
            self._subhalos[_i].extend((_isub, self._halos[_i].sub[_isub]) for _isub in range(length))
        self._subhalos_bak = self._subhalos.copy()

        ###################################################################################
        #most bound particle positions
        """
        _Sub_to_FOF = _halos_tot._fof_properties['FirstSubOfHalo']
        _NsubPerHalo = _halos_tot._fof_properties['NsubPerHalo']
        _FoFParentHalo = _halos_tot._sub_properties['SubParentHalo']
        _NFoF_Tot = len(_NsubPerHalo)
        _survivingfof = (_NsubPerHalo > 0)
        _FOF_with_MBP = np.arange(_NFoF_Tot)[_survivingfof]
        _Sub_to_FOF_list = list(np.array(_Sub_to_FOF[_survivingfof]).astype(int))

        _SubMostBoundID = _halos_tot._sub_properties['SubMostBoundID']
        _FOFMostBoundID = list(np.array(_SubMostBoundID[_Sub_to_FOF_list]).astype(int)) 

        #Number of FoF Halos with a surviving subhalo that has a most bound particle
        _NFoF_MBP = len(_FOFMostBoundID)
        
        _FOF_MBP = pn.array.SimArray( -1. * np.ones([_NFoF_Tot,3]), _halos_tot[0].s['pos'].units)

        _matched = 0
        for nfof, ID in zip(_FOF_with_MBP, _FOFMostBoundID): 
            ID = np.long(ID)    
            print("KEYS:", _halos_tot[nfof].s.all_keys())
            for arr in (_halos_tot[nfof].s, _halos_tot[nfof].g, _halos_tot[nfof].dm):
                if len(_halos_tot[nfof].s)>0 and len(_halos_tot[nfof].g)>0:
                    indx = np.argwhere(arr['ParticleIDs'] == ID)
                    if len(indx) == 1:
                        _FOF_MBP[nfof,:] = arr['pos'][indx[0]]
                        _matched += 1
                        break
                else:
                    continue
        if _matched != _NFoF_MBP:
            print("We only matched", _matched, "of the", _NFoF_MBP, "FoF halos we should have matches for, out of the Total FoF set of", _NFoF_Tot)
        
        self._mbp = _FOF_MBP[FOF_with_MBP]
        """
        ###################################################################################

        self._vars = {} #dictionary to store variables


    def _get_r200(self, idx):
        try:
            _r200 = self._halos[idx].properties['Halo_R_Crit200'].in_units('kpc')
            self._vars["r200_"+str(idx)] = _r200
            return _r200
        except KeyError:
            try:
                return self._halos_bak[idx].properties['Halo_R_Crit200'].in_units('kpc') 
            except KeyError:
                return self._vars["r200_"+str(idx)]

    def _get_mbp(self, idx):
        return self._mbp[idx]

    def _check_backup(self, idx, sub_idx=-1):
        if sub_idx < 0:
            return self._halos[idx]==self._halos_bak[idx]
        else:
            return self._subhalos[idx][sub_idx]==self._subhalos_bak[idx][sub_idx]

    def get_shape(self, idx, sub_idx=-1):
        """
        Returns a pynbody shape. 
        Returns -1 (error) if: 1) The halo has no subhalos
                               2) Either b/a or c/a are 0 or 1 
                                  (since this is unphysical and leads to undefined triaxiality)
        """
        if sub_idx<0:
            with centering_com(self._halos[idx]):
                shape = pn.analysis.halo.halo_shape(self._halos[idx], 
                                                    N=1, bins='lin')
        else:
            if len(self._subhalos[idx])==0:
                warnings.warn('This halo has no subhalos!')
                return -1
            with centering_com(self._subhalos[idx][sub_idx][1]):
                shape = pn.analysis.halo.halo_shape(self._subhalos[idx][sub_idx][1], 
                                                   N=1, bins='lin')
        if shape[1]>1e-5 and shape[2]>1e-5 and shape[1]<.99999 and shape[2]<.99999:
            return shape
        else:
            return -1

    def get_mass200(self, idx):
        try:
            _m200 = self._halos[idx].properties['Halo_M_Crit200'].in_units('Msol')
            self._vars["m200_"+str(idx)] = _m200
            return _m200
        except KeyError:
            try:
                return self._halos_bak[idx].properties['Halo_M_Crit200'].in_units('Msol') 
            except KeyError:
                return self._vars["m200_"+str(idx)]

    def get_number_halos(self):
        return self._N

    def restore(self, idx, sub_idx=-1):
        if idx >= self._N: raise ValueError('The halo index is too large.')
        if type(idx) != int or type(sub_idx) != int: raise TypeError('The indexes must be integers.')
        if sub_idx<0:
            self._halos[idx], self._halos_bak[idx] = self._halos_bak[idx], self._halos[idx]
        else:
            self._subhalos[idx][sub_idx], self._subhalos_bak[idx][sub_idx] = self._subhalos_bak[idx][sub_idx], self._subhalos[idx][sub_idx]

    def filter_(self, halo_idxs=0, sub_idx=-1, filter_str='Sphere_1.2', option=None, sim=True):
        """
        Arguments:

        halo_idxs: indexes of the halos to filter; if more than one, it must be a set
        sub_idx: index of the subhalo. If set to a negative number, the whole halo is filtered.
        filters:
            1. Sphere_radius, where radius is a number relative to r200. Example: Sphere_1
            2. BandPass_prop_min_max. Example: BandPass_y_1 kpc_2 kpc. It must have units.
        option: None, 'half1', 'half2' or 'all' to select halos faster.
                When option!=None 'halos_idxs' are ignored.
                By default only the largest halo is filtered.
        sim: Whether to filter the full simulation or just the specified halo.
             This is relevant because the full simulation also includes particles
             that were not bounded to any halo or subhalo.

        Note: The inputs for each filter are separated using the underscore '_' charachter.
              More filters can be added to 'filter_dict' in the same fashion with variable 
              number of arguments.
        """
        if self._sim==None and sim:
            ValueError('Using the simulation in the filtering requires defining the simulation in the class constructor.')
        if type(halo_idxs) != int and type(halo_idxs) != set:
            raise TypeError('The introduced indexes have the wrong format')
        if type(halo_idxs) is int:
            halo_idxs = {halo_idxs}

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
            print("filter ", i)
            with centering_com(self._halos[i]):
                if sub_idx < 0: #filter the halo
                    if not self._check_backup(i):
                        warnings.warn('This halo backup is already being used!')
                    self._halos_bak[i] = copy.copy(self._halos[i])

                    if sim:
                        self._halos[i] = self._sim[_filter_func(_r200, _split_values, 
                                                            self._halos[i])]
                    else: 
                        self._halos[i] = self._halos[i][_filter_func(_r200, _split_values, 
                                                                 self._halos[i])]
                else: #filter the halo but centered in the subhalo
                    if not self._check_backup(i):
                        warnings.warn('This subhalo backup is already being used!')
                    if self._subhalos[i][sub_idx][0] != sub_idx:
                        raise RuntimeError('Different subhalos are being mixed!')
                        self._subhalos_bak[i][sub_idx] = self._subhalos[i][sub_idx]
                
                    #with centering_com(self._subhalos[i][sub_idx][1]):
                    if sim:
                        _subhalo = self._sim[_filter_func(_r200, _split_values, 
                                                          self._subhalos[i][sub_idx][1])]
                    else:
                        _subhalo = self._halos[i][_filter_func(_r200, _split_values, 
                                                               self._subhalos[i][sub_idx][1])]
                        self._subhalos[i][sub_idx] = (sub_idx,_subhalo)
    
            mem()

    def is_relaxed(self, halo_idx, sub_idx=-1):
        if sub_idx < 0:
            _h = self._halos[halo_idx]
        else:
            if len(self._subhalos[halo_idx])==0:
                warnings.warn('This halo has no subhalos!')
                return -1
            _h = self._subhalos[halo_idx][sub_idx][1]
        with centering_com(_h): #centering the main halo
            print("COM after centering in is_relaxed():", com(_h))
            _com = pn.analysis.halo.center_of_mass(_h)
            _r200 = self._get_r200(halo_idx)
            _dist = np.sqrt( np.power(_h[0]['pos']-_com,2).sum() )
            #The minimum potential point is obtained with:
            #halo[halo_index].subhalo[subhalo_index][0]['pos']
        return bool(_dist < 0.07*_r200)

    def is_resolved(self, halo_idx, sub_idx=-1, intersect_value=1.):
        """
        Return either a boolean or -1, when the density profile is zero somewhere in the middle.
        """
        if sub_idx < 0:
            _h = self._halos[halo_idx]
        else:
            if len(self._subhalos[halo_idx])==0:
                warnings.warn('This halo has no subhalos!')
                return -1
            _h = self._subhalos[halo_idx][sub_idx][1]
        with centering_com(_h): #centering the main halo
            print("Halo is_resolved():", _h)
            _r200 = self._get_r200(halo_idx)
            print("R200 is resolved():", _r200)
            _prof = self.get_profile(halo_idx, sub_idx, component='dm', bins=(2.5,30.,18),
                                     bin_type='log', normalize_x=False)
            _relax = self.relaxation(halo_idx, sub_idx, _prof)
            if len(_relax[0])==1 and len(_relax[1])==1:
                return -1
            _inter = intersect(_relax[1], _relax[0], intersect_value=intersect_value)
            if _inter[1]==True:            
                return _inter[0] < _r200*np.power(10,-1.25)
            else:
                print("Halo number:", halo_idx)
                warning.warn('No intersection was found for this halo. The resolution criteria cannot thus be checked. The halo is assumed to be relaxed.')
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

    def concentration_200(self, idx, sub_idx=-1, 
                          bins=(2.5,30.,18), bin_type='log', normalize_x=False):
        """
        Calculates the concentration at r200 
        """
        if sub_idx<0:
            _current_halo = self._halos[idx] 
        else:
            print('Concentration index:', idx)
            if len(self._subhalos[idx])==0:
                warnings.warn('This halo has no subhalos!')
                return -1
            _current_halo = self._subhalos[idx][sub_idx][1]

        with centering_com(_current_halo):
            print("Halo concentration_200():", _current_halo)
            print("COM after centering in concentration_200():", com(_current_halo))
            #In pynbody, r200 is not defined for subhalos
            #as such, the main halo is always used to retrieve the r200
            _r200 = self._get_r200(idx)
            print("r200:", _r200)
            _profile = self.get_profile(idx, sub_idx, component='dm', 
                                        bins=bins, bin_type=bin_type, 
                                        normalize_x=normalize_x)
            print(_profile['density'])
            if not np.any(_profile['density']):
                warnings.warn('The density array is filled with zeros. The concentration cannot be correctly evaluated.')
                return -1
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
                    bins=(2.5,30,18), bin_type='log', normalize_x=False):
        """
        Obtain the profile of a single halo.
        Arguments:
        'idx': which halo to get the profile from
        'sub_idx': which subhalo to get the profile from
        'component': 'dm', 'stars', 'gas' or 'all'
        'bins': tuple with (start,stop,nbins)
        'bin_type': either 'linear' or 'log'
        'normalize_x': 'None' for no normalization and 'r200' for normalization
        'centering': whether to centre the halo before obtaining the profile or not
        """
        _h = self._halos[idx] if sub_idx<0 else self._subhalos[idx][sub_idx][1]
        with centering_com(_h):
            components = {'dm':_h.dm, 'stars':_h.s, 'gas':_h.g, 'all':_h}
            if component not in components:
                raise ValueError('The specified component does not exist.')
            if type(bins) is not tuple or len(bins) != 3:
                raise TypeError('Please insert bins with the following pattern: (start, stop, nbins).')
        
            if normalize_x:
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


    def relaxation(self, idx, sub_idx=-1, profile=None, component='dm', bins=(2.5,30,18), bin_type='log'):
        """
        Creates a t_relax(r) / t_circ(r200) array for the provided halo.
        If you do not provide a profile, the method creates one for you. 
        The profile is needed to obtain the density as a function of r.
        
        Returns:
        The above array and the corresponding binned radius values
        """
        _h = self._halos[idx]
        _rho_crit_z = pn.array.SimArray(pn.analysis.cosmology.rho_crit(_h, unit="Msol kpc**-3"), 
                                        "Msol kpc**-3")
        if sub_idx >= 0:
            _h = self._subhalos[idx][sub_idx][1]
        with centering_com(_h):
            if profile==None:
                _prof = self.get_profile(idx, sub_idx, component=component, bins=bins, 
                                         bin_type=bin_type, normalize_x=False)
                _density = _prof['density']
                _rbins = _prof['rbins']
                print(_density)
                print(_rbins)
            else:
                _density = profile['density']
                _rbins= profile['rbins']

            _nparts = pn.array.SimArray([len(_h[_h['r']<_rbins[i]]) for i in range(len(_rbins))],'')

        if 0. in _density:
            print("Halo number", idx, "has some zero in the density.")
            zero_idxs = [i for i,item in enumerate(_density) if item==0]

            if zero_idxs[0]==len(_density)-1:
                _density[_density>0]
                _rbins[_density>0]
            else:
                if _density[zero_idxs[0]+1]>0:
                    print(_density)
                    warnings.warn('One of the density values not in the far end of the distribution is zero. The relaxation cannot thus be calculated.')
                    return np.zeros(1), np.zeros(1)
                else: 
                    _density[_density>0]
                    _rbins[_density>0]
        return (np.sqrt(200)*_nparts)/(8*np.log(_nparts))*np.sqrt(_rho_crit_z/_density), _rbins
