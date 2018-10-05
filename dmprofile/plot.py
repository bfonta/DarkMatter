import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import itertools, collections

from pynbody.analysis.theoretical_profiles import NFWprofile
from dmprofile.fit import nfw_fit
from dmprofile.utilities import intersect, is_list_empty

#decorator: applies the same function to all objects in _Plot._obj_list
def do_all_objects(func):
    def wrapper(self, *args, **kwargs):
        kwargs['i_profile'] = 0
        for i in range(self._nrows):
            for j in range(self._ncols):
                if kwargs['i_profile'] < self._N:
                    kwargs['i'] = i
                    kwargs['j'] = j
                    func(self, *args, **kwargs)
                    kwargs['i_profile'] += 1
    return wrapper


class _Plot():
    """
    Important: each individual object contained inside obj_list is meant to be converted into a plot, so in general each object will be a list of values.
    """
    def __init__(self, obj_list, width, height, nrows, ncols, name=''):
        if type(obj_list) is not list:
            raise TypeError
        if type(width) is not int or type(height) is not int or type(nrows) is not int or type(ncols) is not int:
            raise TypeError
        self._obj_list = obj_list
        self._N = len(self._obj_list)
        self._nrows, self._ncols = nrows, ncols
        self._name = name

    def _geometry(self, n):
        """
        Returns number of rows '_r' and columns '_c' that any plot with more than 2 figure has.
        'n' is the number of figures to plot
        """
        _r = np.floor(np.sqrt(n))
        _c = _r
        _r_iter, _c_iter = itertools.cycle(range(2)), itertools.cycle(range(2))
        next(_r_iter)
        while _r*_c<self._N:
            _r += next(_r_iter)
            _c += next(_c_iter)
        return int(_r), int(_c)
        
    def _set_axis_indexes(self, axis_idx):
        """
        Set the axis indexes into 1D or 2D, due to matplotlib axis indexing.
        1D is needed when one plans to plot less than two pictures.
        """
        if self._N>2:
            return axis_idx[0], axis_idx[1]
        else:
            return axis_idx[0]

    def set_title(self, title):
        self._fig.suptitle(title)
        
    def set_name(self, name):
        self._name = name

    def set_axis(self, idx, xlabel, ylabel, xscale='linear', yscale='linear'):
        """
        idx can be either a single value or a tuple, depending on the total number of profiles
        """
        if self._N>2:
            indexes = idx[0], idx[1]
        else:
            indexes = idx[0]
        self._axis[indexes].set_xscale(xscale)
        self._axis[indexes].set_yscale(yscale)
        self._axis[indexes].set_xlabel(xlabel)
        self._axis[indexes].set_ylabel(ylabel)

    def set_axis_all(self, xlabel, ylabel, xscale='linear', yscale='linear'):
        if self._N>2:
            _r, _c = self._geometry(self._N)
            for i in range(_r):
                for j in range(_c):
                    self.set_axis((i,j), xlabel, ylabel, xscale, yscale)
        else:
            for i in range(self._N):
                self.set_axis((i,0), xlabel, ylabel, xscale, yscale)

    def set_all_properties(self, title="", model=None, 
                       xlabel='R [kpc]', ylabel=None, 
                       xscale='linear', yscale='linear'):
        """
        model options: 'density_profile', 'mass_enc_profile', 'mass_profile',
                       'density_enc_profile', 'relaxation', 'concentration_mass'
        """
        self.set_title(title)
        scale_options = ['linear', 'log']
        if xscale not in scale_options or yscale not in scale_options:
            raise ValueError('The selected scale is not supported.')
        model_options = ['density_profile', 
                         'density_enc_profile',
                         'mass_enc_profile', 
                         'mass_profile',
                         'concentration_mass',
                         'relaxation']
        if model not in model_options:
            raise ValueError('The specified model is not currently predefined.')
        if model!='None':
            self.set_axis_all('R [kpc]', r'$\rho$ [M$_{\odot}$ kpc$^{-3}$]', xscale, yscale) if model=='density_profile' else self.set_axis_all('R [kpc]', r'$\rho_{enc}$ [M$_{\odot}$ kpc$^{-3}$]', xscale, yscale) if model=='density_enc_profile' else self.set_axis_all('R [kpc]', r'M$_{enc}$ [M$_{\odot}$]', xscale, yscale) if model=='mass_enc_profile' else self.set_axis_all(r'M$_{200}$ [M$_{\odot}$]', r'c$_{200}$', xscale, yscale) if model=='concentration_mass' else self.set_axis_all('R [kpc]', r't$_{relax}(R)$/t$_{circ}$(r$_{200}$)', xscale, yscale) if model=='relaxation' else self.set_axis_all('R [kpc]', r'M [M$_{\odot}$]', xscale, yscale)
        else:
            for i in range(self._nrows):
                if self._N>2:
                    for j in range(self._ncols):
                        self.set_axis((i,j), xlabel, ylabel, xscale, yscale)
                else:
                    self.set_axis((i,0), xlabel, ylabel, xscale, yscale)
                    
    def savefig(self, name=''):
        if self._name=='':
            try:
                self._fig.savefig(name)
            except ValueError:
                print("A figure cannot be saved without having a name.")
        else:
            self._fig.savefig(self._name)


class Profile(_Plot):
    """
    Plots all kinds of pynbody related profiles.

    Args:
    p: list of pynbody profile objects
    """
    def __init__(self, p, w=10, h=9, nrows=0, ncols=0, name=''):
        super().__init__(p, w, h, nrows, ncols, name)
        self._p = self._obj_list
        if nrows==0 or ncols==0:
            self._nrows, self._ncols = super()._geometry(self._N)
        self._fig, self._axis = plt.subplots(nrows=self._nrows, ncols=self._ncols, figsize=(w,h))

    def plot(self, idx, axis_idx, x_var, y_var):
        """
        Args:
        1. 'idx':  profile index (int)
        2. 'axis_idx': axis index (int [self._N<=2] or tuple [self._N>2])
        3. 'x_var': variable stored in the profile to plot as x variable
            options: radius, mass, mass_enc, density
        4. 'y_var': variable stored in the profile to plot as y variable
            options: radius, mass, mass_enc, density
        """
        def _define_vars(x_var, y_var):
            _translate = {'radius': self._p[idx]['rbins'],
                          'mass': self._p[idx]['mass'],
                          'mass_enc': self._p[idx]['mass_enc'],
                          'density': self._p[idx]['density'],
                          'density_enc': (3*self._p[idx]['mass_enc'])
                          /(4*np.pi*np.power(self._p[idx]['rbins'],3))}
            if x_var not in _translate.keys() or y_var not in _translate.keys():
                raise ValueError
            return _translate[x_var], _translate[y_var]

        _x, _y = _define_vars(x_var, y_var)
        indexes = super()._set_axis_indexes(axis_idx)
        self._axis[indexes].plot(_x, _y)

    @do_all_objects
    def plot_all(self, x_var, y_var, *args, **kwargs):
        self.plot(kwargs['i_profile'], (kwargs['i'], kwargs['j']), x_var, y_var)

    def fit_and_plot(self, idx, axis_idx, function='nfw'):
        functions_implemented = ['nfw']
        if function not in functions_implemented:
            raise ValueError('The specified function has not yet been implemented.')
        
        if function=='nfw':
            _fit = nfw_fit(self._p[idx])
            _func = NFWprofile.profile_functional_static(np.array(self._p[idx]['rbins']), 
                                                         _fit[0], _fit[1])
            indexes = super()._set_axis_indexes(axis_idx)
            self._axis[indexes].plot(np.array(self._p[idx]['rbins']), 
                                     _func, label='NFW fit')
            self._axis[indexes].legend()

    @do_all_objects
    def fit_and_plot_all(self, function='nfw', *args, **kwargs):
        self.fit_and_plot(kwargs['i_profile'], (kwargs['i'], kwargs['j']))

        
    def draw_line(self, axis_idx, value, orientation='v', label='---', color='grey', 
                  linestyle = '--'):
        """
        Arguments:
        'idx': profile index
        'axis_idx' index of the axis where the line will be plotted
        'value': value at which will be plotted; it is related to the 'orientation' option ('v' or 'h')
        """
        indexes = super()._set_axis_indexes(axis_idx)
        _xmin, _xmax = self._axis[indexes].get_xbound()
        _ymin, _ymax = self._axis[indexes].get_ybound()
        if orientation=='v':
            self._axis[indexes].plot([value,value],[_ymin,_ymax], 
                                     linestyle=linestyle, color=color, label=label)
        elif orientation=='h':
            self._axis[indexes].plot([_xmin,_xmax],[value,value],
                                     linestyle=linestyle, color=color, label=label)
        else:
            raise ValueError('The specified orientation is not supported.')
        self._axis[indexes].legend()

class Shape(_Plot):
    """
    Plots all kinds of pynbody related shapes.

    Args:
    s: list of shape arrays for plotting, i.e., list of lists of pynbody shape objects
    extra_var: list of some additional property common to all objects from which the shapes were taken
           for example: radius, mass
    """
    def __init__(self, s, extra_var=[], w=10, h=9, nrows=0, ncols=0, name=''):
        super().__init__(s, w, h, nrows, ncols, name)
        self._s = self._obj_list
        self._extra_var = extra_var
        if nrows==0 or ncols==0:
            self._nrows, self._ncols = super()._geometry(self._N)
        self._fig, self._axis = plt.subplots(nrows=self._nrows, ncols=self._ncols, figsize=(w,h))

    def _shape_array_to_numpy(self, s):
        """
        Converts the shape formatting to allow correct indexing.
        This is needed for plotting.
        """
        _snew = []
        for i in range(len(s)):
            _a, _b, _c, _d, _e, _f = ([] for i in range(6))
            for j in range(len(s[i])):
                _a.append(s[i][j][0][0]) #pos
                _b.append(s[i][j][1][0]) #b/a
                _c.append(s[i][j][2][0]) #c/a
                _d.append(s[i][j][3][0]) #alignment
                _e.append(s[i][j][4][0]) #rotation matrix
                _f.append((1-s[i][j][1][0]**2)/(1-s[i][j][2][0]**2)) if 1-s[i][j][2][0]**2<1 else _f.append(1.1) #triaxiality 
            _snew.append([_a,_b,_c,_d,_e,_f])
        return _snew

    def scatter_plot(self, idx, axis_idx, x_var, y_var, resolved_bools=[], relaxed_bools=[]):
        """
        Args:
        1. 'idx': pynbody shape list object. It refers to the way self._s is ordered
        2. 'axis_idx': axis index (int [self._N<=2] or tuple [self._N>2])
        3. 'x_var': variable stored in the shape to plot as x variable
            options: pos, b/a, c/a, align, rot, triax, m200, r200
        4. 'y_var': variable stored in the profile to plot as y variable
            options: pos, b/a, c/a, align, rot, triax, m200, r200
        """
        handles = []
        indexes = super()._set_axis_indexes(axis_idx)

        def _define_vars(x_var, y_var):
            """
            Defines which arrays are to be plotted according to user input.
            See online documentation of pynbody.analysis.halo.halo_shape() method.
            """
            _s_for_plot = self._shape_array_to_numpy(self._s)
            _match = {'position': 0,
                      'b/a': 1,
                      'c/a': 2,
                      'align': 3,
                      'rot': 4,
                      'triax': 5}
            if x_var not in _match.keys() and y_var not in _match.keys():
                raise ValueError("At least one of the variables being plot must be related to the shape pynbody method.")
            elif x_var in _match.keys() and y_var in _match.keys():
                return _s_for_plot[idx][_match[x_var]], _s_for_plot[idx][_match[y_var]]
            else: #either x or y are not contained inside the self._s object
                if x_var not in _match.keys():
                    return self._extra_var, _s_for_plot[idx][_match[y_var]]
                else:
                    return _s_for_plot[idx][_match[x_var]], self._extra_var

        if not is_list_empty(resolved_bools):
            grey = (0.2,0.2,0.2)
            lightgrey = (0.6,0.6,0.6)
            resolution_colors = [lightgrey if resolved_bools[i]==False else grey for i in range(len(self._s[idx]))]
            handles.append(mpatches.Patch(color=lightgrey, label='unresolved structures'))
            handles.append(mpatches.Patch(color=grey, label='resolved structures'))
        else:
            resolution_colors = [(0.4,0.4,0.4)]*len(self._s[idx])

        if not is_list_empty(relaxed_bools):
            relaxation_markers = ['o' if relaxed_bools[i]==False else 'v' for i in range(len(self._s[idx]))]  
            handles.append(plt.plot([],[], marker='o',ms=8,c='k',
                                    linestyle='None',label='unrelaxed structure')[0])
            handles.append(plt.plot([],[], marker='v',ms=8,c='k',
                                    linestyle='None',label='relaxed structure')[0])
            
            self._axis[indexes].legend(handles=handles)
        else:
            relaxation_markers = ['v']*len(self._s[idx])
        
        _x, _y = _define_vars(x_var, y_var)
        _sorted_lists = sorted(zip(_x,_y), key=lambda x: x[0])
        _x, _y = [[q[i] for q in _sorted_lists] for i in range(2)]
        
        for xp, yp, c, m in zip(_x, _y, resolution_colors, relaxation_markers):
            self._axis[indexes].scatter(xp, yp, c=c, marker=m)


class Concentration(_Plot):
    """
    Plots all kinds of pynbody related concentrations.

    Arguments:
    c: list of lists of concentration values
    extra_var: list of some additional property common to all objects from which the shapes were taken
           for example: radius, mass
    """
    def __init__(self, c, extra_var=[], w=10, h=9, nrows=0, ncols=0, name=''):
        super().__init__(c, w, h, nrows, ncols, name)
        self._c = self._obj_list
        self._extra_var = extra_var
        if nrows==0 or ncols==0:
            self._nrows, self._ncols = super()._geometry(self._N)
        self._fig, self._axis = plt.subplots(nrows=self._nrows, ncols=self._ncols, figsize=(w,h))

    def scatter_plot(self, idx, axis_idx, concentration_axis='y', 
                     resolved_bools=[], relaxed_bools=[]):
        """
        Args:
        'resolved_bools': Resolved halos. The order refers to the concentration array included when calling the class constructor.
        'relaxed_bools': Relaxed halos. The order refers to the concentration array included when calling the class constructor.
        """
        axis_options = ['x','y']
        handles = []
        if concentration_axis not in axis_options:
            raise ValueError('The axis you selected for plotting the concentration is not valid.')
        if idx >= self._N:
            raise ValueError('You cannot plot something that does not exist :)')
        indexes = super()._set_axis_indexes(axis_idx)
        
        if not is_list_empty(resolved_bools):
            grey = (0.2,0.2,0.2)
            lightgrey = (0.6,0.6,0.6)
            resolution_colors = [lightgrey if resolved_bools[i]==False else grey for i in range(len(self._c[idx]))]
            handles.append(mpatches.Patch(color=lightgrey, label='unresolved structures'))
            handles.append(mpatches.Patch(color=grey, label='resolved structures'))
        else:
            resolution_colors = [(0.4,0.4,0.4)]*len(self._c[idx])

        if not is_list_empty(relaxed_bools):
            relaxation_markers = ['o' if relaxed_bools[i]==False else 'v' for i in range(len(self._c[idx]))]  
            handles.append(plt.plot([],[], marker='o',ms=8,c='k',
                                    linestyle='None',label='unrelaxed structure')[0])
            handles.append(plt.plot([],[], marker='v',ms=8,c='k',
                                    linestyle='None',label='relaxed structure')[0])
            
            self._axis[indexes].legend(handles=handles)
        else:
            relaxation_markers = ['v']*len(self._c[idx])
        for xp, yp, c, m in zip(self._c[idx], self._extra_var, resolution_colors, relaxation_markers):
            if concentration_axis=='x':
                self._axis[indexes].scatter(xp, yp, c=c, marker=m)
            else:
                self._axis[indexes].scatter(yp, xp, c=c, marker=m)

class Relaxation(_Plot):
    """
    Args:
    relax: list of pynbody relaxation arrays (obtained for example from dmprofile.halos.relaxation)
    r: corresponding list of radius
    """
    def __init__(self, relax, radius, w=10, h=9, nrows=0, ncols=0, name=''):
        super().__init__(relax, w, h, nrows, ncols, name)
        self._relax = self._obj_list
        self._radius = radius
        if nrows==0 or ncols==0:
            self._nrows, self._ncols = super()._geometry(self._N)
        self._fig, self._axis = plt.subplots(nrows=self._nrows, ncols=self._ncols, figsize=(w,h))

    def plot(self, idx, axis_idx, x=None, y=None):
        """
        Args:
        'idx':  profile index (int)
        'axis_idx': axis index (int [self._N<=2] or tuple [self._N>2])
        'x' and 'y': data arrays
        """
        indexes = super()._set_axis_indexes(axis_idx)
        if x==None and y==None:
            self._axis[indexes].plot(self._radius[idx], self._relax[idx])
        else:
            self._axis[indexes].plot(x, y)

    @do_all_objects
    def plot_all(self, x, y, *args, **kwargs):
        self.plot(kwargs['i_profile'], (kwargs['i'], kwargs['j']))

    def intersect_and_plot(self, idx, axis_idx, intersect_value=1.):
        indexes = super()._set_axis_indexes(axis_idx)
        degree = 5
        self._axis[indexes].plot(self._radius[idx], self._relax[idx], '.') 
        _intersect_x = intersect(x=self._radius[idx], y=self._relax[idx], deg=degree)
        _polyfit = np.polyfit(self._radius[idx], self._relax[idx], deg=degree)
        _func = np.poly1d(_polyfit)
        _x_new = np.linspace(self._radius[idx][0], self._relax[idx][-1], 10000)
        _y_new = _func(_x_new)
        self._axis[indexes].plot(_x_new, _y_new, markersize=1, label='polynomial fit') #draw pseudo-line
        self._axis[indexes].plot([min(self._radius[idx]), max(self._radius[idx])], 
                                 [intersect_value, intersect_value], 
                                 color='r', linestyle='--', 
                                 label=r't$_{relax}(R)$/t$_{circ}$(r$_{200}$) = '+str(intersect_value))
        _half_bin_radius = (self._radius[idx][1]-self._radius[idx][0])/2
        _half_bin_relax = (self._relax[idx][1]-self._relax[idx][0])/2
        _text_idx = len(self._relax[idx])-int(len(self._relax[idx])/4)
        if _intersect_x[1]==True: #if there was an intersection
            self._axis[indexes].plot(_intersect_x[0], intersect_value, 'ko')
            self._axis[indexes].annotate('Intersection at x='+str(np.round_(_intersect_x[0],3)), xy=(0.015,0.75), xycoords='axes fraction')
        self._axis[indexes].set_xlim([self._radius[idx][0] -_half_bin_radius, 
                                      self._radius[idx][-1] + _half_bin_radius])
        self._axis[indexes].set_ylim([self._relax[idx][0] + _half_bin_relax, 
                                      self._relax[idx][-1] - _half_bin_relax])
        self._axis[indexes].legend()

    @do_all_objects
    def intersect_and_plot_all(self, function='nfw', *args, **kwargs):
        self.intersect_and_plot(kwargs['i_profile'], (kwargs['i'], kwargs['j']))
