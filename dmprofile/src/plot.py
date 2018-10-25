import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import itertools, collections
from scipy.optimize import curve_fit

from pynbody.analysis.theoretical_profiles import NFWprofile
from dmprofile.src.utilities import linear_function_generator, remove_low_occupancy_bins
from .fit import nfw_fit
from .utilities import intersect, is_list_empty

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
    Returns number of rows '_r' and columns '_c' that any plot with more than 2 figure has.
    This behaiour has no effect if the user manually specifies the number of rows and columns.
    The variable 'n' is the number of figures to plot.
    """
    def __init__(self, obj_list, width, height, nrows, ncols, name=''):
        if type(obj_list) is not list:
            raise TypeError
        if type(width) is not int or type(height) is not int or type(nrows) is not int or type(ncols) is not int:
            raise TypeError
        self._obj_list = obj_list
        self._N = len(self._obj_list)
        self._width, self._height = width, height
        self._name = name
        
        if nrows!=0 and ncols!=0:
            self._nrows, self._ncols = nrows, ncols
            self._fig, self._axis = plt.subplots(nrows=self._nrows, ncols=self._ncols,
                                                 figsize=(self._width, self._height))
        else:
            _r = np.floor(np.sqrt(self._N))
            _c = _r
            _r_iter, _c_iter = itertools.cycle(range(2)), itertools.cycle(range(2))
            next(_r_iter)
            while _r*_c<self._N:
                _r += next(_r_iter)
                _c += next(_c_iter)
            self._nrows, self._ncols = int(_r), int(_c)
            self._fig, self._axis = plt.subplots(nrows=self._nrows, ncols=self._ncols,
                                                 figsize=(self._width, self._height))

    def _set_axis_indexes(self, axis_idx):
        """
        Set the axis indexes into 1D or 2D, due to matplotlib axis indexing.
        1D is needed when one plans to plot less than two pictures.
        """
        if self._nrows>=2 and self._ncols>=2:
            return axis_idx[0], axis_idx[1]
        elif self._nrows==1:
            return axis_idx[1]
        else: #ncols==1 & nrows>=2
            return axis_idx[0]

    def set_title(self, title):
        self._fig.suptitle(title)
        
    def set_name(self, name):
        self._name = name

    def set_axis(self, axis_idx, xlabel, ylabel, xscale='linear', yscale='linear',
                 xmin=None, xmax=None, ymin=None, ymax=None):
        """
        axis_idx can be either a single value or a tuple, depending on the total number of profiles
        """
        indexes = self._set_axis_indexes(axis_idx)

        if xmin!=None and xmax!=None:
            self._axis[indexes].set_xlim([xmin,xmax])
        if ymin!=None and ymax!=None:
            self._axis[indexes].set_ylim([ymin,ymax])

        if (xmin==None and xmax!=None) or (xmin!=None and xmax==None):
            raise ValueError('Please provide both xmin and xmax options.')
        if (ymin==None and ymax!=None) or (ymin!=None and ymax==None):
            raise ValueError('Please provide both ymin and ymax options.')

        self._axis[indexes].set_xscale(xscale)
        self._axis[indexes].set_yscale(yscale)
        self._axis[indexes].set_xlabel(xlabel)
        self._axis[indexes].set_ylabel(ylabel)

    def _set_axis_3D(self, axis, xlabel='', xscale='linear', ylabel='', yscale='linear',
                     zlabel='', zscale='linear',
                     xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None):

        if (xmin==None and xmax!=None) or (xmin!=None and xmax==None):
            raise ValueError('Please provide both xmin and xmax options.')
        if (ymin==None and ymax!=None) or (ymin!=None and ymax==None):
            raise ValueError('Please provide both ymin and ymax options.')
        if (zmin==None and zmax!=None) or (zmin!=None and zmax==None):
            raise ValueError('Please provide both zmin and zmax options.')

        if xmin!=None and xmax!=None:
            axis.set_xlim([xmin,xmax])
        if ymin!=None and ymax!=None:
            axis.set_ylim([ymin,ymax])        
        if zmin!=None and zmax!=None:
            axis.set_zlim([zmin,zmax])

        axis.set_xlabel(xlabel)
        axis.set_xscale(xscale)
        axis.set_xscale(xscale)
        axis.set_xlabel(xlabel)
        axis.set_yscale(yscale)
        axis.set_yscale(yscale)
        axis.set_zlabel(zlabel)
        axis.set_zscale(zscale)
        axis.set_zscale(zscale)

    def set_axis_all(self, xlabel='', ylabel='', xscale='linear', yscale='linear',
                     xmin=None, xmax=None, ymin=None, ymax=None):
        if self._N>2:
            _r, _c = self._nrows, self._ncols
            for i in range(_r):
                for j in range(_c):
                    self.set_axis((i,j), xlabel, ylabel, xscale, yscale,
                                  xmin, xmax, ymin, ymax)
        else:
            for i in range(self._N):
                self.set_axis((i,0), xlabel, ylabel, xscale, yscale,
                              xmin, xmax, ymin, ymax)

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

    def scatter_plot_3D(self, x, y, z, axis_idx, label='', linestyle='None', color='b',
                        xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None,
                        azim=0, elev=10):
        """
        Note: The 2D index has to be converted into a 1D index whenever the plots is 3D 
              (see matplotlib three-dimensional plotting)
        """
        _axis_idx_3D = axis_idx[1] + self._ncols*axis_idx[0] + 1
        _axis = self._fig_3D.add_subplot(self._nrows, self._ncols, _axis_idx_3D, projection='3d')
        _axis.scatter(x, y, z, zdir='z', marker='.', linestyle=linestyle,
                      label=label, color=color)
        _axis.azim = azim
        _axis.elev = elev
        if label != '': self._axis[indexes].legend()
        self._set_axis_3D(axis=_axis, zlabel='', zscale='linear',
                         xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax)

    def binned_plot(self, axis_idx, x, y, nbins, xerr, yerr, 
                    marker, color, label, shift, fit, xscale,
                    min_bins=3):
        """x
        idx: which of the objects to plot (ordered by the constructor)
        axis_idx: where to plot the plot. Example: (0,1)
        x: unbinned 'x' data array
        y: unbinned 'y' data array
        nbins: number of bins, based on the 'x' variable distribution
        xerr: 'x' data uncertainties array
        yerr: 'y' data uncertainties array
        """
        if type(nbins)!=int:
            raise TypeError('The number of bins has to be an integer value.')
        indexes = self._set_axis_indexes(axis_idx)
        self._axis[indexes].legend().set_visible(False)
        
        x = np.array(x)
        y = np.array(y)
        if xscale=='log':
            bin_edges = np.logspace(np.log10(np.amin(x)), np.log10(np.amax(x)), nbins+1)
        else:
            _hh = np.histogram(x,nbins)
            print("HHHHHHH", _hh[0])
            bin_edges = remove_low_occupancy_bins(_hh, min_bins)
            nbins = len(bin_edges)-1

        bin_edges_2 = np.roll(bin_edges,-1)[:-1]
        bin_edges = bin_edges[:-1]

        y_binned_values = [[] for _ in range(nbins)]

        for ix,iy in zip(x,y):
            for b in range(nbins):
                if ix>=bin_edges[b] and ix<=bin_edges_2[b]:
                    y_binned_values[b].append(iy)

        y_mean = [sum(y_binned_values[i])/len(y_binned_values[i]) if len(y_binned_values[i])!=0 else 0. for i in range(nbins)]

        #fig, ax = plt.subplots(figsize=(10,6))
        bin_centre = bin_edges + (bin_edges_2-bin_edges)/2
        bin_centre = bin_centre + shift

        if xerr!=None:
            xerr = [list(bin_centre-shift-bin_edges), list(bin_centre-shift-bin_edges)]
        yerr = [[0]*nbins, [0]*nbins]

        self._axis[indexes].errorbar(np.array(bin_centre), np.array(y_mean), 
                                     xerr=xerr, yerr=yerr, fmt='none', ecolor=color, capsize=3)
        self._axis[indexes].scatter(np.array(bin_centre), np.array(y_mean),
                                    marker=marker, s=40, color=color, label=label)
        if label != '': 
            self._axis[indexes].legend()

        if fit:
            slope, intercept = curve_fit(lambda _x,_p0,_p1: _p0*_x+_p1, bin_centre, y_mean)[0]
            print("PARAMETERS:", slope, intercept)
            gen_obj = linear_function_generator(slope, intercept, bin_edges[0], bin_edges_2[-1])
            _x_fit, _y_fit = ([] for i in range(2))
            for i in range(1000):
                _next = next(gen_obj)
                _x_fit.append(_next[0])
                _y_fit.append(_next[1])
            self._axis[indexes].plot(np.array(_x_fit), np.array(_y_fit), color=color, linewidth=.6)
            self._axis[indexes].annotate('Slope: ' + str(np.round_(slope,2)) + 
                                         '\nInterception: ' + str(np.round_(intercept,2)), 
                                         xy=(0.02,0.02), xycoords='axes fraction',
                                         fontname='sans-serif', fontstyle='oblique', 
                                         fontweight='semibold')

    def draw_line(self, axis_idx, value, orientation, label, color, linestyle):
        """
        Arguments:
        'idx': profile index
        'axis_idx' index of the axis where the line will be plotted
        'value': value at which will be plotted; it is related to the 'orientation' option ('v' or 'h')
        """
        indexes = self._set_axis_indexes(axis_idx)
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
        super().draw_line(axis_idx=axis_idx, value=value, orientation=orientation, 
                          label=label, color=color, linestyle=linestyle)
        
class Shape(_Plot):
    """
    Plots all kinds of pynbody related shapes.

    Args:
    s: List of shape arrays for plotting, i.e., list of lists of pynbody shape objects
    extra_var: List of lists of some additional properties. Example: radius, mass
               By default, the first list of this list is chosen for all plots
               It can be controlled with the 'extra_idx' parameter, which follows the order
               of the list of shape objects.
    """
    def __init__(self, s, extra_var=[[]], w=10, h=9, nrows=0, ncols=0, name=''):
        super().__init__(s, w, h, nrows, ncols, name)
        self._s = self._obj_list
        self._extra_var = extra_var

    def _shape_array_to_numpy(self, s):
        """
        Converts the shape formatting to allow correct indexing.
        This is needed for plotting.
        """
        _snew = []
        for i in range(len(s)):
            _a, _b, _c, _d, _e, _f = ([] for i in range(6))
            for j in range(len(s[i])):
                if s[i][j][1][0]==1 or s[i][j][2][0]==1:
                    raise ValueError('The triaxiality becomes NaN when either b/a or c/a are equal to one. Please do not introduce those values.')
                if s[i][j][1][0]==0 or s[i][j][2][0]==0:
                    raise ValueError('The triaxiality becomes NaN when either b/a or c/a are equal to zero. Please do not introduce those values.')
                _a.append(s[i][j][0][0]) #pos
                _b.append(s[i][j][1][0]) #b/a
                _c.append(s[i][j][2][0]) #c/a
                _d.append(s[i][j][3][0]) #alignment
                _e.append(s[i][j][4][0]) #rotation matrix
                _f.append((1-s[i][j][1][0]**2)/(1-s[i][j][2][0]**2)) if 1-s[i][j][2][0]**2<1 else _f.append(1.1) #triaxiality 
            _snew.append([_a,_b,_c,_d,_e,_f])
        return _snew

    def _define_vars(self, idx, x_var, y_var, extra_idx=0):
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
                return self._extra_var[extra_idx], _s_for_plot[idx][_match[y_var]]
            else:
                return _s_for_plot[idx][_match[x_var]], self._extra_var[extra_idx]


    def scatter_plot(self, idx, axis_idx, x_var, y_var, extra_idx=0, resolved_bools=[], relaxed_bools=[]):
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

        if not is_list_empty(resolved_bools):
            grey = (0.2,0.2,0.2)
            lightgrey = (0.6,0.6,0.6)
            resolution_colors = [lightgrey if resolved_bools[i]==False else grey for i in range(len(resolved_bools))]
            handles.append(mpatches.Patch(color=lightgrey, label='unresolved structures'))
            handles.append(mpatches.Patch(color=grey, label='resolved structures'))
        else:
            resolution_colors = [(0.4,0.4,0.4)]*len(self._s[idx])

        if not is_list_empty(relaxed_bools):
            relaxation_markers = ['o' if relaxed_bools[i]==False else 'v' for i in range(len(relaxed_bools))]  
            handles.append(plt.plot([],[], marker='o',ms=8,c='k',
                                    linestyle='None',label='unrelaxed structure')[0])
            handles.append(plt.plot([],[], marker='v',ms=8,c='k',
                                    linestyle='None',label='relaxed structure')[0])
            self._axis[indexes].legend(handles=handles)
        else:
            relaxation_markers = ['v']*len(self._s[idx])
        
        _x, _y = self._define_vars(idx, x_var, y_var, extra_idx)
        _sorted_lists = sorted(zip(_x,_y), key=lambda x: x[0])
        _x, _y = [[q[i] for q in _sorted_lists] for i in range(2)]
        
        for xp, yp, c, m in zip(_x, _y, resolution_colors, relaxation_markers):
            self._axis[indexes].scatter(xp, yp, c=c, marker=m)


    def binned_plot(self, idx, axis_idx, nbins, x_var, y_var, extra_idx=0, xerr=None, yerr=None,
                    marker='o', color='b', label='', shift=0., fit=False, xscale='linear'):
        _x, _y = self._define_vars(idx, x_var, y_var, extra_idx)
        super().binned_plot(axis_idx=axis_idx, 
                            x=_x, y=_y, 
                            nbins=nbins, xerr=xerr, yerr=yerr,
                            marker=marker, color=color, label=label, shift=shift, fit=fit, xscale=xscale)

class Concentration(_Plot):
    """
    Plots all kinds of pynbody related concentrations.

    Arguments:
    c: list of lists of concentration values
    extra_var: List of lists of some additional properties. Example: radius, mass
               By default, the first list of this list is chosen for all plots
               It can be controlled with the 'extra_idx' parameter, which follows the order
               of the list of concentration objects.
    """
    def __init__(self, c, extra_var=[[]], w=10, h=9, nrows=0, ncols=0, name=''):
        super().__init__(c, w, h, nrows, ncols, name)
        self._c = self._obj_list
        self._extra_var = extra_var

    def scatter_plot(self, idx, axis_idx, extra_idx=0, concentration_axis='y', 
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
        for xp, yp, c, m in zip(self._c[idx], self._extra_var[extra_idx], resolution_colors, relaxation_markers):
            if concentration_axis=='x':
                self._axis[indexes].scatter(xp, yp, c=c, marker=m)
            else:
                self._axis[indexes].scatter(yp, xp, c=c, marker=m)

    def binned_plot(self, idx, axis_idx, nbins, extra_idx=0, xerr=None, yerr=None,
                    marker='o', color='b', label='', shift=0., fit=False, xscale='linear'):
        _x = np.array(self._extra_var[extra_idx])
        _y = np.array(self._c[idx])
        super().binned_plot(axis_idx=axis_idx, 
                            x=_x, y=_y, 
                            nbins=nbins, xerr=xerr, yerr=yerr,
                            marker=marker, color=color, label=label, shift=shift, fit=fit, xscale=xscale)


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

    def draw_line(self, axis_idx, value, orientation='v', label='---', color='grey', 
                  linestyle = '--'):
        """
        Arguments:
        'idx': profile index
        'axis_idx' index of the axis where the line will be plotted
        'value': value at which will be plotted; it is related to the 'orientation' option ('v' or 'h')
        """
        super().draw_line(axis_idx=axis_idx, value=value, orientation=orientation, 
                          label=label, color=color, linestyle=linestyle)


class Particles(_Plot):
    def __init__(self, pos, w=10, h=9, nrows=0, ncols=0, name=''):
        """
        The 'pos' argument is a list of 2D lists with the 'x', 'y' and 'z' positions of the particles.
        Each list corresponds to a plot.
        
        Example:
        Using a pynbody simulation directly: instance = Particles([sim, sim_2])
        """
        super().__init__(pos, w, h, nrows, ncols, name)
        self._fig_3D = plt.figure(figsize=(self._width, self._height))
        self._pos = self._obj_list

    def scatter_plot_3D(self, idx, axis_idx, label='', linestyle='None', color='b',
                        xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None,
                        azim=0, elev=10):
        super().scatter_plot_3D(x=self._pos[idx]['x'], y=self._pos[idx]['y'], z=self._pos[idx]['z'],
                                axis_idx=axis_idx, label=label, linestyle=linestyle, color=color,
                                xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax,
                                azim=azim, elev=elev)

    def savefig(self, name='', option='2D'):
        """
        overload saving function
        option: '2D' or '3D'
        """
        options = ['2D', '3D']
        if option not in options:
            raise ValueError('Please select a valid option.')

        if option=='3D': self._fig, self._fig_3D = self._fig_3D, self._fig
        super().savefig(name)        
        if option=='3D': self._fig, self._fig_3D = self._fig_3D, self._fig


###Extra plotting utilities###
def plot_from_file(fname, splitter=',', scatter=False, save=''):
    _x, _y = ([] for _ in range(2))
    with open(fname, 'r') as f:
        for values in f:
            values = values.split(splitter)
            if len(values) != 2:
                raise RuntimeError('The file must have two data columns.')
            _x.append(float(values[0]))
            _y.append(float(values[1]))
    if scatter:
        plt.scatter(np.array(_x), np.array(_y))
    else:
        plt.plot(np.array(_x), np.array(_y))
    if save!='':
        plt.savefig(save)
