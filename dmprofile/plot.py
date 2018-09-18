import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import itertools, collections

class _Plot():
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
        print("Nrows ", _r, " Ncols ", _c)
        return int(_r), int(_c)
        
    def _set_axis_indexes(self, axis_idx):
        """
        Set the axis indexes into 1D or 2D, due to matplotlib axis indexing.
        1D is needed when one plans to plot less than two pictures.
        """
        if self._N>2:
            return axis_idx[0], axis_idx[1]
        else:
            return axis_idx

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
            indexes = idx
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
                self.set_axis(i, xlabel, ylabel, xscale, yscale)

    def set_all_properties(self, title="", model=None, 
                       xlabel='R [kpc]', ylabel=None, 
                       xscale='linear', yscale='linear'):
        """
        model options: 'density_profile', 'mass_enc_profile', 'mass_profile'
        """
        self.set_title(title)
        model_options = ['density_profile', 
                         'mass_enc_profile', 
                         'mass_profile']
        if model not in model_options:
            raise ValueError('The specified model is not currently predefined.')
        if model!='None':
            self.set_axis_all('R [kpc]', r'$\rho$ [M$_{\odot}$ kpc$^{-3}$]', xscale, yscale) if model=='density_profile' else self.set_axis_all('R [kpc]', r'M$_{enc}$ [M$_{\odot}$]', xscale, yscale) if model=='mass_enc_profile' else self.set_axis_all('R [kpc]', r'M [M$_{\odot}$]', xscale, yscale)
        else:
            for i in range(self._nrows):
                if self._N>2:
                    for j in range(self._ncols):
                        self.set_axis((i,j), xlabel, ylabel, xscale, yscale)
                else:
                    self.set_axis(i, xlabel, ylabel, xscale, yscale)
    
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
            _translate = {'radius': 'rbins',
                          'mass': 'mass',
                          'mass_enc': 'mass_enc',
                          'density': 'density'}
            if x_var not in _translate.keys() or y_var not in _translate.keys():
                raise ValueError
            print(idx)
            return self._p[idx][_translate[x_var]], self._p[idx][_translate[y_var]]

        _x, _y = _define_vars(x_var, y_var)
        indexes = super()._set_axis_indexes(axis_idx)
        self._axis[indexes].plot(_x, _y)

    def plot_all(self, x_var, y_var):
        i_profile = 0
        for i in range(self._nrows):
            if self._N>2:
                for j in range(self._ncols):
                    if i_profile<self._N:
                        self.plot(i_profile, (i, j), x_var, y_var)
                        i_profile += 1
            else:
                self.plot(i_profile, i, x_var, y_var)
                i_profile += 1


class Shape(_Plot):
    """
    Plots all kinds of pynbody related shapes.

    Args:
    s: list of shape arrays for plotting, i.e., list of lists of pynbody shape objects
    extra: list of some additional property common to all objects from which the shapes were taken
           for example: radius, mass
    """
    def __init__(self, s, extra=[], w=10, h=9, nrows=0, ncols=0, name=''):
        super().__init__(s, w, h, nrows, ncols, name)
        self._s = self._obj_list
        self._extra = extra
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

    def scatter_plot(self, idx, axis_idx, x_var, y_var):
        """
        Args:
        1. 'idx': pynbody shape list object. It refers to the way self._s is ordered
        2. 'axis_idx': axis index (int [self._N<=2] or tuple [self._N>2])
        3. 'x_var': variable stored in the shape to plot as x variable
            options: pos, b/a, c/a, align, rot, triax, m200, r200
        4. 'y_var': variable stored in the profile to plot as y variable
            options: pos, b/a, c/a, align, rot, triax, m200, r200
        """
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
                    return self._extra, _s_for_plot[idx][_match[y_var]]
                else:
                    return _s_for_plot[idx][_match[x_var]], self._extra

        _x, _y = _define_vars(x_var, y_var)
        _sorted_lists = sorted(zip(_x,_y), key=lambda x: x[0])
        _x, _y = [[q[i] for q in _sorted_lists] for i in range(2)]
        indexes = super()._set_axis_indexes(axis_idx)
        self._axis[indexes].scatter(_x, _y)
