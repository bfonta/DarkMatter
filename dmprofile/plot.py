import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import itertools, collections

class Profile():
    """
    Plots all kinds of pynbody related profiles.

    Args:
    p: list of pynbody profile objects
    """
    def __init__(self, p, w=10, h=9, nrows=0, ncols=0, name=''):
        if type(w) is not int or type(h) is not int or type(nrows) is not int or type(ncols) is not int:
            raise TypeError
        if type(p) is not list:
            raise TypeError
        self.p = p
        self.N = len(self.p)
        self.name = name
        self.nrows = nrows
        self.ncols = ncols

        def geometry():
            r = np.floor(np.sqrt(self.N))
            c = r
            r_iter, c_iter = itertools.cycle(range(2)), itertools.cycle(range(2))
            next(r_iter)
            while r*c<self.N:
                r += next(r_iter)
                c += next(c_iter)
                print(r,c)
            print("Nrows ", r, " Ncols ", c)
            return int(r), int(c)

        if self.nrows==0 or self.ncols==0:
            self.nrows, self.ncols = geometry()
        self.fig, self.axis = plt.subplots(nrows=self.nrows, ncols=self.ncols, figsize=(w,h))

    def set_title(self, title):
        self.fig.suptitle(title)
        
    def set_name(self, name):
        self.name = name

    def set_axis(self, idx, xlabel, ylabel, xscale='linear', yscale='linear'):
        """
        idx can be either a single value or a tuple, depending on the total number of profiles
        """
        if self.N>2:
            indexes = idx[0], idx[1]
        else:
            indexes = idx
        self.axis[indexes].set_xscale(xscale)
        self.axis[indexes].set_yscale(yscale)
        self.axis[indexes].set_xlabel(xlabel)
        self.axis[indexes].set_ylabel(ylabel)

    def set_properties(self, title, xlabel, ylabel, xscale='linear', yscale='linear'):
        self.set_title(title)
        for i in range(self.nrows):
            if self.N>2:
                for j in range(self.ncols):
                    self.set_axis((i,j), xlabel, ylabel, xscale, yscale)
            else:
                self.set_axis(i, xlabel, ylabel, xscale, yscale)
    def plot(self, idx, axis_idx, x_var, y_var):
        """
        Args:
        1. 'idx':  profile index (int)
        2. 'axis_idx': axis index (int [self.N<=2] or tuple [self.N>2])
        3. 'x_var': variable stored in the profile to plot as x variable
            options: radius, mass, mass_enc, density
        4. 'y_var': variable stored in the profile to plot as y variable
            options: radius, mass, mass_enc, density
        """
        if self.N>2:
            indexes = axis_idx[0], axis_idx[1]
        else:
            indexes = axis_idx
        print(indexes)
        
        def _define_vars(x_var, y_var):
            _translate = {'radius': 'rbins',
                          'mass': 'mass',
                          'mass_enc': 'mass_enc',
                          'density': 'density'}
            if x_var not in _translate.keys() or y_var not in _translate.keys():
                raise ValueError
            return self.p[idx][_translate[x_var]], self.p[idx][_translate[y_var]]

        _x, _y = _define_vars(x_var, y_var)
        print(_y)
        self.axis[indexes].plot(_x, _y)

    def plot_all(self, x_var, y_var):
        i_profile = 0
        for i in range(self.nrows):
            if self.N>2:
                for j in range(self.ncols):
                    self.plot(i_profile, (i, j), x_var, y_var)
                    i_profile += 1
            else:
                self.plot(i_profile, i, x_var, y_var)
                i_profile += 1

    def savefig(self, name=''):
        if self.name=='':
            try:
                self.fig.savefig(name)
            except ValueError:
                print("A figure cannot be saved without having a name.")
        else:
            self.fig.savefig(self.name)
