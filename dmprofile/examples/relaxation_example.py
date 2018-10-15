from context import src
from src.halos import Halos
from src.move import centering_com
from src import plot
from src.utilities import intersect

import matplotlib
import matplotlib.pyplot as plt

import pynbody as pn
import pynbody.plot.sph as sph
import pynbody.units as u

import os
import numpy as np
from scipy.optimize import curve_fit


BIN_NUMBER = 20

path_first = '/fred/oz071/aduffy/Smaug/DMONLY_L010N0256/data'
path1 = os.path.join(path_first, 'subhalos_103/subhalo_103')
path2 = os.path.join(path_first, 'snapshot_103/snap_103')

h = Halos(path1, path2, N=10, min_size=300)
N = h.get_number_halos()
h.filter_(halo_idxs=0, sub_idx=0, filter_str='Sphere_1.2', option='all', sim=True)

profiles, relax, rbins = ([] for i in range(3))
for i in range(2,5):
    profiles.append(h.get_profile(i,component='dm',bins=(3.,30,BIN_NUMBER),bin_type='log',normalize=False))
    _rel = h.relaxation(i, profiles[-1])
    relax.append(_rel[0])
    rbins.append(_rel[1])

p = plot.Profile(profiles, name="figs/Density.png")
p.set_all_properties(model='density_profile', xscale='log', yscale='log')
p.plot_all("radius", "density")
p.fit_and_plot_all('nfw')
p.savefig()

r = plot.Relaxation(relax, rbins, h=16, w=16)
r.intersect_and_plot(0, (0,0), intersect_value=1.)
r.intersect_and_plot(1, (0,1), intersect_value=1.)
r.intersect_and_plot(2, (1,1), intersect_value=1.)
r.set_all_properties(model='relaxation', xscale='log', yscale='log')
r.savefig('figs/Relaxation.png') 


if __name__ == 'main':
    unittest.main()
