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


HALO_NUMBER = 497 #after this the main halo has no subhalos
BIN_NUMBER = 40

DataFolder = "/fred/oz071/balves/"
SubhalosFolder = "Test_NOSN_NOZCOOL_L010N0128/data/subhalos_103/subhalo_103"
SnapshotFolder = "Test_NOSN_NOZCOOL_L010N0128/data/snapshot_103/snap_103.hdf5"

h = Halos(os.path.join(DataFolder,SubhalosFolder), HALO_NUMBER)

profiles, relax = ([] for i in range(2))
for i in range(2,5):
    profiles.append(h.get_profile(i,component='dm',bins=(10,30,BIN_NUMBER),bin_type='log',normalize=False))
    relax.append(h.relaxation(i, profiles[-1]))

p = plot.Profile(profiles, name="../../figs/Density.png")
p.set_all_properties(model='density_profile', xscale='log', yscale='log')
p.plot_all("radius", "density") 
p.fit_and_plot_all('nfw')                                                                                     

r = plot.Relaxation(relax, [profiles[i]['rbins'] for i in range(3)], h=16, w=16)
r.intersect_and_plot(0, (0,0), intersect_value=1.)
r.intersect_and_plot(1, (0,1), intersect_value=1.)
r.intersect_and_plot(2, (0,2), intersect_value=1.)
r.set_all_properties(model='relaxation', xscale='log', yscale='log')
r.savefig('Relaxation.png')
#p.draw_line(0, (0,0), intersect(prof1['rbins'], relax1, deg=7))
#p.draw_line(1, (1,0), intersect(prof2['rbins'], relax2, deg=7))
p.savefig()     

if __name__ == 'main':
    unittest.main()
