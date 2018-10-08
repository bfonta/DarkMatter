from context import src
import unittest
from src.halos import Halos
from src.move import centering_com

import os
import numpy as np
from scipy.optimize import curve_fit

import pynbody as pn
import pynbody.plot.sph as sph
import pynbody.units as u

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

HALO_NUMBER = 497 #after this the main halo has no subhalos
BIN_NUMBER = 100

DataFolder = "/fred/oz071/balves/"
SubhalosFolder = "Test_NOSN_NOZCOOL_L010N0128/data/subhalos_103/subhalo_103"
SnapshotFolder = "Test_NOSN_NOZCOOL_L010N0128/data/snapshot_103/snap_103.hdf5"

h = Halos(os.path.join(DataFolder,SubhalosFolder), HALO_NUMBER)

halo_idx = 0
halo = h.get_halo(halo_idx)
centering_com(halo)
subhalo = h.get_subhalo(halo_idx)
h.filter_(halo_idx, sub_idx=0, filter_str='Sphere_1')
halo_f = h.get_subhalo(halo_idx)
#prof1 = h.get_profile(2, component='dm', bins=(2.5,30,BIN_NUMBER), bin_type='log', normalize=False)
h.filter_(halo_idx, sub_idx=0, filter_str='Sphere_0.4')
halo_f2 = h.get_subhalo(halo_idx)
#prof2 = h.get_profile(2, component='dm', bins=(10,20,30), bin_type='log', normalize=False)

#plotting
x, y, z = halo_f['x'], halo_f['y'], halo_f['z']
xx, yy, zz = halo_f2['x'], halo_f2['y'], halo_f2['z']
x_s, y_s, z_s = subhalo['x'], subhalo['y'], subhalo['z']
x_h, y_h, z_h = halo['x'], halo['y'], halo['z']

fig, axis = plt.subplots(nrows=2, ncols=2)
axis[0,0].plot(x, y, marker='.', linestyle='None', label='Sphere filter', color='b')
axis[0,0].legend()
axis[0,0].set_xlim([-35,35])
axis[0,0].set_ylim([-35,35])
axis[0,1].plot(x_s, y_s, marker='.', linestyle='None', label='Main subhalo', color='g')
axis[0,1].legend()
axis[0,1].set_xlim([-35,35])
axis[0,1].set_ylim([-35,35])
axis[1,0].plot(x_h, y_h, marker='.', linestyle='None', label='Halo', color='orange')
axis[1,0].legend()
axis[1,0].set_xlim([-35,35])
axis[1,0].set_ylim([-35,35])
axis[1,1].plot(xx, yy, marker='.', linestyle='None', label='Smaller sphere', color='brown')
axis[1,1].legend()
axis[1,1].set_xlim([-35,35])
axis[1,1].set_ylim([-35,35])
plt.savefig("figs/Filtering_xy.png")

fig, axis = plt.subplots(nrows=2, ncols=2)
axis[0,0].plot(x, z, marker='.', linestyle='None', label='Sphere filter', color='b')
axis[0,0].legend()
axis[0,0].set_xlim([-50,50])
axis[0,0].set_ylim([-50,50])
axis[0,1].plot(x_s, z_s, marker='.', linestyle='None', label='Main subhalo', color='g')
axis[0,1].legend()
axis[0,1].set_xlim([-50,50])
axis[0,1].set_ylim([-50,50])
axis[1,0].plot(x_h, z_h, marker='.', linestyle='None', label='Halo', color='orange')
axis[1,0].legend()
axis[1,0].set_xlim([-50,50])
axis[1,0].set_ylim([-50,50])
axis[1,1].plot(xx, zz, marker='.', linestyle='None', label='Smaller sphere', color='brown')
axis[1,1].legend()
axis[1,1].set_xlim([-50,50])
axis[1,1].set_ylim([-50,50])
plt.savefig("figs/Filtering_xz.png")

fig = plt.figure(figsize=(13,11))
axis = fig.add_subplot(2, 2, 1, projection='3d')
axis.scatter(x, y, z, zdir='z', label="Sphere filter", color='b')
axis.legend()
axis.set_xlim([-40,40])
axis.set_ylim([-40,40])
axis.set_zlim([-40,40])
axis = fig.add_subplot(2, 2, 2, projection='3d')
axis.scatter(x_s, y_s, z_s, zdir='z', label="Main subhalo", color='g')
axis.legend()
axis.set_xlim([-40,40])
axis.set_ylim([-40,40])
axis.set_zlim([-40,40])
axis = fig.add_subplot(2, 2, 3, projection='3d')
axis.scatter(x_h, y_h, z_h, zdir='z', label="Halo", color='orange')
axis.legend()
axis.set_xlim([-40,40])
axis.set_ylim([-40,40])
axis.set_zlim([-40,40])
axis = fig.add_subplot(2, 2, 4, projection='3d')
axis.scatter(xx, yy, zz, zdir='z', label="Smaller sphere", color='brown')
axis.legend()
axis.set_xlim([-40,40])
axis.set_ylim([-40,40])
axis.set_zlim([-40,40])
plt.savefig("figs/3D_2.png")

if __name__ == 'main':
    unittest.main()
