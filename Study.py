import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pynbody as pn
import pynbody.plot.sph as sph
import pynbody.units as u
import numpy as np
from scipy.optimize import curve_fit

import os
from dmprofile.halos import Halos
from dmprofile.move import centering_com
from dmprofile import plot

print("################################################################")
print("################Code is now running#############################")
print("################################################################")

HALO_NUMBER = 500

DataFolder = "/fred/oz071/balves/"
SubhalosFolder = "Test_NOSN_NOZCOOL_L010N0128/data/subhalos_103/subhalo_103"
SnapshotFolder = "Test_NOSN_NOZCOOL_L010N0128/data/snapshot_103/snap_103"

h = Halos(os.path.join(DataFolder,SubhalosFolder), HALO_NUMBER)

#Concentration
c = h.concentration_200(2)
print(c)

#Density
""""
profiles = []
for i in range(2):
    h0 = h.get_halo(2)
    with Centering(h0).com():
        profiles.append(pn.analysis.profile.Profile(h0.dm, ndim=3, nbins=250, min=1, max=35))

p = plot.Profile(profiles, name="Density.png")
p.set_all_properties(title='Density', model='density_profile', yscale='log')
p.plot_all("radius", "density")
p.fit_and_plot_all('nfw')
p.savefig()
"""

#Shape
"""
s, M200 = ([] for i in range(2))
halos = h.get_halos()
for i_halo in range(HALO_NUMBER):
    with centering_com(halos[i_halo]):
        r200 = halos[i_halo].properties['Halo_R_Crit200'].in_units('kpc')
        m200 = halos[i_halo].properties['Halo_M_Crit200'].in_units('Msol')
        if i_halo%100==0: print("Halo number: ", i_halo)
        s_tmp = pn.analysis.halo.halo_shape(halos[i_halo], N=1, rout=r200, bins='lin')
        if s_tmp[1]>1e-6 and s_tmp[2]>1e-6: #avoid infinities
            M200.append(np.log10(m200))
            s.append(s_tmp)


shape = plot.Shape([s, s, s], M200, name="Shape.png", w=16, h=13)
shape.set_axis((0,0), xlabel=r'$\log_{10}(M)$ [M$_{\odot}$]', ylabel='b/a',
               xscale='log', yscale='linear')
shape.set_axis((0,1), xlabel=r'$\log_{10}(M)$ [M$_{\odot}$]', ylabel='c/a',
               xscale='log', yscale='linear') 
shape.set_axis((1,1), xlabel=r'$\log_{10}(M)$ [M$_{\odot}$]', ylabel='triaxiality',
               xscale='log', yscale='linear') 
shape.scatter_plot(0, axis_idx=(0,0), x_var="mass", y_var="b/a")
shape.scatter_plot(1, axis_idx=(0,1), x_var="mass", y_var="c/a")
shape.scatter_plot(2, axis_idx=(1,1), x_var="mass", y_var="triax")
shape.savefig()
"""
