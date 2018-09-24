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
BIN_NUMBER = 100

DataFolder = "/fred/oz071/balves/"
SubhalosFolder = "Test_NOSN_NOZCOOL_L010N0128/data/subhalos_103/subhalo_103"
SnapshotFolder = "Test_NOSN_NOZCOOL_L010N0128/data/snapshot_103/snap_103"

h = Halos(os.path.join(DataFolder,SubhalosFolder), HALO_NUMBER)
halos = h.get_halos()

#Concentration
"""
c, M200 = ([] for i in range(2))
for i in range(2,22):
    M200.append(halos[i].properties['Halo_M_Crit200'].in_units('Msol'))
    c.append(h.concentration_200(i))
c_obj = plot.Concentration([c, c], extra_var=M200, name='Concentration.png') 
c_obj.set_all_properties(model='concentration_mass')
c_obj.scatter_plot(0, (0,0))
c_obj.scatter_plot(1, (1,0))
c_obj.savefig()
"""

#Density
prof1 = h.get_profile(0, 'dm', bins=(3.5,25.,BIN_NUMBER), bin_type='log', normalize=False)
prof2 = h.get_profile(0, 'dm', bins=(3.5,25.,BIN_NUMBER))            
p = plot.Profile([prof1,prof2], name="Density.png")
p.set_all_properties(model='density_profile')#, yscale='log')
p.plot_all("radius", "density")
p.fit_and_plot_all('nfw')
p.savefig()

relax = h.relaxation(0,prof2)
r = plot.Relaxation([relax,relax],
                    [prof2['rbins'],prof2['rbins']])

r.intersect_and_plot(0, (0,0), intersect_value=50)
r.plot(1, (1,0))
r.savefig('Relaxation.png')

#This illustrates a manual log binning procedure,
#np.array([4*10**(i*0.078) for i in range(16)])

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
