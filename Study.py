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
profiles = []
h2 = h.get_halo(0)
with centering_com(h2):
    r200 = h2.properties['Halo_R_Crit200'].in_units('kpc')
    print(r200)
    profiles.append(pn.analysis.profile.Profile(h2.dm, ndim=3, 
                                                nbins=140, min=0.1, max=25))
    profiles.append(pn.analysis.profile.Profile(h2.dm, ndim=3, 
                                                calc_x=lambda x: x['r'], 
                                                nbins=140, min=0.1, max=25))
    profiles.append(pn.analysis.profile.Profile(h2.dm, ndim=3, 
                                                calc_x=lambda x: x['r']/r200, 
                                                nbins=140, min=0.1, max=1.3))
    profiles.append(pn.analysis.profile.Profile(h2.dm, ndim=3, 
                                                calc_x=lambda x: np.log10(x['r']/r200),
                                                nbins=140,min=-2,max=2))

    #This illustrates a manual log binning procedure,
    #since the one in pynbody does not seem to work.
    prof = pn.analysis.profile.Profile(h2.dm, ndim=3,
                                       bins=np.array([4*10**(i*0.078) for i in range(16)]))
    plt.plot(np.log10(prof['rbins']/r200),prof['density'])
    plt.savefig("QQQQQ.png")
    
p = plot.Profile(profiles, name="Density.png")
p.set_all_properties(model='density_profile', yscale='log')
p.plot_all("radius", "density")
#p.fit_and_plot_all('nfw')
p.savefig()

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
