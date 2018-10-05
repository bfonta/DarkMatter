#QUESTIONS:
#the scale radius of the subhalos is the same as the scale radius of the halos
#do I center on the halo or instead on the main subhalo?

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import pynbody as pn
import pynbody.plot.sph as sph
import pynbody.units as u
import numpy as np
from scipy.optimize import curve_fit

import os
from dmprofile.halos import Halos
from dmprofile.move import centering_com
from dmprofile import plot
from dmprofile.utilities import intersect

print("################################################################")
print("################Code is now running#############################")
print("################################################################")

HALO_NUMBER = 497 #after this the main halo has no subhalos
BIN_NUMBER = 500

DataFolder = "/fred/oz071/balves/"
SubhalosFolder = "Test_NOSN_NOZCOOL_L010N0128/data/subhalos_103/subhalo_103"
SnapshotFolder = "Test_NOSN_NOZCOOL_L010N0128/data/snapshot_103/snap_103.hdf5"

h = Halos(os.path.join(DataFolder,SubhalosFolder), HALO_NUMBER)

#Filtering
_idx = 0
r200 = h.get_halo(_idx).properties['Halo_R_Crit200'].in_units('kpc')
halo = h.get_halo(_idx)
centering_com(halo)
subhalo = h.get_subhalo(_idx)
print(len(subhalo))
h.filter_(_idx, 0, filter_str='Sphere_1.2')
subhalo = h.get_subhalo(_idx)
print(len(subhalo))
prof1 = h.get_profile(_idx, 0, 'dm', bins=(1,45,BIN_NUMBER), bin_type='linear', normalize=False)
h.filter_(_idx, 0, filter_str='Sphere_1')
subhalo_f = h.get_subhalo(_idx)
print(len(subhalo_f))
prof2 = h.get_profile(_idx, 0, 'dm', bins=(1,45,BIN_NUMBER), bin_type='linear', normalize=False)

rhoz = pn.array.SimArray(pn.analysis.cosmology.rho_crit(halo, unit="Msol kpc**-3"), "Msol kpc**-3")
profiles = [prof1,prof2]

p = plot.Profile(profiles, name="Density_dm.png")
p.set_all_properties(model='density_enc_profile', xscale='linear', yscale='log')
p.plot_all("radius", "density_enc")
p.draw_line((0,0), rhoz*200, 'h', label=r'200$\rho_{crit}$', color='red')
p.draw_line((0,0), r200, 'v', label=r'r$_{200}$')
p.draw_line((1,0), rhoz*200, 'h',label=r'200$\rho_{crit}$', color='red')
p.draw_line((1,0), r200, 'v',label=r'r$_{200}$')
#p.fit_and_plot_all('nfw')
p.savefig()

"""
x, y, z = halo['x'], halo['y'], halo['z']
xx, yy, zz = subhalo_f['x'], subhalo_f['y'], subhalo_f['z']
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
plt.savefig("Filtering_xy.png")

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
plt.savefig("Filtering_xz.png")

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
plt.savefig("3D.png")
"""

#Concentration
"""
c, M200, res, rel = ([] for i in range(4))

exceptions = [1,11,15,16,18,19,20,22,24]
iterable = iter([item for item in range(1,25) if item not in exceptions])

for i in iterable:
    print("NUMBER:_____________________ ", i)
    M200.append(halos[i].properties['Halo_M_Crit200'].in_units('Msol'))
    c.append(h.concentration_200(idx=i, sub=True, sub_idx=0))
    res.append(h.is_resolved(i))
    rel.append(h.is_relaxed(i))
c_obj = plot.Concentration([c, c], extra_var=M200, name='Concentration.png') 
c_obj.set_all_properties(model='concentration_mass')
c_obj.scatter_plot(0, (0,0), resolved_bools=res, relaxed_bools=rel)
c_obj.scatter_plot(1, (1,0))
c_obj.savefig()
"""

#Density
"""
profiles, relax = ([] for i in range(2))
for i in range(2,5):
    profiles.append(h.get_profile(i, 'dm', bins=(2.5,30,BIN_NUMBER), bin_type='log', normalize=False))
    relax.append(h.relaxation(i, profiles[-1]))
p = plot.Profile(profiles, name="Density.png")
p.set_all_properties(model='density_profile', xscale='log', yscale='log')
p.plot_all("radius", "density")
p.fit_and_plot_all('nfw')
                 
r = plot.Relaxation(relax, [profiles[i]['rbins'] for i in range(8)], h=16, w=16)
r.intersect_and_plot(0, (0,0), intersect_value=1.)
r.intersect_and_plot(1, (0,1), intersect_value=1.)
r.intersect_and_plot(2, (0,2), intersect_value=1.)
r.set_all_properties(model='relaxation', xscale='log', yscale='log')
r.savefig('Relaxation.png')

#p.draw_line(0, (0,0), intersect(prof1['rbins'], relax1, deg=7))
#p.draw_line(1, (1,0), intersect(prof2['rbins'], relax2, deg=7))
p.savefig()
"""

#This illustrates a manual log binning procedure,
#np.array([4*10**(i*0.078) for i in range(16)])

#Shape
"""
s, M200, rel, res = ([] for i in range(4))
exceptions = [1,11,15,16,18,19,20,22,24]
iterable = iter([item for item in range(1,25) if item not in exceptions])

for i_halo in iterable:
    with centering_com(halos[i_halo]):
        r200 = halos[i_halo].properties['Halo_R_Crit200'].in_units('kpc')
        m200 = halos[i_halo].properties['Halo_M_Crit200'].in_units('Msol')
        res.append(h.is_resolved(i_halo))
        rel.append(h.is_relaxed(i_halo))
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
shape.scatter_plot(0, axis_idx=(0,0), x_var="mass", y_var="b/a", 
                   resolved_bools=res, relaxed_bools=rel)
shape.scatter_plot(1, axis_idx=(0,1), x_var="mass", y_var="c/a",
                   resolved_bools=res, relaxed_bools=rel)
shape.scatter_plot(2, axis_idx=(1,1), x_var="mass", y_var="triax",
                   resolved_bools=res, relaxed_bools=rel)
shape.savefig()
"""
