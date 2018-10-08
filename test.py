#QUESTIONS:
#the scale radius of the subhalos is the same as the scale radius of the halos
#do I center on the halo or instead on the main subhalo?

import os
import pynbody as pn
import dmprofile
from dmprofile.src.halos import Halos
from dmprofile.src import plot
from dmprofile.src.utilities import rho_crit
from dmprofile.src.move import centering_com

print("################################################################")
print("################Code is now running#############################")
print("################################################################")

HALO_NUMBER = 497 #after this the main halo has no subhalos
BIN_NUMBER = 100

DataFolder = "/fred/oz071/balves/"
SubhalosFolder = "Test_NOSN_NOZCOOL_L010N0128/data/subhalos_103/subhalo_103"
SnapshotFolder = "Test_NOSN_NOZCOOL_L010N0128/data/snapshot_103/snap_103.hdf5"

sim = pn.load(os.path.join(DataFolder,SnapshotFolder))
sim['eps'] = 200.*pn.units.pc
sim.physical_units()

h = Halos(os.path.join(DataFolder,SubhalosFolder), HALO_NUMBER)
halo = h.get_halo(0)
r200 = halo.properties['Halo_R_Crit200'].in_units('kpc')
center = pn.analysis.halo.center_of_mass(halo)
sim_f = sim[pn.filt.Sphere(r200, center)]
centering_com(sim_f)
center2 = pn.analysis.halo.center_of_mass(sim_f)
sim_f2 = sim_f[pn.filt.Sphere(r200/3, center2)]

plot = plot.Particles([sim_f, sim_f2], w=6, h=8)
plot.scatter_plot_3D(0, (0,0), zmin=-35, zmax=35, xmin=-35, xmax=35, ymin=-35, ymax=35)
plot.scatter_plot_3D(1, (1,0), zmin=-35, zmax=35, xmin=-35, xmax=35, ymin=-35, ymax=35)
plot.savefig('figs/3D.png', '3D')
