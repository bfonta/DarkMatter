#QUESTIONS:
#the scale radius of the subhalos is the same as the scale radius of the halos
#do I center on the halo or instead on the main subhalo?

import os
import dmprofile
from dmprofile.src.halos import Halos
from dmprofile.src import plot
from dmprofile.src.utilities import rho_crit

print("################################################################")
print("################Code is now running#############################")
print("################################################################")

HALO_NUMBER = 497 #after this the main halo has no subhalos
BIN_NUMBER = 100

DataFolder = "/fred/oz071/balves/"
SubhalosFolder = "Test_NOSN_NOZCOOL_L010N0128/data/subhalos_103/subhalo_103"
SnapshotFolder = "Test_NOSN_NOZCOOL_L010N0128/data/snapshot_103/snap_103.hdf5"

h = Halos(os.path.join(DataFolder,SubhalosFolder), HALO_NUMBER)
h.filter_(halo_idx, filter_str='Sphere_1')

prof1 = h.get_profile(0, component='dm', bins=(2.5,30,BIN_NUMBER), bin_type='log', normalize=False)


"""
p = plot.Profile(profiles, name="figs/test.png", h=9, w=8)
p.set_all_properties(model='density_profile', xscale='log', yscale='log')
p.plot_all("radius", "density")
#p.draw_line((0,0), rhoz*200, 'h', label=r'$\rho_{crit}$', color='red')
#p.draw_line((0,0), r200_1, 'v', label=r'r$_{200}$')
#p.draw_line((1,0), rhoz*200, 'h',label=r'$\rho_{crit}$', color='red')
#p.draw_line((1,0), r200_2, 'v',label=r'r$_{200}$')
p.fit_and_plot_all('nfw')
p.savefig()
"""
