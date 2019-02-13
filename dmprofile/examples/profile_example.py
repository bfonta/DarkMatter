from context import src
from src.halos import Halos
from src import plot
from src.move import centering_com
from src.utilities import intersect, write_to_file as wf, binning_using_binwidth
from src.parser import add_args

import matplotlib
import matplotlib.pyplot as plt

import pynbody as pn
import pynbody.plot.sph as sph
import pynbody.units as u

import os
import argparse
import numpy as np
from scipy.optimize import curve_fit

#Run with: python dmprofile/examples/profile_example.py --sym_types TYPE1 --sim_sizes 128

FLAGS, _ = add_args(argparse.ArgumentParser())
print("Parsed arguments:")
for k,v in FLAGS.__dict__.items():
    print('{}: {}'.format(k,v))

st = FLAGS.sim_types
addition = '' 
if FLAGS.sim_sizes==['128']: addition='.hdf5'
path_first = ['/fred/oz071/aduffy/Smaug/'+st[i]+'_L010N0'+FLAGS.sim_sizes[0] +'/data' for i in range(len(st))]
path1 = [os.path.join(path_first[i], 'subhalos_103/subhalo_103') for i in range(len(st))]
path2 = [os.path.join(path_first[i], 'snapshot_103/snap_103'+addition) for i in range(len(st))]

h = [Halos(path1[i], min_size=FLAGS.sim_min_particle_number) for i in range(len(st))]
N = [h[i].get_number_halos() for i in range(len(st))]

prof = []

bins1 = binning_using_binwidth(-1.25, 0, 0.078)
bins2 = binning_using_binwidth(-1.25, 0.2, 0.078)
halo_idx = 0
with centering_com(h[0].get_halo(halo_idx)):
    prof0 = h[0].get_profile(halo_idx, component='dm', bins=bins1,
                                bin_type='custom', normalize_x=True)
with centering_com(h[1].get_halo(halo_idx)):
    prof1 = h[1].get_profile(halo_idx, component='dm', bins=bins2,
                                bin_type='custom', normalize_x=True)

#prof0 will have index=0, prof1 will have index=1 and so on
p = plot.Profile([prof0, prof1], name='figs/DensityProfile.png')
p.set_all_properties(model='density_profile')

#current options are: radius, mass, enclosed mass, density and enclosed_density
p.plot_all(x_var='radius', y_var='density')

#fit all the profiles (NFW is currently the only option but the code can be easily expanded)
#p.fit_and_plot_all(function='nfw')

p.savefig()


if __name__ == 'main':
    unittest.main()
