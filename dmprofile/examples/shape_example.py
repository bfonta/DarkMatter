from context import src
from src.halos import Halos
from src.move import centering_com
from src import plot
from src.utilities import intersect
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

FLAGS, _ = add_args(argparse.ArgumentParser())
print("Parsed arguments:")
for k,v in FLAGS.__dict__.items():
    print('{}: {}'.format(k,v))

st = FLAGS.sim_types
addition = '' 
if FLAGS.sim_size=='128': additon='.hdf5'
path_first = ['/fred/oz071/aduffy/Smaug/'+st[i]+'_L010N0'+FLAGS.sim_size+'/data' for i in range(len(st))]
path1 = [os.path.join(path_first[i], 'subhalos_103/subhalo_103') for i in range(len(st))]
path2 = [os.path.join(path_first[i], 'snapshot_103/snap_103'+addition) for i in range(len(st))]

h = [Halos(path1[i], min_size=FLAGS.sim_min_particle_number) for i in range(len(st))]
N = [h[i].get_number_halos() for i in range(len(st))]

for isim in range(len(st)):
    M200_shape, s = ([] for i in range(2))
    for i in range(N[isim]):
        with centering_com(h[isim].get_halo(i)):
            print(h[isim].get_halo(i))
            isres = h[isim].is_resolved(i, sub_idx=0)        
            isrel = h[isim].is_relaxed(i, sub_idx=0)
            relax_tmp = h[isim].concentration_200(idx=i, sub_idx=0)
            s_tmp = h[isim].get_shape(i, 0)
            if s_tmp!=-1 and isres!=-1 and relax_tmp!=-1 and isrel!=-1:
                M200_shape.append(np.log10(h[isim].get_mass200(i)))
                s.append(s_tmp)
    wf('file.txt', s, M200_shape, mode='shape')

"""
shape = plot.Shape([s, s, s], M200, name="../../figs/Shape.png", w=16, h=13)
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

if __name__ == 'main':
    unittest.main()
