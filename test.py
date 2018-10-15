#QUESTIONS:
#the scale radius of the subhalos is the same as the scale radius of the halos
#do I center on the halo or instead on the main subhalo?
#check halos.is_resolved(). What happens when no intersection is found?

import os
import sys
import glob
import gzip
import argparse
import pynbody as pn

import numpy as np

import dmprofile
from dmprofile.src.halos import Halos
from dmprofile.src import plot
from dmprofile.src.utilities import rho_crit, intersect
from dmprofile.src.move import centering_com
from dmprofile.src.parser import parser

print("################################################################")
print("################Code is now running#############################")
print("################################################################")

BIN_NUMBER = 20

path_first = '/fred/oz071/aduffy/Smaug/DMONLY_L010N0256/data'
path1 = os.path.join(path_first, 'subhalos_103/subhalo_103')
path2 = os.path.join(path_first, 'snapshot_103/snap_103')

h = Halos(path1, path2, N=10, min_size=300)
N = h.get_number_halos()
h.filter_(halo_idxs=0, sub_idx=0, filter_str='Sphere_1.2', option='all', sim=True)

"""
s, M200, rel, res = ([] for i in range(4))
#exceptions = [1,2,15,61,170,173]
exceptions = [1,2,15,61,170,173]
ids = [i for i in range(N) if i not in exceptions]
for i_halo in ids:
    with centering_com(h.get_halo(i_halo)):
        m200 = h.get_mass200(i_halo)
        res.append(h.is_resolved(i_halo))
        rel.append(h.is_relaxed(i_halo))
        s_tmp = h.get_shape(i_halo)
        if s_tmp[1]>1e-6 and s_tmp[2]>1e-6: #avoid infinities
            M200.append(np.log10(m200))
            s.append(s_tmp)
shape = plot.Shape([s, s, s], M200, name="figs/Shape.png", w=11, h=10)
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
shape.binned_plot(0, axis_idx=(1,0), nbins=5, x_var="mass", y_var="triax",
                  color='b') 
shape.savefig()
"""

profiles, relax, rbins = ([] for i in range(3))
for i in range(2,5):
    profiles.append(h.get_profile(i,component='dm',bins=(2.,15,BIN_NUMBER),bin_type='log',normalize=False))
    _rel = h.relaxation(i, profiles[-1])
    relax.append(_rel[0])
    rbins.append(_rel[1])

p = plot.Profile(profiles, name="figs/Density.png")
p.set_all_properties(model='density_profile', xscale='log', yscale='log')
p.plot_all("radius", "density")
p.fit_and_plot_all('nfw')
p.savefig()

r = plot.Relaxation(relax, rbins, h=16, w=16)
r.intersect_and_plot(0, (0,0), intersect_value=1.)
r.intersect_and_plot(1, (0,1), intersect_value=1.)
r.intersect_and_plot(2, (1,1), intersect_value=1.)
r.set_all_properties(model='relaxation', xscale='log', yscale='log')
r.savefig('figs/Relaxation.png')
