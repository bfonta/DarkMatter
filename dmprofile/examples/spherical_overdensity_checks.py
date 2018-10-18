"""
Compare results with a "spherical overdensity" filter relative to results without that filter.
"""

from context import src

import os
import sys
import glob
import gzip
import argparse
import pynbody as pn

import dmprofile
from dmprofile.src.halos import Halos
from dmprofile.src import plot
from dmprofile.src.utilities import rho_crit, intersect
from dmprofile.src.move import centering_com
from dmprofile.src.parser import parser
from dmprofile.src.utilities import write_to_file as wf

path_first = '/fred/oz071/aduffy/Smaug/WTHERM_LateRe_L010N0512/data'
path1 = os.path.join(path_first, 'subhalos_103/subhalo_103')
path2 = os.path.join(path_first, 'snapshot_103/snap_103')

h1 = Halos(path1, path2, N=22, min_size=300)
N1 = h1.get_number_halos()
h1.filter_(halo_idxs=0, sub_idx=0, filter_str='Sphere_1.2', option='all', sim=True)

s1, par1a, par2a, M200, rel, res = ([] for i in range(6))
exceptions = []
ids = [i for i in range(N1) if i not in exceptions]
for i_halo in ids:
    _halo = h1.get_subhalo(i_halo)
    print("BYE!",_halo)
    with centering_com(_halo):
        m200 = h1.get_mass200(i_halo)
        res.append(h1.is_resolved(i_halo, intersect_value=0.7))
        rel.append(h1.is_relaxed(i_halo, 0))
        s_tmp = h1.get_shape(i_halo, 0)
        print(i_halo)
        print(s_tmp)
        if s_tmp[1]>1e-6 and s_tmp[2]>1e-6: #avoid infinities
            M200.append(np.log10(m200))
            par1a.append(s_tmp[1][0])
            par2a.append(s_tmp[2][0])
            s1.append(s_tmp)

shape = plot.Shape([s1, s1, s1], M200, name="figs/Shape_with_filter.png", w=11, h=10)
shape.set_axis((0,0), xlabel=r'$\log_{10}(M)$ [M$_{\odot}$]', ylabel='b/a',
               xscale='log', yscale='linear')
shape.set_axis((0,1), xlabel=r'$\log_{10}(M)$ [M$_{\odot}$]', ylabel='c/a',
               xscale='log', yscale='linear')
shape.set_axis((1,0), xlabel=r'$\log_{10}(M)$ [M$_{\odot}$]', ylabel='triaxiality',
               xscale='log', yscale='linear')
shape.set_axis((1,1), xlabel=r'$\log_{10}(M)$ [M$_{\odot}$]', ylabel='triaxiality',
               xscale='log', yscale='linear')

shape.binned_plot(0, axis_idx=(1,1), nbins=5, x_var="mass", y_var="triax",
                  color='b', label='Blue') 
shape.binned_plot(1, axis_idx=(1,1), nbins=5, x_var="mass", y_var="triax",
                  color='r', shift=0.04, label='Red') 
shape.binned_plot(2, axis_idx=(1,1), nbins=5, x_var="mass", y_var="triax",
                  color='g', shift=0.08, label='Green') 
shape.scatter_plot(0, axis_idx=(0,0), x_var="mass", y_var="b/a",
                   resolved_bools=res, relaxed_bools=rel)
shape.scatter_plot(1, axis_idx=(0,1), x_var="mass", y_var="c/a",
                   resolved_bools=res, relaxed_bools=rel)
shape.scatter_plot(2, axis_idx=(1,0), x_var="mass", y_var="triax",
                   resolved_bools=res, relaxed_bools=rel)
shape.savefig()

par3a = (1-np.array(par1a))**2 / (1-np.array(par2a))**2

h2 = Halos(path1, path2, N=22, min_size=300)
N2 = h2.get_number_halos()

s2, par1b, par2b, M200, rel, res = ([] for i in range(6))

exceptions = []
ids = [i for i in range(N2) if i not in exceptions]
for i_halo in ids:
    with centering_com(h2.get_halo(i_halo)):
        m200 = h2.get_mass200(i_halo)
        res.append(h2.is_resolved(i_halo, intersect_value=0.7))
        rel.append(h2.is_relaxed(i_halo, 0))
        s_tmp = h2.get_shape(i_halo, 0)
        print(i_halo)
        print(s_tmp)
        if s_tmp[1]>1e-6 and s_tmp[2]>1e-6: #avoid infinities
            M200.append(np.log10(m200))
            par1b.append(s_tmp[1][0])
            par2b.append(s_tmp[2][0])
            s2.append(s_tmp)

shape = plot.Shape([s2, s2, s2], M200, name="figs/Shape.png", w=11, h=10)
shape.set_axis((0,0), xlabel=r'$\log_{10}(M)$ [M$_{\odot}$]', ylabel='b/a',
               xscale='log', yscale='linear')
shape.set_axis((0,1), xlabel=r'$\log_{10}(M)$ [M$_{\odot}$]', ylabel='c/a',
               xscale='log', yscale='linear')
shape.set_axis((1,0), xlabel=r'$\log_{10}(M)$ [M$_{\odot}$]', ylabel='triaxiality',
               xscale='log', yscale='linear')
shape.set_axis((1,1), xlabel=r'$\log_{10}(M)$ [M$_{\odot}$]', ylabel='triaxiality',
               xscale='log', yscale='linear')

shape.binned_plot(0, axis_idx=(1,1), nbins=5, x_var="mass", y_var="triax",
                  color='b', label='Blue') 
shape.binned_plot(1, axis_idx=(1,1), nbins=5, x_var="mass", y_var="triax",
                  color='r', shift=0.04, label='Red') 
shape.binned_plot(2, axis_idx=(1,1), nbins=5, x_var="mass", y_var="triax",
                  color='g', shift=0.08, label='Green') 
shape.scatter_plot(0, axis_idx=(0,0), x_var="mass", y_var="b/a",
                   resolved_bools=res, relaxed_bools=rel)
shape.scatter_plot(1, axis_idx=(0,1), x_var="mass", y_var="c/a",
                   resolved_bools=res, relaxed_bools=rel)
shape.scatter_plot(2, axis_idx=(1,0), x_var="mass", y_var="triax",
                   resolved_bools=res, relaxed_bools=rel)
shape.savefig()

par3b = (1-np.array(par1b))**2 / (1-np.array(par2b))**2

wf("data/shape_double_1.txt",np.array(par1a),np.array(par1b))
wf("data/shape_double_2.txt",np.array(par2a),np.array(par2b))
wf("data/shape_double_3.txt",np.array(par3a),np.array(par3b))

"""
c1, M200_1, res, rel = ([] for i in range(4))
exceptions = [2,18]
ids = [i for i in range(N1) if i not in exceptions]
for i in ids:
    print("HALO:", i)
    with centering_com(h1.get_halo(i)):
        M200_1.append(h1.get_mass200(i))
        c1.append(h1.concentration_200(idx=i, sub_idx=0))
        res.append(h1.is_resolved(i))
        rel.append(h1.is_relaxed(i))

c_obj = plot.Concentration([c1, c1], extra_var=M200_1, name='figs/Concentration_with_filter.png')
c_obj.set_all_properties(model='concentration_mass')
c_obj.scatter_plot(0, (0,0), resolved_bools=res, relaxed_bools=rel)
c_obj.scatter_plot(1, (1,0), resolved_bools=res, relaxed_bools=rel)
c_obj.savefig()
"""

#del h1

"""
h2 = Halos(path1, path2, N=22, min_size=300)
N2 = h2.get_number_halos()

c2, M200_2, res, rel = ([] for i in range(4))
exceptions = [2,18]
ids = [i for i in range(N2) if i not in exceptions]
for i in ids:
    print("HALO:", i)
    with centering_com(h2.get_halo(i)):
        M200_2.append(h2.get_mass200(i))
        c2.append(h2.concentration_200(idx=i, sub_idx=0))
        res.append(h2.is_resolved(i))
        rel.append(h2.is_relaxed(i))

c_obj = plot.Concentration([c2, c2], extra_var=M200_2, name='figs/Concentration.png')
c_obj.set_all_properties(model='concentration_mass')
c_obj.scatter_plot(0, (0,0), resolved_bools=res, relaxed_bools=rel)
c_obj.scatter_plot(1, (1,0), resolved_bools=res, relaxed_bools=rel)
c_obj.savefig()

for i in range(len(M200_1)):
    print(M200_1[i], M200_2[i])
print("Before assert")
assert M200_1==M200_2
print("After assert")

wf("data/concentration_double.txt",np.array(c1),np.array(c2))
"""

if __name__ == 'main':
    unittest.main()
