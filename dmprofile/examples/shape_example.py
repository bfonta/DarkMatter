from context import src
from src.halos import Halos
from src.move import centering_com
from src import plot
from src.utilities import intersect

import matplotlib
import matplotlib.pyplot as plt

import pynbody as pn
import pynbody.plot.sph as sph
import pynbody.units as u

import os
import numpy as np
from scipy.optimize import curve_fit

HALO_NUMBER = 497 #after this the main halo has no subhalos
BIN_NUMBER = 100

DataFolder = "/fred/oz071/balves/"
SubhalosFolder = "Test_NOSN_NOZCOOL_L010N0128/data/subhalos_103/subhalo_103"
SnapshotFolder = "Test_NOSN_NOZCOOL_L010N0128/data/snapshot_103/snap_103.hdf5"

h = Halos(os.path.join(DataFolder,SubhalosFolder), HALO_NUMBER)
halos = h.get_halos()

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

if __name__ == 'main':
    unittest.main()
