#QUESTIONS:
#the scale radius of the subhalos is the same as the scale radius of the halos
#do I center on the halo or instead on the main subhalo?
#check halos.is_resolved(). What happens when no intersection is found?

import os
import sys
import glob
import gzip
import argparse

import numpy as np

import dmprofile
from dmprofile.src.halos import Halos
from dmprofile.src import plot
from dmprofile.src.utilities import rho_crit, intersect
from dmprofile.src.move import centering_com, centering_mbp
from dmprofile.src.parser import add_args
from dmprofile.src.utilities import write_to_file as wf

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
    c, M200, M200_shape, res, rel, s = ([] for i in range(6))
    #2,6,14,25,27,46,54,91,116,117,122]
    for i in range(N[isim]):
        print("SIM", st[isim], "  HALO:", i)
        with centering_com(h[isim].get_halo(i)):
            print(h[isim].get_halo(i))
            isres = h[isim].is_resolved(i, sub_idx=0)        
            isrel = h[isim].is_relaxed(i, sub_idx=0)
            relax_tmp = h[isim].concentration_200(idx=i, sub_idx=0)
            s_tmp = h[isim].get_shape(i, 0)
            if s_tmp!=-1 and isres!=-1 and relax_tmp!=-1 and isrel!=-1:
                M200.append(h[isim].get_mass200(i))
                res.append(isres)     
                rel.append(isrel)
                c.append(relax_tmp)
                M200_shape.append(np.log10(h[isim].get_mass200(i)))
                s.append(s_tmp)
    wf('data/Concentration_'+st[isim]+'_'+str(FLAGS.sim_min_particle_number)+
       '_'+FLAGS.sim_size+'.txt', c, M200, res, rel)
    wf('data/Shape_'+st[isim]+'_'+str(FLAGS.sim_min_particle_number)+
       '_'+FLAGS.sim_size+'.txt', s, M200_shape, mode='shape')
