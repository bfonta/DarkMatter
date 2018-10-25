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
from dmprofile.src.move import centering_com, centering_mbp
from dmprofile.src.parser import parser
from dmprofile.src.utilities import write_to_file as wf

print("################################################################")
print("################Code is now running#############################")
print("################################################################")

size = '256'
addition = ''
if size=='128': additon='.hdf5'
path_first = '/fred/oz071/aduffy/Smaug/WTHERM_LateRe_L010N0'+size+'/data'
path1 = os.path.join(path_first, 'subhalos_103/subhalo_103')
path2 = os.path.join(path_first, 'snapshot_103/snap_103'+addition)

min_size1=1000
h1 = Halos(path1, min_size=min_size1)
N1 = h1.get_number_halos()

c1, M200_1, M200_shape_1, res1, rel1, s1 = ([] for i in range(6))
#2,6,14,25,27,46,54,91,116,117,122]
for i in range(N1):
    print("HALO:", i)
    with centering_com(h1.get_halo(i)):
        print(h1.get_halo(i))
        isres = h1.is_resolved(i, sub_idx=0)        
        isrel = h1.is_relaxed(i, sub_idx=0)
        relax_tmp = h1.concentration_200(idx=i, sub_idx=0)
        s_tmp = h1.get_shape(i, 0)
        if s_tmp!=-1 and isres!=-1 and relax_tmp!=-1 and isrel!=-1:
            M200_1.append(h1.get_mass200(i))
            res1.append(isres)     
            rel1.append(isrel)
            c1.append(relax_tmp)
            M200_shape_1.append(np.log10(h1.get_mass200(i)))
            s1.append(s_tmp)

#########
min_size2=500
h2 = Halos(path1, min_size=min_size2)
N2 = h2.get_number_halos()

c2, M200_2, M200_shape_2, res2, rel2, s2 = ([] for i in range(6))
for i in range(N2):
    print("HALO:", i)
    with centering_com(h2.get_halo(i)):
        isres = h2.is_resolved(i, sub_idx=0)
        isrel = h2.is_relaxed(i, sub_idx=0)
        relax_tmp = h2.concentration_200(idx=i, sub_idx=0)
        s_tmp = h2.get_shape(i, 0)
        if isres!=-1 and relax_tmp!=-1 and isrel!=-1 and s_tmp!=-1:
            M200_2.append(h2.get_mass200(i))
            res2.append(isres)     
            rel2.append(isrel)
            c2.append(relax_tmp)
            M200_shape_2.append(np.log10(h2.get_mass200(i)))
            s2.append(s_tmp)

file_w1 = len(s1)
file_w2 = len(s2)
wf('data/Concentration_'+str(min_size1)+'_'+size+'.txt', c1, M200_1, res1, rel1)
wf('data/Concentration_'+str(min_size2)+'_'+size+'.txt', c2, M200_2, res2, rel2)
wf('data/Shape_'+str(min_size1)+'_'+size+'.txt', s1, M200_shape_1, mode='shape')
wf('data/Shape_'+str(min_size2)+'_'+size+'.txt', s2, M200_shape_2, mode='shape')
