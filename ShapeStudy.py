import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pynbody as pn
import pynbody.plot.sph as sph
import pynbody.units as u
import numpy as np
from scipy.optimize import curve_fit

import os

print("################################################################")
print("################Code is now running#############################")
print("################################################################")

DataFolder = "/fred/oz071/balves/"
SnapshotFolder = "Test_NOSN_NOZCOOL_L010N0128/data/subhalos_103/subhalo_103"
SubhalosFolder = "Test_NOSN_NOZCOOL_L010N0128/data/snapshot_103/snap_103"

s = pn.load(os.path.join(DataFolder,SnapshotFolder))
s['eps'] = 200.*pn.units.pc
s.physical_units()

HALO_NUMBER = 500
assert HALO_NUMBER <= len(s.halos())

halos = s.halos()[:HALO_NUMBER]

def centering(sim):
    com = pn.analysis.halo.center_of_mass(sim)
    return pn.transformation.inverse_translate(sim,com)

with centering(halos[0]):    
    p = pn.analysis.profile.Profile(halos[0], ndim=3)
    plt.yscale('log')
    plt.xscale('log')
    plt.plot(p['rbins'], p['mass_enc'], '--', color='brown')
    plt.xlabel('$R [kpc]$')
    plt.ylabel(r'M$_{enc}$ [M$_{\odot}$]')
    plt.savefig("Mass.png")


M200, ba, ca, triax = ([] for i in range(4))
for i_halo in range(HALO_NUMBER):
    with centering(halos[i_halo]):
        r200 = halos[i_halo].properties['Halo_R_Crit200'].in_units('kpc')
        m200 = halos[i_halo].properties['Halo_M_Crit200'].in_units('Msol')
        shape = pn.analysis.halo.halo_shape(halos[i_halo], N=1, rout=r200, bins='lin')
        if shape[1]>1e-10 and shape[2]>1e-10: #avoid infinities
            M200.append(m200)
            ba.append(shape[1])
            ca.append(shape[2])
            triax.append((1-shape[1]**2)/(1-shape[2]**2))

print(len(M200), len(ba), len(ca))
ba = np.array(ba)
ca = np.array(ca)
triax = np.array(triax)

#plotting
f, axs = plt.subplots(3,1,figsize=(9,10))
axs[0].scatter(np.log10(M200), ba)
axs[0].set_xlabel(r'$\log_{10}$(M$_{200}$ [M$_{\odot}$])')
axs[0].set_ylabel(r'$\frac{b}{a}$')
axs[1].scatter(np.log10(M200), ca)
axs[1].set_xlabel(r'$\log_{10}$(M$_{200}$ [M$_{\odot}$])')
axs[1].set_ylabel(r'$\frac{c}{a}$')
axs[2].scatter(np.log10(M200), triax)
axs[2].set_xlabel(r'$\log_{10}$(M$_{200}$ [M$_{\odot}$])')
axs[2].set_ylabel(r'Triaxiality: $\frac{a²-b²}{a²-c²}$')
f.savefig("Shape.png")
