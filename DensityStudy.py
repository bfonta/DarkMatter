import matplotlib
matplotlib.use('Agg')
from matplotlib import pylab as plt

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
SnapshotFolder = "Test_NOSN_NOZCOOL_L010N0128/data/snapshot_103/snap_103.hdf5"
SubhalosFolder = "Test_NOSN_NOZCOOL_L010N0128/data/subhalos_103/subhalo_103"

s = pn.load(os.path.join(DataFolder,SubhalosFolder))
s['eps'] = 200.*pn.units.pc
s.physical_units()

HALO_NUMBER = 5
MIN_R, MAX_R = 0, 35
COMPONENTS_NUMBER = 4 #total, gas, stars and dark matter
assert HALO_NUMBER < len(s.halos())
halos = s.halos()[:HALO_NUMBER]
var = 0


r200 = halos[var].properties['Halo_R_Crit200'].in_units('kpc')
with pn.analysis.halo.center(halos[var], mode='com'):   
    def xvar():
        return lambda x: x['r']
    
    p = []
    halo_component = [halos[var], halos[var].g, halos[var].s, halos[var].d]
    for i in range(COMPONENTS_NUMBER): 
        p.append(pn.analysis.profile.Profile(halo_component[i], nbins=250, max=MAX_R, ndim=3))

#plotting

WIDTH, HEIGHT = 11, 10
labels = ["total", "gas", "stars", "dark matter"]

f, axs = plt.subplots(2,2,figsize=(WIDTH,HEIGHT))
f.suptitle("Enclosed mass as a function of radial distance")
for i in range(COMPONENTS_NUMBER):
    idx1 = int(np.floor(i/2))
    idx2 = int(i%2)
    axs[idx1,idx2].set_yscale('log')
    axs[idx1,idx2].set_xscale('log')
    axs[idx1,idx2].plot(p[i]['rbins'], p[i]['mass_enc'], '--', color='brown', label=labels[i])
    axs[idx1,idx2].set_xlabel('$R [kpc]$')
    axs[idx1,idx2].set_ylabel(r'M$_{enc}$ [M$_{\odot}$]')
    axs[idx1,idx2].legend()
f.savefig("EnclosedMass")

#define density
def density(profile, max=MAX_R, min=MIN_R):
    rbins = profile["rbins"][(profile["rbins"]>min) & (profile["rbins"]<max)]
    rbins_shift = np.roll(profile["rbins"][(profile["rbins"]>min) & (profile["rbins"]<max)], 1)
    rbins_shift[0] = 0.
    return profile["mass"][(profile["rbins"]>min) & (profile["rbins"]<max)] / (1.333*np.pi*(np.power(rbins,3)-np.power(rbins_shift,3)))

#fitting
def nfw(x, d_char, rs):
    """nfw profile with two independent parameters: the scale radius and the characteristic density"""
    return d_char / ((x/rs)*np.power(1+x/rs,2))
def nfw_rsquared(x, d_char, scale_radius):
    return d_char*x**2 / ((x/scale_radius)*np.power(1+x/scale_radius,2))


MIN_FIT, MAX_FIT = 0, 30
params, covar = ([] for _ in range(2))
params2, covar2 = ([] for _ in range(2))
p0 = [[3.5e+08, 30],[6.4e+8, 30],[4.1e+8, 30],[7.1e+08, 30]]
p0_2 = [[7.6e+07, 6.4e+02],[7.6e+08, 6.4e+02],[7.6e+08, 6.4e+02],[7.6e+08, 6.4e+02]]
for i in range(COMPONENTS_NUMBER):
    fit = curve_fit(nfw, p[i]['rbins'][p[i]['rbins']<MAX_FIT], density(p[i], MAX_FIT), bounds=([1e6,0], [1e12,1e3]), p0=p0[i], maxfev=2000000)
    params.append(fit[0])
    covar.append(fit[1])
    
    fit2 = curve_fit(nfw_rsquared, p[i]['rbins'][p[i]['rbins']<MAX_FIT], density(p[i], MAX_FIT)*p[i]['rbins'][p[i]['rbins']<MAX_FIT]**2, p0_2[i],maxfev=2000000)
    params2.append(fit2[0])
    covar2.append(fit2[1])
    
print("Fit parameters: ", params)
print("Squared fit parameters: ", params)

#plotting
#######################################
f, axs = plt.subplots(2,2,figsize=(WIDTH,HEIGHT))
f.suptitle(r"$\rho(r)$")
for i in range(COMPONENTS_NUMBER):
    idx1 = int(np.floor(i/2))
    idx2 = int(i%2)
    axs[idx1,idx2].set_yscale('log')
    #axs[idx1,idx2].set_xscale('log')
    axs[idx1,idx2].plot(p[i]['rbins'], density(p[i]), '--', color='black', label=labels[i])
    axs[idx1,idx2].plot(p[i]['rbins'][p[i]['rbins']<MAX_FIT], nfw(p[i]['rbins'][p[i]['rbins']<MAX_FIT], *params[i]), color='red', label=labels[i]+' fit')
    axs[idx1,idx2].set_xlabel('$R [kpc]$')
    axs[idx1,idx2].set_ylabel(r'$\rho$ [M$_{\odot}$ kpc$^{-3}$]')
    axs[idx1,idx2].legend()
f.savefig("Density")

#######################################
f, axs = plt.subplots(2,2,figsize=(WIDTH,HEIGHT))
f.suptitle(r"$\rho(r)$")
for i in range(COMPONENTS_NUMBER):
    idx1 = int(np.floor(i/2))
    idx2 = int(i%2)
    axs[idx1,idx2].set_yscale('log')
    axs[idx1,idx2].set_xscale('log')
    axs[idx1,idx2].plot(p[i]['rbins'], density(p[i])*p[i]['rbins']**2, '--', color='blue', label=labels[i])
    axs[idx1,idx2].plot(p[i]['rbins'], nfw_rsquared(p[i]['rbins'], *params2[i]), color='red', label=labels[i])
    axs[idx1,idx2].set_xlabel('$R [kpc]$')
    axs[idx1,idx2].set_ylabel(r'$\rho$ $r^{2}$ [M$_{\odot}$ kpc$^{-1}$]')
    axs[idx1,idx2].legend()
f.savefig("Density_squared.png")



s = pn.load(os.path.join(DataFolder,SnapshotFolder))
rho_crit_z = pn.array.SimArray(pn.analysis.cosmology.rho_crit(s, unit="Msol kpc**-3"), "Msol kpc**-3")
v=200
rv = np.array(halos[var].properties['Halo_R_Crit200'].in_units('kpc a h**-1'))
print("Critical density: ", rho_crit_z, rho_crit_z.units, "Virial radius ", rv)

def nfw_rs_only(x, rs, rv, v):
    """
    NFW profile with one independent parameter: the scale radius.
    The profile is normalized to the critical density at the time of the simulation.
    """
    def g(rs):
        return 1 / (np.log(1+rv/rs)-(rv/rs)/(1+rv/rs))
    def denominator(x, rs): return (x/rs)*np.power(1+x/rs,2)
    return 200*np.power(rv/rs,3)*g(rs) / (3*denominator(x,rs))

params, covar = ([] for _ in range(2))
p0 = [[0.5],[1.16],[0.7],[30.]]
for i in range(COMPONENTS_NUMBER):
    fit = curve_fit(lambda x,rs: nfw_rs_only(x,rs,rv,v), p[i]['rbins'][p[i]['rbins']<MAX_FIT], 
                    density(p[i],MAX_FIT)/rho_crit_z, p0=p0[i], bounds=(0,45), maxfev=2000000)
    params.append(fit[0])
    covar.append(fit[1])
    print("Scale radius of the", labels[i], "distribution: ", fit[0][0])

#plotting
f, axs = plt.subplots(2,2,figsize=(WIDTH,HEIGHT))
f.suptitle(r"One parameter dependent $\rho(r)$")
for i in range(COMPONENTS_NUMBER):
    idx1 = int(np.floor(i/2))
    idx2 = int(i%2)
    axs[idx1,idx2].set_yscale('log')
    #axs[idx1,idx2].set_xscale('log')
    axs[idx1,idx2].plot(p[i]['rbins'], density(p[i])/rho_crit_z, '--', color='black', label=labels[i])
    axs[idx1,idx2].plot(p[i]['rbins'][p[i]['rbins']<MAX_FIT], nfw_rs_only(p[i]['rbins'][p[i]['rbins']<MAX_FIT], *params[i], rv=rv, v=v), color='red', label=labels[i]+' fit')
    axs[idx1,idx2].set_xlabel('$R [kpc]$')
    axs[idx1,idx2].set_ylabel(r'$\rho$ / $\rho_{crit}^{0}$')
    axs[idx1,idx2].legend()
f.savefig("Density_complex")
