import argparse
import numpy as np
import dmprofile
from dmprofile.src import plot
from dmprofile.src.utilities import read_from_file as rf
from dmprofile.src.parser import add_args

FLAGS, _ = add_args(argparse.ArgumentParser())
print("Parsed arguments:")
for k,v in FLAGS.__dict__.items():
    print('{}: {}'.format(k,v))

st = FLAGS.sim_types
min_size = FLAGS.sim_min_particle_number
sim_sizes = FLAGS.sim_sizes[0]
redshift_dict = {5: '103', 6: '080', 7: '065', 8: '054', 9: '045'}
rshift = redshift_dict[FLAGS.redshift]

_c = [rf('data/Concentration_'+st[i]+'_'+str(min_size)+'_'+sim_sizes+'_redshift'+rshift+'.txt') 
      for i in range(len(st))]

c, M200, res, rel = ([[] for _ in range(len(st))] for _ in range(4))
for isim in range(len(st)):
    c[isim] = _c[isim][0]
    M200[isim] = np.log10(_c[isim][1])
    #res[isim] = _c[isim][2]
    #rel[isim] = _c[isim][3]

c_obj = plot.Concentration(c, extra_var=M200, nrows=1, ncols=1, 
                           name='figs/Concentration'+str(FLAGS.redshift)+'.png', h=8, w=7)
c_obj.set_all_properties(model='concentration_mass_log')

names_dict = {'DMONLY': 'DMONLY', 
              'REF_LateRe': 'ZC_WSNe_Kinetic_LateRe',
              'REF_NOZCOOL': 'PrimC_WSNe_Kinetic_EarlyRe',
              'REF_NOZCOOL_LateRe': 'PrimC_WSNe_Kinetic_LateRe',
              'WTHERM_LateRe': 'ZC_SSNe_Thermal_LateRe'}
colors_dict = {'DMONLY': 'black', 
              'REF_LateRe': 'green',
              'REF_NOZCOOL': 'red',
              'REF_NOZCOOL_LateRe': 'darkorange',
              'WTHERM_LateRe': 'blueviolet'}
for isim in range(len(st)):
    c_obj.binned_plot(isim, axis_idx=(0,0), nbins=8, min_bins=3, rg = (7.,11.),
                      label = names_dict[st[isim]], 
                      color = colors_dict[st[isim]],
                      n_resampling=10000,
                      xscale='linear', xlim=[6.8,11.2],
                      extra_idx=isim, fit=False)
c_obj.savefig()
