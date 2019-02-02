import argparse
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
    M200[isim] = _c[isim][1]
    #res[isim] = _c[isim][2]
    #rel[isim] = _c[isim][3]

c_obj = plot.Concentration(c, extra_var=M200, nrows=1, ncols=1, 
                           name='figs/Concentration'+str(FLAGS.redshift)+'.png', h=8, w=7)
c_obj.set_all_properties(model='concentration_mass')

for isim in range(len(st)):
    c_obj.binned_plot(isim, axis_idx=(0,0), nbins=5, 
                      label='> '+str(min_size)+' '+st[isim]+' '+sim_sizes, min_bins=10,
                      extra_idx=isim, xscale='log', xerr=0, fit=False)
c_obj.savefig()
