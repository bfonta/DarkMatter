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
sim_size = FLAGS.sim_size

_c = [rf('data/Concentration_'+st[i]+'_'+str(min_size)+'_'+sim_size+'.txt') for i in range(len(st))]
_s = [rf('data/Shape_'+st[i]+'_'+str(min_size)+'_'+sim_size+'.txt',mode='shape') for i in range(len(st))]

c, M200, res, rel, s, M200_shape = ([] for i in range(6))
for isim in range(len(st)):
    c.append(_c[isim][0])
    M200.append(_c[isim][1])
    res.append(_c[isim][2])
    rel.append(_c[isim][3])
    s.append(_s[isim][0])
    M200_shape.append(_s[isim][1])


shape = plot.Shape([s[i] for i in range(len(st))], extra_var=[M200_shape[i] for i in range(len(st))], 
                   nrows=1, ncols=3, name="figs/Shape.png", h=8, w=17)
shape.set_axis((0,0), xlabel=r'$\log_{10}(M)$ [M$_{\odot}$]', ylabel='b/a',
               xscale='linear', yscale='linear')
shape.set_axis((0,1), xlabel=r'$\log_{10}(M)$ [M$_{\odot}$]', ylabel='c/a',
               xscale='linear', yscale='linear')
shape.set_axis((0,2), xlabel=r'$\log_{10}(M)$ [M$_{\odot}$]', ylabel='triaxiality',
               xscale='linear', yscale='linear')

for isim in range(len(st)):
    fit=False
    if isim==len(st)-1: fit=True
    shape.binned_plot(isim, axis_idx=(0,0), nbins=5, x_var="mass", y_var="b/a", color='blue', 
                      label='> '+str(min_size)+' '+st[isim]+' '+sim_size, 
                      extra_idx=0, xscale='linear', xerr=0, ylim=[.55,.85],fit=fit)

for isim in range(len(st)):
    fit=False
    if isim==len(st)-1: fit=True
    shape.binned_plot(isim, axis_idx=(0,1), nbins=5, x_var="mass", y_var="c/a", color='blue', 
                      label='> '+str(min_size)+' '+st[0]+' '+sim_size, 
                      extra_idx=0, xscale='linear', xerr=0, ylim=[.4,.7], fit=fit)

for isim in range(len(st)):
    fit=False
    if isim==len(st)-1: fit=True
    shape.binned_plot(isim, axis_idx=(0,2), nbins=5, x_var="mass", y_var="triax", color='blue', 
                      label='> '+str(min_size)+' '+st[0]+' '+sim_size, 
                      extra_idx=0, xscale='linear', xerr=0, ylim=[.53,.9], fit=fit)

shape.savefig()




"""
shape.scatter_plot(0, axis_idx=(0,0), x_var="mass", y_var="b/a",
                   resolved_bools=res1, relaxed_bools=rel1, extra_idx=0)
shape.scatter_plot(0, axis_idx=(0,1), x_var="mass", y_var="c/a",
                   resolved_bools=res1, relaxed_bools=rel1, extra_idx=0)
shape.scatter_plot(0, axis_idx=(0,2), x_var="mass", y_var="triax",
                   resolved_bools=res1, relaxed_bools=rel1, extra_idx=0)
shape.scatter_plot(1, axis_idx=(1,0), x_var="mass", y_var="b/a",
                   resolved_bools=res2, relaxed_bools=rel2, extra_idx=1)
shape.scatter_plot(1, axis_idx=(1,1), x_var="mass", y_var="c/a",
                   resolved_bools=res2, relaxed_bools=rel2, extra_idx=1)
shape.scatter_plot(1, axis_idx=(1,2), x_var="mass", y_var="triax",
                   resolved_bools=res2, relaxed_bools=rel2, extra_idx=1)
"""

"""c_obj = plot.Concentration([c1, c2, c3, c4], extra_var=[M200_1, M200_2, M200_3, M200_4],
                           nrows=1, ncols=1, name='figs/Concentration.png')
c_obj.set_all_properties(model='concentration_mass')
c_obj.binned_plot(0, (0,0), 10, extra_idx=0, color='green')
c_obj.binned_plot(1, (0,0), 10, extra_idx=1, shift=0.1e11, color='orange')
c_obj.binned_plot(2, (0,0), 10, extra_idx=2, shift=0.2e11, color='blue')
c_obj.binned_plot(3, (0,0), 10, extra_idx=3, shift=0.2e11, color='red')

c_obj.scatter_plot(0, (0,0), extra_idx=0, resolved_bools=res1, relaxed_bools=rel1)
c_obj.scatter_plot(1, (1,0), extra_idx=1, resolved_bools=res2, relaxed_bools=rel2)
c_obj.savefig()
"""
