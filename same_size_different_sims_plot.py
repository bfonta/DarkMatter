import numpy as np
import dmprofile
from dmprofile.src import plot
from dmprofile.src.utilities import read_from_file as rf

size = '128'
st = ['DMONLY', 'REF_LateRe', 'REF_NOZCOOL', 'REF_NOZCOOL_LateRe', 'WTHERM_LateRe']
min_size=300

_c = [rf('data/Concentration_'+st[i]+'_'+str(min_size)+'_'+size+'.txt') for i in range(len(st))]
_s = [rf('data/Shape_'+st[i]+'_'+str(min_size)+'_'+size+'.txt', mode='shape') for i in range(len(st))]

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

shape.binned_plot(0, axis_idx=(0,0), nbins=5, x_var="mass", y_var="b/a", color='blue', 
                  label='> '+str(min_size)+' '+st[0]+' '+size, 
                  extra_idx=0, xscale='linear', xerr=0, fit=True)
shape.binned_plot(1, axis_idx=(0,0), nbins=5, x_var="mass", y_var="b/a", color='red',
                  label='> '+str(min_size)+' '+st[1]+' '+size, 
                  extra_idx=1, xscale='linear', xerr=0)
shape.binned_plot(2, axis_idx=(0,0), nbins=5, x_var="mass", y_var="b/a", color='green',
                  label='> '+str(min_size)+' '+st[2]+' '+size, 
                  extra_idx=2, xscale='linear', xerr=0)
shape.binned_plot(3, axis_idx=(0,0), nbins=5, x_var="mass", y_var="b/a", color='orange',
                  label='> '+str(min_size)+' '+st[3]+' '+size, 
                  extra_idx=3, xscale='linear', xerr=0)
shape.binned_plot(4, axis_idx=(0,0), nbins=5, x_var="mass", y_var="b/a", color='grey',
                  label='> '+str(min_size)+' '+st[4]+' '+size, 
                  extra_idx=4, xscale='linear', xerr=0, ylim=[.55,.85])

shape.binned_plot(0, axis_idx=(0,1), nbins=5, x_var="mass", y_var="c/a", color='blue', 
                  label='> '+str(min_size)+' '+st[0]+' '+size, 
                  extra_idx=0, xscale='linear', xerr=0, fit=True)
shape.binned_plot(1, axis_idx=(0,1), nbins=5, x_var="mass", y_var="c/a", color='red',
                  label='> '+str(min_size)+' '+st[1]+' '+size, 
                  extra_idx=1, xscale='linear', xerr=0)
shape.binned_plot(2, axis_idx=(0,1), nbins=5, x_var="mass", y_var="c/a", color='green',
                  label='> '+str(min_size)+' '+st[2]+' '+size, 
                  extra_idx=2, xscale='linear', xerr=0)
shape.binned_plot(3, axis_idx=(0,1), nbins=5, x_var="mass", y_var="c/a", color='orange',
                  label='> '+str(min_size)+' '+st[3]+' '+size, 
                  extra_idx=3, xscale='linear', xerr=0)
shape.binned_plot(4, axis_idx=(0,1), nbins=5, x_var="mass", y_var="c/a", color='grey',
                  label='> '+str(min_size)+' '+st[4]+' '+size, 
                  extra_idx=4, xscale='linear', xerr=0, ylim=[.4,.7])

shape.binned_plot(0, axis_idx=(0,2), nbins=5, x_var="mass", y_var="triax", color='blue', 
                  label='> '+str(min_size)+' '+st[0]+' '+size, 
                  extra_idx=0, xscale='linear', xerr=0, fit=True)
shape.binned_plot(1, axis_idx=(0,2), nbins=5, x_var="mass", y_var="triax", color='red',
                  label='> '+str(min_size)+' '+st[1]+' '+size, 
                  extra_idx=1, xscale='linear', xerr=0)
shape.binned_plot(2, axis_idx=(0,2), nbins=5, x_var="mass", y_var="triax", color='green',
                  label='> '+str(min_size)+' '+st[2]+' '+size, 
                  extra_idx=2, xscale='linear', xerr=0)
shape.binned_plot(3, axis_idx=(0,2), nbins=5, x_var="mass", y_var="triax", color='orange',
                  label='> '+str(min_size)+' '+st[3]+' '+size, 
                  extra_idx=3, xscale='linear', xerr=0)
shape.binned_plot(4, axis_idx=(0,2), nbins=5, x_var="mass", y_var="triax", color='grey',
                  label='> '+str(min_size)+' '+st[4]+' '+size, 
                  extra_idx=4, xscale='linear', xerr=0, ylim=[.53,.9])

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
