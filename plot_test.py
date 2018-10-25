import numpy as np
import dmprofile
from dmprofile.src import plot
from dmprofile.src.utilities import read_from_file as rf

min_size1=1000
min_size2=500
c1, M200_1, res1, rel1 = rf('data/Concentration_'+str(min_size1)+'_128.txt')
c2, M200_2, res2, rel2 = rf('data/Concentration_'+str(min_size2)+'_128.txt')
c3, M200_3, res3, rel3 = rf('data/Concentration_'+str(min_size1)+'_256.txt')
c4, M400_4, res4, rel4 = rf('data/Concentration_'+str(min_size2)+'_256.txt')

s1, M200_shape_1 = rf('data/Shape_'+str(min_size1)+'_128.txt', mode='shape')
s2, M200_shape_2 = rf('data/Shape_'+str(min_size2)+'_128.txt', mode='shape')
s3, M200_shape_3 = rf('data/Shape_'+str(min_size1)+'_256.txt', mode='shape')
s4, M200_shape_4 = rf('data/Shape_'+str(min_size2)+'_256.txt', mode='shape')

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

shape = plot.Shape([s1, s2, s3, s4], extra_var=[M200_shape_1, M200_shape_2, M200_shape_3, M200_shape_4], 
                    nrows=1, ncols=3, name="figs/Shape.png", h=8, w=17)
shape.set_axis((0,0), xlabel=r'$\log_{10}(M)$ [M$_{\odot}$]', ylabel='b/a',
               xscale='linear', yscale='linear')
shape.set_axis((0,1), xlabel=r'$\log_{10}(M)$ [M$_{\odot}$]', ylabel='c/a',
               xscale='linear', yscale='linear')
shape.set_axis((0,2), xlabel=r'$\log_{10}(M)$ [M$_{\odot}$]', ylabel='triaxiality',
               xscale='linear', yscale='linear')

shape.binned_plot(0, axis_idx=(0,0), nbins=5, x_var="mass", y_var="b/a", color='blue', 
                  label='> '+str(min_size1)+' (128)', extra_idx=0, xscale='linear', xerr=0)
shape.binned_plot(1, axis_idx=(0,0), nbins=5, x_var="mass", y_var="b/a", color='red', 
                  label='> '+str(min_size2)+' (128)', extra_idx=1, xscale='linear', xerr=0)
shape.binned_plot(2, axis_idx=(0,0), nbins=5, x_var="mass", y_var="b/a", color='green', 
                  label='> '+str(min_size1)+' (256)', extra_idx=2, xscale='linear', xerr=0)
shape.binned_plot(3, axis_idx=(0,0), nbins=5, x_var="mass", y_var="b/a", color='orange', 
                  label='> '+str(min_size2)+' (256)', extra_idx=3, xscale='linear', fit=True, xerr=0)

shape.binned_plot(0, axis_idx=(0,1), nbins=5, x_var="mass", y_var="c/a", color='blue', 
                  label='> '+str(min_size1)+' (128)', extra_idx=0, xscale='linear', xerr=0)
shape.binned_plot(1, axis_idx=(0,1), nbins=5, x_var="mass", y_var="c/a", color='red',
                  label='> '+str(min_size2)+' (128)', extra_idx=1, xscale='linear', xerr=0)
shape.binned_plot(2, axis_idx=(0,1), nbins=5, x_var="mass", y_var="c/a", color='green',
                  label='> '+str(min_size1)+' (256)', extra_idx=2, xscale='linear', xerr=0)
shape.binned_plot(3, axis_idx=(0,1), nbins=5, x_var="mass", y_var="c/a", color='orange',
                  label='> '+str(min_size2)+' (256)', extra_idx=3, xscale='linear',fit=True, xerr=0)

shape.binned_plot(0, axis_idx=(0,2), nbins=5, x_var="mass", y_var="triax", color='blue', 
                  label='> '+str(min_size1)+' (128)', extra_idx=0, xscale='linear', xerr=0)
shape.binned_plot(1, axis_idx=(0,2), nbins=5, x_var="mass", y_var="triax", color='red',
                  label='> '+str(min_size2)+' (128)', extra_idx=1, xscale='linear', xerr=0)
shape.binned_plot(2, axis_idx=(0,2), nbins=5, x_var="mass", y_var="triax", color='green',
                  label='> '+str(min_size1)+' (256)', extra_idx=2, xscale='linear', xerr=0)
shape.binned_plot(3, axis_idx=(0,2), nbins=5, x_var="mass", y_var="triax", color='orange',
                  label='> '+str(min_size2)+' (256)', extra_idx=3, xscale='linear',fit=True, xerr=0)
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

