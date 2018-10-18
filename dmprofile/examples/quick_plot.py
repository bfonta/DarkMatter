"""
Basic plots using text files. The files have two columns of values (x and y data)
"""
from context import src
from src.plot import plot_from_file as pf

s1 = "data/shape_double_1"
s2 = "data/shape_double_2"
s3 = "data/shape_double_3"

pf(s1+'.txt', scatter=True)
pf(s2+'.txt', scatter=True)
pf(s3+'.txt', scatter=True, save='shape_double.png')
