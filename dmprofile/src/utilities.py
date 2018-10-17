import numpy as np
import random
import pynbody as pn

def rho_crit(sim):
    return pn.array.SimArray(pn.analysis.cosmology.rho_crit(sim, unit="Msol kpc**-3"), "Msol kpc**-3")

def intersect(x=[], y=[], intersect_value=1., deg=1):
    """
    Calculates de intersection between two binned distributions.
    Returns 
    1. The abcissa of the input data point closest to the intersection.
    2. A boolean that says whether there was at least one intersection or not.
    """
    if len(x)<2 or len(y)<2:
        raise ValueError('Please provide the x and y data.')
    _polyfit = np.polyfit(x, y, deg=deg)
    _func = np.poly1d(_polyfit)
    _x2 = np.linspace(x[0], x[-1], 10000)
    _y2 = _func(_x2)
    _idx = np.argwhere(np.diff(np.sign(_y2 - intersect_value))).flatten()
    if len(_idx)==0:
        return -1, False
    _idx = int(_idx[-1]) #consider the last intersection index only
    if abs(_y2[_idx]-intersect_value) < abs(_y2[_idx+1]-intersect_value):
        return _x2[_idx], True
    else:
        return _x2[_idx+1], True

def is_list_empty(l):
    if isinstance(l, list): #it is a list
        return all( map(is_list_empty, l)) #all([]) is True
    return False

def linear_function_generator(p1, p2, left, right):
    while True:
        r = random.uniform(left,right)
        yield r, r*p1 + p2

def memory_usage_resource():
    import os
    import psutil
    _process = psutil.Process(os.getpid())
    print(_process.memory_info().rss / float(2**20), "MiB being used.")

def write_to_file(fname, arr1, arr2):
    with open(fname, 'w') as f:
        for x, y in zip(arr1,arr2):
            f.write('{},{}\n'.format(x,y))
