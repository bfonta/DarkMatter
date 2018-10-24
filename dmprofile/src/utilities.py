import numpy as np
import random
import pynbody as pn
import h5py

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

def remove_low_occupancy_bins(h, nbins_min):
    """
    Returns an array representing the edges of the bins of the 'h' numpy.histogram
    where the bins containing less than 'nbins_min' samples were removed.
    It only removes bins at the edges. If a bins has less than 'nbins_min' samples but lies between
    bins that have 'nbins_min' samples or more it will not be removed.
    If many consecutive bins have less than 'nbins_min' samples and are at the edge, 
    all will be removed.
    """
    bin_samples = h[0]
    bin_edges = h[1]
    for i in reversed(range(len(bin_samples))):
        flag = False
        if bin_samples[i]<3:
            flag = True
            bin_edges = bin_edges[:-1]
        if flag==False:
            break
    for i in range(len(bin_samples)):
        flag = False
        if bin_samples[i]<3:
            flag = True
            bin_edges = bin_edges[1:]
        if flag==False:
            break
    return bin_edges

def memory_usage_resource():
    import os
    import psutil
    _process = psutil.Process(os.getpid())
    print(_process.memory_info().rss / float(2**20), "MiB being used.")

def write_to_file(fname, *args, mode='standard'):
    """
    Write to files a length variable group of arrays.
    Works both with numpy arrays and lists.
    """
    extension = {'txt'}
    ext = fname.split('.')
    if len(ext)==1:
        raise ValueError('Please add an extension to the name of the file.')
    if ext[1] not in extension:
        raise ValueError('The file cannot be saved with that extension.')
    modes = {'standard','shape'}
    if mode not in modes:
        raise ValueError('The specified mode is not supported.')
    if len(args)<2:
        raise ValueError('At least to arrays must be written to the file.')
    if not type(fname)==str:
        raise TypeError('The first argument must be a string, and represents the name of the file.')
    for arg in args:
        if not isinstance(arg, (list, np.ndarray)):
            raise TypeError('Except for the file name, all other arguments must be numpy arrays.')

    if mode=='shape':
        with open(fname, 'w') as f:    
            for x, *data in zip(*args):
                blanks = '{},'*len(args)+'\n'
                f.write('{},{},'.format(x[1][0],x[2][0]))
                f.write('{}\n'.format(data[0]))

    else:
        with open(fname, 'w') as f:    
            for x, *data in zip(*args):
                blanks = '{},'*len(args)+'\n'
                f.write('{},'.format(x))
                for i in range(len(args)-1):
                    if i==len(args)-2:
                        f.write('{}\n'.format(data[i]))
                    else:
                        f.write('{},'.format(data[i]))

def read_from_file(fname, splitter=',', mode='standard'):
    """
    File types: 'txt'

    Modes: 'standard', 'shape'
    The 'standard' mode assumes that the data is stored in a plain format, like:
    1.1,4.3,0.5,...
    2.1,5.9,2.4,...
    3.4,10.2,8.2,...
    ...

    The 'shape' mode was created to deal with the way the plotting of the pynbody shape
    is done in this package ('dmprofile'). The shape, when stored, is a pynbody array 
    which includes the position, b/a, c/a, alignment and rotation matrix. As such, it 
    cannot be read with the 'standard' mode. The method return the shape plus an additional
    variable, such as the mass or the radius, for instance. 
    """
    extension = {'txt'}
    ext = fname.split('.')
    if len(ext)==1:
        raise ValueError('Please add an extension to the name of the file.')
    if ext[1] not in extension:
        raise ValueError('The file cannot be saved with that extension.')
    modes = {'standard', 'shape'}
    if mode not in modes:
        raise ValueError('The specified mode is not supported.')

    if mode=='shape':
        ba, ca, extra = ([] for i in range(3))
        with open(fname, 'r') as f:
            for values in f:
                values = values.replace('\n','').split(splitter)
                ba.append(float(values[0]))
                ca.append(float(values[1]))
                extra.append(float(values[2]))
        assert len(ba)==len(ca)
        s = []
        for i in range(len(ba)):
            triax = (1-ba[i])**2 / (1-ca[i])**2
            s.append([[0.],[ba[i]],[ca[i]],[0.],[0.],[triax]])
        return s, extra

    else:
        with open(fname) as f:
            for values in f:
                values = values.split(splitter)
                ncols = len(values)
                break
            data = [[] for i in range(ncols)]
            with open(fname) as f:
                for values in f:
                    values = values.replace('\n','').split(splitter)
                    for i,val in enumerate(values):
                        if val=='False':
                            data[i].append(False)
                        elif val=='True':
                            data[i].append(True)
                        else:
                            data[i].append(float(val))
        return [np.array(data[i]) for i in range(ncols)]
