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

def write_to_file(fname, *args, mode='normal', hdf5_shape=(0,0)):
    """
    Write to files a length variable group of arrays.
    Works both with numpy arrays and lists.
    """
    modes = ['normal','hdf5']
    if mode not in modes:
        raise ValueError('The specified mode is not supported.')
    if mode=='hdf5' and hdf5_shape==(0,0,0):
        raise ValueError('Please specify a valid shape for the hdf5 dataset.')
    if mode=='hdf5' and fname[-5:]!='.hdf5':
        raise ValueError('Please specify an appropriate name for a hdf5 file.')
    if len(hdf5_shape)!=2:
        raise TypeError('Currently only 2D hdf5 datasets can be stored.')
    if len(args)<2:
        raise ValueError('At least to arrays must be written to the file.')
    if not type(fname)==str:
        raise TypeError('The first argument must be a string, and represents the name of the file.')
    for arg in args:
        if not isinstance(arg, (list, np.ndarray)):
            raise TypeError('Except for the file name, all other arguments must be numpy arrays.')

    if mode=='hdf5':
        with h5py.File(fname, 'w') as f:
            dset = f.create_dataset('shape_dataset', hdf5_shape)
            counter = 0
            for x, *data in zip(*args):
                print("TYPE::", x)
                print("ID:", counter)
                dset[counter,0]=x[1][0]
                dset[counter,1]=x[2][0]
                dset[counter,2]=data[0]
                counter += 1
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

def read_from_file(fname, splitter=',', mode='normal'):
    modes = ['normal','hdf5']
    if mode not in modes:
        raise ValueError('The specified mode is not supported.')

    if mode=='hdf5':
        ba, ca, extra = ([] for i in range(3))
        with h5py.File(fname, 'r') as f:
            dset = f['shape_dataset']
            counter = 0
            for i in range(dset.shape[0]):
                ba.append(dset[i,0])
                ca.append(dset[i,1])
                extra.append(dset[i,2])
        return ba, ca, extra

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
