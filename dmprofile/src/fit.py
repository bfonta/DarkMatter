import numpy as np
from pynbody.analysis.theoretical_profiles import NFWprofile

def nfw_fit(prof):
    """
    NFW fit. Returns the characteristic density and the scale radius.
    """
    _fit = NFWprofile.fit(np.array(prof['rbins']),
                          np.array(prof['density']))
    print("FIT: ", _fit)
    return _fit[0][0], _fit[0][1]
