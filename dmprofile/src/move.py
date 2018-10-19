import pynbody as pn

def centering_com(sim):
    _com = pn.analysis.halo.center_of_mass(sim)
    return pn.transformation.inverse_translate(sim,_com)

def centering_mbp(halos, idx):
    """
    Center the 'sim' object using the most bound particle (mbp).
    'halos' must be an instance of the 'Halos' class
    'idx' is the index of the halos one wishes to center.
    """
    #if not isinstance(halos, Halos):
    #    raise TypeError('Centering with the most bound particle has to be done with a dmrpfile.src.halos.Halos() object.')
    _mbp = halos._get_mbp(idx)
    print(_mbp)
    #return pn.transformation.inverse_translate(sim,_mbp)
