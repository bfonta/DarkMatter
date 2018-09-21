import pynbody as pn


def centering_com(sim):
    _com = pn.analysis.halo.center_of_mass(sim)
    return pn.transformation.inverse_translate(sim,_com)
   
