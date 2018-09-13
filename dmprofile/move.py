import pynbody as pn

class Centering():
    """
    Deals with centering operations aplied to a simulation object
    """
    def __init__(self, sim):
        self.sim = sim

    def com(self):
        _com = pn.analysis.halo.center_of_mass(self.sim)
        return pn.transformation.inverse_translate(self.sim,_com)
   
