import chianti
import chianti.util as util
import chianti.constants as const
import matplotlib.pyplot as pl
#from sources import *
class pne:
    ''' classs top model a planetary nebula'''
    def __init__(self,  source, distance, temperature, density, elements = ['h','he']):
        self.Source = source
        self.Distance = distance
        self.Temperature = temperature
        self.Density = density
        self.IonStages = []
        for element in elements:
            z = util.el2z(element)
            for stage in range(1, z+2):
                gname = util.zion2name(z, stage)
                print ' gname = ', gname
                self.IonStages.append(chianti.core.ion(gname, temperature, density))
        for astage in self.IonStages:
            if hasattr(astage, 'Photox'):
                print astage.IonStr, astage.Photox.keys()
                energy = const.ryd2erg*astage.Photox['energy']
                source.incident(distance, energy)
                pl.loglog(energy, source.Incident['photons'])


