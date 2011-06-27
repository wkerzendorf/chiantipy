class blackStar:
    '''temperature in K, radius is the stellar radius in cm '''
    def __init__(self, temperature, radius):
        self.temperature = temperature
        self.radius = radius
    def incident(self, distance, energy):
        ''' distance in cm and energy in erg'''
        bb = pi*(self.radius/distance)**2*blackbody(temperature, energy)
        self.Incident = bb
    #
    # ---------------------------------------------------------------------
    #
def blackbody(temperature, variable, hnu=1):
    ''' to calculate the black body photon distribution as a function of energy in erg (hnu = 1) or as a function
    of wavelength (Angstroms) (hnu=0)
    photons cm^-2 s^-1 str^-1 ergs^-1'''
    if hnu:
        energy = variable
        bb =(2./(const.planck*(const.hc**2)))*energy**2/(np.exp(energy/(const.boltzmann*temperature)) - 1.)
        return {'photons':bb, 'temperature':temperature, 'energy':energy}
    else:
        wvl = 1.e-8*variable
        bb = ((2.*const.pi*const.light)/wvl**4)/(np.exp(const.hc/(wvl*const.boltzmann*temperature)) - 1.)
        return {'photons':bb, 'temperature':temperature, 'wvl':wvl}
