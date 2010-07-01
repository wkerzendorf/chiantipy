''' functions needed for multiprocessing module mspectrum'''
def doFf(inputs):
    ''' helper for freefree'''
    ionS = inputs[0]
    temperature = inputs[1]
    wavelength = inputs[2]
    cont = chianti.core.continuum(ionS, temperature)
    cont.freeFree(wavelength)
    return cont
    #
    # ----------------------------------------------
    #
def doFb(inputs):
    ''' helper for freebound'''
    ionS = inputs[0]
    temperature = inputs[1]
    wavelength = inputs[2]
    cont = chianti.core.continuum(ionS, temperature)
    cont.freeBound(wavelength)
    return cont
    #
    # ----------------------------------------------
    #
def doIon(inputs):
    ''' helper for ion'''
    ionS = inputs[0]
    temperature = inputs[1]
    density = inputs[2]
    wavelength = inputs[3]
    wvlRange = [wavelength.min(), wavelength.max()]
    filter = inputs[4]
    thisIon = chianti.core.ion(ionS, temperature, density)
    thisIon.intensity(wvlRange = wvlRange)
    thisIon.spectrum(wavelength,  filter=filter)
    if (thisIon.Z - thisIon.Ion) in [0, 1]:
        thisIon.twoPhoton(wavelength)
    return thisIon
    #
    # ----------------------------------------------
    #
def doDiel(inputs):
    ''' helper for dielectronic ions'''
    ionS = inputs[0]
    temperature = inputs[1]
    density = inputs[2]
    wavelength = inputs[3]
    wvlRange = [wavelength.min(), wavelength.max()]
    filter = inputs[4]
    thisIon = chianti.core.ion(ionS, temperature, density)
    thisIon.intensity(wvlRange = wvlRange)
    thisIon.spectrum(wavelength, filter=filter)
    return thisIon
