''' functions needed for multiprocessing module mspectrum'''
import chianti
def doFfQ(inQ, outQ):
    ''' helper for freefree'''
    for inputs in iter(inQ.get, 'STOP'):
        ionS = inputs[0]
        temperature = inputs[1]
        wavelength = inputs[2]
        cont = chianti.core.continuum(ionS, temperature)
        cont.freeFree(wavelength)
        outQ.put(cont.FreeFree)
    return
    #
    # ----------------------------------------------
    #
def doFbQ(inQ, outQ):
    ''' helper for freefree'''
    for inputs in iter(inQ.get, 'STOP'):
        ionS = inputs[0]
        temperature = inputs[1]
        wavelength = inputs[2]
        cont = chianti.core.continuum(ionS, temperature)
        cont.freeBound(wavelength)
        outQ.put(cont.FreeBound)
    return
    #
    # ----------------------------------------------
    #
def doIonQ(inQueue, outQueue):
    ''' helper for ion, also does two-photon'''
    for inputs in iter(inQueue.get, 'STOP'):
        ionS = inputs[0]
        temperature = inputs[1]
        density = inputs[2]
        wavelength = inputs[3]
        wvlRange = [wavelength.min(), wavelength.max()]
        filter = inputs[4]
        thisIon = chianti.core.ion(ionS, temperature, density)
        thisIon.spectrum(wavelength,  filter=filter)
        outList = [ionS, thisIon.Spectrum]
        if not thisIon.Dielectronic:
            if (thisIon.Z - thisIon.Ion) in [0, 1]:
                thisIon.twoPhoton(wavelength)
                outList.append(thisIon.TwoPhoton)
        outQueue.put(outList)
    return
