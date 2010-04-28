''' line profile filters from creating synthetic spectra

Copyright 2009, 2010 Kenneth P. Dere

This software is distributed under the terms of the GNU General Public License
that is found in the LICENSE file

'''
import numpy as np
import types
def gaussianR(wvl,wvl0, factor=None):
    '''a gaussian filter where factor is the resolving power, so that the gaussian width (standard deviation)
    is given by wvl0/factor'''
    if type(factor) != types.NoneType:
        std = wvl0/factor
    else:
        print ' the resolving power of the gaussianR filter is undefined'
        return None
    wvl = np.asarray(wvl, 'float64')
    dwvl = wvl - np.roll(wvl, 1)
    dwvl[0] = dwvl[1]
    return np.exp(-((wvl - wvl0)/(2.*std))**2)/(dwvl*np.sqrt(2.*np.pi)*std)

def gaussian(wvl,wvl0, factor=None):
    '''a gaussian filter where factor is the gaussian width (standard deviation)'''
    if type(factor) != types.NoneType:
        std = factor
    else:
        print ' the width of the gaussian filter is undefined'
        return None
    wvl = np.asarray(wvl, 'float64')
    dwvl = wvl - np.roll(wvl, 1)
    dwvl[0] = dwvl[1]
    return np.exp(-((wvl - wvl0)/(2.*std))**2)/(dwvl*np.sqrt(2.*np.pi)*std)

def boxcar(wvl, wvl0, factor=None):
    ''' box-car filter, factor is the full width of the box filter'''
    wvl = np.asarray(wvl, 'float64')
    dwvl = wvl - np.roll(wvl, 1)
    dwvl[0] = dwvl[1]
    one = np.ones_like(wvl)
    zed = np.zeros_like(wvl)
    if type(factor) != types.NoneType:
        # width must be at least equal to the wavelength step
        width = max(factor, dwvl.min())
    else:
        print ' the width of the box filter is undefined'
        return None
    good1 = (wvl > wvl0 - width/2.)
    good2 = (wvl < wvl0 + width/2.)
    realgood = np.logical_and(good1, good2)
    return np.where(realgood, one, zed)/(dwvl*width)

def lorentz(wvl, wvl0, factor=None):
    '''the lorentz profile with the exception that all factors are in wavelength units
    rather than frequency as the lorentz profile is usually defined
    factor is the value of the so-called constant gamma'''
    if type(factor) != types.NoneType:
        gamma = factor
    else:
        print ' the factor gamma of the lorentz filter is undefined'
        return None
    wvl = np.asarray(wvl, 'float64')
    dwvl = wvl - np.roll(wvl, 1)
    dwvl[0] = dwvl[1]
    return (gamma/(2.*np.pi)**2)/(dwvl(wvl - wvl0)**2 + (gamma/(4.*np.pi))**2)




