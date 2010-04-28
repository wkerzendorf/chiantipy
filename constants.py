'''A set of physical constants.

Copyright 2009, 2010 Kenneth P. Dere

This software is distributed under the terms of the GNU General Public License
that is found in the LICENSE file


Most are from http://physics.nist.gov/cuu - the NIST Reference on
Constants, Units and Uncertainty'''

import numpy as np
planck = 6.6260693e-27   #erg s
planckEv = 4.13566743e-15  # ev s
light = 29979245800.  # cm/s
ev2Ang = 12.39841875e+3
ev2Erg = 1.602176487e-12
pi = 3.1415926535897931
boltzmann = 1.3806504e-16  # cgs
boltzmannEv = 8.617343e-5
invCm2Ev = 1./8.06554465e+3
ryd2Ev = 13.6056923
ryd2erg = 2.17987197e-11  #ergs
fine = 7.2973525376e-3  # fine structure constant ~ 1./137
emass = 9.10938215e-28  #  electron mass in g
bohr = 0.52917720859e-8  # bohr radius in cm
#
# derived constants
std2fwhm = 2.*np.sqrt(2.*np.log(2.))
#
invCm2Erg = planck*light
#
boltzmannRyd = boltzmannEv/ryd2Ev
# collision produces the 8.63e-6 factor
collision = planck**2/((2.*pi*emass)**1.5*np.sqrt(boltzmann))
#
freeFree = 1.e+8*(light/(3.*emass))*(fine*planck/pi)**3*np.sqrt((2.*pi)/(3.*emass*boltzmann))
#
sutherland = (2./(3.*np.sqrt(3.)))*np.sqrt(pi/(2.*boltzmann*emass**3))*(planck*fine/pi)**3
#
freeBound = 1.e+8*(8.*fine*(planck**3))/(3.*np.sqrt(3.)*np.pi*(emass**4)*light)*(emass/(2.*np.pi*boltzmann))**1.5
#
verner = (1.e-8/(planck*light**3*emass**3))*(emass/(2.*pi*boltzmann))**1.5
#
freeFreeLoss = (8./3.)*np.sqrt(pi*boltzmann/(6.*emass**3))*(planck/pi)**2*fine**3
#
freeBoundLoss = ((16.*fine*(planck**2))/(3.*pi*np.sqrt(3.)*(emass**3)*(light**2)))*np.sqrt(emass/(2.*pi*boltzmann))

El=['h','he','li','be','b','c','n','o','f','ne','na', \
    'mg','al','si','p','s','cl','ar','k','ca','sc','ti', \
    'v','cr','mn','fe','co','ni','cu','zn',\
    'ga','ge','as','se','br','kr']
Ionstage=['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII', \
    'XIV','XV','XVI','XVII','XVIII','XIX','XX','XXI',' XXII','XXIII','XXIV', \
    'XXV','XXVI','XXVII','XXVIII','XXIX','XXX','XXXI','XXXII','XXXIII','XXXIV', \
    'XXXV','XXXVI','XXXVII']
#
