'''Contains the main classes for ChiantiPy users.

Copyright 2009, 2010 Kenneth P. Dere

This software is distributed under the terms of the GNU General Public License
that is found in the LICENSE file


'''
import os
import types
#from ConfigParser import *
import numpy as np
from scipy import interpolate
#
# the following is necessary to make chiantipy non interactive for the web
try:
    chInteractive = int(os.environ['CHIANTIPY_INTERACTIVE'])
except:
    chInteractive = 1
#
if chInteractive:
    import pylab as pl
else:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as pl
try:
    from matplotlib.delaunay.triangulate import Triangulation
except:
    from scikits.delaunay.triangulate import Triangulation
import chianti
from chianti import util
import chianti.constants as const
import chianti.filters as chfilters
#
xuvtop=os.environ['XUVTOP']


#if chInteractive:
#    import matplotlib as pl
#else:
#    import pylab as pl
#
ip = util.ipRead()
MasterList = util.masterListRead()
#
defaults = util.defaultsRead(verbose = chInteractive)
#
#if defaults['gui'] == 'qt':
#   import chianti.gui_qt.gui as gui
#elif defaults['gui'] == 'wx':
#    import chianti.gui_wx.gui as gui
#else:
#    print ' unknown gui - ',defaults['gui']
#    print ' the full functionality of the chiant.core.ion class will not be available'
    #
if pl.rcParams['backend'].lower() == 'qt4agg':
    import chianti.gui_qt.gui as gui
elif pl.rcParams['backend'].lower() == 'wxagg':
    import chianti.gui_wx.gui as gui
elif pl.rcParams['backend'].lower() == 'gtkagg':
    import chianti.gui_cl.gui as gui
elif pl.rcParams['backend'].lower() == 'agg':
    import chianti.gui_cl.gui as gui
elif pl.rcParams['backend'].lower() == 'agg':
    import chianti.gui_cl.gui as gui
elif pl.rcParams['backend'].lower() == 'macosx ':
    import chianti.gui_cl.gui as gui
else:
    print ' - Warning - '
    print ' - in order to use the various gui dialogs, the matlpotlib/pylab backend needs'
    print ' - to be either Qt4Agg or WXAgg - '
    print ' - in order to use the command line dialogs, the matlpotlib/pylab backend needs'
    print ' - to be GTKAgg or MacOSX - '
    print ' - current backend is ',pl.rcParams['backend']
    print ' - the full functionality of the chianti.core.ion class may not be available'
    print ' - it would probably be better to set your matplotlib backend to either'
    print ' - Qt4Agg, WXAgg, GTKAgg, or MacOSX'
    print ' - using the command line dialogs for now but there could be problems -'
    import chianti.gui_cl.gui as gui
    #
class continuum:
    '''The top level class for continuum calculations.

    includes methods for the calculation of the free-free and free-bound continua.'''
    def __init__(self, ionStr,  temperature=None,  density=None):
        nameDict = util.convertName(ionStr)
        self.Z = nameDict['Z']
        self.Ion = nameDict['Ion']
        self.IonStr = ionStr
        self.Dielectronic = 0
        self.Defaults=defaults
        self.AbundanceName = defaults['abundfile']
        self.IoneqName = defaults['ioneqfile']
        #
        #  ip in eV, reading Ip of next lower level, needed for freeBound

        if self.Ion > 1:
            self.Ip=ip[self.Z-1, self.Ion-2]
        else:
            print ' in continuum, trying to use the neutral ion'
            self.FreeBound={}
            return
        #
        if type(temperature) != types.NoneType:
            self.Temperature = np.array(temperature,'float64')
        #
        if type(density) != types.NoneType:
            self.Density = np.asarray(density,'float64')
        #
        #----------------------------------------------------------------------------------------
        #
    def freeBoundEmiss(self, wvl, verner=1):
        '''Calculates the free-bound (radiative recombination) continuum emissivity of an ion.

        Uses the Gaunt factors of Karzas, W.J, Latter, R, 1961, ApJS, 6, 167
        for recombination to the ground level, the photoionization cross sections of
        Verner and Yakovlev, 1995, A&ASS, 109, 125
        are used to develop the free-bound cross section
        provides emissivity = ergs cm^-2 s^-1 str^-1 Angstrom ^-1 for an individual ion
        does not include the elemental abundance or ionization fraction
        the specified ion is the target ion'''
        #
        wvl = np.asarray(wvl, 'float64')
        try:
            temperature = self.Temperature
        except:
            print ' temperature undefined'
            return
        #
        try:
            fblvl = self.Fblvl
        except:
            fblvlname = util.zion2filename(self.Z,self.Ion-1)+'.fblvl'
            self.Fblvl = util.fblvlRead(fblvlname)
            fblvl = self.Fblvl
        #  need some data for the recombining ion
        try:
            rFblvl = self.rFblvl
        except:
            if self.Ion == self.Z+1:
                # then we looking for the bare ion
                rFblvl = {'mult':[1., 1.]}
            else:
                rfblvlname = util.zion2filename(self.Z,self.Ion)+'.fblvl'  # previously self.Ion)
                self.rFblvl = util.fblvlRead(rfblvlname)
                rFblvl = self.rFblvl
        #  6/9/2010 the recombining iion is the present ion
        #
        # for the ionization potential, Ip, must use that of the recombined ion
        ipcm = self.Ip/const.invCm2Ev
        #
        nlvls = len(fblvl['lvl'])
        # pqn = principle quantum no. n
        pqn = fblvl['pqn']
        # l is angular moment quantum no. L
        l = fblvl['l']
        # energy level in inverse cm
        ecm = fblvl['ecm']
        for i in range(nlvls):
            print ' lvl ecm wvl = ',i, ecm[i], 1.e+8/(ipcm-ecm[i])
        # statistical weigths/multiplicities
        mult = fblvl['mult']
        multr = rFblvl['mult']
        #
        #
        ipcm = self.Ip/const.invCm2Ev
        wecm = ipcm-ecm
        #
        # get Karzas-Latter Gaunt factors
        try:
            klgfb = self.Klgfb
        except:
            self.Klgfb = util.klgfbRead()
            klgfb = self.Klgfb
        #
        nWvl = wvl.size
        nTemp = temperature.size
        #
        if verner:
            lvl1 = 1
        else:
            lvl1 = 0
            #
        if (nTemp > 1) and (nWvl > 1):
            mask = np.zeros((nlvls,nTemp,nWvl),'Bool')
            fbrate = np.zeros((nlvls,nTemp,nWvl),'float64')
            expf = np.zeros((nlvls,nTemp,nWvl),'float64')
            ratg = np.zeros((nlvls),'float64')
            fbRate = np.zeros((nTemp,nWvl),'float64')
            if verner:
                self.vernerCross(wvl)
                vCross = self.VernerCross
            #
            ratg[0] = float(mult[0])/float(multr[0])
            ipLvlEv = self.Ip - const.invCm2Ev*ecm[0]
            ipLvlErg = const.ev2Erg*ipLvlEv
            for itemp in range(nTemp):
                mask[0,itemp] = 1.e+8/wvl < (ipcm - ecm[0])
                expf[0,itemp] = np.exp((ipLvlErg - 1.e+8*const.planck*const.light/wvl)/(const.boltzmann*temperature[itemp]))
                fbrate[0,itemp] = (const.planck*const.light/(1.e-8*wvl))**5*const.verner*ratg[0]*expf[0,itemp]*vCross/temperature[itemp]**1.5
            for ilvl in range(lvl1,nlvls):
                ipLvlEv = self.Ip - const.invCm2Ev*ecm[ilvl]
                ipLvlErg = const.ev2Erg*ipLvlEv
                scaledE = np.log(const.ev2Ang/(ipLvlEv*wvl))
                thisGf = klgfb['klgfb'][pqn[ilvl]-1, l[ilvl]]
                spl = interpolate.splrep(klgfb['pe'], thisGf)
                gf = np.exp(interpolate.splev(scaledE, spl))
                ratg[ilvl] = float(mult[ilvl])/float(multr[0]) # ratio of statistical weights
            #
                for itemp in range(nTemp):
                    expf[ilvl] = np.exp((ipLvlErg - 1.e+8*const.planck*const.light/wvl)/(const.boltzmann*temperature[itemp]))
                    ipLvlErg = const.ev2Erg*ipLvlEv
                    expf[ilvl,itemp] = np.exp((ipLvlErg - 1.e+8*const.planck*const.light/wvl)/(const.boltzmann*temperature[itemp]))
                    mask[ilvl,itemp] = 1.e+8/wvl < (ipcm - ecm[ilvl])
                    fbrate[ilvl,itemp] = const.freeBound*ratg[ilvl]*(ipLvlErg**2/float(pqn[ilvl]))*gf*expf[ilvl,itemp]/(temperature[itemp]**1.5*(wvl)**2)
            fbrma = np.ma.array(fbrate)
            fbrma.mask =  mask
            fbrma.fill_value = 0.
            fbRate = (fbrma).sum(axis=0)
            fbRate.fill_value = 0.
            self.FreeBoundEmiss = {'emiss':fbRate.data, 'temperature':temperature,'wvl':wvl}
            #
        elif (nTemp == 1) and (nWvl > 1):
            mask = np.zeros((nlvls,nWvl),'Bool')
            fbrate = np.zeros((nlvls,nWvl),'float64')
            expf = np.zeros((nlvls,nWvl),'float64')
            ratg = np.zeros((nlvls),'float64')
            if verner:
                self.vernerCross(wvl)
                vCross = self.VernerCross
                #
                mask[0] = 1.e+8/wvl < (ipcm - ecm[0])
                ratg[0] = float(mult[0])/float(multr[0])
                ipLvlEv = self.Ip - const.invCm2Ev*ecm[0]
                print ' verner ipLvlEv =', ipLvlEv
                ipLvlErg = const.ev2Erg*ipLvlEv
                expf[0] = np.exp((ipLvlErg - 1.e+8*const.planck*const.light/wvl)/(const.boltzmann*temperature))
                fbrate[0] = (const.planck*const.light/(1.e-8*wvl))**5*const.verner*ratg[0]*vCross/temperature**1.5
#                print ' expf verner = ', expf[0]
#                print ' fbrate verner = ', fbrate[0]
            for ilvl in range(lvl1,nlvls):
                # scaled energy is relative to the ionization potential of each individual level
                ipLvlEv = self.Ip - const.invCm2Ev*ecm[ilvl]
#                print ' ilvl, ipLvlEv = ', ilvl, ipLvlEv
                scaledE = np.log(const.ev2Ang/(ipLvlEv*wvl))
                thisGf = klgfb['klgfb'][pqn[ilvl]-1, l[ilvl]]
                spl = interpolate.splrep(klgfb['pe'], thisGf)
                gf = np.exp(interpolate.splev(scaledE, spl))
                mask[ilvl] = 1.e+8/wvl < (ipcm - ecm[ilvl])
                ratg[ilvl] = float(mult[ilvl])/float(multr[0]) # ratio of statistical weights
                ipLvlErg = const.ev2Erg*ipLvlEv
                expf[ilvl] = np.exp((ipLvlErg - 1.e+8*const.planck*const.light/wvl)/(const.boltzmann*temperature))
                fbrate[ilvl] = const.freeBound*ratg[ilvl]*(ipLvlErg**2/float(pqn[ilvl]))*gf/(temperature**1.5*(wvl)**2)
#                print ' ilvl, scaledE = ', ilvl, scaledE
#                print ' ilvl, thisGf = ', ilvl, thisGf
#                print ' ilvl, gf = ', ilvl, gf
#                print ' ilvl ratg = ', ilvl, ratg[ilvl]
#                print ' ilvl mask = ', ilvl, mask[ilvl]
#                print ' ilvl expf = ', ilvl, expf[ilvl]
#                print ' fbrate = ', ilvl,  fbrate[ilvl]
            fbrma = np.ma.array(fbrate)
            fbrma.mask =  mask
            fbrma.fill_value = 0.
            # factor of 1.e-8 converts to Angstrom^-1, otherwise it would be cm^-1
            fbRate = (expf*fbrma).sum(axis=0)
            fbRate.fill_value = 0.
            self.FreeBoundEmiss = {'emiss':fbRate.data, 'temperature':temperature,'wvl':wvl}
        else:
            mask = np.zeros((nlvls,nTemp),'Bool')
            fbrate = np.zeros((nlvls,nTemp),'float64')
            expf = np.zeros((nlvls,nTemp),'float64')
            ratg = np.zeros((nlvls),'float64')
            if verner:
                self.vernerCross(wvl)
                vCross = self.VernerCross
            #
            mask[0] = 1.e+8/wvl < (ipcm - ecm[0])
            ratg[0] = float(mult[0])/float(multr[0])
            ipLvlEv = self.Ip - const.invCm2Ev*ecm[0]
            ipLvlErg = const.ev2Erg*ipLvlEv
            expf[0] = np.exp((ipLvlErg - 1.e+8*const.planck*const.light/wvl)/(const.boltzmann*temperature))
            fbrate[0] = (const.planck*const.light/(1.e-8*wvl))**5*const.verner*ratg[0]*expf[0]*vCross/temperature**1.5
            for ilvl in range(lvl1,nlvls):
                # scaled energy is relative to the ionization potential of each individual level
                ipLvlEv = self.Ip - const.invCm2Ev*ecm[ilvl]
                scaledE = np.log(const.ev2Ang/(ipLvlEv*wvl))
                thisGf = klgfb['klgfb'][pqn[ilvl]-1, l[ilvl]]
                spl = interpolate.splrep(klgfb['pe'], thisGf)
                gf = np.exp(interpolate.splev(scaledE, spl))
                mask[ilvl] = 1.e+8/wvl < (ipcm - ecm[ilvl])
                ratg[ilvl] = float(mult[ilvl])/float(multr[0]) # ratio of statistical weights
                ipLvlErg = const.ev2Erg*ipLvlEv
                expf[ilvl] = np.exp((ipLvlErg - 1.e+8*const.planck*const.light/wvl)/(const.boltzmann*temperature))
                fbrate[ilvl] = const.freeBound*ratg[ilvl]*(ipLvlErg**2/float(pqn[ilvl]))*expf[ilvl]*gf/(temperature**1.5*(wvl)**2)
            fbrma = np.ma.array(fbrate)
            fbrma.mask =  mask
            fbrma.fill_value = 0.
            # factor of 1.e-8 converts to Angstrom^-1, otherwise it would be cm^-1
            fbRate = (fbrma).sum(axis=0)
            fbRate.fill_value = 0.
            self.FreeBoundEmiss = {'emiss':fbRate.data, 'temperature':temperature,'wvl':wvl}
            #
            # ----------------------------------------------------------------------------
            #
    def freeBound(self, wvl, verner=1):
        '''to calculate the free-bound (radiative recombination) continuum rate coefficient of an ion, where
        the ion is taken to be the recombined iion,
        including the elemental abundance and the ionization equilibrium population
        uses the Gaunt factors of Karzas, W.J, Latter, R, 1961, ApJS, 6, 167
        for recombination to the ground level, the photoionization cross sections of Verner and Yakovlev, 1995, A&ASS, 109, 125
        are used to develop the free-bound cross section
        includes the elemental abundance and the ionization fraction
        provides emissivity = ergs cm^-2 s^-1 str^-1 Angstrom ^-1'''
        wvl = np.asarray(wvl, 'float64')
        #
        #  ip in eV, for freebound
        if self.Ion > 1 :
            self.Ip=ip[self.Z-1, self.Ion-2]
        else:
            print ' in freeBound, trying to use the neutral ion as the recombining ion'
            self.FreeBound={}
            return
        try:
            temperature = self.Temperature
        except:
            print ' temperature undefined'
            return
        # the recombined ion contains that data for fblvl
        try:
            fblvl = self.Fblvl
        except:
            fblvlname = util.zion2filename(self.Z,self.Ion-1)+'.fblvl'
            self.Fblvl = util.fblvlRead(fblvlname)
            fblvl = self.Fblvl
            # in case there is no fblvl file
            if 'errorMessage' in fblvl.keys():
                self.FreeBound = fblvl
                return
        #  need data for the current/recombining ion
        try:
            rFblvl = self.rFblvl
        except:
            if self.Ion+1 == self.Z:
                # this is a bare ion
                rFblvl = {'mult':[1., 1.]}
            else:
                rfblvlname = util.zion2filename(self.Z,self.Ion)+'.fblvl'
                self.rFblvl = util.fblvlRead(fblvlname)
                rFblvl = self.rFblvl
            if 'errorMessage' in rFblvl.keys():
                self.FreeBound = fblvl
                return
        try:
            gIoneq = self.IoneqOne
        except:
            self.ioneqOne()
            gIoneq = self.IoneqOne
        #
        try:
            abund = self.Abundance
        except:
            self.AbundanceAll = util.abundanceRead(abundancename = self.AbundanceName)
            self.Abundance = self.AbundanceAll['abundance'][self.Z-1]
            abund = self.Abundance
        #
        ipcm = self.Ip/const.invCm2Ev
        iperg = self.Ip*const.ev2Erg
        #
        nlvls = len(fblvl['lvl'])
        # pqn = principle quantum no. n
        pqn = np.asarray(fblvl['pqn'], 'float64')
        # l is angular moment quantum no. L
        l = fblvl['l']
        # energy level in inverse cm
        ecm = fblvl['ecm']
        eerg = np.asarray(ecm, 'float64')*const.invCm2Erg
        # get log of photon energy relative to the ionization potential
        mult = fblvl['mult']
        multr = rFblvl['mult']
#       ipcm = self.Ip/const.invCm2Ev
        #
        #
        wecm=1.e+8/(ipcm-ecm)
        #
        # get karzas-latter Gaunt factors
        try:
            klgfb = self.Klgfb
        except:
            self.Klgfb = util.klgfbRead()
            klgfb = self.Klgfb
        #
        nWvl = wvl.size
        nTemp = temperature.size
    # statistical weigths/multiplicities
        mult = fblvl['mult']
        multr = rFblvl['mult']
        #
        if verner:
            lvl1 = 1
        else:
            lvl1 = 0
            #
        nWvl = wvl.size
        nTemp = temperature.size
        #
        if (nTemp > 1) and (nWvl > 1):
            mask = np.zeros((nlvls,nTemp,nWvl),'Bool')
            fbrate = np.zeros((nlvls,nTemp,nWvl),'float64')
            expf = np.zeros((nlvls,nTemp,nWvl),'float64')
            ratg = np.zeros((nlvls),'float64')
            fbRate = np.zeros((nTemp,nWvl),'float64')
            if verner:
                self.vernerCross(wvl)
                vCross = self.VernerCross
                ratg[0] = float(mult[0])/float(multr[0])
                ipLvlEv = self.Ip - const.invCm2Ev*ecm[0]
                ipLvlErg = const.ev2Erg*ipLvlEv
                ipLvlErg = iperg - eerg[0]
                for itemp in range(nTemp):
                    mask[0,itemp] = 1.e+8/wvl < (ipcm - ecm[0])
                    expf[0,itemp] = np.exp((ipLvlErg - 1.e+8*const.planck*const.light/wvl)/(const.boltzmann*temperature[itemp]))
                    fbrate[0,itemp] = (const.planck*const.light/(1.e-8*wvl))**5*const.verner*gIoneq[itemp]*ratg[0]*expf[0,itemp]*vCross/temperature[itemp]**1.5
            for ilvl in range(lvl1,nlvls):
                ipLvlEv = self.Ip - const.invCm2Ev*ecm[ilvl]
                ipLvlErg = const.ev2Erg*ipLvlEv
                ipLvlErg = iperg - eerg[0]
                scaledE = np.log(const.ev2Ang/(ipLvlEv*wvl))
                thisGf = klgfb['klgfb'][pqn[ilvl]-1, l[ilvl]]
                spl = interpolate.splrep(klgfb['pe'], thisGf)
                gf = np.exp(interpolate.splev(scaledE, spl))
                ratg[ilvl] = float(mult[ilvl])/float(multr[0]) # ratio of statistical weights
            #
                for itemp in range(nTemp):
                    expf[ilvl] = np.exp((ipLvlErg - 1.e+8*const.planck*const.light/wvl)/(const.boltzmann*temperature[itemp]))
                    ipLvlErg = const.ev2Erg*ipLvlEv
                    expf[ilvl,itemp] = np.exp((ipLvlErg - 1.e+8*const.planck*const.light/wvl)/(const.boltzmann*temperature[itemp]))
                    mask[ilvl,itemp] = 1.e+8/wvl < (ipcm - ecm[ilvl])
                    fbrate[ilvl,itemp] = const.freeBound*gIoneq[itemp]*ratg[ilvl]*(ipLvlErg**2/float(pqn[ilvl]))*gf*expf[ilvl,itemp]/(temperature[itemp]**1.5*(wvl)**2)
            fbrma = np.ma.array(fbrate)
            fbrma.mask =  mask
            fbrma.fill_value = 0.
            fbRate = (fbrma).sum(axis=0)
            fbRate.fill_value = 0.
            self.FreeBound = {'rate':abund*fbRate.data, 'temperature':temperature,'wvl':wvl}
            #
        elif (nTemp == 1) and (nWvl > 1):
            mask = np.zeros((nlvls,nWvl),'Bool')
            fbrate = np.zeros((nlvls,nWvl),'float64')
            expf = np.zeros((nlvls,nWvl),'float64')
            ratg = np.zeros((nlvls),'float64')
            if verner:
                self.vernerCross(wvl)
                vCross = self.VernerCross
                # mask is true for bad values
                mask[0] = 1.e+8/wvl < (ipcm - ecm[0])
                ratg[0] = float(mult[0])/float(multr[0])
                ipLvlEv = self.Ip - const.invCm2Ev*ecm[0]
                ipLvlErg = const.ev2Erg*ipLvlEv
#               ipLvlErg = iperg - eerg[0]
#               print ' verner ipLvlErg = ', ipLvlErg,  ipLvlErg/(const.boltzmann*temperature)
                expf[0] = np.exp((ipLvlErg - 1.e+8*const.planck*const.light/wvl)/(const.boltzmann*temperature))
                fbrate[0] = (const.planck*const.light/(1.e-8*wvl))**5*const.verner*ratg[0]*expf[0]*vCross/temperature**1.5
            #
            for ilvl in range(lvl1,nlvls):
            # scaled energy is relative to the ionization potential of each individual level
                ipLvlEv = self.Ip - const.invCm2Ev*ecm[ilvl]
                scaledE = np.log(const.ev2Ang/(ipLvlEv*wvl))
                thisGf = klgfb['klgfb'][pqn[ilvl]-1, l[ilvl]]
                spl = interpolate.splrep(klgfb['pe'], thisGf)
                gf = np.exp(interpolate.splev(scaledE, spl))
                mask[ilvl] = 1.e+8/wvl < (ipcm - ecm[ilvl])
                ratg[ilvl] = float(mult[ilvl])/float(multr[0]) # ratio of statistical weights
                ipLvlErg = const.ev2Erg*ipLvlEv
#               ipLvlErg = iperg - eerg[ilvl]
#               print ' ilvl epLvlErg = ', ipLvlErg,  ipLvlErg/(const.boltzmann*temperature)
                expf[ilvl] = np.exp((ipLvlErg - 1.e+8*const.planck*const.light/wvl)/(const.boltzmann*temperature))
                fbrate[ilvl] = const.freeBound*ratg[ilvl]*(ipLvlErg**2/float(pqn[ilvl]))*expf[ilvl]*gf/(temperature**1.5*(wvl)**2)
#               fbrate[ilvl] = const.freeBound*ratg[ilvl]*(eerg[ilvl]**2/float(pqn[ilvl]))*expf[ilvl]*gf/(temperature**1.5*(wvl)**2)
            fbrma = np.ma.array(fbrate)
            fbrma.mask =  mask
            fbrma.fill_value = 0.
            # factor of 1.e-8 converts to Angstrom^-1, otherwise it would be cm^-1
            fbRate = abund*gIoneq*fbrma.sum(axis=0)
            fbRate.fill_value = 0.
            self.FreeBound = {'rate':fbRate.data, 'temperature':temperature,'wvl':wvl}
        #elif (nTemp > 1) and (nWvl == 1):
        else:
            mask = np.zeros((nlvls,nTemp),'Bool')
            fbrate = np.zeros((nlvls,nTemp),'float64')
            expf = np.zeros((nlvls,nTemp),'float64')
            ratg = np.zeros((nlvls),'float64')
            if verner:
                self.vernerCross(wvl)
                vCross = self.VernerCross
            mask[0] = 1.e+8/wvl < (ipcm - ecm[0])
            ratg[0] = float(mult[0])/float(multr[0])
            ipLvlEv = self.Ip - const.invCm2Ev*ecm[0]
            ipLvlErg = const.ev2Erg*ipLvlEv
            expf[0] = np.exp((ipLvlErg - 1.e+8*const.planck*const.light/wvl)/(const.boltzmann*temperature))
            fbrate[0] = (const.planck*const.light/(1.e-8*wvl))**5*const.verner*ratg[0]*expf[0]*vCross/temperature**1.5
            for ilvl in range(lvl1,nlvls):
                # scaled energy is relative to the ionization potential of each individual level
                ipLvlEv = self.Ip - const.invCm2Ev*ecm[ilvl]
                scaledE = np.log(const.ev2Ang/(ipLvlEv*wvl))
                thisGf = klgfb['klgfb'][pqn[ilvl]-1, l[ilvl]]
                spl = interpolate.splrep(klgfb['pe'], thisGf)
                gf = np.exp(interpolate.splev(scaledE, spl))
                mask[ilvl] = 1.e+8/wvl < (ipcm - ecm[ilvl])
                ratg[ilvl] = float(mult[ilvl])/float(multr[0]) # ratio of statistical weights
                ipLvlErg = const.ev2Erg*ipLvlEv
                expf[ilvl] = np.exp((ipLvlErg - 1.e+8*const.planck*const.light/wvl)/(const.boltzmann*temperature))
                fbrate[ilvl] = const.freeBound*ratg[ilvl]*(ipLvlErg**2/float(pqn[ilvl]))*expf[ilvl]*gf/(temperature**1.5*(wvl)**2)
            fbrma = np.ma.array(fbrate)
            fbrma.mask =  mask
            fbrma.fill_value = 0.
            # factor of 1.e-8 converts to Angstrom^-1, otherwise it would be cm^-1
            fbRate = abund*gIoneq*(fbrma).sum(axis=0)
            fbRate.fill_value = 0.
            self.FreeBound = {'rate':fbRate.data, 'temperature':temperature,'wvl':wvl}
        #
        # ----------------------------------------------------------------------------------------
        #
            #
    def freeBoundLoss(self):
        '''to calculate the free-bound (radiative recombination) energy loss rate coefficient of an ion,
        the ion is taken to be the recombined iion,
        including the elemental abundance and the ionization equilibrium population
        uses the Gaunt factors of Karzas, W.J, Latter, R, 1961, ApJS, 6, 167
        provides rate = ergs cm^-2 s^-1 '''
        try:
            temperature = self.Temperature
        except:
            print ' temperature undefined'
            return
        #
        if self.Ion > 1:
            self.Ip=ip[self.Z-1, self.Ion-1]
        else:
            print ' in freeBound, trying to use the neutral ion as the recombining ion'
            self.FreeBound={}
            return
        #
        try:
            fblvl = self.Fblvl
        except:
            fblvlname = util.zion2filename(self.Z,self.Ion-1)+'.fblvl'
            self.Fblvl = util.fblvlRead(fblvlname)
            fblvl = self.Fblvl
    #  need some data for the recombining/target ion
        try:
            rFblvl = self.rFblvl
        except:
            if self.Ion == self.Z+1:
                # then we looking for the bare ion
                rFblvl = {'mult':[1., 1.]}
            else:
                rfblvlname = util.zion2filename(self.Z,self.Ion)+'.fblvl'
                self.rFblvl = util.fblvlRead(rfblvlname)
                rFblvl = self.rFblvl
        try:
            gIoneq = self.IoneqOne
        except:
            self.ioneqOne()
            gIoneq = self.IoneqOne
        #
        try:
            abund = self.Abundance
        except:
            self.AbundanceAll = util.abundanceRead(abundancename = self.AbundanceName)
            self.Abundance = self.AbundanceAll['abundance'][self.Z-1]
            abund = self.Abundance
            #
            #ipcm = self.Ip/const.invCm2Ev
            #
            nlvls = len(fblvl['lvl'])
            # pqn = principle quantum no. n
            pqn = fblvl['pqn']
            # l is angular moment quantum no. L
            l = fblvl['l']
            # energy level in inverse cm
            ecm = fblvl['ecm']
            # get log of photon energy relative to the ionization potential
            mult = fblvl['mult']
            multr = rFblvl['mult']
            #
            #
            #wecm=1.e+8/(ipcm-ecm)
            #
            # get karzas-latter Gaunt factors
            try:
                klgfb = self.Klgfb
            except:
                self.Klgfb = util.klgfbRead()
                klgfb = self.Klgfb
            #
            nTemp = temperature.size
        # statistical weigths/multiplicities
            mult = fblvl['mult']
            multr = rFblvl['mult']
            #
            #
            #wecm=1.e+8/(ipcm-ecm)
            #
        fbrate = np.zeros((nlvls,nTemp),'float64')
        ratg = np.zeros((nlvls),'float64')
        for ilvl in range(nlvls):
            # scaled energy is relative to the ionization potential of each individual level
            # will add the average energy of a free electron to this to get typical photon energy to
            # evaluate the gaunt factor
            hnuEv = 1.5*const.boltzmann*temperature/const.ev2Erg
            ipLvlEv = self.Ip - const.invCm2Ev*ecm[ilvl]
            scaledE = np.log(hnuEv/ipLvlEv)
            thisGf = klgfb['klgfb'][pqn[ilvl]-1, l[ilvl]]
            spl = interpolate.splrep(klgfb['pe'], thisGf)
            gf = np.exp(interpolate.splev(scaledE, spl))
            ratg[ilvl] = float(mult[ilvl])/float(multr[0]) # ratio of statistical weights
            ipLvlErg = const.ev2Erg*ipLvlEv
            fbrate[ilvl] = ratg[ilvl]*(ipLvlErg**2/float(pqn[ilvl]))*gf/np.sqrt(temperature)
        fbRate = abund*gIoneq*const.freeBoundLoss*(fbrate.sum(axis=0))
        self.FreeBoundLoss = {'rate':fbRate, 'temperature':temperature}
        #
        # ----------------------------------------------------------------------------------------
        #
    def vernerCross(self,wvl):
        '''Calculates the photoionization cross section.

        The data are from Verner and Yakovlev
        the cross section refers to the next lower ionization stage'''
        try:
            vernerDat = self.VernerDat
        except:
            self.VernerDat = util.vernerRead()
            vernerDat = self.VernerDat
        z = self.Z
        stage = self.Ion
        ip = self.Ip
        ipcm = self.Ip/const.invCm2Ev
        ecm = self.Fblvl['ecm']
        #
        en = const.ev2Ang/wvl
        y = en/vernerDat['e0'][z,stage-1]
        fy= vernerDat['sig0'][z,stage-1]*((y - 1.)**2 + vernerDat['yw'][z,stage-1]**2) * y**(-5.5 - vernerDat['l'][z,stage-1] + 0.5*vernerDat['p'][z,stage-1]) * (1. + np.sqrt(y/vernerDat['ya'][z,stage-1]))**(-vernerDat['p'][z,stage-1])
#       mask = en < vernerDat['eth'][z,stage]
        # will use Chianti values for energy of ground level
        mask = (1.e+8/wvl) < (ipcm - ecm[0])
        vCross = np.ma.array(fy)
        vCross.mask = mask
        vCross.fill_value = 0.
        # cross-section will be output in cm^2
        self.VernerCross = vCross*1.e-18
        #
        # ----------------------------------------------------------------------------------------
        #
    def freeFree(self, wvl):
        '''Calculates the free-free emission for a single ion.

        Includes elemental abundance and ionization equilibrium population.
        Uses Itoh where valid and Sutherland elsewhere'''
        if self.Ion == 1:
            self.FreeFree = {'errorMessage':' freefree is not produced by neutrals'}
        else:
            wvl = np.asarray(wvl, 'float64')
            ffs = self.sutherland(wvl)
            #
            ffi = self.itoh(wvl)
            ff = ffs['suthFf']
            if not 'errorMessage' in ffi.keys():
                iff = ffi['itohFf']
                itohMask = np.logical_not(iff.mask)
                ff[itohMask] = iff[itohMask]
        #
        try:
            gIoneq = self.IoneqOne
        except:
            self.ioneqOne()
            gIoneq = self.IoneqOne
        if type(gIoneq) == types.FloatType:
            # only one temperature specified
            if gIoneq == 0.:
                ffRate = np.zeros(wvl.size)
                self.FreeFree = {'rate':ffRate, 'temperature':self.Temperature,'wvl':wvl}
                return
        else:
            if gIoneq.sum() == 0.:
                ffRate = np.zeros((self.Temperature.size, wvl.size), 'float64')
                self.FreeFree = {'rate':ffRate, 'temperature':self.Temperature,'wvl':wvl}
                return
#       print ' gIoneq = ', gIoneq
        if wvl.size > 1:
            gIoneq = gIoneq.repeat(wvl.size).reshape(self.Temperature.size,wvl.size)
            #
            try:
                abund = self.Abundance
            except:
                self.AbundanceAll = util.abundanceRead(abundancename = self.AbundanceName)
                self.Abundance = self.AbundanceAll['abundance'][self.Z-1]
                abund = self.Abundance
                #
            ffRate = (const.freeFree*(self.Z)**2*abund*gIoneq*ff).squeeze()
            self.FreeFree = {'rate':ffRate, 'temperature':self.Temperature,'wvl':wvl}
        #
        # ----------------------------------------------------------------------------------------
        #
    def freeFreeEmiss(self, wvl):
        '''Calculates the free-free emissivity for a single ion.
        does not include element abundance or ionization fraction
        Uses Itoh where valid and Sutherland elsewhere'''
        if self.Ion == 1:
            self.FreeFreeEmiss = {'errorMessage':' freefree is not produced by neutrals'}
        else:
            wvl = np.asarray(wvl, 'float64')
            ffs = self.sutherland(wvl)
            ff = ffs['suthFf']
            ffi = self.itoh(wvl)
            if 'errorMessage' not in ffi.keys():
                iff = ffi['itohFf']
                itohMask = np.logical_not(iff.mask)
                ff[itohMask] = iff[itohMask]
            #
            ffRate = (const.freeFree*(self.Z)**2*ff).squeeze()
            self.FreeFree = {'rate':ffRate, 'temperature':self.Temperature,'wvl':wvl}
        #
        # ----------------------------------------------------------------------------------------
        #
    def freeFreeLoss(self):
        '''Calculates the total free-free emission for a single ion.

        Includes elemental abundance and ionization equilibrium population.

        Uses Itoh where valid and Sutherland elsewhere'''
        #
        if self.Ion == 1:
            self.FreeFree = {'errorMessage':' freefree is not produced by neutrals'}
        else:
            try:
                temperature = self.Temperature
            except:
                print ' temperature undefined'
                return
            try:
                gffint = self.Gffint['gffint']
                g2 = self.Gffint['g2']
            except:
                self.Gffint = self.gffintRead()
                gffint = self.Gffint['gffint']
                g2 = self.Gffint['g2']
            #
            gamma2 = self.Ip*const.ev2Erg/(const.boltzmann*temperature)
            #
            spl = interpolate.splrep(g2, gffint)
            gff = interpolate.splev(np.log(gamma2), spl)
            #
            try:
                gIoneq = self.IoneqOne
            except:
                self.ioneqOne()
                gIoneq = self.IoneqOne
            #
            try:
                abund = self.Abundance
            except:
                self.AbundanceAll = util.abundanceRead(abundancename = self.AbundanceName)
                self.Abundance = self.AbundanceAll['abundance'][self.Z-1]
                abund = self.Abundance
                #
            ffRate = const.freeFreeloss*(self.Z)**2*abund*gIoneq*gff*np.sqrt(temperature)
            self.FreeFreeLoss = {'rate':ffRate, 'temperature':temperature}
        #
        # ----------------------------------------------------------------------------------------
        #
    def klgfbInterp(self, wvl, n, l):
        '''A Python version of the CHIANTI IDL procedure karzas_xs.

        Interpolates free-bound gaunt factor of Karzas and Latter, (1961, Astrophysical Journal
        Supplement Series, 6, 167) as a function of wavelength (wvl).'''
        try:
            klgfb = self.Klgfb
        except:
            self.Klgfb = util.klgfbRead()
            klgfb = self.Klgfb
        # get log of photon energy relative to the ionization potential
        sclE = np.log(self.Ip/(wvl*const.ev2ang))
        thisGf = klgfb['klgfb'][n-1, l]
        spl = interpolate.splrep(klgfb['pe'], thisGf)
        gf = interpolate.splev(sclE, spl)
        return gf

        #
        # ----------------------------------------------------------------------------------------
        #
    def itoh(self, wvl):
        '''Calculates free-free emission with the free-free gaunt factors of Itoh et al. (ApJS 128, 125, 2000).

        the relativistic values are valid for 6. < log10(T) < 8.5 and -4. < log10(u) < 1.'''
        wvl = np.array(wvl, 'float64')
        try:
            itohCoef = self.ItohCoef
        except:
            self.ItohCoef = util.itohRead()['itohCoef'][self.Z].reshape(11, 11)
            itohCoef = self.ItohCoef
        try:
            t = (np.log10(self.Temperature) -7.25)/1.25
        except:
            errorMessage = ' temperature undefined in continuum.itoh'
            print(errorMessage)
            return {'errorMessage':errorMessage}
        if type(self.Temperature) == types.FloatType:
            nTemp = 1
        else:
            nTemp = self.Temperature.size
        #
        nWvl = wvl.size
        #
        if (nTemp > 1) and (nWvl > 1):
            u = (const.planck*const.light*1.e+8/const.boltzmann)*np.outer(1./self.Temperature, 1./wvl )
            lU = (np.log10(u) + 1.5)/2.5
            lT = (np.log10(self.Temperature) -7.25)/1.25
            g = np.zeros((nTemp, nWvl), 'float64')
            rad = np.ma.zeros((nTemp, nWvl), 'float64')
            for itemp in range(nTemp):
                for j in range(11):
                    for i in range(11):
                        g[itemp]+= itohCoef[i,j]*(lT[itemp]**i)*(lU[itemp]**j)
                rad[itemp] = (1./wvl)**2*g[itemp]*np.exp(-u[itemp])/np.sqrt(self.Temperature[itemp])
            tArray = np.zeros((1, len(self.Temperature)), 'float64')
            tArray[0] = self.Temperature
            t2Array = np.repeat(tArray.transpose(), len(wvl), axis=1)
            nonValidT1 = np.log10(t2Array) < 6.
            nonValidT2 = np.log10(t2Array) > 8.5
            nonValidT = np.logical_or(nonValidT1, nonValidT2)
            nonValidU1 = np.log10(u) < -4.
            nonValidU2 = np.log10(u) > 1.
            nonValidU = np.logical_or(nonValidU1, nonValidU2)
            nonValid = np.logical_or(nonValidT, nonValidU)
            rad.mask = nonValid
            rad.set_fill_value(0.)
            g=np.ma.array(g, mask=nonValid, fill_value=0.)
            return {'itohGff':g, 'itohFf':rad}
        elif (nTemp > 1) and (nWvl == 1):
            u = (const.planck*const.light*1.e+8/const.boltzmann)/(self.Temperature*wvl )
            lU = (np.log10(u) + 1.5)/2.5
            lT = (np.log10(self.Temperature) -7.25)/1.25
            g = np.zeros((nTemp), 'float64')
            rad = np.ma.zeros((nTemp), 'float64')
            for itemp in range(nTemp):
                for j in range(11):
                    for i in range(11):
                        g[itemp]+= itohCoef[i,j]*(lT[itemp]**i)*(lU[itemp]**j)
                rad[itemp] = (1./wvl)**2*g[itemp]*np.exp(-u[itemp])/np.sqrt(self.Temperature[itemp])
            nonValidT1 = np.log10(self.Temperature) < 6.
            nonValidT2 = np.log10(self.Temperature) > 8.5
            nonValidT = np.logical_or(nonValidT1, nonValidT2)
            nonValidU1 = np.log10(u) < -4.
            nonValidU2 = np.log10(u) > 1.
            nonValidU = np.logical_or(nonValidU1, nonValidU2)
            nonValid = np.logical_or(nonValidT, nonValidU)
            rad.mask = nonValid
            rad.set_fill_value(0.)
            g=np.ma.array(g, mask=nonValid, fill_value=0.)
            return {'itohGff':g, 'itohFf':rad}
        elif (nTemp == 1) and (nWvl > 1):
            if (np.log10(self.Temperature) < 6.) or (np.log10(self.Temperature > 8.5)):
                errorMessage ='invalid temperature in continuum.itoh()'
                return {'errorMessage':errorMessage}
            else:
                u = (const.planck*const.light*1.e+8/const.boltzmann)/(self.Temperature*wvl )
                lU = (np.log10(u) + 1.5)/2.5
                lT = (np.log10(self.Temperature) -7.25)/1.25
                g = np.zeros(nWvl, 'float64')
                rad = np.ma.zeros((nWvl), 'float64')
                for j in range(11):
                    for i in range(11):
                        g+= itohCoef[i,j]*(lT**i)*(lU**j)
                rad = np.ma.array((1./wvl)**2*g*np.exp(-u)/np.sqrt(self.Temperature), 'Float64')
                nonValidU1 = np.log10(u) < -4.
                nonValidU2 = np.log10(u) > 1.
                nonValidU = np.logical_or(nonValidU1, nonValidU2)
                nonValid = nonValidU
                rad.mask = nonValid
                rad.set_fill_value(0.)
                g=np.ma.array(g, mask=nonValid, fill_value=0.)
                return {'itohGff':g, 'itohFf':rad}
        elif (nTemp == 1) and (nWvl == 1):
            u = (const.planck*const.light*1.e+8/const.boltzmann)/(self.Temperature*wvl )
            if (np.log10(self.Temperature) < 6.) or (np.log10(self.Temperature > 8.5)):
                errorMessage ='invalid temperature in continuum.itoh()'
                return {'errorMessage':errorMessage}
            elif (np.log10(u) < -4.) or (np.log10(u) > 8.5):
                errorMessage ='invalid temperature/wavelength in continuum.itoh()'
                return {'errorMessage':errorMessage}
            else:
                u = (const.planck*const.light*1.e+8/const.boltzmann)/(self.Temperature*wvl )
                lU = (np.log10(u) + 1.5)/2.5
                lT = (np.log10(self.Temperature) -7.25)/1.25
                g = np.zeros(nWvl, 'float64')
                rad = np.ma.zeros((nWvl), 'float64')
                for j in range(11):
                    for i in range(11):
                        g+= itohCoef[i,j]*(lT**i)*(lU**j)
                rad = np.ma.array((1./wvl)**2*g*np.exp(-u)/np.sqrt(self.Temperature), 'Float64')
                nonValidU1 = np.log10(u) < -4.
                nonValidU2 = np.log10(u) > 1.
                nonValidU = np.logical_or(nonValidU1, nonValidU2)
                nonValid = nonValidU
                rad.mask = nonValid
                rad.set_fill_value(0.)
                return {'itohGff':g, 'itohFf':rad}
        #
        # - - - - - - - - - - - - - - - - - - - - - - -
        #
    def sutherland(self, wvl):
        '''Calculates the free-free continuum using the free-free gaunt factor calculations of Sutherland, 1998, MNRAS, 300, 321.

        the wavelengths (wvl) will be sorted, first thing'''
        #
        wvl = np.array(wvl, 'float64')
        nWvl = wvl.size
#       factor = 5.44436e-39
        try:
            temperature = self.Temperature
        except:
            errorMessage = ' temperature undefined in continuum.sutherland'
            print(errorMessage)
            return {'errorMessage':errorMessage}
        #  read in the gaunt factors, if necessary and get interpolator
        try:
            gffInterpolator = self.GffInterpolator
        except:
            self.Gff = util.gffRead()
            gff = self.Gff
            iu=(np.log10(gff['u1d']) + 4.)*10.
            ig=(np.log10(gff['g21d']) + 4.)*5.
            gaunt = gff['gff']
            #tr=Triangulation(iu.flatten(),ig.flatten())
            tr=Triangulation(iu,ig)
            self.GffInterpolator = tr.nn_interpolator(gaunt.flatten())
            gffInterpolator = self.GffInterpolator
    #
        gga = np.array((float(self.Z)**2*const.ryd2erg/const.boltzmann)*(1./temperature),'float64')
        nonValidGg1 = np.log10(gga) < -4.
        nonValidGg2 = np.log10(gga) > 4.
        nonValidGg = np.logical_or(nonValidGg1, nonValidGg2)
        ggOut = np.ma.array(gga, mask = nonValidGg, fill_value=True)
        iGg = np.ma.array((np.log10(gga) + 4.)*5., mask=nonValidGg,  fill_value=0.)
        #
        if nonValidGg.sum():
            errorMessage = 'no valid temperatures in continuum.sutherland'
            print(errorMessage)
            return {'errorMessage':errorMessage}
        else:
                #iUu = np.ma.array((np.log10(uu) + 4.)*10., mask=nonValidUu,  fill_value=0.)
        #iGg = (np.log10(gg) + 4.)*5.
        #print ' iGg.shape = ',iGg, iGg.shape
            #
            nWvl = wvl.size
            nTemp = temperature.size
            #
            if (nTemp > 1) and (nWvl > 1):
                ff = np.ma.zeros((nWvl, nTemp), 'float64')
                gffOut1 = np.ma.zeros((nWvl, nTemp), 'float64')
                gffOutMask = np.zeros((nWvl, nTemp), 'Bool')
                uuOut = np.zeros((nWvl, nTemp), 'float64')
                for iwvl in range(nWvl):
                    uu = ((const.planck*const.light*1.e+8/const.boltzmann)/(wvl[iwvl]*temperature))  #.flatten()
                    nonValidUu1 = np.log10(uu) < -4.
                    nonValidUu2 = np.log10(uu) > 4.
                    nonValidUu = np.logical_or(nonValidUu1, nonValidUu2)
                    gffOutMask[iwvl] = nonValidUu
                    uuOut[iwvl] = np.ma.array(uu, mask=nonValidUu, fill_value=True)
                    iUu = np.ma.array((np.log10(uu) + 4.)*10., mask=nonValidUu,  fill_value=0.)
                    gffOut1[iwvl] = gffInterpolator(iUu, iGg)
                    wvlt = 1./(wvl[iwvl]**2*np.sqrt(temperature))  # was sortedTemperature
                    ff[iwvl] = (np.exp(-uuOut[iwvl])*gffOut1[iwvl]*wvlt)
                gffOut1.mask = gffOutMask
                gffOut1.set_fill_value(0.)
                gffOut = gffOut1.transpose()
                ff.mask = gffOutMask
                ff.set_fill_value(0.)
                return {'suthFf':ff.transpose(), 'suthGff':gffOut}
            #
            if (nTemp == 1) and (nWvl > 1):
                uu = ((const.planck*const.light*1.e+8/const.boltzmann)/(wvl*temperature))  # .flatten()
                nonValidUu1 = np.log10(uu) < -4.
                nonValidUu2 = np.log10(uu) > 4.
                nonValidUu = np.logical_or(nonValidUu1, nonValidUu2)
                gffOutMask = nonValidUu
                iUu = (np.log10(uu) + 4.)*10.
                gffOut1 = gffInterpolator(iUu, iGg.repeat(nWvl))
                wvlt = 1./(wvl**2*np.sqrt(temperature))
                ff = np.ma.array(np.exp(-uu)*gffOut1*wvlt)
                ff.mask=gffOutMask
                ff.set_fill_value(0.)
                gffOut = np.ma.array(gffOut1, mask=gffOutMask, fill_value=0.)
                return {'suthFf':ff, 'suthGff':gffOut, 'iUu':iUu, 'gffOut1':gffOut1, 'wvlt':wvlt,  'iGg':iGg.repeat(nWvl), 'gffInterpolator':gffInterpolator}
        #elif (nTemp > 1) and (nWvl == 1):
            else:
                #print ' igg.shape = ',iGg.shape
                #gffOut1 = np.ma.zeros((nTemp), 'float64')
                #gffOutMask = np.zeros((nTemp), 'Bool')
                #uuOut = np.zeros((nTemp), 'float64')
                #
                uu = (const.planck*const.light*1.e+8/const.boltzmann) /(wvl*temperature).flatten()
                nonValidUu1 = np.log10(uu) < -4.
                nonValidUu2 = np.log10(uu) > 4.
                nonValidUu = np.logical_or(nonValidUu1, nonValidUu2)
                gffOutMask = nonValidUu
                uuOut = np.ma.array(uu, mask=nonValidUu, fill_value=True)
                #iUu = np.ma.array((np.log10(uu) + 4.)*10., mask=nonValidUu,  fill_value=0.)
                iUu = (np.log10(uu) + 4.)*10.
                #print ' iUu.shape = ',iUu.shape
                gffOut1 = gffInterpolator(iUu, iGg.flatten())
                #
                wvlt = 1./(wvl**2*np.sqrt(temperature))
                ff1 = np.exp(-uuOut)*gffOut1*wvlt
                ff = np.ma.array(ff1, mask=gffOutMask, fill_value=0.)
                gffOut = np.ma.array(gffOut1, mask=gffOutMask, fill_value=0.)
        return {'suthFf':ff, 'suthGff':gffOut}
            #elif (nTemp == 1) and (nWvl == 1):
        ##iGg = (np.log10(gg) + 4.)*5.
        #uu = (const.planck*const.light*1.e+8/const.boltzmann)/(wvl*temperature)
        #nonValidUu1 = np.log10(uu) < -4.
        #nonValidUu2 = np.log10(uu) > 4.
        #nonValidUu = np.logical_or(nonValidUu1, nonValidUu2)
        #if not nonValidUu:
            #iUu = (np.log10(uu) + 4.)*10.
            #gffOut = gffInterpolator(iUu, iGg)
            #wvlt = 1./(wvl**2*np.sqrt(temperature))
            #ff = np.exp(-uu)*gffOut*wvlt
            #return {'suthFf':ff, 'suthGff':gffOut}
        #else:
            #errorMessage = 'invalid Temperature/Wavlength in continuum.sutherland'
            #print(errorMessage)
            #return {'errorMessage':errorMessage}
        #
        # ----------------------------------------------------------------------------------------
        #
    def ioneqOne(self):
        '''Determine the ionization equilibrium for the selected ion as a function of temperature.'''
        #
        try:
            temperature = self.Temperature
        except:
            print ' self.Temperature undefined in ioneqOne'
            return
        #
        try:
            ioneqAll = self.IoneqAll
        except:
            self.IoneqAll = util.ioneqRead(ioneqname = self.Defaults['ioneqfile'])
            ioneqAll=self.IoneqAll
        #
        ioneqTemperature = ioneqAll['ioneqTemperature']
        Z=self.Z
        Ion=self.Ion
        Dielectronic=self.Dielectronic
        #
        thisIoneq=ioneqAll['ioneqAll'][Z-1,Ion-1].squeeze()
#        thisIoneq = self.Ioneq
        gioneq = thisIoneq > 0.
        y2 = interpolate.splrep(np.log(ioneqTemperature[gioneq]),np.log(thisIoneq[gioneq]),s=0)
        goodt1 = self.Temperature >= ioneqTemperature[gioneq].min()
        goodt2 = self.Temperature <= ioneqTemperature[gioneq].max()
        goodt = np.logical_and(goodt1,goodt2)
        #
        if goodt.sum() > 0:
            gIoneq=interpolate.splev(np.log(self.Temperature),y2)   #,der=0)
            gIoneq=np.exp(gIoneq)
        else:
            gIoneq=0.
        #
        self.IoneqOne=gIoneq
        return  # gIoneq

        #
        # -------------------------------------------------------------------------------------
        #
class ion:
    '''The top level class for performing spectral calculations for an ion in the CHIANTI database.

    ionStr is a string corresponding such as 'c_5' that corresponds to the C VI ion.
    temperature in Kelvin
    density in cm^-3
    radTemperature, the radiation black-body temperature in Kelvin
    rPlot, the distance from the center of the star in stellar radii
    '''
    def __init__(self,ionStr,temperature=None,density=None,pDensity='default', radTemperature=0,rPhot=1.,verbose=0, setup=True):
        #
        #
        self.__version__ = chianti.__version__
        self.IonStr=ionStr
        self.Defaults=defaults
        self.AbundanceName = defaults['abundfile']
        self.IoneqName = defaults['ioneqfile']
        #
        self.Z=util.convertName(ionStr)['Z']
        self.Ion=util.convertName(ionStr)['Ion']
        self.Dielectronic=util.convertName(ionStr)['Dielectronic']
        self.Spectroscopic=util.zion2spectroscopic(self.Z,self.Ion)
        self.FileName=util.zion2filename(self.Z, self.Ion,dielectronic=self.Dielectronic )
        #
        self.RadTemperature = radTemperature
        self.RPhot = rPhot
        #
        #  ip in eV, but don't read for bare ions
        if self.Ion <= self.Z:
            self.Ip=ip[self.Z-1, self.Ion-1-self.Dielectronic]
        #
        if type(temperature) != types.NoneType:
            self.Temperature = np.array(temperature,'float64')
        #
        #
        self.AbundanceAll = util.abundanceRead(abundancename = self.AbundanceName)
        self.Abundance = self.AbundanceAll['abundance'][self.Z-1]
        #
        self.IoneqAll = util.ioneqRead(ioneqname = self.Defaults['ioneqfile'])
        self.ioneqOne()
        #
        #  this needs to go after setting temperature and reading ionization equilibria
        if pDensity == 'default':
            self.p2eRatio()
        #
        if type(density) != types.NoneType:
            self.Density = np.asarray(density,'float64')
            #
            if pDensity == 'default':
                self.PDensity = self.ProtonDensityRatio*self.Density
            else:
                self.PDensity = pDensity
        if setup:
            #
            # read in all data if in masterlist
            #  if not, there should still be ionization and recombination rates
            if self.IonStr in MasterList:
                self.Elvlc = util.elvlcRead(self.IonStr, verbose=verbose)
                self.Wgfa = util.wgfaRead(self.IonStr)
                self.Nwgfa=len(self.Wgfa['lvl1'])
                self.Splups = util.splupsRead(self.IonStr)
                self.Nsplups=len(self.Splups['lvl1'])
                #
                # need to determine the number of levels that can be populated
                nlvlElvlc = len(self.Elvlc['lvl'])
                nlvlWgfa = max(self.Wgfa['lvl2'])
                nlvlSplups = max(self.Splups['lvl2'])
                self.Nlvls = min([nlvlElvlc, nlvlWgfa, nlvlSplups])
#               self.Cilvl = util.cireclvlRead(self.IonStr, 'cilvl')
                #  .cisplups may not exist
                #
                self.CiSplups = util.splupsRead(self.IonStr,ci=1)
                if type(self.CiSplups) != types.NoneType:
                    self.Ncisplups=len(self.CiSplups["lvl1"])
                else:
                    self.Ncisplups = -1
                #  .reclvl file may not exist
                self.Reclvl = util.cireclvlRead(self.IonStr, 'reclvl')
                if type(self.Reclvl) != types.NoneType:
                    self.Nreclvl = len(self.Reclvl['lvl1'])
                else:
                    self.Nreclvl = -1
                #
                #  psplups file may not exist
                self.Psplups = util.splupsRead(self.IonStr, prot=True)
                if type(self.Psplups) != types.NoneType:
                    self.Npsplups=len(self.Psplups["lvl1"])
                else:
                    self.Npsplups = -1
        #
        # ------------------------------------------------------------------------------
        #
    def diCross(self, energy=None, verbose=False):
        '''Calculate the direct ionization cross section.

        Given as a function of the incident electron energy in eV, puts values into DiCross'''
        iso=self.Z -self.Ion + 1
        if type(energy) == types.NoneType:
            btenergy=0.1*np.arange(10)
            btenergy[0]=0.01
            dum=np.ones(len(btenergy))
            [energy, dum]=util.descale_bti(btenergy, dum, 2., self.Ip)
        energy=np.asarray(energy, 'float64')
        #
        if iso == 1 and self.Z >= 6:
            #  hydrogenic sequence
            ryd=27.2113845/2.
            u=energy/self.Ip
            ev1ryd=self.Ip/ryd
            a0=0.5291772108e-8
            a_bohr=const.pi*a0**2   # area of bohr orbit
            if self.Z >= 20:
                ff = (140.+(self.Z/20.)**3.2)/141.
            else:
                ff = 1.
            qr = util.qrp(self.Z,u)*ff
            bb = 1.  # hydrogenic
            qh = bb*a_bohr*qr/ev1ryd**2
            self.DiCross = {'energy':energy, 'cross':qh}
        elif iso == 2 and self.Z >= 10:
            #  use
            ryd=27.2113845/2.
            u=energy/self.Ip
            ev1ryd=self.Ip/ryd
            a0=0.5291772108e-8
            a_bohr=const.pi*a0**2   # area of bohr orbit
            if self.Z >= 20:
                ff=(140.+(self.Z/20.)**3.2)/141.
            else:
                ff=1.
            qr=util.qrp(self.Z,u)*ff
            bb=2.  # helium-like
            qh=bb*a_bohr*qr/ev1ryd**2
            self.DiCross={'energy':energy, 'cross':qh}
        else:
            try:
                diparams = self.DiParams
            except:
                self.DiParams = util.diRead(self.IonStr)
                diparmas = self.DiParams
            cross=np.zeros(len(energy), 'float32')

            for ifac in range(self.DiParams['info']['nfac']):
                # prob. better to do this with masked arrays
                goode=energy > self.DiParams['ev1'][ifac]
                if goode.sum() > 0:
                    dum=np.ones(len(energy))
                    btenergy, btdum=util.scale_bti(energy[goode],dum[goode], self.DiParams['btf'][ifac], self.DiParams['ev1'][ifac] )
                    # these interpolations were made with the scipy routine used here
                    y2=interpolate.splrep(self.DiParams['xsplom'][ifac], self.DiParams['ysplom'][ifac], s=0.)
                    btcross=interpolate.splev(btenergy, y2, der=0.)
                    energy1, cross1=util.descale_bti(btenergy, btcross,self.DiParams['btf'][ifac], self.DiParams['ev1'][ifac] )
                    offset=len(energy)-goode.sum()
                    if verbose:
                        pl.plot(self.DiParams['xsplom'][ifac], self.DiParams['ysplom'][ifac])
                        pl.plot(btenergy, btcross)
                    if offset > 0:
                        seq=[np.zeros(offset, 'float32'), cross1]
                        cross1=np.hstack(seq)
                    cross+=cross1*1.e-14
            self.DiCross={'energy':energy, 'cross':cross}
        #
        #-----------------------------------------------------------
        #
    def diRate(self, temperature=None):
        '''Calculate the direct ionization rate coefficient as a function of temperature (K).'''
        try:
            DiParams = self.DiParams
        except:
            DiParams = util.diRead(self.IonStr)
        #
        if type(temperature) == types.NoneType:
            try:
                temperature=self.Temperature
            except:
                print ' temperature is not defined'
                return
        elif type(temperature) != numpy.ndarray:
            temperature = np.array(temperature,'float64')
        #
        #   gauss laguerre n=12
        #
        ngl=12
        xgl=np.asarray([0.115722117358021,0.611757484515131,1.512610269776419,2.833751337743509
            ,4.599227639418353,6.844525453115181,9.621316842456871,13.006054993306350
            ,17.116855187462260,22.151090379396983,28.487967250983992,37.099121044466926], 'float64')


        wgl=np.asarray([2.647313710554435e-01,3.777592758731382e-01,2.440820113198774e-01,9.044922221168074e-02
            ,2.010238115463406e-02,2.663973541865321e-03,2.032315926629993e-04,8.365055856819753e-06
            ,1.668493876540914e-07,1.342391030515027e-09,3.061601635035012e-12,8.148077467426124e-16], 'float64')
        #
        alpha=5.287e+13
        tev=const.boltzmannEv*temperature
        #
        ntemp=temperature.size
        #
        #
        if ntemp == 1:
            x0=self.Ip/tev  # Ip in eV
            beta=np.sqrt(const.boltzmann*temperature)
            egl=self.Ip+xgl*tev
            self.diCross(energy=egl)
            crossgl=self.DiCross['cross']
            term1=wgl*xgl*crossgl
            term2=wgl*crossgl
            newcross=alpha*beta*np.exp(-x0)*(term1.sum()+x0*term2.sum())
            rate=newcross
        else:
            rate=np.zeros(ntemp, 'float64')
            for itemp in range(ntemp):
                x0=self.Ip/tev[itemp]  # Ip in eV
                beta=np.sqrt(const.boltzmann*temperature[itemp])
                egl=self.Ip+xgl*tev[itemp]
                self.diCross(energy=egl)
                crossgl=self.DiCross['cross']
                term1=wgl*xgl*crossgl
                term2=wgl*crossgl
                newcross=alpha*beta*np.exp(-x0)*(term1.sum()+x0*term2.sum())
                rate[itemp]=newcross
        self.DiRate={'temperature':temperature, 'rate':rate}
            #
            #-----------------------------------------------------------
            #
    def eaDescale(self,temperature=None):
        #
        """Calculates the effective collision strengths (upsilon) for excitation-autoionization as a function of temperature."""
        #
        #  xt=kt/de
        #
        #  need to make sure elvl is >0, except for ground level
        #
        try:
            eaparams=self.EaParams
        except:
            self.EaParams = util.eaRead(self.IonStr)
            eaparams=self.EaParams
        #
        if type(temperature) == types.NoneType:
            try:
                temperature=self.Temperature
            except:
                print ' temperature is not defined'
                return
        ntemp=temperature.size
        nsplups=len(eaparams['de'])
        if ntemp > 1:
            ups=np.zeros((nsplups,ntemp),"Float32")
        else:
            ups=np.zeros(nsplups,"Float32")
        #
        for isplups in range(0,nsplups):
            l1=self.EaParams["lvl1"][isplups]-1
            l2=self.EaParams["lvl2"][isplups]-1
            ttype=self.EaParams["ttype"][isplups]
            cups=self.EaParams["cups"][isplups]
            nspl=self.EaParams["nspl"][isplups]
            de=self.EaParams["de"][isplups]
            dx=1./(float(nspl)-1.)
##                print self.EaParams["EaParams"][l1,l2]
            splups=self.EaParams["splups"][isplups,0:nspl]
            kte=const.boltzmannEv*temperature/(const.ryd2Ev*de)
            #
            if ttype ==1:
                st=1.-np.log(cups)/np.log(kte+cups)
                xs=dx*np.arange(nspl)
                y2=interpolate.splrep(xs,splups,s=0)  #allow smoothing,s=0)
                sups=interpolate.splev(st,y2,der=0)
                ups[isplups]=sups*np.log(kte+np.exp(1.))
            #
            if ttype == 2:
                st=kte/(kte+cups)
                xs=dx*np.arange(nspl)
                y2=interpolate.splrep(xs,splups,s=0)
                sups=interpolate.splev(st,y2,der=0)
                ups[isplups]=sups
            #
            if ttype == 3:
                st=kte/(kte+cups)
                xs=dx*np.arange(nspl)
                y2=interpolate.splrep(xs,splups,s=0)
                sups=interpolate.splev(st,y2,der=0)
                ups[isplups]=sups/(kte+1.)
            #
            if ttype == 4:
                st=1.-np.log(cups)/np.log(kte+cups)
                xs=dx*np.arange(nspl)
                y2=interpolate.splrep(xs,splups,s=0)
                sups=interpolate.splev(st,y2,der=0)
                ups[isplups]=sups*np.log(kte+cups)
            #
            if ttype == 5:
                st=kte/(kte+cups)
                xs=dx*np.arange(nspl)
                y2=interpolate.splrep(xs,splups,s=0)  #allow smoothing,s=0)
                sups=interpolate.splev(st,y2,der=0)
                ups[isplups]=sups/(kte+0.)
            #
            #
            elif ttype > 5:  print ' t_type ne 1,2,3,4,5=',ttype,l1,l2
        #
        #
        ups=np.where(ups > 0.,ups,0.)
        #
        self.EaParams['ups']=ups
        return ups
        #
        # -------------------------------------------------------------------------------------
        #
    def eaCross(self, energy=None, verbose=False):
        '''Provide the excitation-autoionization cross section.

        Energy is given in eV.'''
        # get neaev from diparams file
        #
        try:
            diparams=self.DiParams
        except:
            self.diRead()
        #
        if self.DiParams['info']['neaev'] == 0:
#            print ' no EA rates'
            return
        else:
            if type(energy) == types.NoneType:
                energy=self.Ip*10.**(0.05*arange(31))
            try:
                easplom=self.Easplom
            except:
                self.splomRead()
                easplom=self.Easplom
            #
#        q=np.ma.array(u, 'float32', mask=bu, fill_value=0.)
            #  need to replicate neaev
            ntrans=len(easplom['de'])
            nsplom=easplom['splom'].shape[1]
            x=0.25*arange(nsplom)
            eaev=self.DiParams['eaev']
            if len(eaev) ==1:
                for itrans in range(ntrans):
                    eaev.append(eaev[0])

            for itrans in range(ntrans):
                x0=const.ev*eaparams['de'][itrans]
                good = energy > x0

                thisenergy=np.ma.array(energy, 'float32', mask=good, fill_value=0.)

                earate+=eaev[iups]*8.63e-6*eaparams['ups'][iups]*np.exp(-x0)/(np.sqrt(temperature))
            self.EaRate={'rate':earate, 'temperature':temperature}
            return
        #
        # -------------------------------------------------------------------------------------
        #
    def eaRate(self, temperature=None):
        '''Calculate the excitation-autoionization rate coefficient.'''
        # get neaev from diparams file
        #
        try:
            diparams=self.DiParams
        except:
            self.DiParams = util.diRead(self.IonStr)
        #
        if self.DiParams['info']['neaev'] == 0:
#            print ' no EA rates'
            return
        else:
            if type(temperature) == types.NoneType:
                try:
                    temperature=self.Temperature
                except:
                    bte=0.1*np.arange(10)
                    bte[0]=0.01
                    dum=np.ones(10, 'float32')
                    [temperature, dum]=util.descale_bt(bte, dum, self.EaParams['cups'][0], self.DiParams['de'][0])
                    self.Temperature=temperature
            try:
                eaparams=self.EaParams
            except:
                self.eaParams = util.eaRead(self.IonStr)
                self.eaDescale(temperature=temperature)
                eaparams=self.EaParams
            #
            #  need to replicate neaev
            nups=len(eaparams['de'])
            tev=const.boltzmannEv*temperature
            earate=np.zeros(temperature.size, 'float32')
            eaev=self.DiParams['eaev']
            if len(eaev) ==1:
                for iups in range(nups):
                    eaev.append(eaev[0])

            for iups in range(nups):
                x0=const.ryd2Ev*eaparams['de'][iups]/tev
                earate+=eaev[iups]*8.63e-6*eaparams['ups'][iups]*np.exp(-x0)/(np.sqrt(temperature))
            self.EaRate={'rate':earate, 'temperature':temperature}
            return
        #
        # -------------------------------------------------------------------------------------
        #
    def ionizRate(self, temperature=None):
        '''Provides the total ionization rate.

        Calls diRate and eaRate.'''
        if self.Z < self.Ion:
#            print ' this is a bare nucleus and has no ionization rate'
            self.IonizRate=None
            return
        self.diRate(temperature=temperature)
        self.eaRate(temperature=temperature)
        if self.DiParams['info']['neaev'] == 0:
            ionizrate=self.DiRate['rate']
        else:
            ionizrate=self.DiRate['rate']+self.EaRate['rate']
        self.IonizRate={'rate':ionizrate, 'temperature':self.DiRate['temperature']}
        return
        #
        # -------------------------------------------------------------------------------------
        #
        #
        # -------------------------------------------------------------------------------------
        #
    def rrRate(self, temperature=None):
        '''Provide the radiative recombination rate coefficient as a function of temperature (K).'''
        try:
            rrparams=self.RrParams
        except:
            self.RrParams = util.rrRead(self.IonStr)
            rrparams=self.RrParams
        #
        if type(temperature) == types.NoneType:
            try:
                temperature=self.Temperature
            except:
                print ' temperature is not defined'
                return
#        print ' rr params type = ', rrparams['rrtype']
        #
        if rrparams['rrtype'] == 1:
            a=rrparams['params'][3]
            b=rrparams['params'][4]
            t0=rrparams['params'][5]
            t1=rrparams['params'][6]
            rate=a/(np.sqrt(temperature/t0)*(1.+np.sqrt(temperature/t0))**(1.-b)*(1.+np.sqrt(temperature/t1))**(1.+b))
            self.RrRate={'temperature':temperature, 'rate':rate}
        elif rrparams['rrtype'] == 2:
            a=rrparams['params'][3]
            b=rrparams['params'][4]
            t0=rrparams['params'][5]
            t1=rrparams['params'][6]
            c=rrparams['params'][7]
            t2=rrparams['params'][8]
            b+=c*np.exp(-t2/temperature)
            rate=a/(np.sqrt(temperature/t0)*(1.+np.sqrt(temperature/t0))**(1.-b)*(1.+np.sqrt(temperature/t1))**(1.+b))
            self.RrRate={'temperature':temperature, 'rate':rate}
        elif rrparams['rrtype'] == 3:
            a=rrparams['params'][2]
            b=rrparams['params'][3]
            rate=a/(temperature/1.e+4)**b
            self.RrRate={'temperature':temperature, 'rate':rate}
        #
        # -------------------------------------------------------------------------------------
        #
        #
        # -------------------------------------------------------------------------------------
        #
    def drRate(self, temperature=None):
        '''Provide the dielectronic recombination rate coefficient as a function of temperature (K).'''
        try:
            drparams=self.DrParams
        except:
            self.DrParams = util.drRead(self.IonStr)
            drparams=self.DrParams
        #
        if type(temperature) == types.NoneType:
            try:
                temperature=self.Temperature
            except:
                print ' temperature is not defined'
                return
#        print ' dr params type = ', drparams['drtype']
        #
        if type(drparams) == types.NoneType:
            self.DrRate=None
        else:
            if drparams['drtype'] == 1:
                # badnell type
                drenergy=drparams['eparams']
                drcoef=drparams['cparams']
                gcoef = drenergy > 0.
                ncoef=gcoef.sum()
#                print ' ncoef = ', gcoef.sum()
                rate=np.zeros(temperature.size, 'float32')
                for icoef in range(ncoef):
                    rate+=drcoef[icoef]*np.exp(-drenergy[icoef]/temperature)
                rate=rate/temperature**1.5
                self.DrRate={'temperature':temperature, 'rate':rate}
            elif drparams['drtype'] == 2:
                # shull type
                params = drparams['params']
                adi = params[0]
                bdi = params[1]
                t0 = params[2]
                t1 = params[3]
                rate=adi*np.exp(-t0/temperature)*(1.+bdi*np.exp(-t1/temperature))/temperature**1.5
                self.DrRate={'temperature':temperature, 'rate':rate}
        #
        # -------------------------------------------------------------------------------------
        #
    def reclvlDescale(self, temperature=None):
        '''Interpolate and extrapolate reclvl rates.

        Used in level population calculations.'''
        if type(temperature) == types.NoneType:
            try:
                temperature=self.Temperature
            except:
                print ' temperature is not defined'
                self.ReclvlRate=None
        try:
            reclvl = self.Reclvl
            if reclvl == types.NoneType:
                self.ReclvlRate = None
                return
        except:
#           print ' reading reclvl file'
            reclvl = util.cireclvlRead(self.IonStr, 'reclvl')
            if type(reclvl) == types.NoneType:
                self.ReclvlRate = None
                return
        #
        #  the rates and temperatures in reclvl are not all the same
        #
        ntemp = temperature.size
        if ntemp == 1:
            recRate = np.zeros(( len(reclvl['lvl1'])), 'float64')
            if temperature < reclvl['temperature'].min():
                self.ReclvlRate = None
            elif temperature > reclvl['temperature'].max():
                for itrans in range(len(reclvl['lvl1'])):
#                   lvl2 = self.Reclvl['lvl2'][itrans]
                    nrecTemp = reclvl['ntemp'][itrans]
                    recRate[itrans] = self.Reclvl['rate'][itrans,nrecTemp-1]*(reclvl['temperature'][itrans, nrecTemp-1]/temperature)
            else:
                for itrans in range(len(self.Reclvl['lvl1'])):
                    lvl2 = self.Reclvl['lvl2'][itrans]
                    nrecTemp = self.Reclvl['ntemp'][itrans]
#                   print ' nrecTemp = ', nrecTemp
                    y2 = interpolate.splrep(np.log(self.Reclvl['temperature'][itrans, :nrecTemp]), np.log(self.Reclvl['rate'][itrans, :nrecTemp]))
                    rec = np.exp(interpolate.splev(np.log(temperature),y2))
                    recRate[itrans] = rec.squeeze()
        else:
            # ntemp > 1
            recRate = np.zeros(( len(reclvl['lvl1']), temperature.size), 'float64')
            #
            for itrans in range(len(reclvl['lvl1'])):
                lvl2 = self.Reclvl['lvl2'][itrans]
                nrecTemp = self.Reclvl['ntemp'][itrans]
                y2 = interpolate.splrep(np.log(self.Reclvl['temperature'][itrans, :nrecTemp]), np.log(self.Reclvl['rate'][itrans, :nrecTemp]))
                goodLow = temperature < self.Reclvl['temperature'][itrans].min()
                if goodLow.sum() >0:
#                   print ' number of low temperatures  = ', goodLow.sum()
                    lowT = temperature[goodLow]
                good1 = temperature >= self.Reclvl['temperature'][itrans].min()
                good2 = temperature <= self.Reclvl['temperature'][itrans].max()
                realgood = np.logical_and(good1,good2)
                if realgood.sum() > 0:
#                   print ' number of mid temperatures  = ', realgood.sum()
                    midT = temperature[realgood]
                goodHigh = temperature > self.Reclvl['temperature'][itrans].max()
                if goodHigh.sum() > 0:
#                   print ' number of high temperatures  = ', goodHigh.sum()
                    highT = temperature[goodHigh]
                lvl2 = self.Reclvl['lvl2'][itrans]
                nrecTemp = self.Reclvl['ntemp'][itrans]
                newRec = np.zeros(ntemp, 'float64')
                index = 0
                if goodLow.sum() == 1:
                    lowRec = np.exp(interpolate.splev(np.log(lowT),y2))
                    newRec[index] = lowRec
                    index += 1
                elif goodLow.sum() > 1:
                    lowRec = np.exp(interpolate.splev(np.log(lowT),y2))
                    for idx in range(goodLow.sum()):
                        newRec[index] = lowRec[idx]
                        index += 1
                if realgood.sum() == 1:
                    midRec = np.exp(interpolate.splev(np.log(midT),y2))
                    newRec[index] = midRec
                    index += 1
                elif realgood.sum() > 1:
                    midRec = np.exp(interpolate.splev(np.log(midT),y2))
                    for idx in range(realgood.sum()):
                        newRec[index] = midRec[idx]
                        index += 1
                if goodHigh.sum() == 1:
                    highRec = np.exp(interpolate.splev(np.log(highT),y2))
                    newRec[index] = highRec
                    index += 1
                elif goodHigh.sum() > 1:
                    highRec = np.exp(interpolate.splev(np.log(highT),y2))
                    for idx in range(goodHigh.sum()):
#                       print ' index, idx = ', index,  idx
                        newRec[index] = highRec[idx]
                        index += 1
        self.ReclvlRate = {'rate':recRate, 'lvl2':reclvl['lvl2']}
        #
        # -------------------------------------------------------------------------------------
        #
    def recombRate(self, temperature=None):
        '''Provides the total recombination rate coefficient.

        Calls diRate and eaRate'''
        if self.Ion == 1:
#            print ' this is a neutral and has no recombination rate'
            self.RecombRate=None
            return
        #
        if type(temperature) == types.NoneType:
            try:
                temperature=self.Temperature
            except:
                print ' temperature is not defined'
                self.RecombRate=None
        self.rrRate(temperature=temperature)
        self.drRate(temperature=temperature)
        if type(self.DrRate) == types.NoneType:
            rate=self.RrRate['rate']
        else:
            rate=self.RrRate['rate']+self.DrRate['rate']
        self.RecombRate={'rate':rate, 'temperature':temperature}
        return
        #
        # -------------------------------------------------------------------------------------
        #
    def p2eRatio(self):
        '''Calculates the proton density to electron density ratio.

        Uses the abundance and ionization equilibrium.'''
        try:
            ab=self.Abundance
        except:
            abundName = self.Defaults['abundfile']
            util.abundanceRead(abundancename = abundName)
        try:
            ioneq=self.Ioneq
        except:
            ioneqname = self.Defaults['ioneqfile']
            self.IoneqAll = util.ioneqRead(ioneqname = ioneqname)
        #
        try:
            temperature=self.Temperature
        except:
            temperature = self.IoneqAll['ioneqTemperature']
#                temperature=self.IoneqTemperature
        else:  temperature=np.asarray(temperature,'float32')
        #
        nTemp=temperature.size
        nEl=len(self.AbundanceAll['abundance'])
        #
        eDensity=np.zeros(nTemp,'float32')
        pDensity=np.zeros(nTemp,'float32')
        ionDensity=0.
#        zDensity=np.zeros(nTemp, 'float32')
        #
        #
        #  only hydrogen contributes to the proton density
        anEl = 0
        ion = 1
        good = self.IoneqAll['ioneqAll'][anEl,ion]  > 0.
        y2 = interpolate.splrep(np.log(self.IoneqAll['ioneqTemperature'][good]),np.log(self.IoneqAll['ioneqAll'][anEl,ion,good]),s=0)
        bad1 = np.log(temperature) < np.log(self.IoneqAll['ioneqTemperature'][good].min())
        bad2 = np.log(temperature) > np.log(self.IoneqAll['ioneqTemperature'][good].max())
        bad=np.logical_or(bad1,bad2)
        goodt=np.logical_not(bad)
        thisIoneq=np.where(goodt,10.**interpolate.splev(np.log(temperature),y2,der=0),0.)
        pDensity+=self.AbundanceAll['abundance'][anEl]*thisIoneq
#        ionDensity=self.Abundance[anEl]
        #
        # all the rest do contribute to the electron and ion densities
        El=[iEl for iEl in range(50) if  self.AbundanceAll['abundance'][iEl] > 0.]
        for anEl in El:
            ionDensity+=self.AbundanceAll['abundance'][anEl]
            for ion in range(1,anEl+2):
                good = self.IoneqAll['ioneqAll'][anEl,ion]  > 0.
                y2 = interpolate.splrep(np.log(self.IoneqAll['ioneqTemperature'][good]),np.log(self.IoneqAll['ioneqAll'][anEl,ion,good]),s=0)
                bad1 = np.log(temperature) < np.log(self.IoneqAll['ioneqTemperature'][good].min())
                bad2 = np.log(temperature) > np.log(self.IoneqAll['ioneqTemperature'][good].max())
                bad = np.logical_or(bad1,bad2)
                goodt = np.logical_not(bad)
#                good1 = temperature >= self.IoneqAll['ioneqTemperature'][good].min()
#                good2 = temperature <= self.IoneqAll['ioneqTemperature'][good].min()
#                goodt = np.logical_and(good1, good2)
#                thisIoneq=np.where(goodt,10.**interpolate.splev(np.log(temperature),y2,der=0),0.)
                # for temperatures outside the range of the ioneq, the default is set to 1, since
                #  this is the denominator
                thisIoneq=np.where(goodt,10.**interpolate.splev(np.log(temperature),y2,der=0),1.)
                eDensity+=float(ion)*self.AbundanceAll['abundance'][anEl]*thisIoneq
        self.ProtonDensityRatio=pDensity/eDensity
        self.EDensity=eDensity
        self.IonDensity=ionDensity
        self.IonDensityRatio=ionDensity/eDensity
#        return # {'PDensityRatio':self.PDensityRatio,'PDensity':pDensity,'EDensity':eDensity}
        #
        # -------------------------------------------------------------------------------------
        #
    def upsilonDescale(self,temperature=None,prot=0, ci=0):
        #
        """Provides the temperatures and effective collision strengths (upsilons)."""
        #
        #  xt=kt/de
        #
        #
        if prot:
            try:
                nsplups=len(self.Psplups["lvl1"])
            except:
                self.Psplups=util.splupsRead(self.IonStr,prot=1)
                if type(self.Psplups) == types.NoneType:
                    self.PUpsilon = None
                    return
                else:
                    nsplups = len(self.CiSplups["lvl1"])
        elif ci:
            try:
                nsplups = len(self.CiSplups["lvl1"])
            except:
                self.CiSplups = util.splupsRead(self.IonStr,ci=1)
                if type(self.CiSplups) == types.NoneType:
                    self.CiUpsilon = None
                    return
                else:
                    nsplups = len(self.CiSplups["lvl1"])
        else:
            try:
                nsplups=len(self.Splups["lvl1"])
            except:
                self.Splups = util.splupsRead(self.IonStr)
                if type(self.Splups) == types.NoneType:
                    self.Upsilon = None
                    return
                else:
                    nsplups = len(self.Splups["lvl1"])
        #
        #
        if temperature == None:
            try:
                temperature=self.Temperature
            except:
                print ' Temperature unknown'
                return
        else:
            self.Temperature=temperature
        #
        try:
            nlvls=len(self.Elvlc["lvl"])
        except:
            self.elvlcRead()
            nlvls=len(self.Elvlc["lvl"])
        #
        #  need to make sure elvl is >0, except for ground level
        eryd=np.asarray(self.Elvlc["eryd"])
        erydth=np.asarray(self.Elvlc["erydth"])
        elvlc=np.where(eryd > 0.,eryd,erydth)
##        de=self.Elvlc["de"]
        temp=np.asarray(temperature)
        ntemp=temp.size
        if ntemp > 1:
            ups=np.zeros((nsplups,ntemp),"Float64")
        else:
            ups=np.zeros(nsplups,"Float64")
        #
        for isplups in range(0,nsplups):
            if prot:
                # for proton rates
                l1=self.Psplups["lvl1"][isplups]-1
                l2=self.Psplups["lvl2"][isplups]-1
                ttype=self.Psplups["ttype"][isplups]
                cups=self.Psplups["cups"][isplups]
                nspl=self.Psplups["nspl"][isplups]
                ttype=6
                dx=1./(float(nspl)-1.)
                splups=self.Psplups["splups"][isplups,0:nspl]
                de=elvlc[l2]-elvlc[l1]
#                de=self.Psplups['de'][isplups]  # these are generally 0.
                kte=temp/(de*1.57888e+5)
            elif ci:
                # for proton rates
                l1 = self.CiSplups["lvl1"][isplups]-1
                l2 = self.CiSplups["lvl2"][isplups]-1
                ttype = self.CiSplups["ttype"][isplups]
                cups = self.CiSplups["cups"][isplups]
                nspl = self.CiSplups["nspl"][isplups]
                ttype = self.CiSplups["ttype"][isplups]
                dx = 1./(float(nspl)-1.)
                splups = self.CiSplups["splups"][isplups,0:nspl]
                de=self.CiSplups['de'][isplups]
                kte = temp/(de*1.57888e+5)
            else:
                l1=self.Splups["lvl1"][isplups]-1
                l2=self.Splups["lvl2"][isplups]-1
                ttype=self.Splups["ttype"][isplups]
                cups=self.Splups["cups"][isplups]
                nspl=self.Splups["nspl"][isplups]
                dx=1./(float(nspl)-1.)
##                print self.Splups["splups"][l1,l2]
                splups=self.Splups["splups"][isplups,0:nspl]
#                de=elvlc[l2]-elvlc[l1]
                de=self.Splups['de'][isplups]
                kte=temp/(de*1.57888e+5)
            #
            der=0
            if ttype ==1:
                st=1.-np.log(cups)/np.log(kte+cups)
                xs=dx*np.arange(nspl)
                y2=interpolate.splrep(xs,splups,s=0)  #allow smoothing,s=0)
                sups=interpolate.splev(st,y2,der=der)
#                sups=interpolate.spline(xs, splups, st)
                ups[isplups]=sups*np.log(kte+np.exp(1.))
            #
            if ttype == 2:
                st=kte/(kte+cups)
                xs=dx*np.arange(nspl)
                y2=interpolate.splrep(xs,splups,s=0)  #allow smoothing,s=0)
                sups=interpolate.splev(st,y2,der=der)
                ups[isplups]=sups
            #
            if ttype == 3:
                st=kte/(kte+cups)
                xs=dx*np.arange(nspl)
                y2=interpolate.splrep(xs,splups,s=0)  #allow smoothing,s=0)
                sups=interpolate.splev(st,y2,der=der)
                ups[isplups]=sups/(kte+1.)
            #
            if ttype == 4:
                st=1.-np.log(cups)/np.log(kte+cups)
                xs=dx*np.arange(nspl)
                y2=interpolate.splrep(xs,splups,s=0)  #allow smoothing,s=0)
                sups=interpolate.splev(st,y2,der=der)
                ups[isplups]=sups*np.log(kte+cups)
            #
            if ttype == 5:
                st=kte/(kte+cups)
                xs=dx*np.arange(nspl)
                y2=interpolate.splrep(xs,splups,s=0)  #allow smoothing,s=0)
                sups=interpolate.splev(st,y2,der=der)
                ups[isplups]=sups/(kte+0.)
            #
            #  descale proton values
            if ttype == 6:
                st=kte/(kte+cups)
                xs=dx*np.arange(nspl)
                y2=interpolate.splrep(xs,splups,s=0)  #allow smoothing,s=0)
                sups=interpolate.splev(st,y2,der=der)
                ups[isplups] = sups
#                ups[isplups]=10.**sups
            #
            elif ttype > 6:  print ' t_type ne 1,2,3,4,5=',ttype,l1,l2
        #
        #
        ups=np.where(ups > 0.,ups,0.)
        #
        if prot == 1:
            self.PUpsilon = ups
        elif ci == 1:
            self.CiUpsilon = ups
        else:
            self.Upsilon = ups
        #
        # -------------------------------------------------------------------------
        #
    def spectrum(self,wavelength, filter=(chfilters.gaussianR,1000.)):
        '''Calculates the line emission spectrum for the specified ion.

        Convolves the results of intensity to make them look like an observed spectrum
        the default filter is the gaussianR filter with a resolving power of 1000.  Other choices
        include chianti.filters.box and chianti.filters.gaussian.  When using the box filter,
        the width should equal the wavelength interval to keep the units of the continuum and line
        spectrum the same.
        Note:  scipy.ndimage.filters also includes a range of filters.'''
        aspectrum = np.zeros_like(wavelength)
        nTemp = self.Temperature.size
        nDens = self.Density.size
        useFilter = filter[0]
        useFactor= filter[1]
        try:
            intensity = self.Intensity
        except:
            self.intensity()
            intensity = self.Intensity
        #
        if (nTemp == 1) and (nDens == 1):
            aspectrum = np.zeros_like(wavelength)
            for iwvl, wvlCalc in enumerate(intensity['wvl']):
                aspectrum += useFilter(wavelength, wvlCalc, factor=useFactor)*intensity['intensity'][iwvl]
        else:
            aspectrum = np.zeros((nTemp, wavelength.size), 'float64')
            for itemp in xrange(nTemp):
                for iwvl, wvlCalc in enumerate(self.Intensity['wvl']):
                    aspectrum[itemp] += useFilter(wavelength, wvlCalc, factor=useFactor)*self.Intensity['intensity'][itemp, iwvl]

        self.Spectrum = {'intensity':aspectrum,  'wvl':wavelength, 'filter':useFilter.__name__, 'filterWidth':useFactor}
        #
        # -------------------------------------------------------------------------------------
        #
    def populate(self,temperature=None,density=None,pDensity=None, popCorrect=1, radTemperature=0,rPhot=1.):
        """Calculate level populations for specified ion."""
        #
        nlvls=self.Nlvls
        nwgfa=self.Nwgfa
        nsplups=self.Nsplups
        npsplups=self.Npsplups
        #
        if type(temperature) == types.NoneType:
            try:
                temperature=self.Temperature
            except:
                print ' no temperature values have been set'
                return
        else:
            self.Temperature = np.asarray(temperature)
            temperature = self.Temperature
        #
        if type(density) == types.NoneType:
            try:
                density=self.Density
            except:
                print ' no density values have been set'
                return
        else:
            self.Density = np.asarray(density)
            density = self.Density
        #
        if type(pDensity) == types.NoneType:
            self.p2eRatio()
            protonDensity = self.ProtonDensityRatio*self.Density
        #
        elif type(pDensity) == types.StringType:
            # the only string it can be is "default"
            self.p2eRatio()
            protonDensity = self.ProtonDensityRatio*self.Density
        else:
            self.PDensity = np.asarray(pDensity)
            if self.PDensity.size == 1:
                self.PDensity = self.PDensity.repeat(self.Density.size)
            protonDensity = self.PDensity
        #
        if radTemperature:
            self.RadTemperature = radTemperature
            self.RPhot = rPhot
        #
        if popCorrect and (not self.Dielectronic):
            self.upsilonDescale(ci=1)
            if type(self.CiUpsilon) != types.NoneType:
                ci = 1
                cisplups = self.CiSplups
                ciupsilon = self.CiUpsilon
                self.recombRate()
                lowers = util.zion2name(self.Z, self.Ion-1)
                # get the lower ionization stage
                lower = ion(lowers, temperature=self.Temperature)
                lower.ionizRate()
                # need to get multiplicity of lower ionization stage
                lowMult = lower.Elvlc['mult']
            else:
                ci = 0
            try:
                if self.Nreclvl > 0:
                    self.reclvlDescale()
                    rec = 1
                    self.ionizRate()
                    #  get the higher ionization stage
                    highers = convertname(self.Z, self.Ion+1)
                    higher = ion(highers, temperature=self.Temperature)
                    higher.recombRate()
                else:
                    rec = 0
            except:
                self.Reclvl = util.cireclvlRead(self.IonStr,'reclvl' )
                reclvl = self.Reclvl
                self.reclvlDescale()
                if type(self.Reclvl) != types.NoneType:
                    rec = 1
                    self.ionizRate()
                    #  get the higher ionization stage
                    highers = util.zion2name(self.Z, self.Ion+1)
#                   print ' highers = ', highers
                    higher = ion(highers, temperature=self.Temperature)
                    higher.recombRate()
                else:
                    rec = 0
        else:
            ci = 0
            rec = 0
        #
#       print ' dir(lower) = ', dir(lower)
#       print ' dir(higher) = ', dir(higher)
        #
        rad=np.zeros((nlvls+ci+rec,nlvls+ci+rec),"float64")  #  the populating matrix for radiative transitions
        #
        #
        for iwgfa in range(0,nwgfa):
            l1 = self.Wgfa["lvl1"][iwgfa]-1
            l2 = self.Wgfa["lvl2"][iwgfa]-1
            rad[l1+ci,l2+ci] += self.Wgfa["avalue"][iwgfa]
            rad[l2+ci,l2+ci] -= self.Wgfa["avalue"][iwgfa]
            if self.RadTemperature:
                if not self.RPhot:
                    dilute = 0.5
                else:
                    dilute = util.dilution(self.RPhot)
                de = const.invCm2Erg*(self.Elvlc['ecm'][l2] - self.Elvlc['ecm'][l1])
                dekt = de/(const.boltzmann*self.RadTemperature)
                factor = dilute*(self.Elvlc['mult'][l2]/self.Elvlc['mult'][l1])/(np.exp(dekt)-1.)
                rad[l2+ci,l1+ci] += self.Wgfa["avalue"][iwgfa]*factor
                rad[l1+ci,l1+ci] -= self.Wgfa["avalue"][iwgfa]*factor

        #
        self.rad=rad
        #
        self.upsilonDescale(temperature=temperature)
        ups=self.Upsilon
        #
        if npsplups >0:
            self.upsilonDescale(temperature=temperature,prot=1)
            pups=self.PUpsilon
        #
        temp=temperature
        ntemp=temp.size
        #
        cc=const.collision*self.Density
        ndens=cc.size
        if npsplups > 0:
            cp=const.collision*protonDensity
        if ntemp > 1 and ndens >1 and ntemp != ndens:
            print ' unless temperature or density are single values'
            print ' the number of temperatures values must match the '
            print ' the number of density values'
            return
        #
        #
        #  note:  the observed energy levels are used to calculate level populations
        #         this is different from the IDL version
        #
        #  the upsilons are derived using the theoretical energies as in the IDL version
        #
        # get corrections for recombination and excitation
        #
        #
        #  first, for ntemp=ndens=1
        if ndens==1 and ntemp==1:
#           print ' ndens, ntemp = 1'
#            pop=np.zeros((nlvls),"float64")
            popmat=np.copy(rad)
            for isplups in range(0,nsplups):
                l1=self.Splups["lvl1"][isplups]-1
                l2=self.Splups["lvl2"][isplups]-1
                if self.Dielectronic:
                    de=np.abs((self.Elvlc["eryd"][l2]-self.Ip/const.ryd2Ev)-self.Elvlc["eryd"][l1])
                else:
                    de=np.abs(self.Elvlc["eryd"][l2]-self.Elvlc["eryd"][l1])
                ekt=(de*const.ryd2erg)/(const.boltzmann*temp)
                fmult1=float(self.Elvlc["mult"][l1])
                fmult2=float(self.Elvlc["mult"][l2])
                popmat[l1+ci,l2+ci]+=cc*ups[isplups]/(fmult2*np.sqrt(temp))
                popmat[l2+ci,l1+ci]+=cc*ups[isplups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
                popmat[l1+ci,l1+ci]-=cc*ups[isplups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
                popmat[l2+ci,l2+ci]-=cc*ups[isplups]/(fmult2*np.sqrt(temp))
            for isplups in range(0,npsplups):
                l1=self.Psplups["lvl1"][isplups]-1
                l2=self.Psplups["lvl2"][isplups]-1
                # for proton excitation, the levels are all below the ionization potential
                de=np.abs(self.Elvlc["eryd"][l2]-self.Elvlc["eryd"][l1])
                ekt=(de*1.57888e+5)/temp
                fmult1=float(self.Elvlc["mult"][l1])
                fmult2=float(self.Elvlc["mult"][l2])
                popmat[l1+ci,l2+ci]+=cp*pups[isplups]/(fmult2*np.sqrt(temp))
                popmat[l2+ci,l1+ci]+=cp*pups[isplups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
                popmat[l1+ci,l1+ci]-=cp*pups[isplups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
                popmat[l2+ci,l2+ci]-=cp*pups[isplups]/(fmult2*np.sqrt(temp))
            # now include ionization rate from
            if ci:
#               print ' ci = ', ci
                popmat[1, 0] += self.Density*lower.IonizRate['rate']
                popmat[0, 0] -= self.Density*lower.IonizRate['rate']
                popmat[0, 1] += self.Density*self.RecombRate['rate']
                popmat[1, 1] -= self.Density*self.RecombRate['rate']
                #
                # the ciRate can be computed for all temperatures
                #
                for itrans in range(len(cisplups['lvl1'])):
                    lvl1 = cisplups['lvl1'][itrans]
                    lvl2 = cisplups['lvl2'][itrans]
                    de = cisplups['de'][itrans]
                    ekt = (de*1.57888e+5)/temperature
                    mult = lowMult[lvl1-1]
                    cirate = const.collision*self.CiUpsilon[itrans]*np.exp(-ekt)/(np.sqrt(temperature)*mult)
                    # this is kind of double booking the ionization rate components
                    popmat[lvl2, lvl1-1] += self.Density*cirate
                    popmat[lvl1-1, lvl1-1] -= self.Density*cirate
            if rec:
#               print ' rec = ', rec
                popmat[-1,  ci] += self.Density*self.IonizRate['rate']
                popmat[ci, ci] -= self.Density*self.IonizRate['rate']
                popmat[ci, -1] += self.Density*higher.RecombRate['rate']
                popmat[-1, -1] -= self.Density*higher.RecombRate['rate']
                #
                for itrans in range(len(reclvl['lvl1'])):
                    lvl1 = reclvl['lvl1'][itrans]
                    lvl2 = reclvl['lvl2'][itrans]
                    popmat[lvl2+ci+1, -1] += self.Density*self.ReclvlRate['rate'][itrans]
                    popmat[-1, -1] -= self.Density*self.ReclvlRate['rate'][itrans]
            # normalize to unity
            norm=np.ones(nlvls+ci+rec,'float64')
            if ci:
                norm[0] = 0.
            if rec:
                norm[-1] = 0.
            popmat[nlvls+ci+rec-1]=norm
            b=np.zeros(nlvls+ci+rec,'float64')
            b[nlvls+ci+rec-1]=1.
            if rec:
                pop = np.linalg.solve(popmat,b)[ci:-rec]
            else:
                pop = np.linalg.solve(popmat,b)[ci:]
        #   next, in case of a single density value
        elif ndens == 1:
            pop=np.zeros((ntemp,nlvls),"float64")
            for itemp in range(0,ntemp):
                popmat=np.copy(rad)
                for isplups in range(0,nsplups):
                    l1=self.Splups["lvl1"][isplups]-1
                    l2=self.Splups["lvl2"][isplups]-1
                    if self.Dielectronic:
                        de=np.abs((self.Elvlc["eryd"][l2]-self.Ip/const.ryd2Ev)-self.Elvlc["eryd"][l1])
                    else:
                        de=np.abs(self.Elvlc["eryd"][l2]-self.Elvlc["eryd"][l1])
                    ekt=(de*1.57888e+5)/temp[itemp]
                    fmult1=float(self.Elvlc["mult"][l1])
                    fmult2=float(self.Elvlc["mult"][l2])
                    popmat[l1+ci,l2+ci]+=cc*ups[isplups, itemp]/(fmult2*np.sqrt(temp[itemp]))
                    popmat[l2+ci,l1+ci]+=cc*ups[isplups, itemp]*np.exp(-ekt)/(fmult1*np.sqrt(temp[itemp]))
                    popmat[l1+ci,l1+ci]-=cc*ups[isplups, itemp]*np.exp(-ekt)/(fmult1*np.sqrt(temp[itemp]))
                    popmat[l2+ci,l2+ci]-=cc*ups[isplups, itemp]/(fmult2*np.sqrt(temp[itemp]))
                for isplups in range(0,npsplups):
                    l1=self.Psplups["lvl1"][isplups]-1
                    l2=self.Psplups["lvl2"][isplups]-1
                    # for proton excitation, the levels are all below the ionization potential
                    de=np.abs(self.Elvlc["eryd"][l2]-self.Elvlc["eryd"][l1])
                    ekt=(de*1.57888e+5)/temp[itemp]
                    fmult1=float(self.Elvlc["mult"][l1])
                    fmult2=float(self.Elvlc["mult"][l2])
                    popmat[l1+ci,l2+ci]+=cp[itemp]*pups[isplups, itemp]/(fmult2*np.sqrt(temp[itemp]))
                    popmat[l2+ci,l1+ci]+=cp[itemp]*pups[isplups, itemp]*np.exp(-ekt)/(fmult1*np.sqrt(temp[itemp]))
                    popmat[l1+ci,l1+ci]-=cp[itemp]*pups[isplups, itemp]*np.exp(-ekt)/(fmult1*np.sqrt(temp[itemp]))
                    popmat[l2+ci,l2+ci]-=cp[itemp]*pups[isplups, itemp]/(fmult2*np.sqrt(temp[itemp]))
                # now include ionization rate from
                if ci:
#                   print ' ci = ', ci
                    popmat[1, 0] += self.Density*lower.IonizRate['rate'][itemp]
                    popmat[0, 0] -= self.Density*lower.IonizRate['rate'][itemp]
                    popmat[0, 1] += self.Density*self.RecombRate['rate'][itemp]
                    popmat[1, 1] -= self.Density*self.RecombRate['rate'][itemp]
                    #
                    # the ciRate can be computed for all temperatures
                    #
                    for itrans in range(len(cisplups['lvl1'])):
                        lvl1 = cisplups['lvl1'][itrans]
                        lvl2 = cisplups['lvl2'][itrans]
                        de = cisplups['de'][itrans]
                        ekt = (de*1.57888e+5)/temperature
                        mult = lowMult[lvl1-1]
                        cirate = const.collision*self.CiUpsilon[itrans]*np.exp(-ekt)/(np.sqrt(temp[itemp])*mult)
                        # this is kind of double booking the ionization rate components
                        popmat[lvl2, lvl1-1] += self.Density*cirate[itemp]
                        popmat[lvl1-1, lvl1-1] -= self.Density*cirate[itemp]
                if rec:
#                   print ' rec = ', rec
                    popmat[-1,  ci] += self.Density*self.IonizRate['rate'][itemp]
                    popmat[ci, ci] -= self.Density*self.IonizRate['rate'][itemp]
                    popmat[ci, -1] += self.Density*higher.RecombRate['rate'][itemp]
                    popmat[-1, -1] -= self.Density*higher.RecombRate['rate'][itemp]
                    #
                    for itrans in range(len(reclvl['lvl1'])):
                        lvl1 = reclvl['lvl1'][itrans]
                        lvl2 = reclvl['lvl2'][itrans]
                        popmat[lvl2+ci+1, -1] += self.Density*self.ReclvlRate['rate'][itrans, itemp]
                        popmat[-1, -1] -= self.Density*self.ReclvlRate['rate'][itrans, itemp]
                # normalize to unity
                norm=np.ones(nlvls+ci+rec,'float64')
                if ci:
                    norm[0] = 0.
                if rec:
                    norm[-1] = 0.
                popmat[nlvls+ci+rec-1]=norm
                b=np.zeros(nlvls+ci+rec,'float64')
                b[nlvls+ci+rec-1]=1.
                thispop=np.linalg.solve(popmat,b)
                if rec:
                    pop[itemp] = thispop[ci:-rec]
                else:
                    pop[itemp] = thispop[ci:]
            #
        elif ntemp == 1:
            pop=np.zeros((ndens,nlvls),"float64")
            for idens in range(0,ndens):
                popmat=np.copy(rad)
                for isplups in range(0,nsplups):
                    l1=self.Splups["lvl1"][isplups]-1
                    l2=self.Splups["lvl2"][isplups]-1
                    if self.Dielectronic:
                        de=np.abs((self.Elvlc["eryd"][l2]-self.Ip/const.ryd2Ev)-self.Elvlc["eryd"][l1])
                    else:
                        de=np.abs(self.Elvlc["eryd"][l2]-self.Elvlc["eryd"][l1])
                    ekt=(de*1.57888e+5)/temp
                    fmult1=float(self.Elvlc["mult"][l1])
                    fmult2=float(self.Elvlc["mult"][l2])
                    popmat[l1+ci,l2+ci]+=cc[idens]*ups[isplups]/(fmult2*np.sqrt(temp))
                    popmat[l2+ci,l1+ci]+=cc[idens]*ups[isplups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
                    popmat[l1+ci,l1+ci]-=cc[idens]*ups[isplups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
                    popmat[l2+ci,l2+ci]-=cc[idens]*ups[isplups]/(fmult2*np.sqrt(temp))
                for isplups in range(0,npsplups):
                    l1=self.Psplups["lvl1"][isplups]-1
                    l2=self.Psplups["lvl2"][isplups]-1
                    # for proton excitation, the levels are all below the ionization potential
                    de=np.abs(self.Elvlc["eryd"][l2]-self.Elvlc["eryd"][l1])
                    ekt=(de*1.57888e+5)/temp
                    fmult1=float(self.Elvlc["mult"][l1])
                    fmult2=float(self.Elvlc["mult"][l2])
                    popmat[l1+ci,l2+ci]+=cp[idens]*pups[isplups]/(fmult2*np.sqrt(temp))
                    popmat[l2+ci,l1+ci]+=cp[idens]*pups[isplups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
                    popmat[l1+ci,l1+ci]-=cp[idens]*pups[isplups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
                    popmat[l2+ci,l2+ci]-=cp[idens]*pups[isplups]/(fmult2*np.sqrt(temp))
                # now include ionization rate from
                if ci:
#                   print ' ci = ', ci
                    popmat[1, 0] += self.Density[idens]*lower.IonizRate['rate']
                    popmat[0, 0] -= self.Density[idens]*lower.IonizRate['rate']
                    popmat[0, 1] += self.Density[idens]*self.RecombRate['rate']
                    popmat[1, 1] -= self.Density[idens]*self.RecombRate['rate']
                    #
                    #
                    for itrans in range(len(cisplups['lvl1'])):
                        lvl1 = cisplups['lvl1'][itrans]
                        lvl2 = cisplups['lvl2'][itrans]
                        de = cisplups['de'][itrans]
                        ekt = (de*1.57888e+5)/temperature
                        mult = lowMult[lvl1-1]
                        cirate = const.collision*self.CiUpsilon[itrans]*np.exp(-ekt)/(np.sqrt(temp)*mult)
                        # this is kind of double booking the ionization rate components
                        popmat[lvl2, lvl1-1] += self.Density[idens]*cirate
                        popmat[lvl1-1, lvl1-1] -= self.Density[idens]*cirate
                if rec:
#                   print ' rec = ', rec
                    popmat[-1,  ci] += self.Density[idens]*self.IonizRate['rate']
                    popmat[ci, ci] -= self.Density[idens]*self.IonizRate['rate']
                    popmat[ci, -1] += self.Density[idens]*higher.RecombRate['rate']
                    popmat[-1, -1] -= self.Density[idens]*higher.RecombRate['rate']
                    #
                    for itrans in range(len(reclvl['lvl1'])):
                        lvl1 = reclvl['lvl1'][itrans]
                        lvl2 = reclvl['lvl2'][itrans]
                        popmat[lvl2+ci+1, -1] += self.Density[idens]*self.ReclvlRate['rate'][itrans]
                        popmat[-1, -1] -= self.Density[idens]*self.ReclvlRate['rate'][itrans]
                # normalize to unity
                norm=np.ones(nlvls+ci+rec,'float64')
                if ci:
                    norm[0] = 0.
                if rec:
                    norm[-1] = 0.
                popmat[nlvls+ci+rec-1]=norm
                b=np.zeros(nlvls+ci+rec,'float64')
                b[nlvls+ci+rec-1]=1.
                thispop=np.linalg.solve(popmat,b)
                if rec:
                    pop[idens] = thispop[ci:-rec]
                else:
                    pop[idens] = thispop[ci:]
                #
        elif ntemp>1  and ntemp==ndens:
            pop=np.zeros((ntemp,nlvls),"float64")
            for itemp in range(0,ntemp):
                temp=self.Temperature[itemp]
                popmat=np.copy(rad)
                for isplups in range(0,nsplups):
                    l1=self.Splups["lvl1"][isplups]-1
                    l2=self.Splups["lvl2"][isplups]-1
                    if self.Dielectronic:
                        de=np.abs((self.Elvlc["eryd"][l2]-self.Ip/const.ryd2Ev)-self.Elvlc["eryd"][l1])
                    else:
                        de=np.abs(self.Elvlc["eryd"][l2]-self.Elvlc["eryd"][l1])
                    ekt=(de*1.57888e+5)/temp
                    fmult1=float(self.Elvlc["mult"][l1])
                    fmult2=float(self.Elvlc["mult"][l2])
                    popmat[l1+ci,l2+ci]+=cc[itemp]*ups[isplups,itemp]/(fmult2*np.sqrt(temp))
                    popmat[l2+ci,l1+ci]+=cc[itemp]*ups[isplups,itemp]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
                    popmat[l1+ci,l1+ci]-=cc[itemp]*ups[isplups,itemp]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
                    popmat[l2+ci,l2+ci]-=cc[itemp]*ups[isplups,itemp]/(fmult2*np.sqrt(temp))
                # proton rates
                for isplups in range(0,npsplups):
                    l1=self.Psplups["lvl1"][isplups]-1
                    l2=self.Psplups["lvl2"][isplups]-1
                    # for proton excitation, the levels are all below the ionization potential
                    de=np.abs(self.Elvlc["eryd"][l2]-self.Elvlc["eryd"][l1])
                    ekt=(de*1.57888e+5)/temp
                    fmult1=float(self.Elvlc["mult"][l1])
                    fmult2=float(self.Elvlc["mult"][l2])
                    popmat[l1+ci,l2+ci]+=cp[itemp]*pups[isplups,itemp]/(fmult2*np.sqrt(temp))
                    popmat[l2+ci,l1+ci]+=cp[itemp]*pups[isplups,itemp]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
                    popmat[l1+ci,l1+ci]-=cp[itemp]*pups[isplups,itemp]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
                    popmat[l2+ci,l2+ci]-=cp[itemp]*pups[isplups,itemp]/(fmult2*np.sqrt(temp))
                # now include ionization rate from
                if ci:
#                   print ' ci = ', ci
                    popmat[1, 0] += self.Density[itemp]*lower.IonizRate['rate'][itemp]
                    popmat[0, 0] -= self.Density[itemp]*lower.IonizRate['rate'][itemp]
                    popmat[0, 1] += self.Density[itemp]*self.RecombRate['rate'][itemp]
                    popmat[1, 1] -= self.Density[itemp]*self.RecombRate['rate'][itemp]
                    #
                    # the ciRate can be computed for all temperatures
                    #
                    for itrans in range(len(cisplups['lvl1'])):
                        lvl1 = cisplups['lvl1'][itrans]
                        lvl2 = cisplups['lvl2'][itrans]
                        de = cisplups['de'][itrans]
                        ekt = (de*1.57888e+5)/temperature
                        mult = lowMult[lvl1-1]
                        cirate = const.collision*self.CiUpsilon[itrans]*np.exp(-ekt)/(np.sqrt(temp)*mult)
                        # this is kind of double booking the ionization rate components
                        popmat[lvl2, lvl1-1] += self.Density[itemp]*cirate[itemp]
                        popmat[lvl1-1, lvl1-1] -= self.Density[itemp]*cirate[itemp]
                if rec:
#                   print ' rec = ', rec
                    popmat[-1,  ci] += self.Density[itemp]*self.IonizRate['rate'][itemp]
                    popmat[ci, ci] -= self.Density[itemp]*self.IonizRate['rate'][itemp]
                    popmat[ci, -1] += self.Density[itemp]*higher.RecombRate['rate'][itemp]
                    popmat[-1, -1] -= self.Density[itemp]*higher.RecombRate['rate'][itemp]
                    #
                    for itrans in range(len(reclvl['lvl1'])):
                        lvl1 = reclvl['lvl1'][itrans]
                        lvl2 = reclvl['lvl2'][itrans]
                        popmat[lvl2+ci+1, -1] += self.Density[itemp]*self.ReclvlRate['rate'][itrans, itemp]
                        popmat[-1, -1] -= self.Density[itemp]*self.ReclvlRate['rate'][itrans, itemp]
                # normalize to unity
                norm=np.ones(nlvls+ci+rec,'float64')
                if ci:
                    norm[0] = 0.
                if rec:
                    norm[-1] = 0.
                popmat[nlvls+ci+rec-1]=norm
                b=np.zeros(nlvls+ci+rec,'float64')
                b[nlvls+ci+rec-1]=1.
                thispop=np.linalg.solve(popmat,b)
                if rec:
                    pop[itemp] = thispop[ci:-rec]
                else:
                    pop[itemp] = thispop[ci:]
            #
        pop=np.where(pop >0., pop,0.)
        self.Population={"temperature":temperature,"density":density,"population":pop, "protonDensity":protonDensity, "ci":ci, "rec":rec}
        #
        return
        #
        # -------------------------------------------------------------------------------------
        #
    def popPlot(self,top=10, saveFile=0):
        """Plots populations vs temperature or density.

        top specifies the number of the most highly populated levels to plot."""
        #self.Population={"temperature":temperature,"density":density,"population":pop}
        try:
            temperature=self.Population["temperature"]
            density=self.Population["density"]
            pop=self.Population["population"]
        except:
            self.populate()
            temperature=self.Population["temperature"]
            density=self.Population["density"]
            pop=self.Population["population"]
        #
        #
        # find the top most populated levels
        #
        lvl=self.Elvlc["lvl"]
        nlvls=len(lvl)
        maxpop=np.zeros(nlvls,'Float32')
        for ilvl in range(nlvls):
            maxpop[ilvl]=pop[:,ilvl].max()
        #
        lvlsort=np.take(lvl,np.argsort(maxpop))
        toplvl=lvlsort[-top:]
        #
#        temp=np.asarray(temperature,'Float32')
        ntemp=temperature.size
        #
        ndens=density.size
        #
        ylabel='Population'
        title=self.Spectroscopic
        #
        pl.figure()
        #
        if chInteractive:
            pl.ion()
        else:
            pl.ioff()
        #
        #
        fontsize=14
        #  first, for ntemp=ndens=1
        if ndens==1 and ntemp==1:
            if chInteractive:
                print ' only a single temperature and density'
            else:
                self.Message = ' only a single temperature and density'
            return
        elif ndens == 1:
            for lvl in toplvl:
                pl.loglog(temperature,pop[:,lvl-1])
                skip=3
                start=divmod(lvl,ntemp)[1]
                for itemp in range(start,ntemp,ntemp/skip):
                    pl.text(temperature[itemp],pop[itemp,lvl-1],str(lvl))
            xlabel='Temperature (K)'
            pl.xlabel(xlabel,fontsize=fontsize)
            pl.ylabel(ylabel,fontsize=fontsize)
            dstr=' -  Density = %10.2e (cm$^{-3}$)' % density
            pl.title(title+dstr,fontsize=fontsize)
            pl.xlim(temperature.min(),temperature.max())
            yl=pl.ylim()
            pl.ylim(yl[0],1.2)
        elif ntemp == 1:
            xlabel=r'Electron Density (cm$^{-3}$)'
            for lvl in toplvl:
                pl.loglog(density,pop[:,lvl-1])
                skip=min(3, ndens)
                start=divmod(lvl,ndens)[1]
                for idens in range(start,ndens,ndens/skip):
                    pl.text(density[idens],pop[idens,lvl-1],str(lvl))
            pl.xlabel(xlabel,fontsize=fontsize)
            pl.ylabel(ylabel,fontsize=fontsize)
            tstr=' -  T = %10.2e (K)' % temperature
            pl.title(title+tstr,fontsize=fontsize)
            pl.xlim(density[density.nonzero()].min(),density.max())
            yl=pl.ylim()
            pl.ylim(yl[0],1.2)
        else:
#            pl.figure()
            ax = pl.subplot(111)
            for lvl in toplvl:
                pl.loglog(temperature,pop[:,lvl-1])
                skip = min(3, ntemp)
                start=divmod(lvl,ntemp)[1]
                for itemp in range(start,ntemp,ntemp/skip):
                    pl.text(temperature[itemp],pop[itemp,lvl-1],str(lvl))
            xlabel='Temperature (K)'
            pl.xlabel(xlabel,fontsize=fontsize)
            pl.ylabel(ylabel,fontsize=fontsize)
#            pl.title(title,fontsize=fontsize)
            pl.xlim(temperature.min(),temperature.max())
            yl=pl.ylim()
#            pl.ylim(yl[0],1.2)
            pl.text(0.1, 0.5,title, horizontalalignment='center', verticalalignment='center', fontsize=fontsize,  transform = ax.transAxes)
            #
            ax2 = pl.twiny()
            xlabel=r'Electron Density (cm$^{-3}$)'
            pl.xlabel(xlabel, fontsize=fontsize)
            pl.loglog(density,pop[:,toplvl[0]], visible=False)
            ax2.xaxis.tick_top()
#            pl.figure()
#            for lvl in toplvl:
#                pl.loglog(density,pop[:,lvl-1])
#                skip = min(3, ntemp)
#                start=divmod(lvl,ndens)[1]
#                for idens in range(start,ndens,ndens/skip):
#                    pl.text(density[idens],pop[idens,lvl-1],str(lvl))
#            xlabel=r'Electron Density (cm$^{-3}$)'
#            pl.xlabel(xlabel,fontsize=fontsize)
#            pl.ylabel(ylabel,fontsize=fontsize)
#            pl.title(title,fontsize=fontsize)
#            pl.xlim(density.min(),density.max())
#            yl=pl.ylim()
            pl.ylim(yl[0],1.2)
        if saveFile:
            pl.savefig(saveFile)
        return
        #
        # -------------------------------------------------------------------------------------
        #
    def emiss(self,temperature=None,density=None,pDensity=None,  wvlRange = None):
        #
        """Calculate and the emissivities for lines of the specified ion.

        wvlRange can be set to limit the calculation to a particular wavelength range

        units:  ergs cm^-3 s^-1 str^-1

        Does not include elemental abundance or ionization fraction

        Wavelengths are sorted """
        #
        #
        doPopulate=False
        try:
            pop=self.Population['population']
        except:
            doPopulate=True
        #
        if temperature != None:
            self.Temperature=np.asarray(temperature,'float32')
            doPopulate=True
        if density != None:
            self.Density=np.asarray(density,'float32')
            doPopulate=True
        if pDensity != None:
            self.PDensity=pDensity
            doPopulate=True
        #
#       nlvls=len(self.Elvlc['lvl'])
##        good=self.Wgfa['avalue'] > 0.
        # using [:] to make a copy things don't change elsewhere
        wvl = self.Wgfa["wvl"][:]
        l1 = self.Wgfa['lvl1'][:]
        l2 = self.Wgfa["lvl2"][:]
        avalue = self.Wgfa["avalue"][:]
        #
        # make sure there are lines in the wavelength range, if specified
        #wvl = [wvl for wvl in mg2.Wgfa['wvl'] if wvl > 400. and wvl < 2000.]
        if type(wvlRange) != types.NoneType:
            l1 = [self.Wgfa['lvl1'][idx] for idx, wvl in enumerate(self.Wgfa['wvl']) if wvl > wvlRange[0] and  wvl < wvlRange[1]]
            l2 = [self.Wgfa['lvl2'][idx] for idx, wvl in enumerate(self.Wgfa['wvl']) if wvl > wvlRange[0] and  wvl < wvlRange[1]]
            avalue = [self.Wgfa['avalue'][idx] for idx, wvl in enumerate(self.Wgfa['wvl']) if wvl > wvlRange[0] and  wvl < wvlRange[1]]
        # this must be at the end
            wvl = [wvl for wvl in self.Wgfa['wvl']if wvl > wvlRange[0] and  wvl < wvlRange[1]]
        #
        # two-photon decays have wvl=0 and nonzero avalues
        zed = wvl.count(0.)
        for one in range(zed):
            idx = wvl.index(0.)
            bad = wvl.pop(idx)
            bad = avalue.pop(idx)
            bad = l1.pop(idx)
            bad = l2.pop(idx)
        wvl=np.abs(np.asarray(wvl))
        nwvl=len(wvl)
        #
        if nwvl == 0:
            self.Emiss = {'errorMessage':self.Spectroscopic+' no lines in this wavelength range'}
            return
        #
        if doPopulate:
            # new values of temperature or density
            self.populate()
            pop=self.Population["population"]
        try:
            ntempden,nlvls=pop.shape
            em=np.zeros((nwvl, ntempden),'Float32')
        except:
            nlvls=len(pop)
            ntempden=1
            em=np.zeros(nwvl,'Float32')
        #
        #
        plotLabels={}
        #
        if self.Defaults['wavelength'] == 'angstrom':
            plotLabels["xLabel"]="Angstroms"
        elif self.Defaults['wavelength'] == 'nm':
            plotLabels["xLabel"]="nanometers"
        elif self.Defaults['wavelength'] == 'kev':
            plotLabels["xLabel"] = "kev"
        #
        if self.Defaults['flux'] == 'energy':
            factor=const.planck*const.light/(4.*const.pi*1.e-8*wvl)
            plotLabels["yLabel"]="ergs cm^-3 s^-1"
        elif self.Defaults['flux'] == 'photon':
            factor=np.ones((nwvl),'Float32')/(4.*const.pi)
            plotLabels["yLabel"]="photons cm^-3 s^-1"
        #
        if ntempden > 1:
            for itempden in range(ntempden):
                for iwvl in range(nwvl):
                    p = pop[itempden,l2[iwvl]-1]
                    em[iwvl, itempden] = factor[iwvl]*p*avalue[iwvl]
            if self.Defaults['wavelength'] == 'kev':
                wvl = const.ev2Ang/np.asarray(wvl)
            elif self.Defaults['wavelength'] == 'nm':
                wvl = wvl/10.
            em = em.take(wvl.argsort(),axis=0)
            wvl.sort()
        else:
            for iwvl in range(0,nwvl):
                p=pop[l2[iwvl]-1]
                em[iwvl]=factor[iwvl]*p*avalue[iwvl]
            if self.Defaults['wavelength'] == 'kev':
                wvlE=const.ev2Ang/np.asarray(wvl)
            elif self.Defaults['wavelength'] == 'nm':
                wvl=wvl/10.
            em=em.take(wvl.argsort())
            wvl.sort()
        self.Emiss = {"wvl":wvl, "emiss":em, "plotLabels":plotLabels}
        return
        #
        # ---------------------------------------------------------------------------
        #
    def emissPlot(self, index=None,  wvlRange=None,  top=10, linLog='lin', relative=0,  verbose=0, saveFile=0 ):
        '''Plot the emissivities.

        wvlRange, a 2 element tuple, list or array determines the wavelength range

        Top specifies to plot only the top strongest lines, default = 10

        linLog specifies a linear or log plot, want either lin or log, default = lin

        normalize = 1 specifies whether to normalize to strongest line, default = 0'''
        #
        title=self.Spectroscopic
        #
        doEmiss=False
        try:
            em = self.Emiss
        except:
            try:
                self.emiss()
                em = self.Emiss
            except:
                print ' emissivities not calculated and emiss() is unable to calculate them'
                print ' perhaps the temperature and/or density are not set'
                return
        emiss = em['emiss']
        wvl = em['wvl']
        temperature = self.Temperature
        density = self.Density
        #
        ndens = density.size
        ntemp = temperature.size
        #
        if ndens == 1 and ntemp > 1:
            if type(index) == types.NoneType:
                index = ntemp/2
            if chInteractive:
                print 'using index = %5i specifying temperature =  %10.2e'%(index, temperature[index])
            else:
                self.Message = 'using index = %5i specifying temperature =  %10.2e'%(index, temperature[index])
            emiss=emiss[:, index]
            dstr=' -  Density = %10.2e (cm$^{-3}$)' % density
            tstr=' -  T = %10.2e (K)' % temperature[index]
        elif ndens > 1 and ntemp == 1:
            if type(index) == types.NoneType:
                index = ndens/2
            if chInteractive:
                print 'using index =%5i specifying density = %10.2e'%(index, density[index])
            else:
                self.Message = 'using index =%5i specifying density = %10.2e'%(index, density[index])
            emiss=emiss[:, index]
            dstr=' -  Density = %10.2e (cm$^{-3}$)' % density[index]
            tstr=' -  T = %10.2e (K)' % temperature
        elif ndens > 1 and ntemp > 1:
            if type(index) == types.NoneType:
                index = ntemp/2
            if chInteractive:
                print 'using index = %5i specifying temperature = %10.2e, density =  %10.2e'%(index, temperature[index], density[index])
            else:
                self.Message = 'using index = %5i specifying temperature = %10.2e, density =  %10.2e'%(index, temperature[index], density[index])
            emiss=emiss[:, index]
            dstr=' -  Density = %10.2e (cm$^{-3}$)' % density[index]
            tstr=' -  T = %10.2e (K)' % temperature[index]
        if type(wvlRange) != types.NoneType:
            wvlIndex = util.between(wvl, wvlRange)
        else:
            wvlIndex = range(wvl.size)
        emiss = emiss[wvlIndex]
        wvl = wvl[wvlIndex]
        #
        self.Error = 0
        if wvl.size == 0:
            if chInteractive:
                print 'No lines in this wavelength interval'
            else:
                self.Error = 1
                self.Message = 'No lines in this wavelength interval'
            return
        elif top == 0:
            top = wvl.size
        elif wvl.size > top:
            isrt = np.argsort(emiss)
            wvl = wvl[isrt[-top:]]
            emiss = emiss[isrt[-top:]]
        else:
            top = wvl.size
        # must follow setting top
        #
        pl.figure()
        ylabel = 'Emissivity'
        if relative:
            emiss = emiss/emiss[:top].max()
            ylabel += ' (Relative)'
        #
        xlabel = 'Wavelength ('+self.Defaults['wavelength'] +')'
        #
        ymin = 10.**(np.log10(emiss.min()).round(0))
        #
        if chInteractive:
            pl.ion()
        else:
            pl.ioff()
        #
        for idx in range(top):
            xx=[wvl[idx], wvl[idx]]
            if linLog == 'lin':
                yy=[0., emiss[idx]]
                pl.plot(xx, yy)
            else:
                yy=[ymin/10., emiss[idx]]
                pl.semilogy(xx, yy)
        pl.xlabel(xlabel)
        pl.ylabel(ylabel)
        pl.title(title+tstr+dstr)
        if wvlRange:
            pl.axis([wvlRange[0], wvlRange[1], ymin, emiss.max()])
        if saveFile:
            pl.savefig(saveFile)
        #
        idx = np.argsort(wvl)
        self.Emiss['wvlTop'] = wvl[idx]
        self.Emiss['emissTop'] = emiss[idx]
        #
        # ---------------------------------------------------------------------------
        #
    def intensity(self,  wvlRange = None):
        """Calculate  the intensities for lines of the specified ion.

        wvlRange, a 2 element tuple, list or array determines the wavelength range

        units:  ergs cm^-3 s^-1 str^-1

        includes elemental abundance and ionization fraction."""
        # emiss ={"wvl":wvl, "emiss":em, "plotLabels":plotLabels}
        #
        self.emiss(wvlRange = wvlRange)
        emiss = self.Emiss
        if 'errorMessage'  in emiss.keys():
            self.Intensity = {'errorMessage': self.Spectroscopic+' no lines in this wavelength region'}
            return
        em = emiss['emiss']
        wvl = emiss['wvl']
        try:
            ab=self.Abundance
        except:
            self.Abundance = util.abundanceRead()
            ab=self.Abundance
        try:
            thisIoneq=self.IoneqOne
        except:
            self.ioneqOne()
            thisIoneq=self.IoneqOne
        try:
            nwvl, ntempden = em.shape
            intensity = np.zeros((ntempden, nwvl),'Float64')
            if thisIoneq.size == 1:
                thisIoneq = np.ones(ntempden, 'float64')*thisIoneq
            for it in range(ntempden):
                #  already done in emiss
#                if self.Defaults['flux'] == 'energy':
#                    intensity[it] = (const.planck*const.light*1.e+8/wvl)*ab*thisIoneq[it]*em[:, it]
#                else:
#                    intensity[it] = ab*thisIoneq[it]*em[:, it]
                intensity[it] = ab*thisIoneq[it]*em[:, it]/(4.*const.pi)
        except:
            nwvl=len(em)
            ntempden=1
#            intensity = np.zeros(nwvl,'Float32')
# this already done in emiss
#            if self.Defaults['flux'] == 'energy':
#                intensity = (const.planck*const.light*1.e+8/wvl)*ab*thisIoneq*em
#            else:
#                intensity = ab*thisIoneq*em
            intensity = ab*thisIoneq*em/(4.*const.pi)
        self.Intensity = {'intensity':intensity, 'wvl':wvl}
        #
        # -------------------------------------------------------------------------------------
        #
    def intensityRatio(self,wvlRange=None,top=10,temperature=None,density=None,pDensity=None):
        """Plot the ratio of 2 lines or sums of lines.

        Shown as a function of density and/or temperature.

        A plot of relative emissivities is shown and then a dialog appears for the user to

        choose a set of lines."""
        #
        #        self.Emiss={"temperature":temperature,"density":density,"wvl":wvl,"emiss":em,
        #        "plotLabels":plotLabels}
        #
        doEmiss=False
        try:
            em = self.Emiss
        except:
            doEmiss = True
        #
        if temperature != None:
            self.Temperature=np.asarray(temperature,'float32')
            doEmiss=True
        if density != None:
            self.Density=np.asarray(density,'float32')
            doEmiss=True
        if pDensity:
            self.PDensity=pDensity
            doEmiss=True
        #
        if doEmiss:
            # new values of temperature or density
            self.emiss(temperature=temperature, density=density, pDensity=pDensity)
            em=self.Emiss
        #
        #
        fontsize=14
        #
        temperature = self.Temperature
        density = self.Density
        emiss = em['emiss']
        wvl = em["wvl"]
        plotLabels=em["plotLabels"]
        xLabel=plotLabels["xLabel"]
        yLabel=plotLabels["yLabel"]
        #
        # find which lines are in the wavelength range if it is set
        #
        #
        if not wvlRange:
            igvl=range(len(wvl))
        else:
            igvl=util.between(wvl,wvlRange)
        nlines=len(igvl)
        #
#        print ' nlines = ',nlines
#        print ' iglv = ',igvl
        igvl=np.take(igvl,wvl[igvl].argsort())
        # find the top most intense lines
        #
        if top > nlines:  top=nlines
        maxEmiss=np.zeros(nlines,'Float32')
        for iline in range(nlines):
            maxEmiss[iline]=emiss[igvl[iline]].max()
        for iline in range(nlines):
            if maxEmiss[iline]==maxEmiss.max():
                maxAll=emiss[igvl[iline]]
        line=range(nlines)
        igvlsort=np.take(igvl,np.argsort(maxEmiss))
#        print 'igvlsort = ', igvlsort
        topLines=igvlsort[-top:]
#        print ' topLines = ', topLines
        maxWvl='%5.3f' % wvl[topLines[-1]]
        maxline=topLines[-1]
        #
        topLines=topLines[wvl[topLines].argsort()]
        #
        #
        # need to make sure there are no negative values before plotting
        good = np.where(emiss > 0.)
        emissMin=emiss[good].min()
        bad=np.where(emiss <= 0.)
        emiss[bad]=emissMin
        #
        #
        ntemp=self.Temperature.size
        #
        ndens=self.Density.size
        #
        ylabel='Emissivity relative to '+maxWvl
        title=self.Spectroscopic
        #
        #
        if ndens==1 and ntemp==1:
            print ' only a single temperature and density'
            return
        elif ndens == 1:
            xlabel='Temperature (K)'
            xvalues=self.Temperature
            outTemperature=self.Temperature
            outDensity=np.zeros(ntemp,'Float32')
            outDensity.fill(self.Density)
            desc_str=' at  Density = %10.2e (cm)$^{-3}$' % self.Density
        elif ntemp == 1:
            xvalues=self.Density
            outTemperature=np.zeros(ndens,'Float32')
            outTemperature.fill(self.Temperature)
            outDensity=self.Density
            xlabel=r'$\rm{Electron Density (cm)^{-3}}$'
            desc_str=' at Temp = %10.2e (K)' % self.Temperature
        else:
            outTemperature=self.Temperature
            outDensity=self.Density
            xlabel='Temperature (K)'
            xvalues=self.Temperature
            desc_str=' for variable Density'
        #
        # put all actual plotting here
        #
        if chInteractive:
            pl.ion()
        else:
            pl.ioff()
        #
        #  maxAll is an array
        ymax = np.max(emiss[topLines[0]]/maxAll)
        ymin = ymax
        pl.figure()
        ax = pl.subplot(111)
        nxvalues=len(xvalues)
        for iline in range(top):
            tline=topLines[iline]
            pl.loglog(xvalues,emiss[tline]/maxAll)
            if np.min(emiss[tline]/maxAll) < ymin:
                ymin = np.min(emiss[tline]/maxAll)
            if np.max(emiss[tline]/maxAll) > ymax:
                ymax = np.max(emiss[tline]/maxAll)
            skip=2
            start=divmod(iline,nxvalues)[1]
            for ixvalue in range(start,nxvalues,nxvalues/skip):
                pl.text(xvalues[ixvalue],emiss[tline,ixvalue]/maxAll[ixvalue],str(wvl[tline]))
        pl.xlim(xvalues.min(),xvalues.max())
#        pl.ylim(ymin, ymax)
        pl.xlabel(xlabel,fontsize=fontsize)
        pl.ylabel(ylabel,fontsize=fontsize)
        if ndens == ntemp and ntemp > 1:
            pl.text(0.07, 0.5,title, horizontalalignment='left', verticalalignment='center', fontsize=fontsize,  transform = ax.transAxes)
            #
            ax2 = pl.twiny()
            xlabelDen=r'Electron Density (cm$^{-3}$)'
            pl.xlabel(xlabelDen, fontsize=fontsize)
            pl.loglog(density,emiss[topLines[top-1]]/maxAll, visible=False)
            ax2.xaxis.tick_top()
            pl.ylim(ymin/1.2, 1.2*ymax)
        else:
            pl.ylim(ymin/1.2, 1.2*ymax)
            pl.title(title+desc_str,fontsize=fontsize)
        pl.draw()
        # get line selection
        #
        numden = gui.choice2Dialog(wvl[topLines])
        #
        # num_idx and den_idx are tuples
        #
        num_idx=numden.numIndex
        if len(num_idx) == 0:
            print ' no numerator lines were selected'
            return
        #
        den_idx=numden.denIndex
        if len(den_idx) == 0:
            print ' no denominator lines were selected'
            return
        #
        numEmiss=np.zeros(len(xvalues),'float32')
        for aline in num_idx:
            numEmiss+=emiss[topLines[aline]]
        #
        denEmiss=np.zeros(len(xvalues),'float32')
        for aline in den_idx:
            denEmiss+=emiss[topLines[aline]]
        #
        # plot the desired ratio
        #  maxAll is an array
        pl.figure()
        ax = pl.subplot(111)
        pl.loglog(xvalues,numEmiss/denEmiss)
        pl.xlim(xvalues.min(),xvalues.max())
        pl.xlabel(xlabel,fontsize=fontsize)
        pl.ylabel('Ratio ('+self.Defaults['flux']+')',fontsize=fontsize)
        desc = title + ':'
        for aline in num_idx:
            desc += ' ' + str(wvl[topLines[aline]])
        desc +=' / '
        for aline in den_idx:
            desc += ' ' + str(wvl[topLines[aline]])
        if ndens == ntemp and ntemp > 1:
            pl.text(0.07, 0.5,desc, horizontalalignment='left', verticalalignment='center', fontsize=fontsize,  transform = ax.transAxes)
            #
            ax2 = pl.twiny()
            xlabelDen=r'Electron Density (cm$^{-3}$)'
            pl.xlabel(xlabelDen, fontsize=fontsize)
            pl.loglog(density,numEmiss/denEmiss, visible=False)
            ax2.xaxis.tick_top()
        else:
#            pl.ylim(ymin, ymax)
            pl.title(desc,fontsize=fontsize)
#       desc=title+' '+str(wvl[num_line])+' / '+str(wvl[den_line])+' '+desc_str
#        pl.title(desc, fontsize=fontsize)
#       pl.title(title+' '+str(wvl[num_line])+' / '+str(wvl[den_line])+' '+desc_str,fontsize=fontsize)
#        pl.draw()
#        pl.ioff()
#        pl.show()
        #
        intensityRatioFileName=self.IonStr
        for aline in num_idx:
            intensityRatioFileName+= '_%3i'%(wvl[topLines[aline]])
        intensityRatioFileName+='_2'
        for aline in den_idx:
            intensityRatioFileName+= '_%3i'%(wvl[topLines[aline]])
        intensityRatioFileName+='.rat'
        self.IntensityRatio={'ratio':numEmiss/denEmiss,'desc':desc,
                'temperature':outTemperature,'density':outDensity,'filename':intensityRatioFileName}
        #
        # -------------------------------------------------------------------------------------
        #
    def intensityRatioSave(self,outfile=''):
        '''Save the intensity ratio to a file.

        The intensity ratio as a function to temperature and density is saved to an asciii file.

        Descriptive information is included at the top of the file.'''
        if outfile == '':
            outfile=self.IntensityRatio['filename']
            print ' filename = ',outfile
        temperature=self.IntensityRatio['temperature']
        density=self.IntensityRatio['density']
        ratio=self.IntensityRatio['ratio']
        out=open(outfile,'w')
        nvalues=len(ratio)
        #
        #  need to add 7 lines to maintain IDL like files
        #
        out.write(outfile+'\n')    #1
        out.write(self.IntensityRatio['desc']+'\n') #2
        out.write(' created with ChiantiPy version '+ self.__version__ +'\n')   #3
        out.write(' columns are temperature, density, ratio'+'\n')  #5
        tunit = 'K'
        out.write(' temperature in '+tunit+', electron density in cm^(-3)'+'\n')  #6
        out.write(' ratio given in '+self.Defaults['flux']+'\n')   #4
        out.write(' '+'\n') #7
        for ivalue in range(nvalues):
            s='%12.3e %12.3e  %12.3e ' % (temperature[ivalue],density[ivalue],ratio[ivalue])
            out.write(s+os.linesep)
        out.close()
        #
        # -------------------------------------------------------------------------------------
        #
    def ioneqOne(self):
        '''Provide the ionization equilibrium for the selected ion as a function of temperature.
        returned in self.IoneqOne'''
        #
        try:
            temperature = self.Temperature
        except:
            return
        #
        try:
            ioneqAll = self.IoneqAll
        except:
            self.IoneqAll = util.ioneqRead(ioneqname = self.Defaults['ioneqfile'])
            ioneqAll=self.IoneqAll
        #
        ioneqTemperature = ioneqAll['ioneqTemperature']
        Z=self.Z
        Ion=self.Ion
        Dielectronic=self.Dielectronic
        #
        thisIoneq=ioneqAll['ioneqAll'][Z-1,Ion-1-Dielectronic].squeeze()
#        thisIoneq = self.Ioneq
        gioneq=thisIoneq > 0.
        y2=interpolate.splrep(np.log(ioneqTemperature[gioneq]),np.log(thisIoneq[gioneq]),s=0)
        goodt1=self.Temperature >= ioneqTemperature[gioneq].min()
        goodt2=self.Temperature <= ioneqTemperature[gioneq].max()
        goodt=np.logical_and(goodt1,goodt2)
        #
        if goodt.sum() > 0:
            gIoneq=interpolate.splev(np.log(self.Temperature),y2)   #,der=0)
            gIoneq=np.exp(gIoneq)
        else:
            gIoneq=0.
        #
        self.IoneqOne=gIoneq
        #
        # -------------------------------------------------------------------------------------
        #
    def gofnt(self,wvlRange=0,top=10,temperature=None,density=None,pDensity=None,verbose=0):
        """Calculate the 'so-called' G(T) function.

        Given as a function of both temperature and density.

        Only the top( set by 'top') brightest lines are plotted.
        the G(T) function is returned in a dictionary self.Gofnt"""
        #
        #self.emiss={"wvl":wvl,"emiss":em,"units":units,"plotLabels":plotLabels}
        #
        #
        doEmiss=False
        try:
            em=self.Emiss
        except:
            doEmiss=True
        #
        if temperature != None:
            self.Temperature=np.asarray(temperature,'float32')
            doEmiss=True
        if density != None:
            self.Density=np.asarray(density,'float32')
            doEmiss=True
        if pDensity:
            self.PDensity=pDensity
            doEmiss=True
        #
        if doEmiss:
            # new values of temperature or density
            self.emiss()
            em=self.Emiss
        #
        #
        try:
            ab=self.Abundance
        except:
            self.Abundance = util.abundanceRead()
            ab=self.Abundance
        #
        fontsize=12
        #
        emiss=em["emiss"]
        wvl=em["wvl"]
        temperature=self.Temperature
        density=self.Density
        plotLabels=em["plotLabels"]
        xLabel=plotLabels["xLabel"]
        yLabel=plotLabels["yLabel"]
        #
        # find which lines are in the wavelength range if it is set
        #
        #
        if type(wvlRange) != type(1):
            igvl=util.between(wvl,wvlRange)
        else:
            igvl=range(len(wvl))
        nlines=len(igvl)
        if nlines ==0:
            print ' no lines in selected interval'
            return
        # find the top most intense lines
        #
        if top > nlines:
            top=nlines
        maxEmiss=np.zeros(nlines,'Float32')
        for iline in range(nlines):
            maxEmiss[iline]=emiss[igvl[iline]].max()
        for iline in range(nlines):
            if maxEmiss[iline]>=maxEmiss.max():
                maxAll=emiss[igvl[iline]]
                maxIndex = igvl[iline]
#        print ' maxIndex, maxAll = ', maxIndex,  maxAll
        line=range(nlines)
        igvlsort=np.take(igvl,np.argsort(maxEmiss))
        topLines=igvlsort[-top:]
        maxWvl='%5.3f' % wvl[topLines[-1]]
        maxline=topLines[-1]
        #
        # need to make sure there are no negative values before plotting
        good = np.where(emiss > 0.)
        emissMin=emiss[good].min()
        bad=np.where(emiss <= 0.)
        emiss[bad]=emissMin
        #
        topLines=topLines[wvl[topLines].argsort()]
        #
        #
        ntemp=self.Temperature.size
        #
        ndens=self.Density.size
        #
        ylabel = 'Emissivity relative to '+maxWvl
        title = self.Spectroscopic
        #
        #
        if ndens==1 and ntemp==1:
            print ' only a single temperature and density'
            return
        elif ndens == 1:
            xlabel='Temperature (K)'
            ngofnt = temperature.size
            xvalues=temperature
            outTemperature=temperature
            outDensity=np.zeros(ntemp,'Float32')
            outDensity.fill(density)
            desc_str=' at Density = %10.2e' % density
        elif ntemp == 1:
            xvalues=density
            ngofnt = density.size
            outTemperature=np.zeros(ndens,'Float32')
            outTemperature.fill(temperature)
            outDensity=density
            xlabel=r'$\rm{Electron Density (cm}^{-3}\rm{)}$'
            desc_str=' at Temperature = %10.2e' % temperature
        else:
            outTemperature=temperature
            outDensity=density
            xlabel='Temperature (K)'
            xvalues=temperature
            ngofnt = temperature.size
            desc_str=' for variable Density'
            #
        #
        # put all actual plotting here
        #
        if chInteractive:
            pl.ion()
        else:
            pl.ioff()
        #
        pl.figure()
        ax = pl.subplot(111)
        nxvalues=len(xvalues)
        #  maxAll is an array
        ymax = np.max(1.2*emiss[top-1]/maxAll)
        ymin = ymax
        for iline in range(top):
            tline=topLines[iline]
            pl.loglog(xvalues,emiss[tline]/maxAll)
            if np.min(emiss[tline]/maxAll) < ymin:
                ymin = np.min(emiss[tline]/maxAll)
            skip=2
            start=divmod(iline,nxvalues)[1]
            for ixvalue in range(start,nxvalues,nxvalues/skip):
                pl.text(xvalues[ixvalue],emiss[tline,ixvalue]/maxAll[ixvalue],str(wvl[tline]))
        pl.xlim(xvalues.min(),xvalues.max())
        pl.ylim(ymin, ymax)
#       yl=pl.ylim()
#       pl.ylim(yl[0],1.2)
        pl.xlabel(xlabel,fontsize=fontsize)
        pl.ylabel(ylabel,fontsize=fontsize)
        if ndens == ntemp and ntemp > 1:
            pl.text(0.07, 0.5,title, horizontalalignment='left', verticalalignment='center', fontsize=fontsize,  transform = ax.transAxes)
            #
            ax2 = pl.twiny()
            xlabelDen=r'Electron Density (cm$^{-3}$)'
            pl.xlabel(xlabelDen, fontsize=fontsize)
            pl.loglog(density,emiss[topLines[top-1]]/maxAll, visible=False)
            ax2.xaxis.tick_top()
        else:
            pl.ylim(ymin, ymax)
            pl.title(title+desc_str,fontsize=fontsize)
        pl.draw()
        #
#        print ' topLInes = ', wvl[topLines]
        wvlChoices = []
        for one in wvl[topLines]:
            wvlChoices.append('%12.3f'%(one))
        gline = gui.selectorDialog(wvlChoices,label='Select line(s)')
        gline_idx=gline.selectedIndex
        #
        #
        gAbund=self.Abundance
        #
        if verbose:
            print ' abundance, ioneq = ',gAbund,thisIoneq
        try:
            thisIoneq=self.IoneqOne
        except:
            self.ioneqOne()
        #        gioneq=np.where(thisIoneq > 0.)
        #        y2=interpolate.splrep(np.log(self.IoneqAll['ioneqTemperature'][gioneq]),np.log(thisIoneq[gioneq]),s=0)  #allow smoothing,s=0)
        #        gIoneq=interpolate.splev(np.log(temperature),y2)   #,der=0)
        gIoneq=self.IoneqOne/density
        #
        #
        #
        # plot the desired ratio
        pl.figure()
        g_line= topLines[gline_idx]#  [0]
        ##        print ' g_line = ',g_line
        #
        gofnt=np.zeros(ngofnt,'float32')
        for aline in g_line:
            gofnt+=gAbund*gIoneq*emiss[aline].squeeze()
        self.Gofnt={'temperature':outTemperature,'density':outDensity,'gofnt':gofnt}
        #
        pl.loglog(xvalues,gofnt)
        pl.xlim(xvalues.min(),xvalues.max())
        pl.xlabel(xlabel,fontsize=fontsize)
        pl.ylabel('Gofnt',fontsize=fontsize)
        if ndens == ntemp and ntemp > 1:
            newTitle = title+' '+str(wvl[g_line])+' '+desc_str
            pl.text(0.07, 0.5,newTitle, horizontalalignment='left', verticalalignment='center', fontsize=fontsize,  transform = ax.transAxes)
            #
            ax2 = pl.twiny()
#            xlabel=r'Electron Density (cm$^{-3}$)'
            pl.xlabel(xlabelDen, fontsize=fontsize)
            pl.loglog(density,gofnt, visible=False)
            ax2.xaxis.tick_top()
        else:
            pl.title(title+' '+str(wvl[g_line])+' '+desc_str, fontsize=fontsize)
        #pl.ioff()
        #pl.show()
#        return
        #
        # - - - - - - - - - - - - - - - - - - - - - - -
        #
    def twoPhotonEmiss(self, wvl):
        ''' to calculate the two-photon continuum rate coefficient - only for hydrogen- and helium-like ions'''
        wvl = np.array(wvl, 'float64')
        nWvl = wvl.size
        if self.Z -self.Ion > 1 or self.Dielectronic:
            # this is not a hydrogen-like or helium-like ion
            self.TwoPhoton = {'emiss':np.zeros(nWvl, 'float4'), 'wvl':wvl}
            return
        else:
            try:
                pop = self.Population['population']
                nTempDens = max(self.Temperature.size, self.Density.size)
            except:
                self.populate()
                pop = self.Population['population']
                nTempDens = max(self.Temperature.size, self.Density.size)
            if nTempDens > 1:
                emiss = np.zeros((nTempDens, nWvl), 'float64')
                if self.Density.size == 1:
                    density = np.repeat(self.Density, nTempDens)
                else:
                    density = self.Density
            else:
                emiss = np.zeros(nWvl, 'float64')
                density = self.Density
            if self.Z == self.Ion:
                # H seq
                l1 = 1-1
                l2 = 2 - 1
                wvl0 = 1.e+8/(self.Elvlc['ecm'][l2] - self.Elvlc['ecm'][l1])
                goodWvl = wvl > wvl0
                y = wvl0/wvl[goodWvl]
                dist = util.twophotonHRead()
                avalue = dist['avalue'][self.Z-1]
                asum = dist['asum'][self.Z-1]
                distr1 = interpolate.splrep(dist['y0'], dist['psi0'][self.Z-1], s=0)
                distr = avalue*y*interpolate.splev(y, distr1)/(asum*wvl[goodWvl])
                if self.Defaults['flux'] == 'energy':
                    f = (const.light*const.planck*1.e+8)/wvl[goodWvl]
                else:
                    f=1.
                if nTempDens == 1:
                    emiss[goodWvl] = f*pop[l2]*distr/self.Density
                else:
                    for it in range(nTempDens):
                        emiss[it, goodWvl] = f*pop[it, l2]*distr/self.Density[it]
                self.TwoPhotonEmiss = {'wvl':wvl, 'emiss':emiss}
            else:
                # He seq
                l1 = 1-1
                l2 = 3 - 1
                wvl0 = 1.e+8/(self.Elvlc['ecm'][l2] - self.Elvlc['ecm'][l1])
                goodWvl = wvl > wvl0
                y = wvl0/wvl[goodWvl]
                dist = util.twophotonHeRead()
                avalue = dist['avalue'][self.Z-1]
                distr1 = interpolate.splrep(dist['y0'], dist['psi0'][self.Z-1], s=0)
                distr = avalue*y*interpolate.splev(y, distr1)/wvl[goodWvl]
                if self.Defaults['flux'] == 'energy':
                    f = (const.light*const.planck*1.e+8)/wvl[goodWvl]
                else:
                    f=1.
                if nTempDens == 1:
                    emiss[goodWvl] = f*pop[l2]*distr/self.Density
                else:
                    for it in range(nTempDens):
                        emiss[it, goodWvl] = f*pop[it, l2]*distr/self.Density[it]
                self.TwoPhotonEmiss = {'wvl':wvl, 'emiss':emiss}
        #
        #-----------------------------------------------------------------
        #
    def twoPhoton(self, wvl):
        ''' to calculate the two-photon continuum - only for hydrogen- and helium-like ions
        includes the elemental abundance and the ionization equilibrium'''
        wvl = np.array(wvl, 'float64')
        nWvl = wvl.size
        if self.Z -self.Ion > 1 or self.Dielectronic:
            # this is not a hydrogen-like or helium-like ion
            print ' not doing 2 photon for ', self.Ions
            self.TwoPhoton = {'emiss':np.zeros(nWvl, 'float64'), 'wvl':wvl}
            return
        else:
            try:
                ab=self.Abundance
            except:
                self.Abundance = util.abundanceRead()
                ab=self.Abundance
            try:
                thisIoneq=self.IoneqOne
            except:
                self.ioneqOne()
                thisIoneq=self.IoneqOne
            try:
                pop = self.Population['population']
                nTempDens = max(self.Temperature.size, self.Density.size)
            except:
                self.populate()
                pop = self.Population['population']
                nTempDens = max(self.Temperature.size, self.Density.size)
            if nTempDens > 1:
                rate = np.zeros((nTempDens, nWvl), 'float64')
                if self.Density.size == 1:
                    density = np.repeat(self.Density, nTempDens)
                else:
                    density = self.Density
            else:
                rate = np.zeros(nWvl, 'float64')
                density = self.Density
            if self.Z == self.Ion:
                # H seq
                l1 = 1-1
                l2 = 2 - 1
                wvl0 = 1.e+8/(self.Elvlc['ecm'][l2] - self.Elvlc['ecm'][l1])
                goodWvl = wvl > wvl0
                y = wvl0/wvl[goodWvl]
                dist = util.twophotonHRead()
                avalue = dist['avalue'][self.Z-1]
                asum = dist['asum'][self.Z-1]
                distr1 = interpolate.splrep(dist['y0'], dist['psi0'][self.Z-1], s=0)
                distr = avalue*y*interpolate.splev(y, distr1)/(asum*wvl[goodWvl])
                if self.Defaults['flux'] == 'energy':
                    f = (const.light*const.planck*1.e+8)/(4.*const.pi*wvl[goodWvl])
                else:
                    f=1./(4.*const.pi)
                if nTempDens == 1:
                    rate[goodWvl] = f*pop[l2]*distr*ab*thisIoneq/density
                else:
                   for it in range(nTempDens):
                        rate[it, goodWvl] = f*pop[it, l2]*distr*ab*thisIoneq[it]/density[it]
                self.TwoPhoton = {'wvl':wvl, 'rate':rate}
            else:
                # He seq
                l1 = 1-1
                l2 = 3 - 1
                wvl0 = 1.e+8/(self.Elvlc['ecm'][l2] - self.Elvlc['ecm'][l1])
                goodWvl = wvl > wvl0
                y = wvl0/wvl[goodWvl]
                dist = util.twophotonHeRead()
                avalue = dist['avalue'][self.Z-1]
                distr1 = interpolate.splrep(dist['y0'], dist['psi0'][self.Z-1], s=0)
                distr = avalue*y*interpolate.splev(y, distr1)/wvl[goodWvl]
                if self.Defaults['flux'] == 'energy':
                    f = (const.light*const.planck*1.e+8)/(4.*const.pi*wvl[goodWvl])
                else:
                    f=1./(4.*const.pi)
                if nTempDens == 1:
                    rate[goodWvl] = f*pop[l2]*distr*ab*thisIoneq/density
                else:
                   for it in range(nTempDens):
                        rate[it, goodWvl] = f*pop[it, l2]*distr*ab*thisIoneq[it]/density[it]
                self.TwoPhoton = {'wvl':wvl, 'rate':rate}
        #
        #-----------------------------------------------------------------
        #
    def twoPhotonLoss(self):
        ''' to calculate the two-photon energy loss rate - only for hydrogen- and helium-like ions
        includes the elemental abundance and the ionization equilibrium'''
        if self.Z -self.Ion > 1 or self.Dielectronic:
            # this is not a hydrogen-like or helium-like ion
            print ' not doing 2 photon for ', self.Ions
            self.TwoPhoton = {'emiss':np.zeros(nWvl, 'float64'), 'wvl':wvl}
            return
        else:
            try:
                ab=self.Abundance
            except:
                self.Abundance = util.abundanceRead()
                ab=self.Abundance
            try:
                thisIoneq=self.IoneqOne
            except:
                self.ioneqOne()
                thisIoneq=self.IoneqOne
            try:
                pop = self.Population['population']
                nTempDens = max(self.Temperature.size, self.Density.size)
            except:
                self.populate()
                pop = self.Population['population']
                nTempDens = max(self.Temperature.size, self.Density.size)
            if nTempDens > 1:
                rate = np.zeros((nTempDens), 'float64')
                if self.Density.size == 1:
                    density = np.repeat(self.Density, nTempDens)
                else:
                    density = self.Density
            else:
                density = self.Density
            if self.Z == self.Ion:
                # H seq
                l1 = 1-1
                l2 = 2 - 1
                wvl0 = 1.e+8/(self.Elvlc['ecm'][l2] - self.Elvlc['ecm'][l1])
                dist = util.twophotonHRead()
                avalue = dist['avalue'][self.Z-1]
                f = (avalue*const.light*const.planck*1.e+8)/wvl0
                if nTempDens == 1:
                    rate = f*pop[l2]*ab*thisIoneq/density
                else:
                   for it in range(nTempDens):
                        rate[it] = f*pop[it, l2]*ab*thisIoneq[it]/density[it]
                self.TwoPhotonLoss = {'temperature':self.Temperature,'density':self.Density,'rate':rate}
            else:
                # He seq
                l1 = 1-1
                l2 = 3 - 1
                wvl0 = 1.e+8/(self.Elvlc['ecm'][l2] - self.Elvlc['ecm'][l1])
                dist = util.twophotonHeRead()
                avalue = dist['avalue'][self.Z-1]
                f = (avalue*const.light*const.planck*1.e+8)/wvl0
                if nTempDens == 1:
                    rate = f*pop[l2]*ab*thisIoneq/density
                else:
                   for it in range(nTempDens):
                        rate[it] = f*pop[it, l2]*ab*thisIoneq[it]/density[it]
                self.TwoPhotonLoss = {'temperature':self.Temperature,'density':self.Density,'rate':rate}
        #
        # ----------------------------------------------
        #
class ioneq(ion):
    '''Calculates the ionization equilibrium for an element as a function of temperature.
    The variable z is the atomic number of the element.  Acceptable values are from 1 to 30.'''
    def __init__(self,z, temperature, verbose=False):
#        self.Defaults=defaults
        ionList=[]
        chIons=[]
        self.Z=z
        self.Temperature = np.array(temperature, 'float64')
        for stage in range(1, z+2):
            ionStr=util.zion2name(z, stage)
            ionList.append(ionStr)
            print z, stage, ionStr
            atom=ion(ionStr, temperature = self.Temperature)
            atom.ionizRate()
            atom.recombRate()
            chIons.append(atom)
#        for anIon in chIons:
#            print ' this ion = ', anIon.Ions
#            if type(anIon.IonizRate) != NoneType:
#                pl.loglog(anIon.IonizRate['temperature'], anIon.IonizRate['rate'])
#        #
#        for anIon in chIons:
#            print ' this ion = ',  anIon.Ions
#            if type(anIon.RecombRate) != NoneType:
#                pl.loglog(anIon.RecombRate['temperature'], anIon.RecombRate['rate'])
        #
        ntemp=chIons[0].IonizRate['temperature'].size
        print ' ntemp = ',ntemp
        if ntemp == 1:
            ioneq=np.zeros((z+1), 'float32')
            factor = []
            for anIon in chIons:
                if type(anIon.IonizRate) != types.NoneType and type(anIon.RecombRate) != types.NoneType:
                    rat=anIon.IonizRate['rate']/anIon.RecombRate['rate']
                    factor.append(rat**2 + rat**(-2))
                else:
                    factor.append(0.)
            factor[0]=max(factor)
            factor[-1]=max(factor)
            ionmax=factor.index(min(factor))
#            print ' it, ionmax', it, ionmax
            ioneq[ionmax]=1.
            #
            for iz in range(ionmax+1, z+1):
                ionrate=chIons[iz-1].IonizRate['rate']
                recrate=chIons[iz].RecombRate['rate']
                ioneq[iz]=ionrate*ioneq[iz-1]/recrate
            #
            for iz in range(ionmax-1, -1, -1):
                ionrate=chIons[iz].IonizRate['rate']
                recrate=chIons[iz+1].RecombRate['rate']
                ioneq[iz]=recrate*ioneq[iz+1]/ionrate
            ionsum=ioneq.sum()
#            print ' ionsum = ', ionsum
            ioneq=ioneq/ionsum
            self.Ioneq=ioneq
        #  ntemp >1
        else:
            ioneq=np.zeros((z+1,ntemp ), 'float32')
            for it in range(ntemp):
                factor=[]
                for anIon in chIons:
                    if type(anIon.IonizRate) != types.NoneType and type(anIon.RecombRate) != types.NoneType:
                        rat=anIon.IonizRate['rate'][it]/anIon.RecombRate['rate'][it]
                        factor.append(rat**2 + rat**(-2))
                    else:
                        factor.append(0.)
                factor[0]=max(factor)
                factor[-1]=max(factor)
                ionmax=factor.index(min(factor))
    #            print ' it, ionmax', it, ionmax
                ioneq[ionmax, it]=1.
                #
                for iz in range(ionmax+1, z+1):
                    ionrate=chIons[iz-1].IonizRate['rate'][it]
                    recrate=chIons[iz].RecombRate['rate'][it]
                    ioneq[iz, it]=ionrate*ioneq[iz-1, it]/recrate
                #
                for iz in range(ionmax-1, -1, -1):
                    ionrate=chIons[iz].IonizRate['rate'][it]
                    recrate=chIons[iz+1].RecombRate['rate'][it]
                    ioneq[iz, it]=recrate*ioneq[iz+1, it]/ionrate
                ionsum=ioneq[:, it].sum()
    #            print ' ionsum = ', ionsum
                ioneq[:, it]=ioneq[:, it]/ionsum
            self.Ioneq=ioneq
#
    def plot(self, stages=None, xr=None, yr=None, oplot=False, label=True, title=True,  bw=False):
        '''Plots the ionization equilibria.

        self.plot(xr=None, yr=None, oplot=False)
        stages = sequence of ions to be plotted, neutral == 1, fully stripped == Z+1
        xr = temperature range, yr = ion fraction range

        for overplotting:
        oplot="ioneqfilename" such as 'mazzotta'
        or if oplot=True or oplot=1 and a widget will come up so that a file can be selected.'''
        if bw:
            linestyle=['k-','k--', 'k-.', 'k:']
        else:
            linestyle=['b-','r--', 'g-.', 'm:']
        #
        if type(stages) == types.NoneType:
            stages=range(1, self.Z+2)
        elif min(stages) < 1 or max(stages) > self.Z+1:
            stages=range(1, self.Z+2)  #  spectroscopic notation
        if type(xr) == types.NoneType:
            xr=[self.Temperature.min(), self.Temperature.max()]
        if type(yr) == types.NoneType:
            yr=[0.01, 1.1]
        xyr=list(xr)
        xyr.extend(list(yr))
        #
        iz=stages[0]
        pl.loglog(self.Temperature, self.Ioneq[iz-1])
        if label:
            idx=self.Ioneq[iz-1] == self.Ioneq[iz-1].max()
            if idx.sum() > 1:
                jdx=np.arange(len(idx))
                idx=jdx[idx].max()
            ann=const.Ionstage[iz-1]
            pl.annotate(ann, [self.Temperature[idx], 0.7*self.Ioneq[iz-1, idx]], ha='center')
        for iz in stages[1:]:
            pl.plot(self.Temperature, self.Ioneq[iz-1], linestyle[0])
            if label:
                idx=self.Ioneq[iz-1] == self.Ioneq[iz-1].max()
                if idx.sum() > 1:
                    jdx=np.arange(len(idx))
                    idx=jdx[idx].mean()
                ann=const.Ionstage[iz-1]
                pl.annotate(ann, [self.Temperature[idx], 0.7*self.Ioneq[iz-1, idx]], ha='center')
        pl.xlabel('Temperature (K)')
        pl.ylabel('Ion Fraction')
        atitle='Chianti Ionization Equilibrium for '+El[self.Z-1].capitalize()
        #
        if oplot != False:
            if type(oplot) == BooleanType:
                result=self.ioneqRead(ioneqname='',default=False)
                if result != False:
                    atitle+='  & '+result['ioneqname'].replace('.ioneq', '')
                    atitle+=' '+linestyle[0]
                    for iz in ions:
                        pl.plot(self.IoneqTemperature, self.IoneqAll[self.Z-1, iz-1],linestyle[0], linestyle[1])
            elif type(oplot) == StringType:
                atitle+='  & '+oplot+' '+linestyle[0]
                atitle+=' '+linestyle[0]
                result=self.ioneqRead(ioneqname=oplot,default=False)
                if result != False:
                    for iz in ions:
                        pl.plot(self.IoneqTemperature, self.IoneqAll[self.Z-1, iz-1],linestyle[0], linestyle[1])
            elif type(oplot) == ListType:
                for iplot in range(len(oplot)):
                    result=self.ioneqRead(ioneqname=oplot[iplot],default=False)
                    if result != False:
                        atitle+='  & '+oplot[iplot]+' '+linestyle[iplot%3]
                        for iz in ions:
                            pl.plot(self.IoneqTemperature, self.IoneqAll[self.Z-1, iz-1], linestyle[iplot%4])
            else:
                print ' oplot not understood ', oplot
        if title:
            pl.title(atitle)
        pl.axis(xyr)
    #
    # -------------------------------------------------------------------------
    #
class spectrum:
    '''Calculate the emission spectrum as a function of temperature and density.

    temperature and density can be arrays but, unless the size of either is one (1),
    the two must have the same size

    the returned spectrum will be convolved with a filter of the specified width on the
    specified wavelength array

    the default filter is gaussianR with a resolving power of 100.  Other filters,
    such as gaussian, box and lorentz, are available in chianti.filters.  When using the box filter,
    the width should equal the wavelength interval to keep the units of the continuum and line
    spectrum the same.

    A selection of ions can be make with ionList containing the names of
    the desired lines in Chianti notation, i.e. C VI = c_6

    a minimum abundance can be specified so that the calculation can be speeded up by excluding
    elements with a low abundance.  Setting doContinuum =0 will skip the continuum calculation.

    Setting em will multiply the spectrum at each temperature by the value of em.

    em [for emission measure], can be a float or an array of the same length as the
    temperature/density.'''
    def __init__(self, temperature, density, wavelength, filter=(chfilters.gaussianR, 1000.),  ionList = None, minAbund=0., doContinuum=1, em = None):
        if type(ionList) == types.NoneType:
            masterlist = util.masterListRead()
        else:
            masterlist = ionList
        self.Temperature = np.asarray(temperature, 'float64')
        nTemp = self.Temperature.size
        self.Density = np.asarray(density, 'float64')
        nDen = self.Density.size
        nTempDen = max([nTemp, nDen])
        if type(em) != types.NoneType:
            if type(em) == types.FloatType:
                if nTempDen > 1:
                    em = np.ones_like(self.Temperature)*em
                    nEm = nTempDen
                else:
                    nEm = 1
            else:
                em = np.asarray(em, 'float64')
                nEm = em.size
                if nEm != nTempDen:
                    print ' the emission measure array must be the same size as the temperature/density array'
                    return
        self.AbundanceName = defaults['abundfile']
        self.AbundanceAll = util.abundanceRead(abundancename = self.AbundanceName)
        ionInfo = util.masterListInfo()
        wavelength = np.asarray(wavelength)
        nWvl = wavelength.size
        self.Wavelength = wavelength
        wvlRange = [wavelength.min(), wavelength.max()]
        print ' wavelength range = ', wvlRange
        #
        freeFree = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        freeBound = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        twoPhoton = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        lineSpectrum = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        #
        #
        for iz in range(31):
            abundance = self.AbundanceAll['abundance'][iz-1]
            if abundance >= minAbund:
                print ' %5i %5s abundance = %10.2e '%(iz, const.El[iz-1],  abundance)
                #
                for ionstage in range(1, iz+2):
                    ionS = util.zion2name(iz, ionstage)
#                   print ' ionS = ', ionS
                    masterListTest = ionS in masterlist
                    masterListInfoTest = ionS in ionInfo.keys()
                    if masterListTest or masterListInfoTest:
                        wvlTestMin = self.Wavelength.min() <= ionInfo[ionS]['wmax']
                        wvlTestMax = self.Wavelength.max() >= ionInfo[ionS]['wmin']
                        ioneqTest = (self.Temperature.max() >= ionInfo[ionS]['tmin']) and (self.Temperature.min() <= ionInfo[ionS]['tmax'])
                    # construct similar test for the dielectronic files
                    ionSd = util.zion2name(iz, ionstage, dielectronic=1)
                    masterListTestD = ionSd in masterlist
                    masterListInfoTestD = ionSd in ionInfo.keys()
                    if masterListTestD or masterListInfoTestD:
                        wvlTestMinD = self.Wavelength.min() <= ionInfo[ionSd]['wmax']
                        wvlTestMaxD = self.Wavelength.max() >= ionInfo[ionSd]['wmin']
                        ioneqTestD = (self.Temperature.max() >= ionInfo[ionSd]['tmin']) and (self.Temperature.min() <=ionInfo[ionSd]['tmax'])
                    ionstageTest = ionstage > 1
                    if ionstageTest and ioneqTest and doContinuum:
                        # ionS is the target ion, cannot be the neutral for the continuum
                        print ' calculating continuum for :  ',  ionS
                        cont = chianti.core.continuum(ionS, temperature)
                        cont.freeFree(wavelength)
    #                   print dir(thisIon)
    #                   print ' wvl = ', thisIon.FreeFree['wvl']
                        if nTempDen ==1:
                            freeFree += cont.FreeFree['rate']
                        else:
                            for iTempDen in range(nTempDen):
                                freeFree[iTempDen] += cont.FreeFree['rate'][iTempDen]
                    #
                        cont.freeBound(wavelength)
                        if 'errorMessage' not in cont.FreeBound.keys():
                            #  an fblvl file exists for this ions
                            if nTempDen == 1:
                                freeBound += cont.FreeBound['rate']
                            else:
                                freeBound[iTempDen] += cont.FreeBound['rate'][iTempDen]
                    if masterListTest and wvlTestMin and wvlTestMax and ioneqTest:
                        print ' calculating spectrum for  :  ', ionS
                        thisIon = chianti.core.ion(ionS, temperature, density)
#                       print ' dir = ', dir(thisIon)
                        thisIon.intensity(wvlRange = wvlRange)
                        # check that there are lines in this wavelength range
                        if 'errorMessage' not in  thisIon.Intensity.keys():
                            thisIon.spectrum(wavelength, filter=filter)
#                           intensity = thisIon.Intensity['intensity']
                            if nTempDen == 1:
                                lineSpectrum += thisIon.Spectrum['intensity']
                            else:
                                for iTempDen in range(nTempDen):
                                    lineSpectrum[iTempDen] += thisIon.Spectrum['intensity'][iTempDen]
                        # get 2 photon emission for H and He sequences
                        if (iz - ionstage) in [0, 1]:
                            thisIon.twoPhoton(wavelength)
                            twoPhoton += thisIon.TwoPhoton['rate']
                    # get dielectronic lines
                    if masterListTestD and wvlTestMinD and wvlTestMaxD and ioneqTestD:
                        print ' calculating spectrum for  :  ', ionSd
                        thisIon = chianti.core.ion(ionSd, temperature, density)
#                       print ' dir = ', dir(thisIon)
                        thisIon.intensity(wvlRange = wvlRange)
                        # check that there are lines in this wavelength range - probably not redundant
                        if 'errorMessage' not in  thisIon.Intensity.keys():
                            thisIon.spectrum(wavelength, filter=filter)
                            if nTempDen == 1:
                                lineSpectrum += thisIon.Spectrum['intensity']
                            else:
                                for iTempDen in range(nTempDen):
                                    lineSpectrum[iTempDen] += thisIon.Spectrum['intensity'][iTempDen]
        self.FreeFree = {'wavelength':wavelength, 'intensity':freeFree.squeeze()}
        self.FreeBound = {'wavelength':wavelength, 'intensity':freeBound.squeeze()}
        self.LineSpectrum = {'wavelength':wavelength, 'intensity':lineSpectrum.squeeze()}
        self.TwoPhoton = {'wavelength':wavelength, 'intensity':twoPhoton.squeeze()}
        #
        total = freeFree + freeBound + lineSpectrum + twoPhoton
        if type(em) != types.NoneType:
            if nEm == 1:
                integrated = total*em
            else:
                integrated = np.zeros_like(wavelength)
                for iTempDen in range(nTempDen):
                    integrated += total[iTempDen]*em[iTempDen]
            self.Spectrum ={'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'em':em}
        else:
            self.Spectrum ={'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1]}
    #
    # -------------------------------------------------------------------------
    #
class ionWeb(ion):
    """
    a class that contains methods to be used for 'Chianti on the Web'
    """
    def gofntSelectLines(self,wvlRange=0, top=10,  saveFile=0):
        """Provide a selection of lines for calculating the 'so-called' G(T) function.

        Given as a function of both temperature and density.

        Only the top( set by 'top') brightest lines are plotted."""
        #
        #self.emiss={"wvl":wvl,"emiss":em,"units":units,"plotLabels":plotLabels}
        #
        #
        doEmiss=False
        try:
            em=self.Emiss
        except:
            doEmiss=True
        #
        #
        if doEmiss:
            # new values of temperature or density
            self.emiss()
            em=self.Emiss
        #
        #
        try:
            ab=self.Abundance
        except:
            self.Abundance = util.abundanceRead()
            ab=self.Abundance
        #
        fontsize=12
        #
        emiss=em["emiss"]
        wvl=em["wvl"]
        temperature=self.Temperature
        density=self.Density
        plotLabels=em["plotLabels"]
        xLabel=plotLabels["xLabel"]
        yLabel=plotLabels["yLabel"]
        #
        # find which lines are in the wavelength range if it is set
        #
        #
        if type(wvlRange) != type(1):
            igvl=util.between(wvl,wvlRange)
        else:
            igvl=range(len(wvl))
        nlines=len(igvl)
        if nlines ==0:
            if chInteractive:
                print ' no lines in selected interval'
            else:
                self.message = ' no lines in selected interval'
            return
        # find the top most intense lines
        #
        if (top > nlines) or (top == 0):
            top=nlines
        maxEmiss=np.zeros(nlines,'Float32')
        for iline in range(nlines):
            maxEmiss[iline]=emiss[igvl[iline]].max()
        for iline in range(nlines):
            if maxEmiss[iline]>=maxEmiss.max():
                maxAll=emiss[igvl[iline]]
                maxIndex = igvl[iline]
#        print ' maxIndex, maxAll = ', maxIndex,  maxAll
        line=range(nlines)
        igvlsort=np.take(igvl,np.argsort(maxEmiss))
        topLines=igvlsort[-top:]
        maxWvl='%5.3f' % wvl[topLines[-1]]
        maxline=topLines[-1]
        #
        # need to make sure there are no negative values before plotting
        good = np.where(emiss > 0.)
        emissMin=emiss[good].min()
        bad=np.where(emiss <= 0.)
        emiss[bad]=emissMin
        #
        topLines=topLines[wvl[topLines].argsort()]
        #
        #
        ntemp=self.Temperature.size
        #
        ndens=self.Density.size
        #
        ylabel = 'Emissivity relative to '+maxWvl
        title = self.Spectroscopic
        #
        #  normally, ionWeb is only using in the non-interactive mode
        if chInteractive:
            pl.ion()
        else:
            pl.ioff()
        #
        #
        if ndens==1 and ntemp==1:
            print ' only a single temperature and density'
            return
        elif ndens == 1:
            xlabel='Temperature (K)'
            ngofnt = temperature.size
            xvalues=temperature
            outTemperature=temperature
            outDensity=np.zeros(ntemp,'Float32')
            outDensity.fill(density)
            desc_str=' at Density = %10.2e' % density
        elif ntemp == 1:
            xvalues=density
            ngofnt = density.size
            outTemperature=np.zeros(ndens,'Float32')
            outTemperature.fill(temperature)
            outDensity=density
            xlabel=r'$\rm{Electron Density (cm}^{-3}\rm{)}$'
            desc_str=' at Temperature = %10.2e' % temperature
        else:
            outTemperature=temperature
            outDensity=density
            xlabel='Temperature (K)'
            xvalues=temperature
            ngofnt = temperature.size
            desc_str=' for variable Density'
            #
        #
        # put all actual plotting here
        #
#        pl.ion()
        pl.figure()
        nxvalues=len(xvalues)
        for iline in range(top):
            tline=topLines[iline]
            pl.loglog(xvalues,emiss[tline]/maxAll)
            skip=2
            start=divmod(iline,nxvalues)[1]
            for ixvalue in range(start,nxvalues,nxvalues/skip):
                pl.text(xvalues[ixvalue],emiss[tline,ixvalue]/maxAll[ixvalue],str(wvl[tline]))
        pl.xlim(xvalues.min(),xvalues.max())
#       yl=pl.ylim()
#       pl.ylim(yl[0],1.2)
        pl.xlabel(xlabel,fontsize=fontsize)
        pl.ylabel(ylabel,fontsize=fontsize)
        pl.title(title+desc_str,fontsize=fontsize)
        if saveFile:
            pl.savefig(saveFile)
        else:
            pl.draw()
        #
#        print ' topLInes = ', wvl[topLines]
        wvlChoices = []
        for one in wvl[topLines]:
            wvlChoices.append('%12.3f'%(one))
        self.wvlChoices = wvlChoices
#        gline = gui.selectorDialog(wvlChoices,label='Select line(s)')
#        gline_idx=gline.selectedIndex
        #
        #
        # -------------------------------------------------------------------------------------
        #
    def gofntShow(self,wvlRange=0,top=10,index=0, saveFile=0):
        """Return a plot of the 'so-called' G(T) function fron the selected lines in index

        Given as a function of both temperature and density.

        Only the top( set by 'top') brightest lines are plotted."""
        #
        #self.emiss={"wvl":wvl,"emiss":em,"units":units,"plotLabels":plotLabels}
        #
        #
        doEmiss=False
        try:
            em=self.Emiss
        except:
            doEmiss=True
        #
        #
        if doEmiss:
            # new values of temperature or density
            self.emiss()
            em=self.Emiss
        #
        #
        try:
            ab=self.Abundance
        except:
            self.Abundance = util.abundanceRead()
            ab=self.Abundance
        #
        fontsize=12
        #
        emiss=em["emiss"]
        wvl=em["wvl"]
        temperature=self.Temperature
        density=self.Density
        plotLabels=em["plotLabels"]
        xLabel=plotLabels["xLabel"]
        yLabel=plotLabels["yLabel"]
        #
        # find which lines are in the wavelength range if it is set
        #
        #
        if type(wvlRange) != type(1):
            igvl=util.between(wvl,wvlRange)
        else:
            igvl=range(len(wvl))
        nlines=len(igvl)
        if nlines ==0:
            print ' no lines in selected interval'
            self.Message = ' no lines in selected interval '
            return
        # find the top most intense lines
        #
        if top > nlines:
            top=nlines
        maxEmiss=np.zeros(nlines,'Float32')
        for iline in range(nlines):
            maxEmiss[iline]=emiss[igvl[iline]].max()
        for iline in range(nlines):
            if maxEmiss[iline]>=maxEmiss.max():
                maxAll=emiss[igvl[iline]]
                maxIndex = igvl[iline]
#        print ' maxIndex, maxAll = ', maxIndex,  maxAll
        line=range(nlines)
        igvlsort=np.take(igvl,np.argsort(maxEmiss))
        topLines=igvlsort[-top:]
        maxWvl='%5.3f' % wvl[topLines[-1]]
        maxline=topLines[-1]
        #
        # need to make sure there are no negative values before plotting
        good = np.where(emiss > 0.)
        emissMin=emiss[good].min()
        bad=np.where(emiss <= 0.)
        emiss[bad]=emissMin
        #
        topLines=topLines[wvl[topLines].argsort()]
        #
        #
        ntemp=self.Temperature.size
        #
        ndens=self.Density.size
        #
        ylabel = ' Gofnt '
        title = self.Spectroscopic
        #
        #
        if ndens==1 and ntemp==1:
            print ' only a single temperature and density'
            return
        elif ndens == 1:
            xlabel='Temperature (K)'
            ngofnt = temperature.size
            xvalues=temperature
            outTemperature=temperature
            outDensity=np.zeros(ntemp,'Float32')
            outDensity.fill(density)
            desc_str=' at Density = %10.2e' % density
        elif ntemp == 1:
            xvalues=density
            ngofnt = density.size
            outTemperature=np.zeros(ndens,'Float32')
            outTemperature.fill(temperature)
            outDensity=density
            xlabel=r'$\rm{Electron Density (cm}^{-3}\rm{)}$'
            desc_str=' at Temperature = %10.2e' % temperature
        else:
            outTemperature=temperature
            outDensity=density
            xlabel='Temperature (K)'
            xvalues=temperature
            ngofnt = temperature.size
            desc_str=' for variable Density'
#            #
#        #
#        # put all actual plotting here
#        #
#        pl.ion()
#        pl.figure()
#        nxvalues=len(xvalues)
#        for iline in range(top):
#            tline=topLines[iline]
#            pl.loglog(xvalues,emiss[tline]/maxAll)
#            skip=2
#            start=divmod(iline,nxvalues)[1]
#            for ixvalue in range(start,nxvalues,nxvalues/skip):
#                pl.text(xvalues[ixvalue],emiss[tline,ixvalue]/maxAll[ixvalue],str(wvl[tline]))
#        pl.xlim(xvalues.min(),xvalues.max())
##       yl=pl.ylim()
##       pl.ylim(yl[0],1.2)
#        pl.xlabel(xlabel,fontsize=fontsize)
#        pl.ylabel(ylabel,fontsize=fontsize)
#        pl.title(title+desc_str,fontsize=fontsize)
#        pl.draw()
#        #
##        print ' topLInes = ', wvl[topLines]
#        wvlChoices = []
#        for one in wvl[topLines]:
#            wvlChoices.append('%12.3f'%(one))
#        gline = gui.selectorDialog(wvlChoices,label='Select line(s)')
#        gline_idx=gline.selectedIndex
#        #
        gline_idx = index
        # for now
        ngofnt = 1
        #
        gAbund=self.Abundance
        #
        try:
            thisIoneq=self.IoneqOne
        except:
            self.ioneqOne()
        #        gioneq=np.where(thisIoneq > 0.)
        #        y2=interpolate.splrep(np.log(self.IoneqAll['ioneqTemperature'][gioneq]),np.log(thisIoneq[gioneq]),s=0)  #allow smoothing,s=0)
        #        gIoneq=interpolate.splev(np.log(temperature),y2)   #,der=0)
        gIoneq=self.IoneqOne/density
        #
        #
        #  ionWeb is normally only used in the non-interative mode
        if chInteractive:
            pl.ion()
        else:
            pl.ioff()
        #
        #
        # plot the desired ratio
        pl.figure()
        g_line= topLines[gline_idx]#  [0]
        #print ' g_line = ',g_line
        #
        gofnt=np.zeros(ngofnt,'float32')
        if ngofnt > 1:
            for aline in g_line:
    #        for aline in gline_idx:
                gofnt += gAbund*gIoneq*emiss[aline].squeeze()
        else:
            gofnt = gAbund*gIoneq*emiss[index].squeeze()

        self.Gofnt={'temperature':outTemperature,'density':outDensity,'gofnt':gofnt}
        #
        pl.loglog(xvalues,gofnt)
        pl.xlim(xvalues.min(),xvalues.max())
        pl.xlabel(xlabel,fontsize=fontsize)
        pl.ylabel('Gofnt',fontsize=fontsize)
        pl.title(title+' '+str(wvl[g_line])+' '+desc_str, fontsize=fontsize)
        if saveFile:
            pl.savefig(saveFile)
        else:
            pl.show()
        #pl.ioff()
        #pl.show()
#        return
    def intensityRatioSelectLines(self,wvlRange=0,top=10,  saveFile=0):
        """Provide a selection of lines for calculating the 'so-called' G(T) function.

        Given as a function of both temperature and density.

        Only the top( set by 'top') brightest lines are plotted."""
        #
        #self.emiss={"wvl":wvl,"emiss":em,"units":units,"plotLabels":plotLabels}
        #
        #
        doEmiss=False
        try:
            em=self.Emiss
        except:
            doEmiss=True
        #
        #
        if doEmiss:
            # new values of temperature or density
            self.emiss()
            em=self.Emiss
        #
        #
        try:
            ab=self.Abundance
        except:
            self.Abundance = util.abundanceRead()
            ab=self.Abundance
        #
        fontsize=12
        #
        emiss=em["emiss"]
        wvl=em["wvl"]
        temperature=self.Temperature
        density=self.Density
        plotLabels=em["plotLabels"]
        xLabel=plotLabels["xLabel"]
        yLabel=plotLabels["yLabel"]
        #
        # find which lines are in the wavelength range if it is set
        #
        #
        if not isinstance(wvlRange, int):
            igvl=util.between(wvl,wvlRange)
        else:
            igvl=range(len(wvl))
        nlines=len(igvl)
        if nlines ==0:
            if chInteractive:
                print ' no lines in selected interval'
            else:
                self.message = ' no lines in selected interval'
            return
        # find the top most intense lines
        #
        if (top > nlines) or (top == 0):
            top=nlines
        maxEmiss=np.zeros(nlines,'Float32')
        for iline in range(nlines):
            maxEmiss[iline]=emiss[igvl[iline]].max()
        for iline in range(nlines):
            if maxEmiss[iline]>=maxEmiss.max():
                maxAll=emiss[igvl[iline]]
                maxIndex = igvl[iline]
#        print ' maxIndex, maxAll = ', maxIndex,  maxAll
        line=range(nlines)
        igvlsort=np.take(igvl,np.argsort(maxEmiss))
        topLines=igvlsort[-top:]
        maxWvl='%5.3f' % wvl[topLines[-1]]
        maxline=topLines[-1]
        #
        # need to make sure there are no negative values before plotting
        good = np.where(emiss > 0.)
        emissMin=emiss[good].min()
        bad=np.where(emiss <= 0.)
        emiss[bad]=emissMin
        #
        topLines=topLines[wvl[topLines].argsort()]
        #
        #
        ntemp=self.Temperature.size
        #
        ndens=self.Density.size
        #
        ylabel = 'Emissivity relative to '+maxWvl
        title = self.Spectroscopic
        #
        #
        if ndens==1 and ntemp==1:
            print ' only a single temperature and density'
            return
        elif ndens == 1:
            xlabel='Temperature (K)'
            ngofnt = temperature.size
            xvalues=temperature
            outTemperature=temperature
            outDensity=np.zeros(ntemp,'Float32')
            outDensity.fill(density)
            desc_str=' at Density = %10.2e' % density
        elif ntemp == 1:
            xvalues=density
            ngofnt = density.size
            outTemperature=np.zeros(ndens,'Float32')
            outTemperature.fill(temperature)
            outDensity=density
            xlabel=r'$\rm{Electron Density (cm}^{-3}\rm{)}$'
            desc_str=' at Temperature = %10.2e' % temperature
        else:
            outTemperature=temperature
            outDensity=density
            xlabel='Temperature (K)'
            xvalues=temperature
            ngofnt = temperature.size
            desc_str=' for variable Density'
            #
        #
        # put all actual plotting here
        #
#        pl.ion()
        #  topLines are sorted by wavelength
        ymax = np.max(1.2*emiss[topLines[0]]/maxAll)
        ymin = ymax
        pl.figure()
        ax = pl.subplot(111)
        nxvalues=len(xvalues)
        for iline in range(top):
            tline=topLines[iline]
            pl.loglog(xvalues,emiss[tline]/maxAll)
            if np.min(emiss[tline]/maxAll) < ymin:
                ymin = np.min(emiss[tline]/maxAll)
            if np.max(emiss[tline]/maxAll) > ymax:
                ymax = np.max(emiss[tline]/maxAll)
            skip=2
            start=divmod(iline,nxvalues)[1]
            for ixvalue in range(start,nxvalues,nxvalues/skip):
                pl.text(xvalues[ixvalue],emiss[tline,ixvalue]/maxAll[ixvalue],str(wvl[tline]))
        pl.xlim(xvalues.min(),xvalues.max())
#        print ' ymin, ymax = ', ymin, ymax
#        pl.ylim(ymin, ymax)
#       yl=pl.ylim()
#       pl.ylim(yl[0],1.2)
        pl.xlabel(xlabel,fontsize=fontsize)
        pl.ylabel(ylabel,fontsize=fontsize)
        if ndens == ntemp and ntemp > 1:
            pl.text(0.07, 0.5,title, horizontalalignment='left', verticalalignment='center', fontsize=fontsize,  transform = ax.transAxes)
            #
            ax2 = pl.twiny()
            xlabelDen=r'Electron Density (cm$^{-3}$)'
            pl.xlabel(xlabelDen, fontsize=fontsize)
            pl.loglog(density,emiss[topLines[top-1]]/maxAll, visible=False)
            ax2.xaxis.tick_top()
            pl.ylim(ymin/1.2, 1.2*ymax)
        else:
            pl.ylim(ymin/1.2, 1.2*ymax)
            pl.title(title+desc_str,fontsize=fontsize)
        if saveFile:
            pl.savefig(saveFile)
        else:
            pl.draw()
        #
#        print ' topLInes = ', wvl[topLines]
        wvlChoices = []
        for one in wvl[topLines]:
            wvlChoices.append('%12.3f'%(one))
        self.wvlChoices = wvlChoices
        self.topLines = topLines
        #
        #   -----------------------------------
        #
    def intensityRatioShow(self,numIdx, denIdx, saveFile=0):
        """Plot the ratio of 2 lines or sums of lines.

        Shown as a function of density and/or temperature."""
        #
        #        self.Emiss={"temperature":temperature,"density":density,"wvl":wvl,"emiss":em,
        #        "plotLabels":plotLabels}
        #
        #
        em = self.Emiss
        #
#        doEmiss=False
#        try:
#            em=self.Emiss
#        except:
#            doEmiss=True
#        #
#        #
#        if doEmiss:
#            # new values of temperature or density
#            self.emiss()
#            em=self.Emiss
        #
        #
        try:
            ab=self.Abundance
        except:
            self.Abundance = util.abundanceRead()
            ab=self.Abundance
        emiss = em['emiss']
        wvl = em["wvl"]
        plotLabels=em["plotLabels"]
        xLabel=plotLabels["xLabel"]
        yLabel=plotLabels["yLabel"]
        #
        # find which lines are in the wavelength range if it is set
        #
        #
#        if not wvlRange:
#            igvl=range(len(wvl))
#        else:
#            igvl=util.between(wvl,wvlRange)
#        nlines=len(igvl)
        #
#        print ' nlines = ',nlines
#        print ' iglv = ',igvl
#        igvl=np.take(igvl,wvl[igvl].argsort())
        # find the top most intense lines
        #
#        if (top > nlines) or (top == 0):
#            top=nlines
#        maxEmiss=np.zeros(nlines,'Float32')
#        for iline in range(nlines):
#            maxEmiss[iline]=emiss[igvl[iline]].max()
#        for iline in range(nlines):
#            if maxEmiss[iline]==maxEmiss.max():
#                maxAll=emiss[igvl[iline]]
#        line=range(nlines)
#        igvlsort=np.take(igvl,np.argsort(maxEmiss))
##        print 'igvlsort = ', igvlsort
#        topLines=igvlsort[-top:]
#        print ' topLines = ', topLines
#        topLines=topLines[wvl[topLines].argsort()]
        topLines = self.topLines
        maxWvl='%5.3f' % wvl[topLines[-1]]
        maxline=topLines[-1]
        #
        #
        #
        # need to make sure there are no negative values before plotting
        good = np.where(emiss > 0.)
        emissMin=emiss[good].min()
        bad=np.where(emiss <= 0.)
        emiss[bad]=emissMin
        #
        #
        ntemp=self.Temperature.size
        #
        ndens=self.Density.size
        #
        ylabel='Emissivity relative to '+maxWvl
        title=self.Spectroscopic
        #
        #
        if ndens==1 and ntemp==1:
            print ' only a single temperature and density'
            return
        elif ndens == 1:
            xlabel='Temperature (K)'
            xvalues=self.Temperature
            outTemperature=self.Temperature
            outDensity=np.zeros(ntemp,'Float32')
            outDensity.fill(self.Density)
            desc_str=' at  Density = %10.2e (cm)$^{-3}$' % self.Density
        elif ntemp == 1:
            xvalues=self.Density
            outTemperature=np.zeros(ndens,'Float32')
            outTemperature.fill(self.Temperature)
            outDensity=self.Density
            xlabel=r'$\rm{Electron Density (cm)^{-3}}$'
            desc_str=' at Temp = %10.2e (K)' % self.Temperature
        else:
            outTemperature=self.Temperature
            outDensity=self.Density
            xlabel='Temperature (K)'
            xvalues=self.Temperature
            desc_str=' for variable Density'
            #
        #
        # put all actual plotting here
        #
        #
        # num_idx and den_idx are tuples
        #
        if np.iterable(numIdx):
            num_idx=numIdx
        else:
            num_idx = [numIdx]
        if len(num_idx) == 0:
            if chInteractive:
                print ' no numerator lines were selected'
            else:
                self.Message = ' no numerator lines were selected'
            return
        #
        if np.iterable(denIdx):
            den_idx=denIdx
        else:
            den_idx = [denIdx]
        #
        if len(den_idx) == 0:
            if chInteractive:
                print ' no denominator lines were selected'
            else:
                self.Message = ' no denominator lines were selected'
            return
        #
#       print ' num_idx = ', num_idx
#       print ' toplines[num_idx] = ', topLines[num_idx]
#       num_line= topLines[num_idx]
#       den_line= topLines[den_idx]
#       #
#       print ' num_line = ', num_line
        numEmiss=np.zeros(len(xvalues),'float32')
        for aline in num_idx:
            numEmiss+=emiss[topLines[aline]]
        #
        denEmiss=np.zeros(len(xvalues),'float32')
        for aline in den_idx:
            denEmiss+=emiss[topLines[aline]]
        fontsize = 12
        #
        # plot the desired ratio
        pl.figure()
        ax = pl.subplot(111)
        pl.loglog(xvalues,numEmiss/denEmiss)
        pl.xlim(xvalues.min(),xvalues.max())
        pl.xlabel(xlabel,fontsize=fontsize)
        pl.ylabel('Ratio ('+self.Defaults['flux']+')',fontsize=fontsize)
        desc = title + ':'
        for aline in num_idx:
            desc += ' ' + str(wvl[topLines[aline]])
        desc +=' / '
        for aline in den_idx:
            desc += ' ' + str(wvl[topLines[aline]])
        if ndens == ntemp and ntemp > 1:
            pl.text(0.07, 0.5,desc, horizontalalignment='left', verticalalignment='center', fontsize=fontsize,  transform = ax.transAxes)
            #
            ax2 = pl.twiny()
            xlabelDen=r'Electron Density (cm$^{-3}$)'
            pl.xlabel(xlabelDen, fontsize=fontsize)
            pl.loglog(outDensity,numEmiss/denEmiss, visible=False)
            ax2.xaxis.tick_top()
        else:
#            pl.ylim(ymin, ymax)
            pl.title(desc,fontsize=fontsize)
#       desc=title+' '+str(wvl[num_line])+' / '+str(wvl[den_line])+' '+desc_str
#        pl.title(desc, fontsize=fontsize)
#       pl.title(title+' '+str(wvl[num_line])+' / '+str(wvl[den_line])+' '+desc_str,fontsize=fontsize)
#        pl.draw()
#        pl.ioff()
#        pl.show()
        #
        if saveFile:
            pl.savefig(saveFile)
        intensityRatioFileName=self.IonStr
        for aline in num_idx:
            intensityRatioFileName+= '_%3i'%(wvl[topLines[aline]])
        intensityRatioFileName+='_2'
        for aline in den_idx:
            intensityRatioFileName+= '_%3i'%(wvl[topLines[aline]])
        intensityRatioFileName+='.rat'
        self.IntensityRatio={'ratio':numEmiss/denEmiss,'desc':desc,
                'temperature':outTemperature,'density':outDensity,'filename':intensityRatioFileName}
        #
        # -------------------------------------------------------------------------------------
        #


