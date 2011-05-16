import types
from datetime import datetime
import numpy as np
np.seterr(over='ignore')
import chianti.core
import chianti.data as chdata
import chianti.constants as const
import chianti.filters as chfilters
import chianti.util as util
#
defaults = chdata.Defaults
chInteractive = chdata.chInteractive
if chInteractive:
    import pylab as pl
else:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as pl

class radloss:
    '''Calculate the emission spectrum as a function of temperature and density.

    includes elemental abundances or ionization equilibria

    temperature and density can be arrays but, unless the size of either is one (1),
    the two must have the same size

    the returned spectrum will be convolved with a filter of the specified width on the
    specified wavelength array

    the default filter is gaussianR with a resolving power of 1000.  Other filters,
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
    def __init__(self, temperature, density, ionList = 0, minAbund=0, doContinuum=1, verbose=0, allLines=1):
        t1 = datetime.now()
        masterlist = chdata.MasterList
        # use the ionList but make sure the ions are in the database
        if ionList:
            alist=[]
            for one in ionList:
                if masterlist.count(one):
                    alist.append(one)
                else:
                    if chInteractive and verbose:
                        pstring = ' %s not in CHIANTI database'%(one)
                        print('')
            masterlist = alist
        self.Defaults=defaults
        self.Temperature = np.asarray(temperature, 'float64')
        nTemp = self.Temperature.size
        self.Density = np.asarray(density, 'float64')
        nDen = self.Density.size
        nTempDen = max([nTemp, nDen])
        self.AbundanceName = defaults['abundfile']
        self.AbundanceAll = chdata.AbundanceAll
        abundAll = self.AbundanceAll['abundance']
        nonzed = abundAll > 0.
        minAbundAll = abundAll[nonzed].min()
        if minAbund < minAbundAll:
            minAbund = minAbundAll
        self.minAbund = minAbund
        ionInfo = util.masterListInfo()
        #
        freeFreeLoss = np.zeros((nTempDen), 'float64').squeeze()
        freeBoundLoss = np.zeros((nTempDen), 'float64').squeeze()
        twoPhotonLoss = np.zeros((nTempDen), 'float64').squeeze()
        boundBoundLoss = np.zeros((nTempDen), 'float64').squeeze()
        #
        #
        for iz in range(31):
            abundance = self.AbundanceAll['abundance'][iz-1]
            if abundance >= minAbund:
                if chInteractive:
                    print ' %5i %5s abundance = %10.2e '%(iz, const.El[iz-1],  abundance)
                #
                for ionstage in range(1, iz+2):
                    ionS = util.zion2name(iz, ionstage)
#                   print ' ionS = ', ionS
                    masterListTest = ionS in masterlist
                    masterListInfoTest = ionS in ionInfo.keys()
                    if masterListTest or masterListInfoTest:
                        ioneqTest = (self.Temperature.max() >= ionInfo[ionS]['tmin']) and (self.Temperature.min() <= ionInfo[ionS]['tmax'])
                    # construct similar test for the dielectronic files
                    ionSd = util.zion2name(iz, ionstage, dielectronic=1)
                    masterListTestD = ionSd in masterlist
                    masterListInfoTestD = ionSd in ionInfo.keys()
                    if masterListTestD or masterListInfoTestD:
                        ioneqTestD = (self.Temperature.max() >= ionInfo[ionSd]['tmin']) and (self.Temperature.min() <=ionInfo[ionSd]['tmax'])
                    ionstageTest = ionstage > 1
                    if ionstageTest and ioneqTest and doContinuum:
                        # ionS is the target ion, cannot be the neutral for the continuum
                        if chInteractive:
                            print ' calculating continuum for :  ',  ionS
                        cont = chianti.core.continuum(ionS, temperature)
                        cont.freeFreeLoss()
    #                   print dir(thisIon)
    #                   print ' wvl = ', thisIon.FreeFree['wvl']
#                        if nTempDen ==1:
                        freeFreeLoss += cont.FreeFreeLoss['rate']
#                        else:
#                            for iTempDen in range(nTempDen):
#                                freeFreeLoss[iTempDen] += cont.FreeFreeLoss['rate'][iTempDen]
                    #
                        cont.freeBoundLoss()
                        if 'errorMessage' not in cont.FreeBoundLoss.keys():
                            #  an fblvl file exists for this ions
#                            if nTempDen == 1:
                            freeBoundLoss += cont.FreeBoundLoss['rate']
#                            else:
#                                freeBound[iTempDen] += cont.FreeBound['rate'][iTempDen]
                    if masterListTest and ioneqTest:
                        if chInteractive:
                            print ' calculating spectrum for  :  ', ionS
                        thisIon = chianti.core.ion(ionS, temperature, density)
#                       print ' dir = ', dir(thisIon)
#                        thisIon.emiss(wvlRange = wvlRange, allLines=allLines)
                        thisIon.boundBoundLoss( allLines=allLines)
                        # check that there are lines in this wavelength range
                        if 'errorMessage' not in  thisIon.BoundBoundLoss.keys():
                            thisIon.boundBoundLoss()
#                           intensity = thisIon.Intensity['intensity']
#                            if nTempDen == 1:
                            boundBoundLoss += thisIon.BoundBoundLoss['rate']
#                            else:
#                                for iTempDen in range(nTempDen):
#                                    lineSpectrum[iTempDen] += thisIon.Spectrum['intensity'][iTempDen]
                        # get 2 photon emission for H and He sequences
                        if (iz - ionstage) in [0, 1]:
                            thisIon.twoPhotonLoss()
                            twoPhotonLoss += thisIon.TwoPhotonLoss['rate']
                    # get dielectronic lines
                    if masterListTestD and ioneqTestD:
                        print ' calculating spectrum for  :  ', ionSd
                        thisIon = chianti.core.ion(ionSd, temperature, density)
#                       print ' dir = ', dir(thisIon)
#                       have to do all lines for the dielectronic satellites
#                        thisIon.emiss(allLines=1)
                        thisIon.intensity(allLines=allLines)
                        # check that there are lines in this wavelength range - probably not redundant
                        if 'errorMessage' not in  thisIon.Intensity.keys():
                            thisIon.boundBoundLoss()
#                            if nTempDen == 1:
                            boundBoundLoss += thisIon.BoundBoundLoss['rate']
#                            else:
#                                for iTempDen in range(nTempDen):
#                                    lineSpectrum[iTempDen] += thisIon.Spectrum['intensity'][iTempDen]
        self.FreeFreeLoss = freeFreeLoss
        self.FreeBoundLoss = freeBoundLoss
        self.LineSpectrumLoss = boundBoundLoss
        self.TwoPhotonLoss = twoPhotonLoss
        #
        total = freeFreeLoss + freeBoundLoss + boundBoundLoss + twoPhotonLoss
        t2 = datetime.now()
        dt=t2-t1
        if chInteractive:
            print ' elapsed seconds = ', dt.seconds
        self.RadLoss ={'rate':total, 'temperature':self.Temperature, 'density':self.Density}
    #
    # -------------------------------------------------------------------------
    #
    def lineSpectrumPlot(self, saveFile=0, plotContinuum=0, linLog = 'lin'):
        ''' to plot the spectrum as a function of wavelength'''
        # must follow setting top
        #
        pl.figure()
        ylabel = 'Intensity'
        #
        xlabel = 'Wavelength ('+self.Defaults['wavelength'] +')'
        #
#        ymin = 10.**(np.log10(emiss.min()).round(0))
        #
        if chInteractive:
            pl.ion()
        else:
            pl.ioff()
        #
        pl.plot(self.LineSpectrum['wavelength'], self.LineSpectrum['intensity'])
        pl.xlabel(xlabel)
        pl.ylabel(ylabel)
        if saveFile:
            pl.savefig(saveFile)
