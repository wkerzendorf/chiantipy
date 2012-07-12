import os
import types
import numpy as np
from scipy import interpolate
import time
#
import chianti.data as chdata
import chianti.sources as sources
#chInteractive = chdata.chInteractive
import pylab as pl
#if chInteractive:
#    import pylab as pl
#else:
##    import matplotlib
##    matplotlib.use('Agg')
#    import matplotlib.pyplot as pl
    #  backend is set in matplotlibrc
#    pl.rcParams['backend'] = 'Agg'
#
#if chInteractive:
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
elif pl.rcParams['backend'].lower() == 'macosx':
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
#    #
import chianti.filters as chfilters
import chianti.util as util
import chianti.constants as const
#import chianti.pl as pl
#
#try:
#    chInteractive = int(os.environ['CHIANTIPY_INTERACTIVE'])
#except:
#    chInteractive = 1
#
xuvtop = chdata.xuvtop
#
#ip = chianti.Ip
#MasterList = chianti.MasterList
#Defaults = chianti.Defaults
#AbundanceAll = chianti.AbundanceAll
#IoneqAll = chianti.IoneqAll
class ion:
    '''The top level class for performing spectral calculations for an ion in the CHIANTI database.

    ionStr is a string corresponding such as 'c_5' that corresponds to the C V ion.
    temperature in Kelvin
    eDensity in cm^-3
    radTemperature, the radiation black-body temperature in Kelvin
    rPlot, the distance from the center of the star in stellar radii
    '''
    def __init__(self, ionStr, temperature=None, eDensity=None, pDensity='default', radTemperature=0,  rStar=0, verbose=0, setup=True,  **kwargs):
        '''The top level class for performing spectral calculations for an ion in the CHIANTI database.

        ionStr is a string corresponding such as 'c_5' that corresponds to the C V ion.
        temperature in Kelvin
        eDensity in cm^-3
        radTemperature, the radiation black-body temperature in Kelvin
        rPlot, the distance from the center of the star in stellar radii
        '''
        #
        #
#        self.__version__ = chianti.__version__
        self.IonStr=ionStr
        self.Defaults=chdata.Defaults
        self.AbundanceName = self.Defaults['abundfile']
        self.IoneqName = self.Defaults['ioneqfile']
        MasterList = chdata.MasterList
        #
        self.Z=util.convertName(ionStr)['Z']
        self.Ion=util.convertName(ionStr)['Ion']
        self.Dielectronic=util.convertName(ionStr)['Dielectronic']
        self.Spectroscopic=util.zion2spectroscopic(self.Z,self.Ion)
        self.FileName=util.zion2filename(self.Z, self.Ion,dielectronic=self.Dielectronic )
        #
        self.RadTemperature = radTemperature
        self.RStar = rStar
        #
        #  ip in eV, but don't read for bare ions
        if self.Ion <= self.Z:
            self.Ip=chdata.Ip[self.Z-1, self.Ion-1-self.Dielectronic]
            if self.Dielectronic:
                self.UpperIp=chdata.Ip[self.Z-1, self.Ion-1]
        #
        if type(temperature) != types.NoneType:
            self.Temperature = np.array(temperature,'float64')
        #
        #
#        self.AbundanceAll = util.abundanceRead(abundancename = self.AbundanceName)
        self.AbundanceAll = chdata.AbundanceAll
        self.Abundance = self.AbundanceAll['abundance'][self.Z-1]
        #
#        self.IoneqAll = util.ioneqRead(ioneqname = self.Defaults['ioneqfile'])
#        if chInteractive:
        self.IoneqAll = chdata.IoneqAll
        self.ioneqOne()
        #
        #  this needs to go after setting temperature and reading ionization equilibria
        if pDensity == 'default':
            self.p2eRatio()
        #
        if type(eDensity) != types.NoneType:
            self.EDensity = np.asarray(eDensity,'float64')
            ndens = self.EDensity.size
            ntemp = self.Temperature.size
            tst1 = ndens == ntemp
            tst2 = ntemp > 1
            tst3 = ndens > 1
            #
            if pDensity == 'default':
                if tst1 and tst2 and tst3:
                    self.PDensity = np.zeros((ntemp), 'float64')
                    for itemp in range(ntemp):
                        self.PDensity[itemp] = self.ProtonDensityRatio[itemp]*self.EDensity[itemp]
                elif tst2 and tst3 and not tst1:
                    print ' if both temperature and eDensity are arrays, they must be of the same size'
                    return
                else:
                    self.PDensity = self.ProtonDensityRatio*self.EDensity
            else:
                self.PDensity = pDensity
        if setup:
            #
            # read in all data if in masterlist
            #  if not, there should still be ionization and recombination rates
            #
            if self.IonStr in MasterList:
                self.Elvlc = util.elvlcRead(self.IonStr, verbose=verbose)
                self.Wgfa = util.wgfaRead(self.IonStr)
                self.Nwgfa=len(self.Wgfa['lvl1'])
                nlvlWgfa = max(self.Wgfa['lvl2'])
                nlvlList =[nlvlWgfa]
                fileName = util.ion2filename(self.IonStr)
#                print 'fileName = ', fileName
                splupsfile = fileName + '.splups'
                if os.path.isfile(splupsfile):
                    # happens the case of fe_3 and prob. a few others
                    self.Splups = util.splupsRead(self.IonStr)
                    self.Nsplups=len(self.Splups['lvl1'])
                    nlvlSplups = max(self.Splups['lvl2'])
                    nlvlList.append(nlvlSplups)
                else:
                    self.Nsplups = 0
                    nlvlSplups = 0
##                self.Nlvls = nlvlElvlc
                #
                file = fileName +'.cilvl'
                if os.path.isfile(file):
                    self.Cilvl = util.cireclvlRead(self.IonStr, cilvl=1)
                    self.Ncilvl=len(self.Cilvl['lvl1'])
                    nlvlCilvl = max(self.Cilvl['lvl2'])
                    nlvlList.append(nlvlCilvl)
                else:
                    self.Ncilvl = 0
                #  .reclvl file may not exist
                reclvlfile = fileName +'.reclvl'
                if os.path.isfile(reclvlfile):
                    self.Reclvl = util.cireclvlRead(self.IonStr, reclvl=1)
                    self.Nreclvl = len(self.Reclvl['lvl1'])
                    nlvlReclvl = max(self.Reclvl['lvl2'])
                    nlvlList.append(nlvlReclvl)
                else:
                    self.Nreclvl = 0
                #  .dielsplups file may not exist
                dielsplupsfile = fileName +'.dielsplups'
                if os.path.isfile(dielsplupsfile):
                    self.DielSplups = util.splupsRead(self.IonStr, diel=1)
                    self.Ndielsplups=len(self.DielSplups["lvl1"])
                    nlvlDielSplups = max(self.DielSplups['lvl2'])
                    nlvlList.append(nlvlDielSplups)
                else:
                    self.Ndielsplups = 0
                #
                #  psplups file may not exist
                psplupsfile = fileName +'.psplups'
                if os.path.isfile(psplupsfile):
                    self.Psplups = util.splupsRead(self.IonStr, prot=True)
                    self.Npsplups=len(self.Psplups["lvl1"])
                else:
                    self.Npsplups = 0
                #
                drparamsFile = fileName +'.drparams'
                if os.path.isfile(drparamsFile):
                    self.DrParams = util.drRead(self.IonStr)
                #
                rrparamsFile = fileName +'.rrparams'
                if os.path.isfile(rrparamsFile):
                    self.RrParams = util.rrRead(self.IonStr)

                #  not needed for ion, only phion
#                photoxfile = util.ion2filename(self.IonStr)+'.photox'
#                if os.path.isfile(photoxfile):
#                    self.Photox = util.photoxRead(self.IonStr)
                #
                # need to determine the number of levels that can be populated
                nlvlElvlc = len(self.Elvlc['lvl'])
#                print ' nlvlElvlc = ', nlvlElvlc
#                print ' other nlvls = ',  nlvlList
#                nlvlWgfa = max(self.Wgfa['lvl2'])
                #  elvlc file can have more levels than the rate level files
                self.Nlvls = min([nlvlElvlc, max(nlvlList)])
        #
        # ------------------------------------------------------------------------------
        #
    def diCross(self, energy=None, verbose=False):
        '''Calculate the direct ionization cross section.

        Given as a function of the incident electron energy in eV, puts values into DiCross'''
        iso=self.Z - self.Ion + 1
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
                    y2=interpolate.splrep(self.DiParams['xsplom'][ifac], self.DiParams['ysplom'][ifac], s=0)
                    btcross=interpolate.splev(btenergy, y2, der=0)
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
    def diRate(self):
        '''Calculate the direct ionization rate coefficient as a function of temperature (K).'''
        if hasattr(self, 'DiParams'):
            DiParams = self.DiParams
        else:
            DiParams = util.diRead(self.IonStr)
        #
        if hasattr(self, 'Temperature'):
            temperature=self.Temperature
        else:
            print ' temperature is not defined'
            return
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
    def eaDescale(self):
        """Calculates the effective collision strengths (upsilon) for excitation-autoionization as a function of temperature."""
        #
        #  xt=kt/de
        #
        #  need to make sure elvl is >0, except for ground level
        #
        if hasattr(self, 'EaParams'):
            eaparams=self.EaParams
        else:
            self.EaParams = util.eaRead(self.IonStr)
            eaparams=self.EaParams
        #
        if hasattr(self, 'Temperature'):
            temperature=self.Temperature
        else:
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
    def eaRate(self):
        '''Calculate the excitation-autoionization rate coefficient.'''
        # get neaev from diparams file
        #
        if hasattr(self, 'DiParams'):
            diparams=self.DiParams
        else:
            self.DiParams = util.diRead(self.IonStr)
        #
        if self.DiParams['info']['neaev'] == 0:
#            print ' no EA rates'
            return
        else:
            if hasattr(self, 'Temperature'):
                temperature=self.Temperature
            else:
                bte=0.1*np.arange(10)
                bte[0]=0.01
                dum=np.ones(10, 'float32')
                [temperature, dum]=util.descale_bt(bte, dum, self.EaParams['cups'][0], self.DiParams['de'][0])
                self.Temperature=temperature
            if hasattr(self, 'EaParams'):
                eaparams=self.EaParams
            else:
                self.eaParams = util.eaRead(self.IonStr)
                self.eaDescale()
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
    def ionizRate(self):
        '''Provides the total ionization rate.

        Calls diRate and eaRate.'''
        if self.Z < self.Ion:
#            print ' this is a bare nucleus and has no ionization rate'
            self.IonizRate = {'rate':np.zeros_like(self.Temperature), 'temperature':self.Temperature}
            return
        self.diRate()
        self.eaRate()
        if self.DiParams['info']['neaev'] == 0:
            ionizrate=self.DiRate['rate']
        else:
            ionizrate=self.DiRate['rate']+self.EaRate['rate']
        self.IonizRate = {'rate':ionizrate, 'temperature':self.DiRate['temperature']}
        return
        #
        # -------------------------------------------------------------------------------------
        #
        #
        # -------------------------------------------------------------------------------------
        #
    def rrRate(self):
        '''Provide the radiative recombination rate coefficient as a function of temperature (K).'''
        if hasattr(self, 'Temperature'):
            temperature=self.Temperature
        else:
            print ' temperature is not defined'
            return
        rrparamsfile = util.ion2filename(self.IonStr) + '.rrparams'
        if hasattr(self, 'RrParams'):
            rrparams=self.RrParams
        elif os.path.isfile(rrparamsfile):
            self.RrParams = util.rrRead(self.IonStr)
            rrparams=self.RrParams
        else:
            self.RrRate={'temperature':temperature, 'rate':np.zeros_like(temperature)}
            return
        #
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
        else:
            self.RrRate={'temperature':temperature, 'rate':np.zeros_like(temperature)}

        #
        # -------------------------------------------------------------------------------------
        #
    def drRate(self):
        '''Provide the dielectronic recombination rate coefficient as a function of temperature (K).
        '''
        #
        if hasattr(self, 'Temperature'):
            temperature=self.Temperature
        else:
            print ' temperature is not defined'
            return {}
        drparamsfile = util.ion2filename(self.IonStr) + '.drparams'
        if hasattr(self, 'DrParams'):
            drparams=self.DrParams
        elif os.path.isfile(drparamsfile):
            self.DrParams = util.drRead(self.IonStr)
            drparams=self.DrParams
        else:
            self.DrRate = {'rate':np.zeros_like(temperature), 'temperature':temperature}
            return
        #
        if drparams['drtype'] == 1:
            # badnell type
            drenergy=drparams['eparams']
            drcoef=drparams['cparams']
            gcoef = drenergy > 0.
            ncoef=gcoef.sum()
#            print ' ncoef = ', gcoef.sum()
            rate=np.zeros(temperature.size, 'float64')
            for icoef in range(ncoef):
                rate += drcoef[icoef]*np.exp(-drenergy[icoef]/temperature)
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
    def cireclvlDescale(self, lvlType):
        '''Interpolate and extrapolate cilvl and reclvl rates.
        type must be either 'reclvl', 'cilvl' or 'rrlvl'
        Used in level population calculations.
        '''
        if hasattr(self, 'Temperature'):
            temperature=self.Temperature
        else:
            print ' temperature is not defined'
            return
        lvlfile = util.ion2filename(self.IonStr)+'.' + lvlType
        if lvlType == 'reclvl':
            if hasattr(self, 'Reclvl'):
                lvl = self.Reclvl
            elif os.path.isfile(lvlfile):
#           print ' reading reclvl file'
                self.Reclvl = util.cireclvlRead(self.IonStr, '.'+lvlType)
                lvl = self.Reclvl
            else:
                self.ReclvlRate = {'rate':zeros_like(temperature)}
                return
        elif lvlType == 'cilvl':
            if hasattr(self, 'Cilvl'):
                lvl = self.Cilvl
            elif os.path.isfile(lvlfile):
#           print ' reading reclvl file'
                self.Cilvl = util.cireclvlRead(self.IonStr, '.'+lvlType)
                lvl = self.Cilvl
            else:
                self.CilvlRate = {'rate':zeros_like(temperature)}
                return
        elif lvlType == 'rrlvl':
            if hasattr(self, 'Rrlvl'):
                lvl = self.Rrlvl
            elif os.path.isfile(lvlfile):
#           print ' reading reclvl file'
                self.Rrlvl = util.cireclvlRead(self.IonStr, '.'+lvlType)
                lvl = self.Rrlvl
            else:
                self.RrlvlRate = {'rate':zeros_like(temperature)}
                return
        #
        #  the rates and temperatures in reclvl are not all the same
        #
        ntemp = temperature.size
        nlvl = len(lvl['lvl1'])
        if ntemp == 1:
            rate = np.zeros(( nlvl), 'float64')
            # previous takes care of temperatures below reclvl['temperature'].min()
            if temperature > lvl['temperature'].max():
                # extrapolate as 1/temperature
                for itrans in range(nlvl):
#                   lvl2 = self.Reclvl['lvl2'][itrans]
                    lvlTemp = lvl['ntemp'][itrans]
                    rate[itrans] = lvl['rate'][itrans,-1]*(lvl['temperature'][itrans, -1]/temperature)
            elif temperature < lvl['temperature'].min():
                # rate is already set to zero
                pass
            else:
                for itrans in range(nlvl):
                    lvl2 = lvl['lvl2'][itrans]
                    nrecTemp = lvl['ntemp'][itrans]
                    y2 = interpolate.splrep(np.log(lvl['temperature'][itrans]), np.log(lvl['rate'][itrans]))
                    cirec = np.exp(interpolate.splev(np.log(temperature),y2))
                    rate[itrans] = cirec.squeeze()
        else:
            # ntemp > 1
            rate = np.zeros(( nlvl, temperature.size), 'float64')
            #
            for itrans in range(nlvl):
                lvl2 = lvl['lvl2'][itrans]
                nTemp = lvl['ntemp'][itrans]
                y2 = interpolate.splrep(np.log(lvl['temperature'][itrans, :nTemp]), np.log(lvl['rate'][itrans, :nTemp]))
                goodLow = temperature < lvl['temperature'][itrans].min()
                if goodLow.sum() >0:
#                   print ' number of low temperatures  = ', goodLow.sum()
                    lowT = temperature[goodLow]
                good1 = temperature >= lvl['temperature'][itrans].min()
                good2 = temperature <= lvl['temperature'][itrans].max()
                realgood = np.logical_and(good1,good2)
                if realgood.sum() > 0:
#                   print ' number of mid temperatures  = ', realgood.sum()
                    midT = temperature[realgood]
                goodHigh = temperature > lvl['temperature'][itrans].max()
                if goodHigh.sum() > 0:
#                   print ' number of high temperatures  = ', goodHigh.sum()
                    highT = temperature[goodHigh]
                lvl2 = lvl['lvl2'][itrans]
                nTemp = lvl['ntemp'][itrans]
                newRate = np.zeros(ntemp, 'float64')
                index = 0
                if goodLow.sum() == 1:
#                    lowRec = np.exp(interpolate.splev(np.log(lowT),y2))
#                    newRec[index] = lowRec
                    newRate[index] = 0.
                    index += 1
                elif goodLow.sum() > 1:
#                    lowRec = np.exp(interpolate.splev(np.log(lowT),y2))
                    for idx in range(goodLow.sum()):
#                        newRec[index] = lowRec[idx]
                        newRate[index] = 0.
                        index += 1
                if realgood.sum() == 1:
                    midRec = np.exp(interpolate.splev(np.log(midT),y2))
                    newRAte[index] = midRec
                    index += 1
                elif realgood.sum() > 1:
                    midRec = np.exp(interpolate.splev(np.log(midT),y2))
                    for idx in range(realgood.sum()):
                        newRate[index] = midRec[idx]
                        index += 1
                if goodHigh.sum() == 1:
#                    highRec = np.exp(interpolate.splev(np.log(highT),y2))
#                    newRec[index] = highRec
                    newRate[index] = 0.
                    index += 1
                elif goodHigh.sum() > 1:
#                    highRec = np.exp(interpolate.splev(np.log(highT),y2))
                    for idx in range(goodHigh.sum()):
#                       print ' index, idx = ', index,  idx
#                        newRec[index] = highRec[idx]
                        newRate[index] = 0.
                        index += 1
                rate[itrans] = newRate
        if lvlType == 'reclvl':
            self.ReclvlRate = {'rate':rate, 'lvl1':lvl['lvl1'], 'lvl2':lvl['lvl2'], 'temperature':temperature}
        elif lvlType == 'cilvl':
            self.CilvlRate = {'rate':rate, 'lvl1':lvl['lvl1'], 'lvl2':lvl['lvl2'], 'temperature':temperature}
        elif lvlType == 'rrlvl':
            self.RrlvlRate = {'rate':rate, 'lvl1':lvl['lvl1'], 'lvl2':lvl['lvl2'], 'temperature':temperature}

        #
        # -------------------------------------------------------------------------------------
        #
#    def reclvlDescale(self):
#        '''
#        Interpolate and extrapolate reclvl rates.
#
#        Used in level population calculations.
#        I think this is made redundant by cireclvlDescale
#        '''
#        if hasattr(self, 'Temperature'):
#            temperature=self.Temperature
#        else:
#            print ' temperature is not defined'
#            self.ReclvlRate = None
#            return
#        reclvlfile = util.ion2filename(self.IonStr)+'.reclvl'
#        if hasattr(self, 'Reclvl'):
#            reclvl = self.Reclvl
#        elif os.path.isfile(reclvlfile):
##           print ' reading reclvl file'
#            reclvl = util.cireclvlRead(self.IonStr, 'reclvl')
#        else:
#            self.ReclvlRate = {'rate':zeros_like(temperature)}
#            return
#        #
#        #  the rates and temperatures in reclvl are not all the same
#        #
#        ntemp = temperature.size
#        if ntemp == 1:
#            recRate = np.zeros(( len(reclvl['lvl1'])), 'float64')
#            # previous takes care of temperatures below reclvl['temperature'].min()
#            if temperature > reclvl['temperature'].max():
#                # extrapolate as 1/temperature
#                for itrans in range(len(reclvl['lvl1'])):
##                   lvl2 = self.Reclvl['lvl2'][itrans]
#                    nrecTemp = reclvl['ntemp'][itrans]
#                    recRate[itrans] = self.Reclvl['rate'][itrans,nrecTemp-1]*(reclvl['temperature'][itrans, nrecTemp-1]/temperature)
#            else:
#                for itrans in range(len(self.Reclvl['lvl1'])):
#                    lvl2 = self.Reclvl['lvl2'][itrans]
#                    nrecTemp = self.Reclvl['ntemp'][itrans]
#                    y2 = interpolate.splrep(np.log(self.Reclvl['temperature'][itrans, :nrecTemp]), np.log(self.Reclvl['rate'][itrans, :nrecTemp]))
#                    rec = np.exp(interpolate.splev(np.log(temperature),y2))
#                    recRate[itrans] = rec.squeeze()
#        else:
#            # ntemp > 1
#            recRate = np.zeros(( len(reclvl['lvl1']), temperature.size), 'float64')
#            #
#            for itrans in range(len(reclvl['lvl1'])):
#                lvl2 = self.Reclvl['lvl2'][itrans]
#                nrecTemp = self.Reclvl['ntemp'][itrans]
#                y2 = interpolate.splrep(np.log(self.Reclvl['temperature'][itrans, :nrecTemp]), np.log(self.Reclvl['rate'][itrans, :nrecTemp]))
#                goodLow = temperature < self.Reclvl['temperature'][itrans].min()
#                if goodLow.sum() >0:
##                   print ' number of low temperatures  = ', goodLow.sum()
#                    lowT = temperature[goodLow]
#                good1 = temperature >= self.Reclvl['temperature'][itrans].min()
#                good2 = temperature <= self.Reclvl['temperature'][itrans].max()
#                realgood = np.logical_and(good1,good2)
#                if realgood.sum() > 0:
##                   print ' number of mid temperatures  = ', realgood.sum()
#                    midT = temperature[realgood]
#                goodHigh = temperature > self.Reclvl['temperature'][itrans].max()
#                if goodHigh.sum() > 0:
##                   print ' number of high temperatures  = ', goodHigh.sum()
#                    highT = temperature[goodHigh]
#                lvl2 = self.Reclvl['lvl2'][itrans]
#                nrecTemp = self.Reclvl['ntemp'][itrans]
#                newRec = np.zeros(ntemp, 'float64')
#                index = 0
#                if goodLow.sum() == 1:
##                    lowRec = np.exp(interpolate.splev(np.log(lowT),y2))
##                    newRec[index] = lowRec
#                    newRec[index] = 0.
#                    index += 1
#                elif goodLow.sum() > 1:
##                    lowRec = np.exp(interpolate.splev(np.log(lowT),y2))
#                    for idx in range(goodLow.sum()):
##                        newRec[index] = lowRec[idx]
#                        newRec[index] = 0.
#                        index += 1
#                if realgood.sum() == 1:
#                    midRec = np.exp(interpolate.splev(np.log(midT),y2))
#                    newRec[index] = midRec
#                    index += 1
#                elif realgood.sum() > 1:
#                    midRec = np.exp(interpolate.splev(np.log(midT),y2))
#                    for idx in range(realgood.sum()):
#                        newRec[index] = midRec[idx]
#                        index += 1
#                if goodHigh.sum() == 1:
##                    highRec = np.exp(interpolate.splev(np.log(highT),y2))
##                    newRec[index] = highRec
#                    newRec[index] = 0.
#                    index += 1
#                elif goodHigh.sum() > 1:
##                    highRec = np.exp(interpolate.splev(np.log(highT),y2))
#                    for idx in range(goodHigh.sum()):
##                       print ' index, idx = ', index,  idx
##                        newRec[index] = highRec[idx]
#                        newRec[index] = 0.
#                        index += 1
#                recRate[itrans] = newRec
#        self.ReclvlRate = {'rate':recRate, 'lvl2':reclvl['lvl2'], 'temperature':temperature}
        #
        # -------------------------------------------------------------------------------------
        #
    def drRateLvl(self, verbose=0):
        ''' to calculate the level resolved dielectronic rate from the higher ionization stage to
        the ion of interest
        rates are determined from autoionizing A-values
        the dictionary self.DrRateLvl contains
            rate = the dielectronic rate into an autoionizing level
            effRate = the dielectronic rate into an autoionizing level mutilplied by the branching
                ratio for a stabilizing transition
            totalRate = the sum of all the effRates
        '''
        if not hasattr(self, 'Higher'):
            nameStuff = util.convertName(self.IonStr)
            z = nameStuff['Z']
            stage = nameStuff['Ion']
            higherStr = util.zion2name(z, stage+1)
            self.Higher = ion(higherStr, self.Temperature, self.EDensity)
        #
        # (4 pi a0^2)^(3/2) = 6.6011e-24 (Badnell et al, 2003, A&A 406, 1151
        coef1 = 6.6011e-24*(const.hartree/(2.*const.boltzmann*self.Temperature))**1.5
        coef2 = (const.planck)**3/(2.*const.pi*const.emass*const.boltzmann*self.Temperature)**1.5
        # next from Aped
        coef3 = (4.*const.pi/(const.boltzmannEv*self.Temperature/const.ryd2Ev))**(1.5)*(const.bohr)**3
        coef = coef2
#        print ' coefs = ', coef1, coef2, coef3
        nt = self.Temperature.size
        allRate = []
        effRate = []
        de = []
        lvl2 = []
        totalRate = np.zeros(nt, 'float64')
        lvlformat = '%7i%7i%10.2e%10.2e'
        for i, avalue in enumerate(self.Auto['avalue']):
            gUpper = float(self.Elvlc['mult'][self.Auto['lvl2'][i] -1])
            gLower = float(self.Higher.Elvlc['mult'][self.Auto['lvl1'][i] - 1])
    #        print i, autoa['lvl2'][i], gLower, gUpper, avalue
            ecm2 = self.Elvlc['ecm'][self.Auto['lvl2'][i] -1]
            if ecm2 == 0.:
                ecm2 = self.Elvlc['ecmth'][self.Auto['lvl2'][i] -1]
            de1 = ecm2*const.invCm2Erg - self.Ip*const.ev2Erg
            de.append(de1)
            dekt = de1/(const.boltzmann*self.Temperature)
            expkt = np.exp(-de1/(const.boltzmann*self.Temperature))
            rate = coef*gUpper*expkt*avalue/(2.*gLower)
            branch = self.Wgfa['avalueLvl'][self.Auto['lvl2'][i] -1]/(avalue + self.Wgfa['avalueLvl'][self.Auto['lvl2'][i] -1])
#            lvl2.append(self.Auto['lvl2'][i])
#            print i, self.Auto['lvl2'][i], cnt
#                print i, j, self.Wgfa['lvl1'][idx[-1 ]], self.Wgfa['lvl2'][idx[-1]], self.Auto['lvl2'][i]
    #        print ' lvl2, rate = ', autoa['lvl2'][i], rate
#            lvlstr = lvlformat%(self.Auto['lvl1'][i], self.Auto['lvl2'][i], avalue, branch)
#            if verbose:
#                tstr = lvlstr
#                rstr = lvlstr
#                if nt == 1:
#                    tstr += '%10.2e'%(self.Temperature)
#                    rstr += '%10.2e%10.2e'%(rate, rate*branch)
#                else:
#                    for it in range(nt):
#                        tstr += '%10.2e'%(self.Temperature[it])
#                        rstr += '%10.2e'%(rate[it])
#                print tstr
#                print rstr
            allRate.append(rate)
            effRate.append(rate*branch)
            totalRate += rate*branch
#            outpt.write(tstr +'\n')
#            outpt.write(rstr + '\n')
        self.DrRateLvl = {'rate':allRate, 'effRate':effRate, 'totalRate':totalRate,  'de':de, 'avalue':self.Auto['avalue']}   #, 'lvl1':lvl1, 'lvl2':lvl2} - in self.Auto
        #
        # -------------------------------------------------------------------------------------
        #
    def recombRate(self):
        '''Provides the total recombination rate coefficient.

        Calls diRate and eaRate'''
        #
        if hasattr(self, 'Temperature'):
            temperature=self.Temperature
        else:
            print ' temperature is not defined'
            self.RecombRate = None
        if self.Ion == 1:
#            print ' this is a neutral and has no recombination rate'
            self.RecombRate = {'rate':np.zeros_like(temperature), 'temperature':temperature}
            return
        self.rrRate()
        self.drRate()
        if not hasattr(self, 'DrRate'):
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
        if hasattr(self, 'Abundance'):
            ab=self.Abundance
        else:
            abundName = self.Defaults['abundfile']
            util.abundanceRead(abundancename = abundName)
        if hasattr(self, 'Ioneq'):
            ioneq=self.Ioneq
        else:
            ioneqname = self.Defaults['ioneqfile']
            self.IoneqAll = util.ioneqRead(ioneqname = ioneqname)
        #
        if hasattr(self, 'Temperature'):
            temperature=self.Temperature
        else:
            temperature = self.IoneqAll['ioneqTemperature']
#                temperature=self.IoneqTemperature
#        else:  temperature=np.asarray(temperature,'float32')
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
    def upsilonDescale(self, prot=0, diel=0):
        """Provides the temperatures and effective collision strengths (upsilons)
        set prot for proton rates
        otherwise, ce will be set for electron collision rates"""
        #
        #  xt=kt/de
        #
        #
        if prot:
            ce = 0
            try:
                nsplups=len(self.Psplups["lvl1"])
            except:
                self.Psplups=util.splupsRead(self.IonStr,prot=1)
                if type(self.Psplups) == types.NoneType:
                    self.PUpsilon = None
                    return
                else:
                    nsplups = len(self.Cilvl["lvl1"])
        elif diel:
            ce = 0
            try:
                nsplups = len(self.DielSplups["lvl1"])
            except:
                self.DielSplups = util.splupsRead(self.IonStr,diel=1)
                if type(self.DielSplups) == types.NoneType:
                    self.DielUpsilon = None
                    return
                else:
                    nsplups = len(self.DielSplups["lvl1"])
        else:
            ce=1
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
        if hasattr(self, 'Temperature'):
            temperature=self.Temperature
        else:
            print ' Temperature undefined'
            return
        #
        if hasattr(self, 'Elvlc'):
            nlvls=len(self.Elvlc["lvl"])
        else:
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
            ups = np.zeros((nsplups,ntemp),"Float64")
            exRate = np.zeros((nsplups,ntemp),"Float64")
            dexRate = np.zeros((nsplups,ntemp),"Float64")
        else:
            ups = np.zeros(nsplups,"Float64")
            exRate = np.zeros((nsplups,ntemp),"Float64")
            dexRate = np.zeros((nsplups,ntemp),"Float64")
        #
        for isplups in range(nsplups):
            if prot:
                # for proton rates
                l1=self.Psplups["lvl1"][isplups]-1
                l2=self.Psplups["lvl2"][isplups]-1
                ttype=self.Psplups["ttype"][isplups]
                cups=self.Psplups["cups"][isplups]
                nspl=self.Psplups["nspl"][isplups]
                dx=1./(float(nspl)-1.)
#                splups=self.Psplups["splups"][isplups,0:nspl]
                splups=self.Psplups["splups"][isplups]
                de=elvlc[l2]-elvlc[l1]
#                de=self.Psplups['de'][isplups]  # these are generally 0.
                kte = const.boltzmann*temp/(de*const.ryd2erg)
            elif diel:
                #
                l1 = self.DielSplups["lvl1"][isplups]-1
                l2 = self.DielSplups["lvl2"][isplups]-1
                ttype = self.DielSplups["ttype"][isplups]
                cups = self.DielSplups["cups"][isplups]
                nspl = self.DielSplups["nspl"][isplups]
                ttype = self.DielSplups["ttype"][isplups]
                dx = 1./(float(nspl)-1.)
#                splups = self.DielSplups["splups"][isplups,0:nspl]
                splups = self.DielSplups["splups"][isplups]
                de=self.DielSplups['de'][isplups]
                kte = const.boltzmann*temp/(de*const.ryd2erg)
            else:
                # electron collisional excitation
                l1=self.Splups["lvl1"][isplups]-1
                l2=self.Splups["lvl2"][isplups]-1
                ttype=self.Splups["ttype"][isplups]
                cups=self.Splups["cups"][isplups]
                nspl=self.Splups["nspl"][isplups]
                dx=1./(float(nspl)-1.)
##                print self.Splups["splups"][l1,l2]
#                splups=self.Splups["splups"][isplups,0:nspl]
                splups=self.Splups["splups"][isplups]
#                de=elvlc[l2]-elvlc[l1]
                de=self.Splups['de'][isplups]
                kte = const.boltzmann*temp/(de*const.ryd2erg)
            #
            der=0
            if ttype == 1:
                st=1.-np.log(cups)/np.log(kte+cups)
                xs=dx*np.arange(nspl)
                y2=interpolate.splrep(xs,splups,s=0)
                sups=interpolate.splev(st,y2,der=der)
#                sups=interpolate.spline(xs, splups, st)
                ups[isplups]=sups*np.log(kte+np.exp(1.))
            #
            if ttype == 2:
                st=kte/(kte+cups)
                xs=dx*np.arange(nspl)
                y2=interpolate.splrep(xs,splups,s=0)
                sups=interpolate.splev(st,y2,der=der)
                ups[isplups]=sups
            #
            if ttype == 3:
                st=kte/(kte+cups)
                xs=dx*np.arange(nspl)
                y2=interpolate.splrep(xs,splups,s=0)
                sups=interpolate.splev(st,y2,der=der)
                ups[isplups]=sups/(kte+1.)
            #
            if ttype == 4:
                st=1.-np.log(cups)/np.log(kte+cups)
                xs=dx*np.arange(nspl)
                y2=interpolate.splrep(xs,splups,s=0)
                sups=interpolate.splev(st,y2,der=der)
                ups[isplups]=sups*np.log(kte+cups)
            #
            if ttype == 5:
                # dielectronic rates
                st=kte/(kte+cups)
                xs=dx*np.arange(nspl)
                y2=interpolate.splrep(xs,splups,s=0)
                sups=interpolate.splev(st,y2,der=der)
                ups[isplups]=sups/(kte+0.)
            #
            #  descale proton values
            if ttype == 6:
                st=kte/(kte+cups)
                xs=dx*np.arange(nspl)
                y2=interpolate.splrep(xs,splups,s=0)
                sups=interpolate.splev(st,y2,der=der)
#                ups[isplups] = sups
                ups[isplups]=10.**sups
            #
            elif ttype > 6:  print ' t_type ne 1,2,3,4,5=',ttype,l1,l2
            #
            if ce:
                if self.Dielectronic:
                    # the dielectronic ions will eventually be discontinued
                    de = np.abs((self.Elvlc["eryd"][l2] - self.UpperIp/const.ryd2Ev) - self.Elvlc["eryd"][l1])
                else:
                    de = np.abs(self.Elvlc["eryd"][l2] - self.Elvlc["eryd"][l1])
#                print ' ce lvl1 %5i  lvl2 %5i de %10.2e'%(l1, l2, de)
                ekt = (de*const.ryd2erg)/(const.boltzmann*temp)
                fmult1 = float(self.Elvlc["mult"][l1])
                fmult2 = float(self.Elvlc["mult"][l2])
                dexRate[isplups] = const.collision*ups[isplups]/(fmult2*np.sqrt(temp))
                exRate[isplups] = const.collision*ups[isplups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
            elif diel:
#                print ' diel lvl1 %5i  lvl2 %5i de %10.2e'%(l1, l2, de)
                de = np.abs((self.Elvlc["eryd"][l2] - self.Ip/const.ryd2Ev) - self.Elvlc["eryd"][l1])
                ekt = (de*const.ryd2erg)/(const.boltzmann*temp)
                fmult1 = float(self.Elvlc["mult"][l1])
                fmult2 = float(self.Elvlc["mult"][l2])
                exRate[isplups] = const.collision*ups[isplups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
            elif prot:
                de = np.abs(self.Elvlc["eryd"][l2]-self.Elvlc["eryd"][l1])
                ekt = (de*1.57888e+5)/temp
                fmult1 = float(self.Elvlc["mult"][l1])
                fmult2 = float(self.Elvlc["mult"][l2])
                dexRate[isplups] = const.collision*ups[isplups]/(fmult2*np.sqrt(temp))
                exRate[isplups] = const.collision*ups[isplups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
        #
        ups=np.where(ups > 0.,ups,0.)
        #
        if prot == 1:
            self.PUpsilon = {'upsilon':ups, 'temperature':temperature, 'exRate':exRate, 'dexRate':dexRate}
        elif diel == 1:
            self.DielUpsilon = {'upsilon':ups, 'temperature':temperature, 'exRate':exRate}
        else:
            self.Upsilon = {'upsilon':ups, 'temperature':temperature, 'exRate':exRate, 'dexRate':dexRate}
        #
        # -------------------------------------------------------------------------
        #
    def setup(self, dir=0, verbose=0):
        '''
        if ion is initiated with setup=0, this allows the setup to be done at a later point
        perhaps, more importantly,  by setting dir to a directory cotaining the necessary files
        for a ChiantiPy ion, it allows one to setup an ion with files not in the current
        Chianti directory
        '''
        #
        # read in all data if in masterlist
        #  if not, there should still be ionization and recombination rates
        #
        MasterList = chdata.MasterList
        #
        if dir:
            fileName = os.path.join(dir, self.IonStr)
        else:
            fileName = util.ion2filename(self.IonStr)
        if self.IonStr in MasterList:
            self.Elvlc = util.elvlcRead('', filename=fileName+'.elvlc',  verbose=verbose)
            self.Wgfa = util.wgfaRead('', filename=fileName+'.wgfa')
            self.Nwgfa=len(self.Wgfa['lvl1'])
            nlvlWgfa = max(self.Wgfa['lvl2'])
            nlvlList =[nlvlWgfa]
#                print 'fileName = ', fileName
            splupsfile = fileName + '.splups'
            if os.path.isfile(splupsfile):
                # happens the case of fe_3 and prob. a few others
                self.Splups = util.splupsRead('', filename=fileName+'.splups')
                self.Nsplups=len(self.Splups['lvl1'])
                nlvlSplups = max(self.Splups['lvl2'])
                nlvlList.append(nlvlSplups)
            else:
                self.Nsplups = 0
                nlvlSplups = 0
##                self.Nlvls = nlvlElvlc
            #
            file = fileName +'.cilvl'
            if os.path.isfile(file):
                self.Cilvl = util.cireclvlRead('',filename = fileName, cilvl=1)
                self.Ncilvl=len(self.Cilvl['lvl1'])
                nlvlCilvl = max(self.Cilvl['lvl2'])
                nlvlList.append(nlvlCilvl)
            else:
                self.Ncilvl = 0
            #  .reclvl file may not exist
            reclvlfile = fileName +'.reclvl'
            if os.path.isfile(reclvlfile):
                self.Reclvl = util.cireclvlRead('',filename=fileName, reclvl=1)
                self.Nreclvl = len(self.Reclvl['lvl1'])
                nlvlReclvl = max(self.Reclvl['lvl2'])
                nlvlList.append(nlvlReclvl)
            else:
                self.Nreclvl = 0
            #  .dielsplups file may not exist
            dielsplupsfile = fileName +'.dielsplups'
            if os.path.isfile(dielsplupsfile):
                self.DielSplups = util.splupsRead('', filename=dielsplupsfile, diel=1)
                self.Ndielsplups=len(self.DielSplups["lvl1"])
                nlvlDielSplups = max(self.DielSplups['lvl2'])
                nlvlList.append(nlvlDielSplups)
            else:
                self.Ndielsplups = 0
            #
            #  psplups file may not exist
            psplupsfile = fileName +'.psplups'
            if os.path.isfile(psplupsfile):
                self.Psplups = util.splupsRead('', filename=psplupsfile,  prot=True)
                self.Npsplups=len(self.Psplups["lvl1"])
            else:
                self.Npsplups = 0
            #
            drparamsFile = fileName +'.drparams'
            if os.path.isfile(drparamsFile):
                self.DrParams = util.drRead(self.IonStr)
            #
            rrparamsFile = fileName +'.rrparams'
            if os.path.isfile(rrparamsFile):
                self.RrParams = util.rrRead(self.IonStr)

            #  not needed for ion, only phion
#                photoxfile = util.ion2filename(self.IonStr)+'.photox'
#                if os.path.isfile(photoxfile):
#                    self.Photox = util.photoxRead(self.IonStr)
            #
            # need to determine the number of levels that can be populated
            nlvlElvlc = len(self.Elvlc['lvl'])
#                print ' nlvlElvlc = ', nlvlElvlc
#                print ' other nlvls = ',  nlvlList
#                nlvlWgfa = max(self.Wgfa['lvl2'])
            #  elvlc file can have more levels than the rate level files
            self.Nlvls = min([nlvlElvlc, max(nlvlList)])
        #
        # -------------------------------------------------------------------------
        #
    def setupNew(self, dir=0, verbose=0):
        '''
        if ion is initiated with setup=0, this allows the setup to be done at a later point
        perhaps, more importantly,  by setting dir to a directory cotaining the necessary files
        for a ChiantiPy ion, it allows one to setup an ion with files not in the current
        Chianti directory
        this is a development version for integrating auoionizing levels etc into the main ion
        '''
        #
        # read in all data if in masterlist
        #  if not, there should still be ionization and recombination rates
        #
        MasterList = chdata.MasterList
        #
        if dir:
            fileName = os.path.join(dir, self.IonStr)
        else:
            fileName = util.ion2filename(self.IonStr)
        if self.IonStr in MasterList:
            self.Elvlc = util.elvlcRead('', filename=fileName+'.elvlc',  verbose=verbose)
            self.Wgfa = util.wgfaRead('', filename=fileName+'.wgfa', total=1)
            self.Nwgfa=len(self.Wgfa['lvl1'])
            nlvlWgfa = max(self.Wgfa['lvl2'])
            nlvlList =[nlvlWgfa]
#                print 'fileName = ', fileName
            splupsfile = fileName + '.splups'
            if os.path.isfile(splupsfile):
                # happens the case of fe_3 and prob. a few others
                self.Splups = util.splupsRead('', filename=fileName+'.splups')
                self.Nsplups=len(self.Splups['lvl1'])
                nlvlSplups = max(self.Splups['lvl2'])
                nlvlList.append(nlvlSplups)
            else:
                self.Nsplups = 0
                nlvlSplups = 0
##                self.Nlvls = nlvlElvlc
            #
            file = fileName +'.cilvl'
            if os.path.isfile(file):
                self.Cilvl = util.cireclvlRead('',filename = fileName, cilvl=1)
                self.Ncilvl=len(self.Cilvl['lvl1'])
                nlvlCilvl = max(self.Cilvl['lvl2'])
                nlvlList.append(nlvlCilvl)
            else:
                self.Ncilvl = 0
            #
            #  not longer using the reclvl files
            #  using the new rrlvl - radiaitive recombination rates only
            #  dielectronic rates derived from the autoionization values in the .auto file
            rrlvlfile = fileName +'.rrlvl'
            if os.path.isfile(rrlvlfile):
                self.Rrlvl = util.cireclvlRead('',filename=fileName, rrlvl=1)
                self.Nrrlvl = len(self.Rrlvl['lvl1'])
                nlvlRrlvl = max(self.Rrlvl['lvl2'])
                nlvlList.append(nlvlRrlvl)
            else:
                self.Nrrlvl = 0
            # in setupNew, the dielsplups files are disregarded
            #  .dielsplups file may not exist
#            dielsplupsfile = fileName +'.dielsplups'
#            if os.path.isfile(dielsplupsfile):
#                self.DielSplups = util.splupsRead('', filename=dielsplupsfile, diel=1)
#                self.Ndielsplups=len(self.DielSplups["lvl1"])
#                nlvlDielSplups = max(self.DielSplups['lvl2'])
#                nlvlList.append(nlvlDielSplups)
#            else:
#                self.Ndielsplups = 0
            #
            #  psplups file may not exist
            psplupsfile = fileName +'.psplups'
            if os.path.isfile(psplupsfile):
                self.Psplups = util.splupsRead('', filename=psplupsfile,  prot=True)
                self.Npsplups=len(self.Psplups["lvl1"])
            else:
                self.Npsplups = 0
            #
            drparamsFile = fileName +'.drparams'
            if os.path.isfile(drparamsFile):
                self.DrParams = util.drRead(self.IonStr)
            #
            rrparamsFile = fileName +'.rrparams'
            if os.path.isfile(rrparamsFile):
                self.RrParams = util.rrRead(self.IonStr)
            #
            # get autoionizing A-values
            autoFile = fileName + '.auto'
            if os.path.isfile(autoFile):
                self.Auto = util.wgfaRead('', filename=autoFile)

            #  not needed for ion, only phion
#                photoxfile = util.ion2filename(self.IonStr)+'.photox'
#                if os.path.isfile(photoxfile):
#                    self.Photox = util.photoxRead(self.IonStr)
            #
            # need to determine the number of levels that can be populated
            nlvlElvlc = len(self.Elvlc['lvl'])
#                print ' nlvlElvlc = ', nlvlElvlc
#                print ' other nlvls = ',  nlvlList
#                nlvlWgfa = max(self.Wgfa['lvl2'])
            #  elvlc file can have more levels than the rate level files
            self.Nlvls = min([nlvlElvlc, max(nlvlList)])
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

        includes ionization equilibrium and elemental abundances

        Note:  scipy.ndimage.filters also includes a range of filters.'''
        aspectrum = np.zeros_like(wavelength)
        nTemp = self.Temperature.size
        nDens = self.EDensity.size
        useFilter = filter[0]
        useFactor= filter[1]
        if hasattr(self, 'Intensity'):
            intensity = self.Intensity
        else:
            self.intensity()
            intensity = self.Intensity
        #
        lvl1 = []
        lvl2 = []
        if (nTemp == 1) and (nDens == 1):
            aspectrum = np.zeros_like(wavelength)
            if not 'errorMessage' in self.Intensity.keys():
                idx = util.between(self.Intensity['wvl'], [wavelength.min(), wavelength.max()])
                for iwvl in idx:
                    wvlCalc = self.Intensity['wvl'][iwvl]
                    aspectrum += useFilter(wavelength, wvlCalc, factor=useFactor)*intensity['intensity'][iwvl]
        else:
            nVar = max(nTemp, nDens)
            aspectrum = np.zeros((nVar, wavelength.size), 'float64')
            if not 'errorMessage' in self.Intensity.keys():
                idx = util.between(self.Intensity['wvl'], [wavelength.min(), wavelength.max()])
                for itemp in xrange(nVar):
                    for iwvl in idx:
                        wvlCalc = self.Intensity['wvl'][iwvl]
                        aspectrum[itemp] += useFilter(wavelength, wvlCalc, factor=useFactor)*self.Intensity['intensity'][itemp, iwvl]
    #                    for iwvl, wvlCalc in enumerate(self.Intensity['wvl']):
    #                        aspectrum[itemp] += useFilter(wavelength, wvlCalc, factor=useFactor)*self.Intensity['intensity'][itemp, iwvl]
        self.Spectrum = {'intensity':aspectrum,  'wvl':wavelength, 'filter':useFilter.__name__, 'filterWidth':useFactor}
        #
        # -------------------------------------------------------------------------------------
        #
    def populate(self, popCorrect=1, verbose=0, **kwargs):
        """
        Calculate level populations for specified ion.  This is a new version that will enable the calculation
        of dielectronic satellite lines without resorting to the dielectronic ions, such as c_5d
        possible keyword arguments include temperature, eDensity, pDensity, radTemperature and rStar
        """
        #
        #
        for one in kwargs.keys():
            if one not in chdata.keywordArgs:
                print ' following keyword is not understood - ',  one
        #
        nlvls=self.Nlvls
        nwgfa=self.Nwgfa
        nsplups=self.Nsplups
        npsplups=self.Npsplups
        #
        if kwargs.has_key('temperature'):
            self.Temperature = np.asarray(kwargs['temperature'])
            temperature = self.Temperature
        elif hasattr(self, 'Temperature'):
            temperature=self.Temperature
        else:
                print ' no temperature values have been set'
                return
        #
        if kwargs.has_key('eDensity'):
            self.EDensity = np.asarray(kwargs['eDensity'])
            eDensity = self.EDensity
        elif hasattr(self, 'EDensity'):
            eDensity = self.EDensity
        else:
            print ' no eDensity values have been set'
            return
        #
        if kwargs.has_key('pDensity'):
            if kwargs['pDensity'] == 'default':
                self.p2eRatio()
                protonDensity = self.ProtonDensityRatio*self.EDensity
            else:
                try:
                    self.PDensity = np.asarray(kwargs['pDensity'])
                except:
                    print ' could not interpret value for keyword pDensity'
                    print ' should be either "default" or a number or array'
                    return
        else:
            if hasattr(self, 'PDensity'):
                protonDensity = self.PDensity
            else:
                self.p2eRatio()
                self.PDensity = self.ProtonDensityRatio*self.EDensity
                protonDensity = self.PDensity
                print ' proton density not specified, set to "default"'
        #
        if 'radTemperature' in kwargs.keys() and 'rStar' in kwargs.keys():
            self.RadTemperature = np.asarray(kwargs['radTemperature'])
            radTemperature = np.array(self.RadTemperature)
            self.RStar = np.asarray(kwargs['rStar'])
            rStar = np.asarray(self.RStar)
        elif hasattr(self, 'RadTemperature') and hasattr(self, 'RStar'):
            radTemperature = self.RadTemperature
            rStar = self.RStar
        #
        #
        # the Dielectronic test should eventually go away
        if popCorrect and (not self.Dielectronic):
            if self.Ncilvl:
                ci = 1
                cilvl = self.Cilvl
                if hasattr(self, 'CilvlRate'):
                    cilvlRate = self.CilvlRate
                else:
                    self.cireclvlDescale('cilvl')
                    cilvlRate = self.CilvlRate
                self.recombRate()
                #
                lowers = util.zion2name(self.Z, self.Ion-1)
                # get the lower ionization stage
                lower = ion(lowers, temperature=self.Temperature, eDensity = self.EDensity)
                lower.ionizRate()
                # need to get multiplicity of lower ionization stage
                lowMult = lower.Elvlc['mult']
            else:
                ci = 0
#            try:
            if self.Nreclvl:
                rec = 1
                reclvl = self.Reclvl
                if hasattr(self, 'ReclvlRate'):
                    reclvlRate = self.ReclvlRate
                else:
#                    print ' doing reclvlDescale in populate'
                    self.cireclvlDescale('reclvl')
                    reclvlRate = self.ReclvlRate
            elif self.Ndielsplups:
                self.upsilonDescale(diel=1)
                dielexRate = self.DielUpsilon['exRate']
                rec = 1
            else:
                rec = 0
#            except:
#                self.Reclvl = util.cireclvlRead(self.IonStr,'reclvl' )
#                reclvl = self.Reclvl
#                self.reclvlDescale()
#                if type(self.Reclvl) != types.NoneType:
#                    rec = 1
#                    self.ionizRate()
#                    #  get the higher ionization stage
#                    highers = util.zion2name(self.Z, self.Ion+1)
##                   print ' highers = ', highers
#                    higher = ion(highers, temperature=self.Temperature, eDensity=self.EDensity)
#                    higher.recombRate()
#                else:
#                    rec = 0
            #
        else:
            ci = 0
            rec = 0
        #
        if rec:
            if self.Ndielsplups:
                self.upsilonDescale(diel=1)
                dielexRate = self.DielUpsilon['exRate']
            # get ionization rate of this iion
            self.ionizRate()
            #  get the higher ionization stage
            highers = util.zion2name(self.Z, self.Ion+1)
            higher = ion(highers, temperature=self.Temperature, eDensity=self.EDensity)
            higher.recombRate()
#        print ' nlvls, ci, rec = ', nlvls, ci, rec
        #
        rad=np.zeros((nlvls+ci+rec,nlvls+ci+rec),"float64")  #  the populating matrix for radiative transitions
        #
        #
        for iwgfa in range(nwgfa):
            l1 = self.Wgfa["lvl1"][iwgfa]-1
            l2 = self.Wgfa["lvl2"][iwgfa]-1
            rad[l1+ci,l2+ci] += self.Wgfa["avalue"][iwgfa]
            rad[l2+ci,l2+ci] -= self.Wgfa["avalue"][iwgfa]
            # photo-excitation and stimulated emission
            if self.RadTemperature:
                if not self.RStar:
                    dilute = 0.5
                else:
                    dilute = util.dilute(self.RStar)
                # next - don't include autoionization lines
                if abs(self.Wgfa['wvl'][iwgfa]) > 0.:
                    de = const.invCm2Erg*(self.Elvlc['ecm'][l2] - self.Elvlc['ecm'][l1])
                    dekt = de/(const.boltzmann*self.RadTemperature)
                    # photoexcitation
                    phexFactor = dilute*(self.Elvlc['mult'][l2]/self.Elvlc['mult'][l1])/(np.exp(dekt) -1.)
                    rad[l2+ci,l1+ci] += self.Wgfa["avalue"][iwgfa]*phexFactor
                    rad[l1+ci,l1+ci] -= self.Wgfa["avalue"][iwgfa]*phexFactor
                    # stimulated emission
                    stemFactor = dilute/(1. - np.exp(-dekt))
                    rad[l1+ci,l2+ci] += self.Wgfa["avalue"][iwgfa]*stemFactor
                    rad[l2+ci,l2+ci] -= self.Wgfa["avalue"][iwgfa]*stemFactor

        #
        self.rad=rad
        #
        if self.Nsplups:
            self.upsilonDescale()
            ups = self.Upsilon['upsilon']
            exRate = self.Upsilon['exRate']
            dexRate = self.Upsilon['dexRate']
        #
        if npsplups:
            self.upsilonDescale(prot=1)
#            pups = self.PUpsilon['upsilon']
            pexRate = self.PUpsilon['exRate']
            pdexRate = self.PUpsilon['dexRate']
        #
        temp=temperature
        ntemp=temp.size
        #
        cc=const.collision*self.EDensity
        ndens=cc.size
        if npsplups:
            cp=const.collision*protonDensity
        if ntemp > 1 and ndens >1 and ntemp != ndens:
            print ' unless temperature or eDensity are single values'
            print ' the number of temperatures values must match the '
            print ' the number of eDensity values'
            return
        #
        # get corrections for recombination and excitation
        #
        #
        #  first, for ntemp=ndens=1
        if ndens==1 and ntemp==1:
            popmat=np.copy(rad)
            for isplups in range(0,nsplups):
                l1=self.Splups["lvl1"][isplups]-1
                l2=self.Splups["lvl2"][isplups]-1
                #
                popmat[l1+ci,l2+ci] += self.EDensity*dexRate[isplups]
                popmat[l2+ci,l1+ci] += self.EDensity*exRate[isplups]
                popmat[l1+ci,l1+ci] -= self.EDensity*exRate[isplups]
                popmat[l2+ci,l2+ci] -= self.EDensity*dexRate[isplups]
                #
            for isplups in range(0,npsplups):
                l1=self.Psplups["lvl1"][isplups]-1
                l2=self.Psplups["lvl2"][isplups]-1
                 #
                popmat[l1+ci,l2+ci] += self.PDensity*pdexRate[isplups]
                popmat[l2+ci,l1+ci] += self.PDensity*pexRate[isplups]
                popmat[l1+ci,l1+ci] -= self.PDensity*pexRate[isplups]
                popmat[l2+ci,l2+ci] -= self.PDensity*pdexRate[isplups]
           # now include ionization rate from
            if ci:
#                print ' ci = ', ci
                #
                # the ciRate can be computed for all temperatures
                #
                ciTot = 0.
                for itrans in range(len(cilvl['lvl1'])):
                    lvl1 = cilvl['lvl1'][itrans]-1
                    lvl2 = cilvl['lvl2'][itrans]-1
#                    de = cilvl['de'][itrans]
#                    ekt = (de*1.57888e+5)/temperature
#                    mult = lowMult[lvl1-1]
                    # this is kind of double booking the ionization rate components
                    popmat[lvl2+ci, lvl1] += self.EDensity*self.CilvlRate['rate'][itrans]
                    popmat[lvl1, lvl1] -= self.EDensity*self.CilvlRate['rate'][itrans]
                    ciTot += self.EDensity*self.CilvlRate['rate'][itrans]
                #
                popmat[1, 0] += (self.EDensity*lower.IonizRate['rate'] - ciTot)
                popmat[0, 0] -= (self.EDensity*lower.IonizRate['rate'] - ciTot)
                popmat[0, 1] += self.EDensity*self.RecombRate['rate']
                popmat[1, 1] -= self.EDensity*self.RecombRate['rate']
            if rec:
                #
                if self.Ndielsplups:
                    branch = np.zeros(self.Ndielsplups, 'float64')
                    for isplups in range(0,self.Ndielsplups):
                        l1 = self.DielSplups["lvl1"][isplups]-1 + nlvls
                        l2 = self.DielSplups["lvl2"][isplups]-1
                        auto = rad[l1+ci, l2+ci]
                        avalue = rad[:, l2+ci]
                        good = avalue > 0.
                        avalueTot = avalue[good].sum()
                        branch[isplups] = (avalueTot-auto)/avalueTot
#                        print ' l1 %4i l2 %4i auto %10.2e  avalue %10.2e tot %10.3f'%( l1,  l2,  auto,  avalueTot,  branch[isplups])
                        self.DielUpsilon['branch'] =  branch
                    #
                    dielTot = 0.
#                    print ' Ndielsplups > 0 '
                    for isplups in range(0,self.Ndielsplups):
                        l1 = self.DielSplups["lvl1"][isplups]-1 + nlvls
                        l2 = self.DielSplups["lvl2"][isplups]-1
                         #
#                        print ' l1, l2, dielexRate = ', l1, l2, dielexRate[isplups]
                        popmat[l2+ci,l1+ci] += self.EDensity*dielexRate[isplups]
                        popmat[l1+ci,l1+ci] -= self.EDensity*dielexRate[isplups]
                        #
                        dielTot += self.EDensity*dielexRate[isplups]*branch[isplups]
                else:
                    dielTot = 0.
                #
#                print ' rec, dielTot  = ', rec,  dielTot
                #
                for itrans in range(self.Nreclvl):
#                    lvl1 = reclvl['lvl1'][itrans]-1
                    lvl2 = reclvl['lvl2'][itrans]-1
                    popmat[lvl2+ci, -1] += self.EDensity*reclvlRate['rate'][itrans]
                    popmat[-1, -1] -= self.EDensity*reclvlRate['rate'][itrans]
                if self.Nreclvl:
                    reclvlRateTot = reclvlRate['rate'].sum(axis=0)
                else:
                    reclvlRateTot = 0.

                #
                popmat[-1,  ci] += self.EDensity*self.IonizRate['rate']
                popmat[ci, ci] -= self.EDensity*self.IonizRate['rate']
                # next 2 line take care of overbooking
                popmat[ci, -1] += self.EDensity*(higher.RecombRate['rate']- reclvlRateTot) - dielTot
                popmat[-1, -1] -= self.EDensity*(higher.RecombRate['rate']- reclvlRateTot) - dielTot
#                print ' higher, rec , dieltot = ',  self.EDensity*higher.RecombRate['rate'], self.EDensity*reclvlRate['rate'].sum(axis=0),  dielTot
            # normalize to unity
            norm=np.ones(nlvls+ci+rec,'float64')
            if ci:
                norm[0] = 0.
            if rec:
                norm[nlvls+ci+rec-1] = 0.
            popmat[nlvls+ci+rec-1]=norm
            b=np.zeros(nlvls+ci+rec,'float64')
            b[nlvls+ci+rec-1]=1.
#            print ' norm = ', norm
#            print 'popmat, last line',  popmat[-1]
#            print ' b = ', b
#            popmat[nlvls/2]=norm
#            b=np.zeros(nlvls+ci+rec,'float64')
#            b[nlvls/2]=1.
#            if rec:
#                fullpop = np.linalg.solve(popmat,b)
#                pop = fullpop[ci:ci+nlvls+rec-1]
#            else:
#                fullpop = np.linalg.solve(popmat,b)
#                pop = fullpop[ci:]
#            fullpop = np.linalg.solve(popmat,b)
            try:
                fullpop=np.linalg.solve(popmat,b)
                pop = fullpop[ci:ci+nlvls]
            except np.linalg.LinAlgError:
                pop = np.zeros(nlvls, 'float64')
#                print ' error in matrix inversion, setting populations to zero at T = ', ('%8.2e')%(temperature)
            #
        #   next, in case of a single eDensity value
#            pop = np.linalg.solve(popmat,b)
        elif ndens == 1:
            pop=np.zeros((ntemp, nlvls),"float64")
#            pop=np.zeros((ntemp,ci + nlvls + rec),"float64")
            for itemp in range(0,ntemp):
                popmat=np.copy(rad)
                for isplups in range(0,nsplups):
                    l1=self.Splups["lvl1"][isplups]-1
                    l2=self.Splups["lvl2"][isplups]-1
                    popmat[l1+ci,l2+ci] += self.EDensity*dexRate[isplups, itemp]
                    popmat[l2+ci,l1+ci] += self.EDensity*exRate[isplups, itemp]
                    popmat[l1+ci,l1+ci] -= self.EDensity*exRate[isplups, itemp]
                    popmat[l2+ci,l2+ci] -= self.EDensity*dexRate[isplups, itemp]
                for isplups in range(0,npsplups):
                    l1=self.Psplups["lvl1"][isplups]-1
                    l2=self.Psplups["lvl2"][isplups]-1
                    # for proton excitation, the levels are all below the ionization potential
                     #
                    popmat[l1+ci,l2+ci] += self.PDensity[itemp]*pdexRate[isplups, itemp]
                    popmat[l2+ci,l1+ci] += self.PDensity[itemp]*pexRate[isplups, itemp]
                    popmat[l1+ci,l1+ci] -= self.PDensity[itemp]*pexRate[isplups, itemp]
                    popmat[l2+ci,l2+ci] -= self.PDensity[itemp]*pdexRate[isplups, itemp]
                # now include ionization rate from
                if ci:
#                    print ' ci = ', ci
                    #
                    # the ciRate can be computed for all temperatures
                    #
                    ciTot = 0.
                    for itrans in range(len(cilvl['lvl1'])):
                        lvl1 = cilvl['lvl1'][itrans]-1
                        lvl2 = cilvl['lvl2'][itrans]-1
#                        de = cilvl['de'][itrans]
#                        ekt = (de*1.57888e+5)/temperature
                        mult = lowMult[lvl1-1]
                        # this is kind of double booking the ionization rate components
                        popmat[lvl2+ci, lvl1] += self.EDensity*self.CilvlRate['rate'][itrans, itemp]
                        popmat[lvl1, lvl1] -= self.EDensity*self.CilvlRate['rate'][itrans, itemp]
                        ciTot += self.EDensity*self.CilvlRate['rate'][itrans, itemp]
#                        popmat[lvl2, lvl1-1] += self.EDensity*cirate[itemp]
#                        popmat[lvl1-1, lvl1-1] -= self.EDensity*cirate[itemp]
                    popmat[1, 0] += (self.EDensity*lower.IonizRate['rate'][itemp] - ciTot)
                    popmat[0, 0] -= (self.EDensity*lower.IonizRate['rate'][itemp] - ciTot)
                    popmat[0, 1] += self.EDensity*self.RecombRate['rate'][itemp]
                    popmat[1, 1] -= self.EDensity*self.RecombRate['rate'][itemp]
                if rec:
                #
                    if self.Ndielsplups:
                        branch = np.zeros(self.Ndielsplups, 'float64')
                        for isplups in range(0,self.Ndielsplups):
                            l1 = self.DielSplups["lvl1"][isplups]-1 + nlvls
                            l2 = self.DielSplups["lvl2"][isplups]-1
                            auto = rad[l1+ci, l2+ci]
                            avalue = rad[:, l2+ci]
                            good = avalue > 0.
                            avalueTot = avalue[good].sum()
                            branch[isplups] = (avalueTot-auto)/avalueTot
#                            print ' l1 %4i l2 %4i auto %10.2e  avalue %10.2e tot %10.3f'%( l1,  l2,  auto,  avalueTot,  branch[isplups])
                            self.DielUpsilon['branch'] =  branch
                        #
                        dielTot = 0.
#                        print ' Ndielsplups > 0 '
                        for isplups in range(0,self.Ndielsplups):
                            l1 = self.DielSplups["lvl1"][isplups]-1 + nlvls
                            l2 = self.DielSplups["lvl2"][isplups]-1
                             #
                            popmat[l2+ci,l1+ci] += self.EDensity*dielexRate[isplups, itemp]
                            popmat[l1+ci,l1+ci] -= self.EDensity*dielexRate[isplups, itemp]
                            #
                            dielTot += self.EDensity*dielexRate[isplups, itemp]*branch[isplups]
                    else:
                        dielTot = 0.
                    if self.Nreclvl:
                        recTot = self.ReclvlRate['rate'][:, itemp].sum()
                    else:
                        recTot = 0.
                #
                    popmat[-1,  ci] += self.EDensity*self.IonizRate['rate'][itemp]
                    popmat[ci, ci] -= self.EDensity*self.IonizRate['rate'][itemp]
                    popmat[ci, -1] += self.EDensity*(higher.RecombRate['rate'][itemp]- recTot) - dielTot
                    popmat[-1, -1] -= self.EDensity*(higher.RecombRate['rate'][itemp]- recTot) + dielTot
#                    popmat[ci, -1] += self.EDensity*(higher.RecombRate['rate'][itemp]- self.ReclvlRate['rate'][:, itemp].sum()) - dielTot
#                    popmat[-1, -1] -= self.EDensity*(higher.RecombRate['rate'][itemp]- self.ReclvlRate['rate'][:, itemp].sum()) + dielTot
#                    popmat[ci, -1] += self.EDensity*higher.RecombRate['rate'][itemp]
#                    popmat[-1, -1] -= self.EDensity*higher.RecombRate['rate'][itemp]
                    #
#                    for itrans in range(len(reclvl['lvl1'])):
                    for itrans in range(self.Nreclvl):
                        lvl1 = reclvl['lvl1'][itrans]-1
                        lvl2 = reclvl['lvl2'][itrans]-1
                        popmat[lvl2+ci, -1] += self.EDensity*self.ReclvlRate['rate'][itrans, itemp]
                        popmat[-1, -1] -= self.EDensity*self.ReclvlRate['rate'][itrans, itemp]
                # normalize to unity
                norm=np.ones(nlvls+ci+rec,'float64')
                if ci:
                    norm[0] = 0.
                if rec:
                    norm[-1] = 0.
                popmat[nlvls+ci+rec-1]=norm
                b=np.zeros(nlvls+ci+rec,'float64')
                b[nlvls+ci+rec-1]=1.
                try:
                    thispop=np.linalg.solve(popmat,b)
                    pop[itemp] = thispop[ci:ci+nlvls]
                except np.linalg.LinAlgError:
                    pop[itemp] = np.zeros(nlvls, 'float64')
#                    print ' error in matrix inversion, setting populations to zero at T = ', ('%8.2e')%(temperature[itemp])
            #
        elif ntemp == 1:
#            pop=np.zeros((ndens,nlvls),"float64")
            pop=np.zeros((ndens,nlvls),"float64")
            for idens in range(0,ndens):
                popmat=np.copy(rad)
                for isplups in range(0,nsplups):
                    l1=self.Splups["lvl1"][isplups]-1
                    l2=self.Splups["lvl2"][isplups]-1
#                    if self.Dielectronic:
#                        de=np.abs((self.Elvlc["eryd"][l2]-self.Ip/const.ryd2Ev)-self.Elvlc["eryd"][l1])
#                    else:
#                        de=np.abs(self.Elvlc["eryd"][l2]-self.Elvlc["eryd"][l1])
#                    ekt=(de*1.57888e+5)/temp
#                    fmult1=float(self.Elvlc["mult"][l1])
#                    fmult2=float(self.Elvlc["mult"][l2])
#                    popmat[l1+ci,l2+ci]+=cc[idens]*ups[isplups]/(fmult2*np.sqrt(temp))
#                    popmat[l2+ci,l1+ci]+=cc[idens]*ups[isplups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l1+ci,l1+ci]-=cc[idens]*ups[isplups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l2+ci,l2+ci]-=cc[idens]*ups[isplups]/(fmult2*np.sqrt(temp))
                #
                    popmat[l1+ci,l2+ci] += self.EDensity[idens]*dexRate[isplups]
                    popmat[l2+ci,l1+ci] += self.EDensity[idens]*exRate[isplups]
                    popmat[l1+ci,l1+ci] -= self.EDensity[idens]*exRate[isplups]
                    popmat[l2+ci,l2+ci] -= self.EDensity[idens]*dexRate[isplups]
                #
                for isplups in range(0,npsplups):
                    l1=self.Psplups["lvl1"][isplups]-1
                    l2=self.Psplups["lvl2"][isplups]-1
#                    # for proton excitation, the levels are all below the ionization potential
#                    de=np.abs(self.Elvlc["eryd"][l2]-self.Elvlc["eryd"][l1])
#                    ekt=(de*1.57888e+5)/temp
#                    fmult1=float(self.Elvlc["mult"][l1])
#                    fmult2=float(self.Elvlc["mult"][l2])
#                    popmat[l1+ci,l2+ci]+=cp[idens]*pups[isplups]/(fmult2*np.sqrt(temp))
#                    popmat[l2+ci,l1+ci]+=cp[idens]*pups[isplups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l1+ci,l1+ci]-=cp[idens]*pups[isplups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l2+ci,l2+ci]-=cp[idens]*pups[isplups]/(fmult2*np.sqrt(temp))
                 #
                    popmat[l1+ci,l2+ci] += self.PDensity[idens]*pdexRate[isplups]
                    popmat[l2+ci,l1+ci] += self.PDensity[idens]*pexRate[isplups]
                    popmat[l1+ci,l1+ci] -= self.PDensity[idens]*pexRate[isplups]
                    popmat[l2+ci,l2+ci] -= self.PDensity[idens]*pdexRate[isplups]
                # now include ionization rate from
                if ci:
#                    print ' ci = ', ci
                    #
                    #
                    ciTot = 0.
                    for itrans in range(len(cilvl['lvl1'])):
                        lvl1 = cilvl['lvl1'][itrans] -1
                        lvl2 = cilvl['lvl2'][itrans] -1
#                        de = cilvl['de'][itrans]
#                        ekt = (de*1.57888e+5)/temperature
                        # this is kind of double booking the ionization rate components
#                        popmat[lvl2, lvl1-1] += self.EDensity[idens]*cirate
#                        popmat[lvl1-1, lvl1-1] -= self.EDensity[idens]*cirate
                        popmat[lvl2+ci, lvl1] += self.EDensity[idens]*self.CilvlRate['rate'][itrans]
                        popmat[lvl1, lvl1] -= self.EDensity[idens]*self.CilvlRate['rate'][itrans]
                        ciTot += self.EDensity[idens]*self.CilvlRate['rate'][itrans]
                    popmat[1, 0] += (self.EDensity[idens]*lower.IonizRate['rate'] -ciTot)
                    popmat[0, 0] -= (self.EDensity[idens]*lower.IonizRate['rate'] -ciTot)
                    popmat[0, 1] += self.EDensity[idens]*self.RecombRate['rate']
                    popmat[1, 1] -= self.EDensity[idens]*self.RecombRate['rate']
                if rec:
#                    print ' rec = ', rec
                    if self.Ndielsplups:
                        branch = np.zeros(self.Ndielsplups, 'float64')
                        for isplups in range(0,self.Ndielsplups):
                            l1 = self.DielSplups["lvl1"][isplups]-1 + nlvls
                            l2 = self.DielSplups["lvl2"][isplups]-1
                            auto = rad[l1+ci, l2+ci]
                            avalue = rad[:, l2+ci]
                            good = avalue > 0.
                            avalueTot = avalue[good].sum()
                            branch[isplups] = (avalueTot-auto)/avalueTot
#                            print ' l1 %4i l2 %4i auto %10.2e  avalue %10.2e tot %10.3f'%( l1,  l2,  auto,  avalueTot,  branch[isplups])
                            self.DielUpsilon['branch'] =  branch
                        #
                        dielTot = 0.
#                        print ' Ndielsplups > 0 '
                        for isplups in range(0,self.Ndielsplups):
                            l1 = self.DielSplups["lvl1"][isplups]-1 + nlvls
                            l2 = self.DielSplups["lvl2"][isplups]-1
                             #
                            popmat[l2+ci,l1+ci] += self.EDensity[idens]*dielexRate[isplups]
                            popmat[l1+ci,l1+ci] -= self.EDensity[idens]*dielexRate[isplups]
                            #
                            dielTot += self.EDensity[idens]*dielexRate[isplups]*branch[isplups]
                    else:
                        dielTot = 0.
                    if self.Nreclvl:
#                        print ' ReclvlRate.shape = ', self.ReclvlRate['rate'].shape
                        recTot = self.ReclvlRate['rate'].sum()
                    else:
                        recTot = 0.
                    #
                    popmat[-1,  ci] += self.EDensity[idens]*self.IonizRate['rate']
                    popmat[ci, ci] -= self.EDensity[idens]*self.IonizRate['rate']
                    popmat[ci, -1] += self.EDensity[idens]*(higher.RecombRate['rate']
                        - recTot) - dielTot
                    popmat[-1, -1] -= self.EDensity[idens]*(higher.RecombRate['rate']
                        - recTot) + dielTot
#                    popmat[ci, -1] += self.EDensity[idens]*higher.RecombRate['rate']
#                    popmat[-1, -1] -= self.EDensity[idens]*higher.RecombRate['rate']
                    #
#                    for itrans in range(len(reclvl['lvl1'])):
                    for itrans in range(self.Nreclvl):
                        lvl1 = reclvl['lvl1'][itrans]-1
                        lvl2 = reclvl['lvl2'][itrans]-1
                        popmat[lvl2+ci, -1] += self.EDensity[idens]*self.ReclvlRate['rate'][itrans]
                        popmat[-1, -1] -= self.EDensity[idens]*self.ReclvlRate['rate'][itrans]
                # normalize to unity
                norm=np.ones(nlvls+ci+rec,'float64')
                if ci:
                    norm[0] = 0.
                if rec:
                    norm[-1] = 0.
                popmat[nlvls+ci+rec-1]=norm
                b=np.zeros(nlvls+ci+rec,'float64')
                b[nlvls+ci+rec-1]=1.
                try:
                    thispop=np.linalg.solve(popmat,b)
                    pop[idens] = thispop[ci:ci+nlvls]
                except np.linalg.LinAlgError:
                    pop[idens] = np.zeros(nlvls, 'float64')
#                    print ' error in matrix inversion, setting populations to zero at eDensity = ', ('%8.2e')%(eDensity[idens])
#                thispop=np.linalg.solve(popmat,b)
#                if rec:
#                    pop[idens] = thispop[ci:ci+nlvls+rec-1]
#                else:
#                    pop[idens] = thispop[ci:]
#                pop[idens] = thispop[ci:ci+nlvls]
                #
        elif ntemp>1  and ntemp==ndens:
            pop=np.zeros((ntemp,nlvls),"float64")
#            pop=np.zeros((ntemp,ci+nlvls+rec),"float64")
            for itemp in range(0,ntemp):
                temp=self.Temperature[itemp]
                popmat=np.copy(rad)
                for isplups in range(0,nsplups):
                    l1=self.Splups["lvl1"][isplups]-1
                    l2=self.Splups["lvl2"][isplups]-1
#                    if self.Dielectronic:
#                        de=np.abs((self.Elvlc["eryd"][l2]-self.Ip/const.ryd2Ev)-self.Elvlc["eryd"][l1])
#                    else:
#                        de=np.abs(self.Elvlc["eryd"][l2]-self.Elvlc["eryd"][l1])
#                    ekt=(de*1.57888e+5)/temp
#                    fmult1=float(self.Elvlc["mult"][l1])
#                    fmult2=float(self.Elvlc["mult"][l2])
#                    popmat[l1+ci,l2+ci]+=cc[itemp]*ups[isplups,itemp]/(fmult2*np.sqrt(temp))
#                    popmat[l2+ci,l1+ci]+=cc[itemp]*ups[isplups,itemp]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l1+ci,l1+ci]-=cc[itemp]*ups[isplups,itemp]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l2+ci,l2+ci]-=cc[itemp]*ups[isplups,itemp]/(fmult2*np.sqrt(temp))
                    #
                    popmat[l1+ci,l2+ci] += self.EDensity[itemp]*dexRate[isplups, itemp]
                    popmat[l2+ci,l1+ci] += self.EDensity[itemp]*exRate[isplups, itemp]
                    popmat[l1+ci,l1+ci] -= self.EDensity[itemp]*exRate[isplups, itemp]
                    popmat[l2+ci,l2+ci] -= self.EDensity[itemp]*dexRate[isplups, itemp]
                # proton rates
                for isplups in range(0,npsplups):
                    l1=self.Psplups["lvl1"][isplups]-1
                    l2=self.Psplups["lvl2"][isplups]-1
                    # for proton excitation, the levels are all below the ionization potential
#                    de=np.abs(self.Elvlc["eryd"][l2]-self.Elvlc["eryd"][l1])
#                    ekt=(de*1.57888e+5)/temp
#                    fmult1=float(self.Elvlc["mult"][l1])
#                    fmult2=float(self.Elvlc["mult"][l2])
#                    popmat[l1+ci,l2+ci]+=cp[itemp]*pups[isplups,itemp]/(fmult2*np.sqrt(temp))
#                    popmat[l2+ci,l1+ci]+=cp[itemp]*pups[isplups,itemp]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l1+ci,l1+ci]-=cp[itemp]*pups[isplups,itemp]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l2+ci,l2+ci]-=cp[itemp]*pups[isplups,itemp]/(fmult2*np.sqrt(temp))
                     #
                    popmat[l1+ci,l2+ci] += self.PDensity[itemp]*pdexRate[isplups, itemp]
                    popmat[l2+ci,l1+ci] += self.PDensity[itemp]*pexRate[isplups, itemp]
                    popmat[l1+ci,l1+ci] -= self.PDensity[itemp]*pexRate[isplups, itemp]
                    popmat[l2+ci,l2+ci] -= self.PDensity[itemp]*pdexRate[isplups, itemp]
                # now include ionization rate from
                if ci:
#                   print ' ci = ', ci
                    #
                    # the ciRate can be computed for all temperatures
                    #
                    ciTot = 0.
                    for itrans in range(len(cilvl['lvl1'])):
                        lvl1 = cilvl['lvl1'][itrans] -1
                        lvl2 = cilvl['lvl2'][itrans] -1
                        # this is kind of double booking the ionization rate components
#                        popmat[lvl2, lvl1-1] += self.EDensity[itemp]*cirate[itemp]
#                        popmat[lvl1-1, lvl1-1] -= self.EDensity[itemp]*cirate[itemp]
                        popmat[lvl2+ci, lvl1] += self.EDensity[itemp]*self.CilvlRate['rate'][itrans, itemp]
                        popmat[lvl1, lvl1] -= self.EDensity[itemp]*self.CilvlRate['rate'][itrans, itemp]
                        ciTot += self.EDensity[itemp]*self.CilvlRAte['rate'][itrans, itemp]
                    popmat[1, 0] += (self.EDensity[itemp]*lower.IonizRate['rate'][itemp] - ciTot)
                    popmat[0, 0] -= (self.EDensity[itemp]*lower.IonizRate['rate'][itemp] - ciTot)
                    popmat[0, 1] += self.EDensity[itemp]*self.RecombRate['rate'][itemp]
                    popmat[1, 1] -= self.EDensity[itemp]*self.RecombRate['rate'][itemp]
                if rec:
                #
                    if self.Ndielsplups:
                        branch = np.zeros(self.Ndielsplups, 'float64')
                        for isplups in range(0,self.Ndielsplups):
                            l1 = self.DielSplups["lvl1"][isplups]-1 + nlvls
                            l2 = self.DielSplups["lvl2"][isplups]-1
                            auto = rad[l1+ci, l2+ci]
                            avalue = rad[:, l2+ci]
                            good = avalue > 0.
                            avalueTot = avalue[good].sum()
                            branch[isplups] = (avalueTot-auto)/avalueTot
#                            print ' l1 %4i l2 %4i auto %10.2e  avalue %10.2e tot %10.3f'%( l1,  l2,  auto,  avalueTot,  branch[isplups])
                            self.DielUpsilon['branch'] =  branch
                        #
                        dielTot = 0.
                        for isplups in range(0,self.Ndielsplups):
                            l1 = self.DielSplups["lvl1"][isplups]-1 + nlvls
                            l2 = self.DielSplups["lvl2"][isplups]-1
                             #
                            popmat[l2+ci,l1+ci] += self.EDensity[itemp]*dielexRate[isplups, itemp]
                            popmat[l1+ci,l1+ci] -= self.EDensity[itemp]*dielexRate[isplups, itemp]
                            #
                            dielTot += self.EDensity[itemp]*dielexRate[isplups, itemp]*branch[isplups]
                    else:
                        dielTot = 0.
                    if self.Nreclvl:
                        recTot = self.ReclvlRate['rate'][:, itemp].sum()
                    else:
                        recTot = 0.
                #
#                   print ' rec = ', rec
                    popmat[-1,  ci] += self.EDensity[itemp]*self.IonizRate['rate'][itemp]
                    popmat[ci, ci] -= self.EDensity[itemp]*self.IonizRate['rate'][itemp]
                    popmat[ci, -1] += self.EDensity[itemp]*(higher.RecombRate['rate'][itemp]
                        - recTot) - dielTot
                    popmat[-1, -1] -= self.EDensity[itemp]*(higher.RecombRate['rate'][itemp]
                        - recTot) + dielTot
#                    popmat[ci, -1] += self.EDensity[itemp]*higher.RecombRate['rate'][itemp]
#                    popmat[-1, -1] -= self.EDensity[itemp]*higher.RecombRate['rate'][itemp]
                    #
                    for itrans in range(self.Nreclvl):
                        lvl1 = reclvl['lvl1'][itrans]-1
                        lvl2 = reclvl['lvl2'][itrans]-1
                        popmat[lvl2+ci, -1] += self.EDensity[itemp]*self.ReclvlRate['rate'][itrans, itemp]
                        popmat[-1, -1] -= self.EDensity[itemp]*self.ReclvlRate['rate'][itrans, itemp]
                # normalize to unity
                norm=np.ones(nlvls+ci+rec,'float64')
                if ci:
                    norm[0] = 0.
                if rec:
                    norm[-1] = 0.
                popmat[nlvls+ci+rec-1]=norm
                b=np.zeros(nlvls+ci+rec,'float64')
                b[nlvls+ci+rec-1]=1.
                try:
                    thispop=np.linalg.solve(popmat,b)
                    pop[itemp] = thispop[ci:ci+nlvls]
                except np.linalg.LinAlgError:
                    pop[itemp] = np.zeros(nlvls, 'float64')
#                    print ' error in matrix inversion, setting populations to zero at T = ', ('%8.2e')%(temperature[itemp])
#                thispop=np.linalg.solve(popmat,b)
#                if rec:
#                    pop[itemp] = thispop[ci:ci+nlvls+rec-1]
#                else:
#                    pop[itemp] = thispop[ci:]
#                pop[itemp] = thispop[ci:ci+nlvls]
            #
        pop=np.where(pop >0., pop,0.)
        self.Population={"temperature":temperature,"eDensity":eDensity,"population":pop, "protonDensity":protonDensity, "ci":ci, "rec":rec, 'popmat':popmat}
        #
        return
        #
        # -------------------------------------------------------------------------------------
        #
    def populateNew(self, popCorrect=1, verbose=0, **kwargs):
        """
        Calculate level populations for specified ion.  This is a new version that will enable the calculation
        of dielectronic satellite lines without resorting to the dielectronic ions, such as c_5d
        possible keyword arguments include temperature, eDensity, pDensity, radTemperature and rStar
        this is a developmental method for using the autoionizing A-values to determine level resolved
        dielectronic recombination rates
        """
        #
        #
        for one in kwargs.keys():
            if one not in chdata.keywordArgs:
                print ' following keyword is not understood - ',  one
        #
        nlvls = self.Nlvls
        nwgfa = self.Nwgfa
        nsplups = self.Nsplups
        npsplups = self.Npsplups
        #
        if kwargs.has_key('temperature'):
            self.Temperature = np.asarray(kwargs['temperature'])
            temperature = self.Temperature
        elif hasattr(self, 'Temperature'):
            temperature=self.Temperature
        else:
            print ' no temperature values have been set'
            return
        #
        if kwargs.has_key('eDensity'):
            self.EDensity = np.asarray(kwargs['eDensity'])
            eDensity = self.EDensity
        elif hasattr(self, 'EDensity'):
            eDensity = self.EDensity
        else:
            print ' no eDensity values have been set'
            return
        #
        if kwargs.has_key('pDensity'):
            if kwargs['pDensity'] == 'default':
                self.p2eRatio()
                protonDensity = self.ProtonDensityRatio*self.EDensity
            else:
                try:
                    self.PDensity = np.asarray(kwargs['pDensity'])
                except:
                    print ' could not interpret value for keyword pDensity'
                    print ' should be either "default" or a number or array'
                    return
        else:
            if hasattr(self, 'PDensity'):
                protonDensity = self.PDensity
            else:
                self.p2eRatio()
                self.PDensity = self.ProtonDensityRatio*self.EDensity
                protonDensity = self.PDensity
                print ' proton density not specified, set to "default"'
        #
        if 'radTemperature' in kwargs.keys() and 'rStar' in kwargs.keys():
            self.RadTemperature = np.asarray(kwargs['radTemperature'])
            radTemperature = np.array(self.RadTemperature)
            self.RStar = np.asarray(kwargs['rStar'])
            rStar = np.asarray(self.RStar)
        elif hasattr(self, 'RadTemperature') and hasattr(self, 'RStar'):
            radTemperature = self.RadTemperature
            rStar = self.RStar
        #
        #
        #
        if self.Ncilvl:
            ci = 1
            cilvl = self.Cilvl
            if hasattr(self, 'CilvlRate'):
                cilvlRate = self.CilvlRate
            else:
                self.cireclvlDescale('cilvl')
                cilvlRate = self.CilvlRate
            self.recombRate()
            #
            lowers = util.zion2name(self.Z, self.Ion-1)
            # get the lower ionization stage
            self.Lower = ion(lowers, temperature=self.Temperature, eDensity = self.EDensity)
            self.Lower.ionizRate()
            # need to get multiplicity of lower ionization stage
            lowMult = self.Lower.Elvlc['mult']
        else:
            ci = 0
        #  evetually will be looking for just an rrLvl attribute
        if self.Nrrlvl:
            rec = 1
            rrlvl = self.Rrlvl
            if hasattr(self, 'RrlvlRate'):
                rrlvlRate = self.RrlvlRate
            else:
#                    print ' doing reclvlDescale in populate'
                self.cireclvlDescale('rrlvl')
                rrlvlRate = self.RrlvlRate
        elif hasattr(self, 'Auto'):
            # will calculate the dielectronic rates below
            rec = 1
        else:
            rec = 0
        #
        if rec:
            if hasattr(self, 'Auto'):
                self.drRateLvl()
            # get ionization rate of this ion
            self.ionizRate()
            #  get the higher ionization stage
            highers = util.zion2name(self.Z, self.Ion+1)
            self.Higher = ion(highers, temperature=self.Temperature, eDensity=self.EDensity)
            self.Higher.recombRate()
#        print ' nlvls, ci, rec = ', nlvls, ci, rec
        #
        rad=np.zeros((nlvls+ci+rec,nlvls+ci+rec),"float64")  #  the populating matrix for radiative transitions
        #
        #
        for iwgfa in range(nwgfa):
            l1 = self.Wgfa["lvl1"][iwgfa]-1
            l2 = self.Wgfa["lvl2"][iwgfa]-1
            rad[l1+ci,l2+ci] += self.Wgfa["avalue"][iwgfa]
            rad[l2+ci,l2+ci] -= self.Wgfa["avalue"][iwgfa]
            # photo-excitation and stimulated emission
            if self.RadTemperature:
                if not self.RStar:
                    dilute = 0.5
                else:
                    dilute = util.dilute(self.RStar)
                # next - don't include autoionization lines
                if abs(self.Wgfa['wvl'][iwgfa]) > 0.:
                    de = const.invCm2Erg*(self.Elvlc['ecm'][l2] - self.Elvlc['ecm'][l1])
                    dekt = de/(const.boltzmann*self.RadTemperature)
                    # photoexcitation
                    phexFactor = dilute*(self.Elvlc['mult'][l2]/self.Elvlc['mult'][l1])/(np.exp(dekt) -1.)
                    rad[l2+ci,l1+ci] += self.Wgfa["avalue"][iwgfa]*phexFactor
                    rad[l1+ci,l1+ci] -= self.Wgfa["avalue"][iwgfa]*phexFactor
                    # stimulated emission
                    stemFactor = dilute/(1. - np.exp(-dekt))
                    rad[l1+ci,l2+ci] += self.Wgfa["avalue"][iwgfa]*stemFactor
                    rad[l2+ci,l2+ci] -= self.Wgfa["avalue"][iwgfa]*stemFactor
        if hasattr(self, 'Auto'):
            for iauto, avalue in enumerate(self.Auto['avalue']):
                l1 = self.Auto["lvl1"][iauto] -1
                l2 = self.Auto["lvl2"][iauto] -1
                rad[l1 + ci + nlvls, l2 + ci] += avalue
                rad[l2 + ci,  l2 + ci] -= avalue

        #
        self.rad=rad
        #
        if self.Nsplups:
            self.upsilonDescale()
            ups = self.Upsilon['upsilon']
            exRate = self.Upsilon['exRate']
            dexRate = self.Upsilon['dexRate']
        #
        if npsplups:
            self.upsilonDescale(prot=1)
#            pups = self.PUpsilon['upsilon']
            pexRate = self.PUpsilon['exRate']
            pdexRate = self.PUpsilon['dexRate']
        #
        temp=temperature
        ntemp=temp.size
        #
        cc=const.collision*self.EDensity
        ndens=cc.size
        if npsplups:
            cp=const.collision*protonDensity
        if ntemp > 1 and ndens >1 and ntemp != ndens:
            print ' unless temperature or eDensity are single values'
            print ' the number of temperatures values must match the '
            print ' the number of eDensity values'
            return
        #
        # get corrections for recombination and excitation
        #
        #
        #  first, for ntemp=ndens=1
        if ndens==1 and ntemp==1:
            popmat=np.copy(rad)
            for isplups in range(0,nsplups):
                l1=self.Splups["lvl1"][isplups]-1
                l2=self.Splups["lvl2"][isplups]-1
                #
                popmat[l1+ci,l2+ci] += self.EDensity*dexRate[isplups]
                popmat[l2+ci,l1+ci] += self.EDensity*exRate[isplups]
                popmat[l1+ci,l1+ci] -= self.EDensity*exRate[isplups]
                popmat[l2+ci,l2+ci] -= self.EDensity*dexRate[isplups]
                #
            for isplups in range(0,npsplups):
                l1=self.Psplups["lvl1"][isplups]-1
                l2=self.Psplups["lvl2"][isplups]-1
                 #
                popmat[l1+ci,l2+ci] += self.PDensity*pdexRate[isplups]
                popmat[l2+ci,l1+ci] += self.PDensity*pexRate[isplups]
                popmat[l1+ci,l1+ci] -= self.PDensity*pexRate[isplups]
                popmat[l2+ci,l2+ci] -= self.PDensity*pdexRate[isplups]
           # now include ionization rate from
            if ci:
#                print ' ci = ', ci
                #
                # the ciRate can be computed for all temperatures
                #
                ciTot = 0.
                for itrans in range(len(cilvl['lvl1'])):
                    lvl1 = cilvl['lvl1'][itrans]-1
                    lvl2 = cilvl['lvl2'][itrans]-1
                    # this is kind of double booking the ionization rate components
                    popmat[lvl2+ci, lvl1] += self.EDensity*self.CilvlRate['rate'][itrans]
                    popmat[lvl1, lvl1] -= self.EDensity*self.CilvlRate['rate'][itrans]
                    ciTot += self.EDensity*self.CilvlRate['rate'][itrans]
                #
                popmat[1, 0] += (self.EDensity*self.Lower.IonizRate['rate'] - ciTot)
                popmat[0, 0] -= (self.EDensity*self.Lower.IonizRate['rate'] - ciTot)
                popmat[0, 1] += self.EDensity*self.RecombRate['rate']
                popmat[1, 1] -= self.EDensity*self.RecombRate['rate']
            if rec:
                #
                if hasattr(self, 'DrRateLvl'):
#                    branch = np.zeros(self.Ndielsplups, 'float64')
                    for idr, rate in enumerate(self.DrRateLvl['rate']):
                        l1 = self.Auto["lvl1"][idr] - 1
                        l2 = self.Auto["lvl2"][idr] - 1
                        popmat[l2+ci,l1+ci+nlvls] += self.EDensity*self.DrRateLvl['rate'][idr]
                        popmat[l1+ci,l1+ci] -= self.EDensity*self.DrRateLvl['rate'][idr]
                        #
                    dielTot = self.EDensity*self.DrRateLvl['totalRate']
                else:
                    dielTot = 0.
                #
#                print ' rec, dielTot  = ', rec,  dielTot
                #
                for itrans in range(self.Nrrlvl):
#                    lvl1 = reclvl['lvl1'][itrans]-1
                    lvl2 = rrlvl['lvl2'][itrans]-1
                    popmat[lvl2+ci, -1] += self.EDensity*rrlvlRate['rate'][itrans]
                    popmat[-1, -1] -= self.EDensity*rrlvlRate['rate'][itrans]
                if self.Nrrlvl:
                    rrlvlRateTot = rrlvlRate['rate'].sum(axis=0)
                else:
                    rrlvlRateTot = 0.

                #
                popmat[-1,  ci] += self.EDensity*self.IonizRate['rate']
                popmat[ci, ci] -= self.EDensity*self.IonizRate['rate']
                # next 2 line take care of overbooking
                netRecomb = self.EDensity*(self.Higher.RecombRate['rate']- rrlvlRateTot) - dielTot
                if netRecomb > 0.:
                    popmat[ci, -1] += netRecomb
                    popmat[-1, -1] -= netRecomb
#                print ' higher, rec , dieltot = ',  self.EDensity*higher.RecombRate['rate'], self.EDensity*rrlvlRate['rate'].sum(axis=0),  dielTot
            # normalize to unity
            norm=np.ones(nlvls+ci+rec,'float64')
            if ci:
                norm[0] = 0.
            if rec:
                norm[nlvls+ci+rec-1] = 0.
            popmat[nlvls+ci+rec-1]=norm
            b=np.zeros(nlvls+ci+rec,'float64')
            b[nlvls+ci+rec-1]=1.
#            print ' norm = ', norm
#            print 'popmat, last line',  popmat[-1]
#            print ' b = ', b
#            popmat[nlvls/2]=norm
#            b=np.zeros(nlvls+ci+rec,'float64')
#            b[nlvls/2]=1.
#            if rec:
#                fullpop = np.linalg.solve(popmat,b)
#                pop = fullpop[ci:ci+nlvls+rec-1]
#            else:
#                fullpop = np.linalg.solve(popmat,b)
#                pop = fullpop[ci:]
#            fullpop = np.linalg.solve(popmat,b)
            try:
                fullpop=np.linalg.solve(popmat,b)
                pop = fullpop[ci:ci+nlvls]
            except np.linalg.LinAlgError:
                pop = np.zeros(nlvls, 'float64')
#                print ' error in matrix inversion, setting populations to zero at T = ', ('%8.2e')%(temperature)
            #
        #   next, in case of a single eDensity value
#            pop = np.linalg.solve(popmat,b)
        #
        # ------------------------------------------------------------------------
        elif ndens == 1:
            pop=np.zeros((ntemp, nlvls),"float64")
#            pop=np.zeros((ntemp,ci + nlvls + rec),"float64")
            for itemp in range(ntemp):
                popmat=np.copy(rad)
                for isplups in range(0,nsplups):
                    l1=self.Splups["lvl1"][isplups]-1
                    l2=self.Splups["lvl2"][isplups]-1
                    popmat[l1+ci,l2+ci] += self.EDensity*dexRate[isplups, itemp]
                    popmat[l2+ci,l1+ci] += self.EDensity*exRate[isplups, itemp]
                    popmat[l1+ci,l1+ci] -= self.EDensity*exRate[isplups, itemp]
                    popmat[l2+ci,l2+ci] -= self.EDensity*dexRate[isplups, itemp]
                for isplups in range(0,npsplups):
                    l1=self.Psplups["lvl1"][isplups]-1
                    l2=self.Psplups["lvl2"][isplups]-1
                    # for proton excitation, the levels are all below the ionization potential
                     #
                    popmat[l1+ci,l2+ci] += self.PDensity[itemp]*pdexRate[isplups, itemp]
                    popmat[l2+ci,l1+ci] += self.PDensity[itemp]*pexRate[isplups, itemp]
                    popmat[l1+ci,l1+ci] -= self.PDensity[itemp]*pexRate[isplups, itemp]
                    popmat[l2+ci,l2+ci] -= self.PDensity[itemp]*pdexRate[isplups, itemp]
                # now include ionization rate from
                if ci:
#                    print ' ci = ', ci
                    #
                    # the ciRate can be computed for all temperatures
                    #
                    ciTot = 0.
                    for itrans in range(len(cilvl['lvl1'])):
                        lvl1 = cilvl['lvl1'][itrans]-1
                        lvl2 = cilvl['lvl2'][itrans]-1
                        mult = lowMult[lvl1-1]
                        popmat[lvl2+ci, lvl1] += self.EDensity*self.CilvlRate['rate'][itrans, itemp]
                        popmat[lvl1, lvl1] -= self.EDensity*self.CilvlRate['rate'][itrans, itemp]
                        ciTot += self.EDensity*self.CilvlRate['rate'][itrans, itemp]
                    #
                    popmat[1, 0] += (self.EDensity*lower.IonizRate['rate'][itemp] - ciTot)
                    popmat[0, 0] -= (self.EDensity*lower.IonizRate['rate'][itemp] - ciTot)
                    popmat[0, 1] += self.EDensity*self.RecombRate['rate'][itemp]
                    popmat[1, 1] -= self.EDensity*self.RecombRate['rate'][itemp]
                if rec:
                #
                #
                    if hasattr(self, 'DrRateLvl'):
    #                    branch = np.zeros(self.Ndielsplups, 'float64')
                        for idr, rate in enumerate(self.DrRateLvl['rate']):
                            l1 = self.Auto["lvl1"][idr] - 1
                            l2 = self.Auto["lvl2"][idr] - 1
                            popmat[l2+ci,l1+ci+nlvls] += self.EDensity*self.DrRateLvl['rate'][idr][itemp]
                            popmat[l1+ci,l1+ci] -= self.EDensity*self.DrRateLvl['rate'][idr][itemp]
                            #
                        dielTot = self.EDensity*self.DrRateLvl['totalRate'][itemp]
                    else:
                        dielTot = 0.
                #
                    if self.Nrrlvl:
                        recTot = self.RrlvlRate['rate'][:, itemp].sum()
                    else:
                        recTot = 0.
                #
                    popmat[-1,  ci] += self.EDensity*self.IonizRate['rate'][itemp]
                    popmat[ci, ci] -= self.EDensity*self.IonizRate['rate'][itemp]
                    netRecomb = self.EDensity*(self.Higher.RecombRate['rate'][itemp]- recTot) - dielTot
                    if netRecomb < 0.:
                        popmat[ci, -1] += netRecomb
                        popmat[-1, -1] -= netRecomb
                    #
#                    for itrans in range(len(rrlvl['lvl1'])):
                    for itrans in range(self.Nrrlvl):
                        lvl1 = rrlvl['lvl1'][itrans]-1
                        lvl2 = rrlvl['lvl2'][itrans]-1
                        popmat[lvl2+ci, -1] += self.EDensity*self.RrlvlRate['rate'][itrans, itemp]
                        popmat[-1, -1] -= self.EDensity*self.RrlvlRate['rate'][itrans, itemp]
                # normalize to unity
                norm=np.ones(nlvls+ci+rec,'float64')
                if ci:
                    norm[0] = 0.
                if rec:
                    norm[-1] = 0.
                popmat[nlvls+ci+rec-1]=norm
                b=np.zeros(nlvls+ci+rec,'float64')
                b[nlvls+ci+rec-1]=1.
                try:
                    thispop=np.linalg.solve(popmat,b)
                    pop[itemp] = thispop[ci:ci+nlvls]
                except np.linalg.LinAlgError:
                    pop[itemp] = np.zeros(nlvls, 'float64')
#                    print ' error in matrix inversion, setting populations to zero at T = ', ('%8.2e')%(temperature[itemp])
            #
        elif ntemp == 1:
#            pop=np.zeros((ndens,nlvls),"float64")
            pop=np.zeros((ndens,nlvls),"float64")
            for idens in range(0,ndens):
                popmat=np.copy(rad)
                for isplups in range(0,nsplups):
                    l1=self.Splups["lvl1"][isplups]-1
                    l2=self.Splups["lvl2"][isplups]-1
                #
                    popmat[l1+ci,l2+ci] += self.EDensity[idens]*dexRate[isplups]
                    popmat[l2+ci,l1+ci] += self.EDensity[idens]*exRate[isplups]
                    popmat[l1+ci,l1+ci] -= self.EDensity[idens]*exRate[isplups]
                    popmat[l2+ci,l2+ci] -= self.EDensity[idens]*dexRate[isplups]
                #
                for isplups in range(0,npsplups):
                    l1=self.Psplups["lvl1"][isplups]-1
                    l2=self.Psplups["lvl2"][isplups]-1
                 #
                    popmat[l1+ci,l2+ci] += self.PDensity[idens]*pdexRate[isplups]
                    popmat[l2+ci,l1+ci] += self.PDensity[idens]*pexRate[isplups]
                    popmat[l1+ci,l1+ci] -= self.PDensity[idens]*pexRate[isplups]
                    popmat[l2+ci,l2+ci] -= self.PDensity[idens]*pdexRate[isplups]
                # now include ionization rate from
                if ci:
#                    print ' ci = ', ci
                    #
                    #
                    ciTot = 0.
                    for itrans in range(len(cilvl['lvl1'])):
                        lvl1 = cilvl['lvl1'][itrans] -1
                        lvl2 = cilvl['lvl2'][itrans] -1
                        popmat[lvl2+ci, lvl1] += self.EDensity[idens]*self.CilvlRate['rate'][itrans]
                        popmat[lvl1, lvl1] -= self.EDensity[idens]*self.CilvlRate['rate'][itrans]
                        ciTot += self.EDensity[idens]*self.CilvlRate['rate'][itrans]
                    popmat[1, 0] += (self.EDensity[idens]*lower.IonizRate['rate'] -ciTot)
                    popmat[0, 0] -= (self.EDensity[idens]*lower.IonizRate['rate'] -ciTot)
                    popmat[0, 1] += self.EDensity[idens]*self.RecombRate['rate']
                    popmat[1, 1] -= self.EDensity[idens]*self.RecombRate['rate']
                if rec:
#                    print ' rec = ', rec
                    if hasattr(self, 'DrRateLvl'):
    #                    branch = np.zeros(self.Ndielsplups, 'float64')
                        for idr, rate in enumerate(self.DrRateLvl['rate']):
                            l1 = self.Auto["lvl1"][idr] - 1
                            l2 = self.Auto["lvl2"][idr] - 1
                            popmat[l2+ci,l1+ci+nlvls] += self.EDensity[idens]*self.DrRateLvl['rate'][idr]
                            popmat[l1+ci,l1+ci] -= self.EDensity[idens]*self.DrRateLvl['rate'][idr]
                            #
                        dielTot = self.EDensity[idens]*self.DrRateLvl['totalRate']
                    else:
                        dielTot = 0.
                    if self.Nrrlvl:
#                        print ' RrlvlRate.shape = ', self.RrlvlRate['rate'].shape
                        recTot = self.RrlvlRate['rate'].sum()
                    else:
                        recTot = 0.
                    #
                    popmat[-1,  ci] += self.EDensity[idens]*self.IonizRate['rate']
                    popmat[ci, ci] -= self.EDensity[idens]*self.IonizRate['rate']
                    netRec = self.EDensity[idens]*(self.Higher.RecombRate['rate'] - recTot) - dielTot
                    if netRec > 0.:
                        popmat[ci, -1] += netRec
                        popmat[-1, -1] -= netRec
                    #
#                    for itrans in range(len(rrlvl['lvl1'])):
                    for itrans in range(self.Nrrlvl):
                        lvl1 = rrlvl['lvl1'][itrans]-1
                        lvl2 = rrlvl['lvl2'][itrans]-1
                        popmat[lvl2+ci, -1] += self.EDensity[idens]*self.RrlvlRate['rate'][itrans]
                        popmat[-1, -1] -= self.EDensity[idens]*self.RrlvlRate['rate'][itrans]
                # normalize to unity
                norm=np.ones(nlvls+ci+rec,'float64')
                if ci:
                    norm[0] = 0.
                if rec:
                    norm[-1] = 0.
                popmat[nlvls+ci+rec-1]=norm
                b=np.zeros(nlvls+ci+rec,'float64')
                b[nlvls+ci+rec-1]=1.
                try:
                    thispop=np.linalg.solve(popmat,b)
                    pop[idens] = thispop[ci:ci+nlvls]
                except np.linalg.LinAlgError:
                    pop[idens] = np.zeros(nlvls, 'float64')
#                    print ' error in matrix inversion, setting populations to zero at eDensity = ', ('%8.2e')%(eDensity[idens])
#                thispop=np.linalg.solve(popmat,b)
#                if rec:
#                    pop[idens] = thispop[ci:ci+nlvls+rec-1]
#                else:
#                    pop[idens] = thispop[ci:]
#                pop[idens] = thispop[ci:ci+nlvls]
                #
        elif ntemp>1  and ntemp==ndens:
            pop=np.zeros((ntemp,nlvls),"float64")
#            pop=np.zeros((ntemp,ci+nlvls+rec),"float64")
            for itemp in range(0,ntemp):
                temp=self.Temperature[itemp]
                popmat=np.copy(rad)
                for isplups in range(0,nsplups):
                    l1=self.Splups["lvl1"][isplups]-1
                    l2=self.Splups["lvl2"][isplups]-1
                    #
                    popmat[l1+ci,l2+ci] += self.EDensity[itemp]*dexRate[isplups, itemp]
                    popmat[l2+ci,l1+ci] += self.EDensity[itemp]*exRate[isplups, itemp]
                    popmat[l1+ci,l1+ci] -= self.EDensity[itemp]*exRate[isplups, itemp]
                    popmat[l2+ci,l2+ci] -= self.EDensity[itemp]*dexRate[isplups, itemp]
                # proton rates
                for isplups in range(0,npsplups):
                    l1=self.Psplups["lvl1"][isplups]-1
                    l2=self.Psplups["lvl2"][isplups]-1
                    # for proton excitation, the levels are all below the ionization potential
                     #
                    popmat[l1+ci,l2+ci] += self.PDensity[itemp]*pdexRate[isplups, itemp]
                    popmat[l2+ci,l1+ci] += self.PDensity[itemp]*pexRate[isplups, itemp]
                    popmat[l1+ci,l1+ci] -= self.PDensity[itemp]*pexRate[isplups, itemp]
                    popmat[l2+ci,l2+ci] -= self.PDensity[itemp]*pdexRate[isplups, itemp]
                # now include ionization rate from
                if ci:
#                   print ' ci = ', ci
                    #
                    # the ciRate can be computed for all temperatures
                    #
                    ciTot = 0.
                    for itrans in range(len(cilvl['lvl1'])):
                        lvl1 = cilvl['lvl1'][itrans] -1
                        lvl2 = cilvl['lvl2'][itrans] -1
                        popmat[lvl2+ci, lvl1] += self.EDensity[itemp]*self.CilvlRate['rate'][itrans, itemp]
                        popmat[lvl1, lvl1] -= self.EDensity[itemp]*self.CilvlRate['rate'][itrans, itemp]
                        ciTot += self.EDensity[itemp]*self.CilvlRAte['rate'][itrans, itemp]
                    popmat[1, 0] += (self.EDensity[itemp]*lower.IonizRate['rate'][itemp] - ciTot)
                    popmat[0, 0] -= (self.EDensity[itemp]*lower.IonizRate['rate'][itemp] - ciTot)
                    popmat[0, 1] += self.EDensity[itemp]*self.RecombRate['rate'][itemp]
                    popmat[1, 1] -= self.EDensity[itemp]*self.RecombRate['rate'][itemp]
                if rec:
                #
                    if hasattr(self, 'DrRateLvl'):
    #                    branch = np.zeros(self.Ndielsplups, 'float64')
                        for idr, rate in enumerate(self.DrRateLvl['rate']):
                            l1 = self.Auto["lvl1"][idr] - 1
                            l2 = self.Auto["lvl2"][idr] - 1
                            popmat[l2+ci,l1+ci+nlvls] += self.EDensity[itemp]*self.DrRateLvl['rate'][idr][itemp]
                            popmat[l1+ci,l1+ci] -= self.EDensity[itemp]*self.DrRateLvl['rate'][idr][itemp]
                            #
                        dielTot = self.EDensity[itemp]*self.DrRateLvl['totalRate'][itemp]
                    else:
                        dielTot = 0.
                #
                    if self.Nrrlvl:
                        recTot = self.RrlvlRate['rate'][:, itemp].sum()
                    else:
                        recTot = 0.
                #
#                   print ' rec = ', rec
                    popmat[-1,  ci] += self.EDensity[itemp]*self.IonizRate['rate'][itemp]
                    popmat[ci, ci] -= self.EDensity[itemp]*self.IonizRate['rate'][itemp]
                    netRec = self.EDensity[itemp]*(self.Higher.RecombRate['rate'][itemp] - recTot) - dielTot
                    if netRec > 0.:
                        popmat[ci, -1] += netRec
                        popmat[-1, -1] -= netRec
                    #
                    for itrans in range(self.Nrrlvl):
                        lvl1 = rrlvl['lvl1'][itrans]-1
                        lvl2 = rrlvl['lvl2'][itrans]-1
                        popmat[lvl2+ci, -1] += self.EDensity[itemp]*self.RrlvlRate['rate'][itrans, itemp]
                        popmat[-1, -1] -= self.EDensity[itemp]*self.RrlvlRate['rate'][itrans, itemp]
                # normalize to unity
                norm=np.ones(nlvls+ci+rec,'float64')
                if ci:
                    norm[0] = 0.
                if rec:
                    norm[-1] = 0.
                popmat[nlvls+ci+rec-1]=norm
                b=np.zeros(nlvls+ci+rec,'float64')
                b[nlvls+ci+rec-1]=1.
                try:
                    thispop=np.linalg.solve(popmat,b)
                    pop[itemp] = thispop[ci:ci+nlvls]
                except np.linalg.LinAlgError:
                    pop[itemp] = np.zeros(nlvls, 'float64')
#                    print ' error in matrix inversion, setting populations to zero at T = ', ('%8.2e')%(temperature[itemp])
#                thispop=np.linalg.solve(popmat,b)
#                if rec:
#                    pop[itemp] = thispop[ci:ci+nlvls+rec-1]
#                else:
#                    pop[itemp] = thispop[ci:]
#                pop[itemp] = thispop[ci:ci+nlvls]
            #
        pop=np.where(pop >0., pop,0.)
        self.Population={"temperature":temperature,"eDensity":eDensity,"population":pop, "protonDensity":protonDensity, "ci":ci, "rec":rec, 'popmat':popmat}
        #
        return
        #
        # -------------------------------------------------------------------------------------
        #
    def popPlot(self,top=10, plotFile=0, saveFile=0):
        """Plots populations vs temperature or eDensity.

        top specifies the number of the most highly populated levels to plot."""
        #self.Population={"temperature":temperature,"eDensity":eDensity,"population":pop}
        fontsize=14
        #
        if hasattr(self, 'Population'):
            temperature=self.Population["temperature"]
            eDensity=self.Population["eDensity"]
            pop=self.Population["population"]
        else:
            self.populate()
            temperature=self.Population["temperature"]
            eDensity=self.Population["eDensity"]
            pop=self.Population["population"]
        #
        #  for case of only a single density and temperature
        if len(pop.shape) == 1:
            spop = np.sort(pop)
            idx = np.argsort(pop)
            minPop = spop[-top:].min()/2.
            if top > pop.size:
                top = pop.size
            for itop in range(1, top+1):
                x = [idx[-itop], idx[-itop], idx[-itop]+1, idx[-itop]+1]
                y = [minPop, spop[-itop], spop[-itop], minPop]
                pl.semilogy(x, y, 'k')
            pl.axis([0, max(idx[-top:])+1, minPop, 1.])
            pl.xlabel('Level', fontsize=fontsize)
            pl.ylabel('Population', fontsize=fontsize)
            return
        #
        # find the top most populated levels
        #
        lvl=self.Elvlc["lvl"]
#        nlvls=len(lvl)
        nlvls = self.Nlvls
        if top > nlvls:
            top = nlvls
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
        ndens=eDensity.size
        #
        ylabel='Population'
        title=self.Spectroscopic
        #
        pl.figure()
        #
        pl.ion()
        #
        #
        #  first, for ntemp=ndens=1
        if ndens==1 and ntemp==1:
            print ' only a single temperature and eDensity'
            self.Message = ' only a single temperature and eDensity'
#            if chInteractive:
#                print ' only a single temperature and eDensity'
#            else:
#                self.Message = ' only a single temperature and eDensity'
            return
        elif ndens == 1:
            toppops = np.zeros((top, ntemp), 'float64')
            for ilvl in range(top):
                toppops[ilvl] = pop[:, toplvl[ilvl]-1]
            nonzero = toppops > 0.
            ymin = min(toppops[nonzero])
            for lvl in toplvl:
                # for some low temperature, populations can not be calculated
                good = pop[:, lvl-1] > 0
                pl.loglog(temperature[good],pop[good,lvl-1])
                skip=3
                if good.sum() == ntemp:
                    start=divmod(lvl,ntemp)[1]
                    for itemp in range(start,ntemp,ntemp/skip):
                        pl.text(temperature[itemp],pop[itemp,lvl-1],str(lvl))
                else:
                    newtemp=[]
                    for i, one in enumerate(temperature):
                        if good[i]:
                            newtemp.append(one)
                    start = divmod(lvl, len(newtemp))[1] + ntemp - good.sum()
                    for itemp in range(start,ntemp,ntemp/skip):
                        pl.text(temperature[itemp],pop[itemp,lvl-1],str(lvl))
            xlabel='Temperature (K)'
            pl.xlabel(xlabel,fontsize=fontsize)
            pl.ylabel(ylabel,fontsize=fontsize)
            dstr=' -  Density = %10.2e (cm$^{-3}$)' % eDensity
            pl.title(title+dstr,fontsize=fontsize)
            pl.xlim(temperature.min(),temperature.max())
#            nonzero = pop
#            yl=pl.ylim()
            pl.ylim(ymin,1.2)
        elif ntemp == 1:
            xlabel=r'Electron Density (cm$^{-3}$)'
#            for lvl in toplvl:
#                pl.loglog(eDensity,pop[:,lvl-1])
#                skip=min(3, ndens)
#                start=divmod(lvl,ndens)[1]
#                for idens in range(start,ndens,ndens/skip):
#                    pl.text(eDensity[idens],pop[idens,lvl-1],str(lvl))
            toppops = np.zeros((top, ndens), 'float64')
            for ilvl in range(top):
                toppops[ilvl] = pop[:, toplvl[ilvl]-1]
            nonzero = toppops > 0.
            ymin = min(toppops[nonzero])
            for lvl in toplvl:
                # for some low temperature, populations can not be calculated
                good = pop[:, lvl-1] > 0
                pl.loglog(eDensity[good],pop[good,lvl-1])
                skip=3
                if good.sum() == ndens:
                    start=divmod(lvl,ndens)[1]
                    for idens in range(start,ndens,ndens/skip):
                        pl.text(eDensity[idens],pop[idens,lvl-1],str(lvl))
                else:
                    newdens=[]
                    for i, one in enumerate(eDensity):
                        if good[i]:
                            newdens.append(one)
                    start = divmod(lvl, len(newdens))[1] + ndens - good.sum()
                    for idens in range(start,ndens,ndens/skip):
                        pl.text(eDensity[idens],pop[idens, lvl-1],str(lvl))
            pl.xlabel(xlabel,fontsize=fontsize)
            pl.ylabel(ylabel,fontsize=fontsize)
            tstr=' -  T = %10.2e (K)' % temperature
            pl.title(title+tstr,fontsize=fontsize)
            pl.xlim(eDensity[eDensity.nonzero()].min(),eDensity.max())
            yl=pl.ylim()
            pl.ylim(yl[0],1.2)
        else:
#            pl.figure()
            ax = pl.subplot(111)
#            for lvl in toplvl:
#                pl.loglog(temperature,pop[:,lvl-1])
#                skip = min(3, ntemp)
#                start=divmod(lvl,ntemp)[1]
#                for itemp in range(start,ntemp,ntemp/skip):
#                    pl.text(temperature[itemp],pop[itemp,lvl-1],str(lvl))
            toppops = np.zeros((top, ntemp), 'float64')
            for ilvl in range(top):
                toppops[ilvl] = pop[:, toplvl[ilvl]-1]
            nonzero = toppops > 0.
            ymin = min(toppops[nonzero])
            for lvl in toplvl:
                # for some low temperature, populations can not be calculated
                good = pop[:, lvl-1] > 0
                pl.loglog(temperature[good],pop[good,lvl-1])
                skip=3
                if good.sum() == ntemp:
                    start=divmod(lvl,ntemp)[1]
                    for itemp in range(start,ntemp,ntemp/skip):
                        pl.text(temperature[itemp],pop[itemp,lvl-1],str(lvl))
                else:
                    newtemp=[]
                    for i, one in enumerate(temperature):
                        if good[i]:
                            newtemp.append(one)
                    start = divmod(lvl, len(newtemp))[1] + ntemp - good.sum()
                    for itemp in range(start,ntemp,ntemp/skip):
                        pl.text(temperature[itemp],pop[itemp,lvl-1],str(lvl))
            xlabel='Temperature (K)'
            pl.xlabel(xlabel,fontsize=fontsize)
            pl.ylabel(ylabel,fontsize=fontsize)
#            pl.title(title,fontsize=fontsize)
#            pl.xlim(temperature.min(),temperature.max())
#            yl=pl.ylim()
#            pl.ylim(ymin,1.2)
            pl.axis([temperature.min(),temperature.max(), ymin, 1.2])
            pl.text(0.1, 0.5,title, horizontalalignment='center', verticalalignment='center', fontsize=fontsize,  transform = ax.transAxes)
            #
            ax2 = pl.twiny()
            xlabel=r'Electron Density (cm$^{-3}$)'
            pl.xlabel(xlabel, fontsize=fontsize)
            pl.loglog(eDensity,pop[:,toplvl[0]], visible=False)
            ax2.xaxis.tick_top()
#            pl.figure()
#            for lvl in toplvl:
#                pl.loglog(eDensity,pop[:,lvl-1])
#                skip = min(3, ntemp)
#                start=divmod(lvl,ndens)[1]
#                for idens in range(start,ndens,ndens/skip):
#                    pl.text(eDensity[idens],pop[idens,lvl-1],str(lvl))
#            xlabel=r'Electron Density (cm$^{-3}$)'
#            pl.xlabel(xlabel,fontsize=fontsize)
#            pl.ylabel(ylabel,fontsize=fontsize)
#            pl.title(title,fontsize=fontsize)
#            pl.xlim(eDensity.min(),eDensity.max())
#            yl=pl.ylim()
#            pl.ylim(yl[0],1.2)
        if plotFile:
            pl.savefig(plotFile)
        self.Population['toplvl'] = toplvl
        return
        #
        # -------------------------------------------------------------------------------------
        #
    def emiss(self, wvlRange = 0,  allLines=1):
        """Calculate the emissivities for lines of the specified ion.

        wvlRange can be set to limit the calculation to a particular wavelength range

        units:  ergs cm^-3 s^-1 str^-1

        Does not include elemental abundance or ionization fraction

        Wavelengths are sorted
        set allLines = 1 to include unidentified lines
        """
        #
        #
        if hasattr(self, 'Population'):
            doPopulate=False
            pop=self.Population['population']
        else:
            self.populate()
            pop=self.Population["population"]
        #
        #
#       nlvls=len(self.Elvlc['lvl'])
##        good=self.Wgfa['avalue'] > 0.
        # using [:] to make a copy things don't change elsewhere
        wvl = np.asarray(self.Wgfa["wvl"][:], 'float64')
        if allLines:
            wvl=np.abs(wvl)
        l1 = np.asarray(self.Wgfa['lvl1'][:], 'int64')
        l2 = np.asarray(self.Wgfa["lvl2"][:], 'int64')
        avalue = np.asarray(self.Wgfa["avalue"][:], 'float64')
        #
        # make sure there are lines in the wavelength range, if specified

        if wvlRange:
            realgood = util.between(wvl, wvlRange)
            l1 = l1[realgood]
            l2 = l2[realgood]
            wvl = wvl[realgood]
            avalue = avalue[realgood]
        #
        # two-photon decays have wvl=0 and nonzero avalues
#        zed = wvl.count(0.)
        nonzed = wvl != 0.
        wvl = wvl[nonzed]
        l1 = l1[nonzed]
        l2 = l2[nonzed]
        avalue = avalue[nonzed]
        nwvl=len(wvl)
        #
        if nwvl == 0:
            self.Emiss = {'errorMessage':self.Spectroscopic+' no lines in this wavelength range'}
            return
        #
        #
        try:
            ntempden,nlvls=pop.shape
            em=np.zeros((nwvl, ntempden),'Float64')
            if self.Temperature.size < ntempden:
                temperature = np.repeat(self.Temperature, ntempden)
            else:
                temperature = self.Temperature
            if self.EDensity.size < ntempden:
                eDensity = np.repeat(self.EDensity, ntempden)
            else:
                eDensity = self.EDensity
        except:
            nlvls=len(pop)
            ntempden=1
            em=np.zeros(nwvl,'Float32')
            eDensity = self.EDensity
            temperature = self.Temperature
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
                    em[iwvl, itempden] = factor[iwvl]*p*avalue[iwvl]/eDensity[itempden]
            if self.Defaults['wavelength'] == 'kev':
                wvl = const.ev2Ang/np.asarray(wvl)
            elif self.Defaults['wavelength'] == 'nm':
                wvl = wvl/10.
            em = em.take(wvl.argsort(),axis=0)
            wvl.sort()
            idx = np.argsort(wvl)
            l1 = l1[idx]
            l2 = l2[idx]
        else:
            for iwvl in range(0,nwvl):
                p=pop[l2[iwvl]-1]
                em[iwvl]=factor[iwvl]*p*avalue[iwvl]/eDensity
            if self.Defaults['wavelength'] == 'kev':
                wvlE=const.ev2Ang/np.asarray(wvl)
            elif self.Defaults['wavelength'] == 'nm':
                wvl=wvl/10.
            idx = np.argsort(wvl)
            wvl = wvl[idx]
            em = em[idx]
            l1 = l1[idx]
            l2 = l2[idx]
        lvl1 = l1.tolist()
        lvl2 = l2.tolist()
        self.Emiss = {"wvl":wvl, "emiss":em, "plotLabels":plotLabels, 'lvl1':lvl1, 'lvl2':lvl2}
        return
        #
        # ---------------------------------------------------------------------------
        #
    def emissPlot(self, index=None,  wvlRange=None,  top=10, linLog='lin', relative=0,  verbose=0, plotFile = 0, saveFile=0 ):
        '''Plot the emissivities.

        wvlRange, a 2 element tuple, list or array determines the wavelength range

        Top specifies to plot only the top strongest lines, default = 10

        linLog specifies a linear or log plot, want either lin or log, default = lin

        normalize = 1 specifies whether to normalize to strongest line, default = 0'''
        #
        title=self.Spectroscopic
        #
        doEmiss=False
        if hasattr(self, 'Emiss'):
            em = self.Emiss
        else:
            try:
                self.emiss()
                em = self.Emiss
            except:
                print ' emissivities not calculated and emiss() is unable to calculate them'
                print ' perhaps the temperature and/or eDensity are not set'
                return
        emiss = em['emiss']
        wvl = em['wvl']
        temperature = self.Temperature
        eDensity = self.EDensity
        #
        ndens = eDensity.size
        ntemp = temperature.size
        #
        if ndens == 1 and ntemp == 1:
            dstr = ' -  Density = %10.2e (cm$^{-3}$)' %(eDensity)
            tstr = ' -  T = %10.2e (K)' %(temperature)
        elif ndens == 1 and ntemp > 1:
            if type(index) == types.NoneType:
                index = ntemp/2
                print 'using index = %5i specifying temperature =  %10.2e'%(index, temperature[index])
                self.Message = 'using index = %5i specifying temperature =  %10.2e'%(index, temperature[index])
#            if chInteractive:
#                print 'using index = %5i specifying temperature =  %10.2e'%(index, temperature[index])
#            else:
#                self.Message = 'using index = %5i specifying temperature =  %10.2e'%(index, temperature[index])
            emiss=emiss[:, index]
            dstr=' -  Density = %10.2e (cm$^{-3}$)' % eDensity
            tstr=' -  T = %10.2e (K)' % temperature[index]
        elif ndens > 1 and ntemp == 1:
            if type(index) == types.NoneType:
                index = ndens/2
                print 'using index =%5i specifying eDensity = %10.2e'%(index, eDensity[index])
                self.Message = 'using index =%5i specifying eDensity = %10.2e'%(index, eDensity[index])
#            if chInteractive:
#                print 'using index =%5i specifying eDensity = %10.2e'%(index, eDensity[index])
#            else:
#                self.Message = 'using index =%5i specifying eDensity = %10.2e'%(index, eDensity[index])
            emiss=emiss[:, index]
            dstr=' -  Density = %10.2e (cm$^{-3}$)' % eDensity[index]
            tstr=' -  T = %10.2e (K)' % temperature
        elif ndens > 1 and ntemp > 1:
            if type(index) == types.NoneType:
                index = ntemp/2
                print 'using index = %5i specifying temperature = %10.2e, eDensity =  %10.2e'%(index, temperature[index], eDensity[index])
                self.Message = 'using index = %5i specifying temperature = %10.2e, eDensity =  %10.2e'%(index, temperature[index], eDensity[index])
#             if chInteractive:
#                print 'using index = %5i specifying temperature = %10.2e, eDensity =  %10.2e'%(index, temperature[index], eDensity[index])
#            else:
#                self.Message = 'using index = %5i specifying temperature = %10.2e, eDensity =  %10.2e'%(index, temperature[index], eDensity[index])
            emiss=emiss[:, index]
            dstr=' -  Density = %10.2e (cm$^{-3}$)' % eDensity[index]
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
            print 'No lines in this wavelength interval'
            self.Error = 1
            self.Message = 'No lines in this wavelength interval'
#            if chInteractive:
#                print 'No lines in this wavelength interval'
#            else:
#                self.Error = 1
#                self.Message = 'No lines in this wavelength interval'
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
        ymin = 10.**(np.log10(emiss.min()).round(0)-0.5 )
        #
        pl.ion()
#        if chInteractive:
#            pl.ion()
#        else:
#            pl.ioff()
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
        if plotFile:
            pl.savefig(plotFile)
        #
        idx = np.argsort(wvl)
        self.Emiss['wvlTop'] = wvl[idx]
        self.Emiss['emissTop'] = emiss[idx]
        #
        # ---------------------------------------------------------------------------
        #
    def intensity(self,  wvlRange = None,  allLines=1):
        """Calculate  the intensities for lines of the specified ion.

        wvlRange, a 2 element tuple, list or array determines the wavelength range

        units:  ergs cm^-3 s^-1 str^-1

        includes elemental abundance and ionization fraction."""
        # emiss ={"wvl":wvl, "emiss":em, "plotLabels":plotLabels}
        #
        self.emiss(wvlRange = wvlRange, allLines=allLines)
        emiss = self.Emiss
        if 'errorMessage'  in emiss.keys():
            self.Intensity = {'errorMessage': self.Spectroscopic+' no lines in this wavelength region'}
            return
        em = emiss['emiss']
        wvl = emiss['wvl']
        if hasattr(self, 'Abundance'):
            ab=self.Abundance
        else:
            self.Abundance = util.abundanceRead()
            ab=self.Abundance
        if hasattr(self, 'IoneqOne'):
            thisIoneq=self.IoneqOne
        else:
            self.ioneqOne()
            thisIoneq=self.IoneqOne
        try:
            nwvl, ntempden = em.shape
            intensity = np.zeros((ntempden, nwvl),'Float64')
            if thisIoneq.size == 1:
                thisIoneq = np.ones(ntempden, 'float64')*thisIoneq
            for it in range(ntempden):
                intensity[it] = ab*thisIoneq[it]*em[:, it]
        except:
            nwvl=len(em)
            ntempden=1
            intensity = ab*thisIoneq*em
        self.Intensity = {'intensity':intensity, 'wvl':wvl}
        #
        # -------------------------------------------------------------------------------------
        #
        #
        # ---------------------------------------------------------------------------
        #
    def boundBoundLoss(self,  wvlRange = None,  allLines=1):
        """Calculate  the summed radiative loss rate for all lines of the specified ion.

        wvlRange, a 2 element tuple, list or array determines the wavelength range

        units:  ergs cm^-3 s^-1

        includes elemental abundance and ionization fraction."""
        # emiss ={"wvl":wvl, "emiss":em, "plotLabels":plotLabels}
        #
        self.emiss(wvlRange = wvlRange, allLines=allLines)
        emiss = self.Emiss
        if 'errorMessage'  in emiss.keys():
            self.Intensity = {'errorMessage': self.Spectroscopic+' no lines in this wavelength region'}
            return
        em = emiss['emiss']
        wvl = emiss['wvl']
        if hasattr(self, 'Abundance'):
            ab=self.Abundance
        else:
            self.Abundance = util.abundanceRead()
            ab=self.Abundance
        if hasattr(self, 'IoneqOne'):
            thisIoneq=self.IoneqOne
        else:
            self.ioneqOne()
            thisIoneq=self.IoneqOne
        try:
            nwvl, ntempden = em.shape
            intensity = np.zeros((ntempden, nwvl),'Float64')
            if thisIoneq.size == 1:
                thisIoneq = np.ones(ntempden, 'float64')*thisIoneq
            for it in range(ntempden):
                if self.Defaults['flux'] != 'energy':
                    intensity[it] = 4.*const.pi*(const.planck*const.light*1.e+8/wvl)*ab*thisIoneq[it]*em[:, it]
                else:
                    intensity[it] = 4.*const.pi*ab*thisIoneq[it]*em[:, it]
            loss = intensity.sum(axis=1)
        except:
            nwvl=len(em)
            ntempden=1
            if self.Defaults['flux'] != 'energy':
                intensity = 4.*const.pi*(const.planck*const.light*1.e+8/wvl)*ab*thisIoneq*em
            else:
                intensity = 4.*const.pi*ab*thisIoneq*em
            loss = intensity.sum()
        self.BoundBoundLoss = {'rate':loss, 'wvlRange':wvlRange, 'temperature':self.Temperature, 'eDensity':self.EDensity}
        #
        # -------------------------------------------------------------------------------------
        #
    def intensityRatio(self,wvlRange=None,top=10):
        """Plot the ratio of 2 lines or sums of lines.

        Shown as a function of density and/or temperature.

        A plot of relative emissivities is shown and then a dialog appears for the user to

        choose a set of lines."""
        #
        #        self.Emiss={"temperature":temperature,"density":density,"wvl":wvl,"emiss":em,
        #        "plotLabels":plotLabels}
        #
        if hasattr(self, 'Emiss'):
            doEmiss=False
            em = self.Emiss
        else:
            doEmiss = True
        #
        #
        if doEmiss:
            # new values of temperature or eDensity
            self.emiss()
            em=self.Emiss
        #
        #
        fontsize=14
        #
        temperature = self.Temperature
        eDensity = self.EDensity
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
        ndens=self.EDensity.size
        #
        ylabel='Emissivity relative to '+maxWvl
        title=self.Spectroscopic
        #
        #
        if ndens==1 and ntemp==1:
            print ' only a single temperature and eDensity'
            return
        elif ndens == 1:
            xlabel='Temperature (K)'
            xvalues=self.Temperature
            outTemperature=self.Temperature
            outDensity=np.zeros(ntemp,'Float32')
            outDensity.fill(self.EDensity)
            desc_str=' at  Density = %10.2e (cm)$^{-3}$' % self.EDensity
        elif ntemp == 1:
            xvalues=self.EDensity
            outTemperature=np.zeros(ndens,'Float32')
            outTemperature.fill(self.Temperature)
            outDensity=self.EDensity
            xlabel=r'$\rm{Electron Density (cm)^{-3}}$'
            desc_str=' at Temp = %10.2e (K)' % self.Temperature
        else:
            outTemperature=self.Temperature
            outDensity=self.EDensity
            xlabel='Temperature (K)'
            xvalues=self.Temperature
            desc_str=' for variable Density'
        #
        # put all actual plotting here
        #
        pl.ion()
#        if chInteractive:
#            pl.ion()
#        else:
#            pl.ioff()
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
            pl.loglog(eDensity,emiss[topLines[top-1]]/maxAll, visible=False)
            ax2.xaxis.tick_top()
            pl.ylim(ymin/1.2, 1.2*ymax)
        else:
            pl.ylim(ymin/1.2, 1.2*ymax)
            pl.title(title+desc_str,fontsize=fontsize)
        pl.draw()
        #  need time to let matplotlib finish plotting
        time.sleep(0.5)
        #
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
            pl.loglog(eDensity,numEmiss/denEmiss, visible=False)
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
                'temperature':outTemperature,'eDensity':outDensity,'filename':intensityRatioFileName, 'numIdx':num_idx, 'denIdx':den_idx}
        #
        # -------------------------------------------------------------------------------------
        #
    def intensityRatioSave(self,outFile=''):
        '''Save the intensity ratio to a file.

        The intensity ratio as a function to temperature and eDensity is saved to an asciii file.

        Descriptive information is included at the top of the file.'''
        if outFile == '':
            outfile=self.IntensityRatio['filename']
#            if chInteractive:
            print ' saving ratio to filename = ',outfile
        if hasattr(self, 'IntensityRatio'):
            temperature=self.IntensityRatio['temperature']
            eDensity=self.IntensityRatio['eDensity']
            ratio=self.IntensityRatio['ratio']
            out=open(outFile,'w')
            nvalues=len(ratio)
            #
            #  need to add 7 lines to maintain IDL like files
            #
            out.write(outFile+'\n')    #1
            out.write(self.IntensityRatio['desc']+'\n') #2
            out.write(' created with ChiantiPy version '+ chdata.__version__ +'\n')   #3
            out.write(' columns are temperature, eDensity, ratio'+'\n')  #5
            tunit = 'K'
            out.write(' temperature in '+tunit+', electron eDensity in cm^(-3)'+'\n')  #6
            out.write(' ratio given in '+self.Defaults['flux']+'\n')   #4
            out.write(' '+'\n') #7
            for ivalue in range(nvalues):
                s='%12.3e %12.3e  %12.3e ' % (temperature[ivalue],eDensity[ivalue],ratio[ivalue])
                out.write(s+os.linesep)
            out.close()
        else:
#            if chInteractive:
            print ' in .intensityRatioSave(), no IntensityRatio is found'
        #
        # -------------------------------------------------------------------------------------
        #
    def ioneqOne(self):
        '''Provide the ionization equilibrium for the selected ion as a function of temperature.
        returned in self.IoneqOne'''
        #
        if hasattr(self, 'Temperature'):
            temperature = self.Temperature
        else:
            return
        #
        if hasattr(self, 'IoneqAll'):
            ioneqAll = self.IoneqAll
        else:
            ioneqAll = util.ioneqRead(ioneqname = self.Defaults['ioneqfile'])
#            if chInteractive:
            self.ioneqAll=self.IoneqAll
        #
        ioneqTemperature = ioneqAll['ioneqTemperature']
        Z=self.Z
        Ion=self.Ion
        Dielectronic=self.Dielectronic
        ioneqOne = np.zeros_like(temperature)
        #
        thisIoneq=ioneqAll['ioneqAll'][Z-1,Ion-1 + Dielectronic].squeeze()
        del ioneqAll
#        thisIoneq = self.Ioneq
        gioneq=thisIoneq > 0.
        goodt1=self.Temperature >= ioneqTemperature[gioneq].min()
        goodt2=self.Temperature <= ioneqTemperature[gioneq].max()
        goodt=np.logical_and(goodt1,goodt2)
        y2=interpolate.splrep(np.log(ioneqTemperature[gioneq]),np.log(thisIoneq[gioneq]),s=0)
        #
        if goodt.sum() > 0:
            if self.Temperature.size > 1:
                gIoneq=interpolate.splev(np.log(self.Temperature[goodt]),y2)   #,der=0)
                ioneqOne[goodt] = np.exp(gIoneq)
            else:
                gIoneq=interpolate.splev(np.log(self.Temperature),y2)
                ioneqOne = np.exp(gIoneq)
        else:
            ioneqOne = 0.
        #
        self.IoneqOne = ioneqOne
        #
        # -------------------------------------------------------------------------------------
        #
    def gofnt(self,wvlRange=0,top=10, verbose=0):
        """Calculate the 'so-called' G(T) function.

        Given as a function of both temperature and eDensity.

        Only the top( set by 'top') brightest lines are plotted.
        the G(T) function is returned in a dictionary self.Gofnt"""
        #
        #self.emiss={"wvl":wvl,"emiss":em,"units":units,"plotLabels":plotLabels}
        #
        #
        if hasattr(self, 'Emiss'):
            doEmiss=False
            em=self.Emiss
        else:
            doEmiss=True
        #
        #
        if doEmiss:
            self.emiss()
            em=self.Emiss
        #
        #
        if hasattr(self, 'Abundance'):
            ab=self.Abundance
        else:
            self.Abundance = util.abundanceRead()
            ab=self.Abundance
        #
        fontsize=12
        #
        emiss=em["emiss"]
        wvl=em["wvl"]
        temperature=self.Temperature
        eDensity=self.EDensity
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
        ndens=self.EDensity.size
        #
        ylabel = 'Emissivity relative to '+maxWvl
        title = self.Spectroscopic
        #
        #
        if ndens==1 and ntemp==1:
            print ' only a single temperature and eDensity'
            return
        elif ndens == 1:
            xlabel='Temperature (K)'
            ngofnt = temperature.size
            xvalues=temperature
            outTemperature=temperature
            outDensity=np.zeros(ntemp,'Float32')
            outDensity.fill(eDensity)
            desc_str=' at Density = %10.2e' % eDensity
        elif ntemp == 1:
            xvalues=eDensity
            ngofnt = eDensity.size
            outTemperature=np.zeros(ndens,'Float32')
            outTemperature.fill(temperature)
            outDensity=eDensity
            xlabel=r'$\rm{Electron Density (cm}^{-3}\rm{)}$'
            desc_str=' at Temperature = %10.2e' % temperature
        else:
            outTemperature=temperature
            outDensity=eDensity
            xlabel='Temperature (K)'
            xvalues=temperature
            ngofnt = temperature.size
            desc_str=' for variable Density'
            #
        #
        # put all actual plotting here
        #
        pl.ion()
#        if chInteractive:
#            pl.ion()
#        else:
#            pl.ioff()
        #
        pl.figure()
        ax = pl.subplot(111)
        nxvalues=len(xvalues)
        #  maxAll is an array
#        print ' emiss = ', np.max(emiss[top-1]), np.max(emiss[0])
#        print ' maxAll = ', maxAll
#        ymax = np.max(1.2*emiss[top-1]/maxAll)
        ymax = 1.2
#        print ' ymax = ', ymax
        ymin = ymax  #  np.min(emiss[0]/maxAll)  # was originally  = ymax
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
            pl.loglog(eDensity,emiss[topLines[top-1]]/maxAll, visible=False)
            ax2.xaxis.tick_top()
        else:
            pl.ylim(ymin, ymax)
            pl.title(title+desc_str,fontsize=fontsize)
        pl.draw()
        #
        time.sleep(0.5)
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
        gIoneq=self.IoneqOne/eDensity
        #
        #
        #
        # plot the desired ratio
        pl.figure()
        g_line= topLines[gline_idx]#  [0]
        ##        print ' g_line = ',g_line
        #
        gofnt=np.zeros(ngofnt,'float64')
        for aline in g_line:
            gofnt+=gAbund*gIoneq*emiss[aline].squeeze()
        self.Gofnt={'temperature':outTemperature,'eDensity':outDensity,'gofnt':gofnt, 'index':g_line, 'wvl':wvl[g_line]}
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
            pl.loglog(eDensity,gofnt, visible=False)
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
                nTempDens = max(self.Temperature.size, self.EDensity.size)
            except:
                self.populate()
                pop = self.Population['population']
                nTempDens = max(self.Temperature.size, self.EDensity.size)
            if nTempDens > 1:
                emiss = np.zeros((nTempDens, nWvl), 'float64')
                if self.EDensity.size == 1:
                    eDensity = np.repeat(self.EDensity, nTempDens)
                else:
                    eDensity = self.EDensity
            else:
                emiss = np.zeros(nWvl, 'float64')
                eDensity = self.EDensity
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
                    emiss[goodWvl] = f*pop[l2]*distr/self.EDensity
                else:
                    for it in range(nTempDens):
                        emiss[it, goodWvl] = f*pop[it, l2]*distr/self.EDensity[it]
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
                    emiss[goodWvl] = f*pop[l2]*distr/self.EDensity
                else:
                    for it in range(nTempDens):
                        emiss[it, goodWvl] = f*pop[it, l2]*distr/self.EDensity[it]
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
                nTempDens = max(self.Temperature.size, self.EDensity.size)
            except:
                self.populate()
                pop = self.Population['population']
                nTempDens = max(self.Temperature.size, self.EDensity.size)
            if nTempDens > 1:
                rate = np.zeros((nTempDens, nWvl), 'float64')
                if self.EDensity.size == 1:
                    eDensity = np.repeat(self.EDensity, nTempDens)
                else:
                    eDensity = self.EDensity
                if self.Temperature.size == 1:
                    temperature = np.repeat(self.Temperature, nTempDens)
                    thisIoneq = np.repeat(thisIoneq, nTempDens)
                else:
                    temperature = self.Temperature
            else:
                rate = np.zeros(nWvl, 'float64')
                eDensity = self.EDensity
                temperature = self.Temperature
            if self.Z == self.Ion:
                # H seq
                l1 = 1-1
                l2 = 2 - 1
                wvl0 = 1.e+8/(self.Elvlc['ecm'][l2] - self.Elvlc['ecm'][l1])
                goodWvl = wvl > wvl0
                if goodWvl.sum() > 0:
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
                        rate[goodWvl] = f*pop[l2]*distr*ab*thisIoneq/eDensity
                    else:
                       for it in range(nTempDens):
                            rate[it, goodWvl] = f*pop[it, l2]*distr*ab*thisIoneq[it]/eDensity[it]
                self.TwoPhoton = {'wvl':wvl, 'rate':rate}

            else:
                # He seq
                l1 = 1-1
                l2 = 3 - 1
                wvl0 = 1.e+8/(self.Elvlc['ecm'][l2] - self.Elvlc['ecm'][l1])
                goodWvl = wvl > wvl0
                if goodWvl.sum() > 0:
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
                        rate[goodWvl] = f*pop[l2]*distr*ab*thisIoneq/eDensity
                    else:
                       for it in range(nTempDens):
                            rate[it, goodWvl] = f*pop[it, l2]*distr*ab*thisIoneq[it]/eDensity[it]
                self.TwoPhoton = {'wvl':wvl, 'rate':rate}
        #
        #-----------------------------------------------------------------
        #
    def twoPhotonLoss(self):
        ''' to calculate the two-photon energy loss rate - only for hydrogen- and helium-like ions
        includes the elemental abundance and the ionization equilibrium'''
        if self.Z -self.Ion > 1 or self.Dielectronic:
            # this is not a hydrogen-like or helium-like ion
            nTempDens = max(self.Temperature.size, self.EDensity.size)
#            if nTempDens > 1:
            rate = np.zeros((nTempDens), 'float64')
            self.TwoPhotonLoss = {'rate':rate}
        else:
            if hasattr(self, 'Abundance'):
                ab=self.Abundance
            else:
                self.Abundance = util.abundanceRead()
                ab=self.Abundance
            if hasattr(self, 'IoneqOne'):
                thisIoneq=self.IoneqOne
            else:
                self.ioneqOne()
                thisIoneq=self.IoneqOne
            if hasattr(self, 'Population'):
                pop = self.Population['population']
                nTempDens = max(self.Temperature.size, self.EDensity.size)
            else:
                self.populate()
                pop = self.Population['population']
                nTempDens = max(self.Temperature.size, self.EDensity.size)
            if nTempDens > 1:
                rate = np.zeros((nTempDens), 'float64')
                if self.EDensity.size == 1:
                    eDensity = np.repeat(self.EDensity, nTempDens)
                else:
                    eDensity = self.EDensity
            else:
                eDensity = self.EDensity
            if self.Z == self.Ion:
                # H seq
                l1 = 1 - 1
                l2 = 2 - 1
                wvl0 = 1.e+8/(self.Elvlc['ecm'][l2] - self.Elvlc['ecm'][l1])
                dist = util.twophotonHRead()
                avalue = dist['avalue'][self.Z-1]
                f = (avalue*const.light*const.planck*1.e+8)/wvl0
                if nTempDens == 1:
                    rate = f*pop[l2]*ab*thisIoneq/eDensity
                else:
                   for it in range(nTempDens):
                        rate[it] = f*pop[it, l2]*ab*thisIoneq[it]/eDensity[it]
                self.TwoPhotonLoss = {'temperature':self.Temperature,'eDensity':self.EDensity,'rate':rate}
            else:
                # He seq
                l1 = 1 - 1
                l2 = 3 - 1
                wvl0 = 1.e+8/(self.Elvlc['ecm'][l2] - self.Elvlc['ecm'][l1])
                dist = util.twophotonHeRead()
                avalue = dist['avalue'][self.Z-1]
                f = (avalue*const.light*const.planck*1.e+8)/wvl0
                if nTempDens == 1:
                    rate = f*pop[l2]*ab*thisIoneq/eDensity
                else:
                   for it in range(nTempDens):
                        rate[it] = f*pop[it, l2]*ab*thisIoneq[it]/eDensity[it]
                self.TwoPhotonLoss = {'temperature':self.Temperature,'eDensity':self.EDensity,'rate':rate}
        #
        # ----------------------------------------------
        #
        #
        # ------------------------------------------------------------------------------
        #
class ionWeb(ion):
    """
    a class that contains methods to be used for 'Chianti on the Web'
    """
    def gofntSelectLines(self,wvlRange=0, top=10,  saveFile=0):
        """Provide a selection of lines for calculating the 'so-called' G(T) function.

        Given as a function of both temperature and eDensity.

        Only the top( set by 'top') brightest lines are plotted."""
        #
        #self.emiss={"wvl":wvl,"emiss":em,"units":units,"plotLabels":plotLabels}
        #
        #
        if hasattr(self, 'Emiss'):
            doEmiss=False
            em=self.Emiss
        else:
            doEmiss=True
        #
        #
        if doEmiss:
            # new values of temperature or eDensity
            self.emiss()
            em=self.Emiss
        #
        #
        if hasattr(self, 'Abundance'):
            ab=self.Abundance
        else:
            self.Abundance = util.abundanceRead()
            ab=self.Abundance
        #
        fontsize=12
        #
        emiss=em["emiss"]
        wvl=em["wvl"]
        temperature=self.Temperature
        eDensity=self.EDensity
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
            self.message = ' no lines in selected interval'
#            if chInteractive:
#                print ' no lines in selected interval'
#            else:
#                self.message = ' no lines in selected interval'
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
        ndens=self.EDensity.size
        #
        ylabel = 'Emissivity relative to '+maxWvl
        title = self.Spectroscopic
        #
        #  normally, ionWeb is only using in the non-interactive mode
        pl.ioff()
#        if chInteractive:
#            pl.ion()
#        else:
#            pl.ioff()
        #
        #
        if ndens==1 and ntemp==1:
            print ' only a single temperature and eDensity'
            return
        elif ndens == 1:
            xlabel='Temperature (K)'
            ngofnt = temperature.size
            xvalues=temperature
            outTemperature=temperature
            outDensity=np.zeros(ntemp,'Float32')
            outDensity.fill(eDensity)
            desc_str=' at Density = %10.2e' % eDensity
        elif ntemp == 1:
            xvalues=eDensity
            ngofnt = eDensity.size
            outTemperature=np.zeros(ndens,'Float32')
            outTemperature.fill(temperature)
            outDensity=eDensity
            xlabel=r'$\rm{Electron Density (cm}^{-3}\rm{)}$'
            desc_str=' at Temperature = %10.2e' % temperature
        else:
            outTemperature=temperature
            outDensity=eDensity
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
        self.topLines = topLines
#        gline = gui.selectorDialog(wvlChoices,label='Select line(s)')
#        gline_idx=gline.selectedIndex
        #
        #
        # -------------------------------------------------------------------------------------
        #
    def gofntShow(self, wvlRange=0, top=10, index=0, saveFile=0):
        """Return a plot of the 'so-called' G(T) function fron the selected lines in index

        Given as a function of both temperature and eDensity.

        Only the top( set by 'top') brightest lines are plotted."""
        #
        #self.emiss={"wvl":wvl,"emiss":em,"units":units,"plotLabels":plotLabels}
        #
        #
        if hasattr(self, 'Emiss'):
            doEmiss=False
            em=self.Emiss
        else:
            doEmiss=True
        #
        #
        if doEmiss:
            # new values of temperature or eDensity
            self.emiss()
            em=self.Emiss
        #
        #
        if hasattr(self, 'Abundance'):
            ab=self.Abundance
        else:
            self.Abundance = util.abundanceRead()
            ab=self.Abundance
        #
        fontsize=12
        #
        emiss=em["emiss"]
        wvl=em["wvl"]
        temperature=self.Temperature
        eDensity=self.EDensity
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
        ndens=self.EDensity.size
        #
        ylabel = ' Gofnt '
        title = self.Spectroscopic
        #
        #
        if ndens==1 and ntemp==1:
            print ' only a single temperature and eDensity'
            return
        elif ndens == 1:
            xlabel='Temperature (K)'
            ngofnt = temperature.size
            xvalues=temperature
            outTemperature=temperature
            outDensity=np.zeros(ntemp,'Float32')
            outDensity.fill(eDensity)
            desc_str=' at Density = %10.2e' % eDensity
        elif ntemp == 1:
            xvalues=eDensity
            ngofnt = eDensity.size
            outTemperature=np.zeros(ndens,'Float32')
            outTemperature.fill(temperature)
            outDensity=eDensity
            xlabel=r'$\rm{Electron Density (cm}^{-3}\rm{)}$'
            desc_str=' at Temperature = %10.2e' % temperature
        else:
            outTemperature=temperature
            outDensity=eDensity
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
        nWvl = len(index)
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
        gIoneq=self.IoneqOne/eDensity
        #
        #
        #  ionWeb is normally only used in the non-interative mode
        pl.ioff()
#        if chInteractive:
#            pl.ion()
#        else:
#            pl.ioff()
        #
        #
        # plot the desired ratio
        pl.figure()
#        g_line = gline_idx#  [0]
        #print ' g_line = ',g_line
        #
        if nWvl > 1:
            gofnt=np.zeros((ngofnt) ,'float64')
#            for aline in g_line:
            for aline in gline_idx:
                gofnt += gAbund*gIoneq*emiss[aline].squeeze()
        else:
            gofnt = gAbund*gIoneq*emiss[index].squeeze()

        self.Gofnt={'temperature':outTemperature,'eDensity':outDensity,'gofnt':gofnt, 'index':gline_idx}
        #
        pl.loglog(xvalues,gofnt)
        pl.xlim(xvalues.min(),xvalues.max())
        pl.xlabel(xlabel,fontsize=fontsize)
        pl.ylabel('Gofnt',fontsize=fontsize)
        pl.title(title+' '+str(wvl[index])+' '+desc_str, fontsize=fontsize)
        if saveFile:
            pl.savefig(saveFile)
        else:
            pl.show()
        #pl.ioff()
        #pl.show()
#        return
    def intensityRatioSelectLines(self, wvlRange=0, top=10,  saveFile=0):
        """Provide a selection of lines for calculating the 'so-called' G(T) function.

        Given as a function of both temperature and eDensity.

        Only the top( set by 'top') brightest lines are plotted."""
        #
        #self.emiss={"wvl":wvl,"emiss":em,"units":units,"plotLabels":plotLabels}
        #
        #
        if hasattr(self, 'Emiss'):
            doEmiss=False
            em=self.Emiss
        else:
            doEmiss=True
        #
        #
        if doEmiss:
            # new values of temperature or eDensity
            self.emiss()
            em=self.Emiss
        #
        #
#        try:
#            ab=self.Abundance
#        except:
#            self.Abundance = util.abundanceRead()
#            ab=self.Abundance
        #
        fontsize=12
        #
        emiss=em["emiss"]
        wvl=em["wvl"]
        temperature=self.Temperature
        eDensity=self.EDensity
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
        if nlines < 2:
            print ' less than 2 lines in selected interval'
            self.message = ' less than 2 lines in selected interval'
            self.Error = 1
#            if chInteractive:
#                print ' less than 2 lines in selected interval'
#            else:
#                self.message = ' less than 2 lines in selected interval'
#                self.Error = 1
            return
        self.Error = 0
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
        ndens=self.EDensity.size
        #
        ylabel = 'Emissivity relative to '+maxWvl
        title = self.Spectroscopic
        #
        #
        if ndens==1 and ntemp==1:
            print ' only a single temperature and eDensity'
            return
        elif ndens == 1:
            xlabel='Temperature (K)'
            ngofnt = temperature.size
            xvalues=temperature
            outTemperature=temperature
            outDensity=np.zeros(ntemp,'Float32')
            outDensity.fill(eDensity)
            desc_str=' at Density = %10.2e' % eDensity
        elif ntemp == 1:
            xvalues=eDensity
            ngofnt = eDensity.size
            outTemperature=np.zeros(ndens,'Float32')
            outTemperature.fill(temperature)
            outDensity=eDensity
            xlabel=r'$\rm{Electron Density (cm}^{-3}\rm{)}$'
            desc_str=' at Temperature = %10.2e' % temperature
        else:
            outTemperature=temperature
            outDensity=eDensity
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
            pl.loglog(eDensity,emiss[topLines[top-1]]/maxAll, visible=False)
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
    def intensityRatioShow(self,numIdx, denIdx, plotDir=0, saveDir=0):
        """Plot the ratio of 2 lines or sums of lines.

        Shown as a function of eDensity and/or temperature.

        to save a plot or txt, only the directory name is needed"""
        #
        #        self.Emiss={"temperature":temperature,"eDensity":eDensity,"wvl":wvl,"emiss":em,
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
#            # new values of temperature or eDensity
#            self.emiss()
#            em=self.Emiss
        #
        #
#        try:
#            ab=self.Abundance
#        except:
#            self.Abundance = util.abundanceRead()
#            ab=self.Abundance
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
        ndens=self.EDensity.size
        #
        ylabel='Emissivity relative to '+maxWvl
        title=self.Spectroscopic
        #
        #
        if ndens==1 and ntemp==1:
            print ' only a single temperature and eDensity'
            return
        elif ndens == 1:
            xlabel='Temperature (K)'
            xvalues=self.Temperature
            outTemperature=self.Temperature
            outDensity=np.zeros(ntemp,'Float32')
            outDensity.fill(self.EDensity)
            desc_str=' at  Density = %10.2e (cm)$^{-3}$' % self.EDensity
        elif ntemp == 1:
            xvalues=self.EDensity
            outTemperature=np.zeros(ndens,'Float32')
            outTemperature.fill(self.Temperature)
            outDensity=self.EDensity
            xlabel=r'$\rm{Electron Density (cm)^{-3}}$'
            desc_str=' at Temp = %10.2e (K)' % self.Temperature
        else:
            outTemperature=self.Temperature
            outDensity=self.EDensity
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
            print ' no numerator lines were selected'
            self.Message = ' no numerator lines were selected'
#            if chInteractive:
#                print ' no numerator lines were selected'
#            else:
#                self.Message = ' no numerator lines were selected'
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
        #
        intRatio = numEmiss/denEmiss
        fontsize = 12
        #
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
            pl.loglog(outDensity,intRatio, visible=False)
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
        intensityRatioFileName = self.IonStr
        for aline in num_idx:
            addstr = '%10.3f'%(wvl[topLines[aline]])
            intensityRatioFileName += '_' + addstr.strip()
        intensityRatioFileName+='_2'
        for aline in den_idx:
            addstr = '%10.3f'%(wvl[topLines[aline]])
            intensityRatioFileName += '_' + addstr.strip()
        #
        #  need to so the before the next save statements
        self.IntensityRatio = {'ratio':intRatio,'desc':desc,
                'temperature':outTemperature,'eDensity':outDensity,'filename':intensityRatioFileName}
        #
        if plotDir:
            plotFile = os.path.join(plotDir, intensityRatioFileName+'.png')
            pl.savefig(plotFile)
        #
        if saveDir:
            txtFile = os.path.join(saveDir, intensityRatioFileName+'.txt')
            self.intensityRatioSave(outFile = txtFile)
        #
        # ------------------------------------------------------------------------
        #
class ioneq(ion):
    '''
    Reads, calculates, and/or plots the ionization equilibrium for an element as a function of temperature.
    The variable z is the atomic number of the element.  Acceptable values are from 1 to 30.
    '''
    def __init__(self,z,  verbose=False):
        '''
        a class to calculate and plot ionization equilibria
        '''
        self.Z=z
        #
        # ---------------------------------------------------
        #
    def load(self, ioneqName):
        '''
        read in an existing file ionization equilibrium file
        ioneqName should be something like 'chianti', or 'chianti_version6'
        '''
        ioneqAll = util.ioneqRead(ioneqName)
        self.Temperature = ioneqAll['ioneqTemperature']
        self.Ioneq = ioneqAll['ioneqAll'][self.Z - 1]
        #
        # ---------------------------------------------------
        #
    def calculate(self, temperature):
        self.Temperature = np.array(temperature, 'float64')
        ionList=[]
        chIons=[]
        z = self.Z
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
        if ntemp == 1:
            ioneq=np.zeros((z+1), 'float32')
            factor = []
            for anIon in chIons:
                if hasattr(anIon, 'RecombRate') and hasattr(anIon, 'IonizRate'):
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
    def plot(self, stages=0, xr=0, yr=0, oplot=0, label=1, title=1,  bw=0):
        '''
        Plots the ionization equilibria.

        self.plot(xr=None, yr=None, oplot=False)
        stages = sequence of ions to be plotted, neutral == 1, fully stripped == Z+1
        xr = temperature range, yr = ion fraction range

        for overplotting:
        oplot="ioneqfilename" such as 'mazzotta'
        or if oplot=True or oplot=1 and a widget will come up so that a file can be selected.
        '''
        if bw:
            linestyle=['k-','k--', 'k-.', 'k:']
        else:
            linestyle=['b-','r--', 'g-.', 'm:']
        #
        if not stages:
            stages=range(1, self.Z+2)
        elif min(stages) < 1 or max(stages) > self.Z+1:
            stages=range(1, self.Z+2)  #  spectroscopic notation
        if not xr:
            xr=[self.Temperature.min(), self.Temperature.max()]
        if not yr:
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
        atitle='Chianti Ionization Equilibrium for '+const.El[self.Z-1].capitalize()
        #
        if oplot:
            if type(oplot) == types.BooleanType:
                result=util.ioneqRead(ioneqname='')
                print 'keys = ', result.keys()
                if result != False:
                    atitle+='  & '+result['ioneqname'].replace('.ioneq', '')
                    atitle+=' '+linestyle[0]
                    for iz in stages:
                        pl.plot(result['ioneqTemperature'], result['ioneqAll'][self.Z-1, iz-1],linestyle[1])
            elif type(oplot) == types.StringType:
                atitle+='  & ' + oplot
                result = util.ioneqRead(ioneqname=oplot)
                print 'keys = ', result.keys()
#                print result
                if result != False:
                    for iz in stages:
                        pl.plot(result['ioneqTemperature'], result['ioneqAll'][self.Z-1, iz-1],linestyle[1])
            elif type(oplot) == types.ListType:
                for iplot in range(len(oplot)):
                    result = util.ioneqRead(ioneqname=oplot[iplot])
                    print 'keys = ', result.keys()
                    if result != False:
                        atitle+='  & '+oplot[iplot]+' '+linestyle[iplot%3]
                        for iz in stages:
                            pl.plot(result['ioneqTemperature'], result['ioneqAll'][self.Z-1, iz-1],linestyle[1])
            else:
                print ' oplot not understood ', oplot
        if title:
            pl.title(atitle)
        pl.axis(xyr)
    #
    # -------------------------------------------------------------------------
    #
