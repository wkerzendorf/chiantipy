import numpy as np
#class ioneq(ion):
#    '''Calculates the ionization equilibrium for an element as a function of temperature.
#    The variable z is the atomic number of the element.  Acceptable values are from 1 to 30.'''
#    def __init__(self,z, temperature, verbose=False):
##        self.Defaults=defaults
#        ionList=[]
#        chIons=[]
#        self.Z=z
#        self.Temperature = np.array(temperature, 'float64')
#        for stage in range(1, z+2):
#            ionStr=util.zion2name(z, stage)
#            ionList.append(ionStr)
#            print z, stage, ionStr
#            atom=ion(ionStr, temperature = self.Temperature)
#            atom.ionizRate()
#            atom.recombRate()
#            chIons.append(atom)
##        for anIon in chIons:
##            print ' this ion = ', anIon.Ions
##            if type(anIon.IonizRate) != NoneType:
##                pl.loglog(anIon.IonizRate['temperature'], anIon.IonizRate['rate'])
##        #
##        for anIon in chIons:
##            print ' this ion = ',  anIon.Ions
##            if type(anIon.RecombRate) != NoneType:
##                pl.loglog(anIon.RecombRate['temperature'], anIon.RecombRate['rate'])
#        #
#        ntemp=chIons[0].IonizRate['temperature'].size
#        print ' ntemp = ',ntemp
#        if ntemp == 1:
#            ioneq=np.zeros((z+1), 'float32')
#            factor = []
#            for anIon in chIons:
#                if type(anIon.IonizRate) != types.NoneType and type(anIon.RecombRate) != types.NoneType:
#                    rat=anIon.IonizRate['rate']/anIon.RecombRate['rate']
#                    factor.append(rat**2 + rat**(-2))
#                else:
#                    factor.append(0.)
#            factor[0]=max(factor)
#            factor[-1]=max(factor)
#            ionmax=factor.index(min(factor))
##            print ' it, ionmax', it, ionmax
#            ioneq[ionmax]=1.
#            #
#            for iz in range(ionmax+1, z+1):
#                ionrate=chIons[iz-1].IonizRate['rate']
#                recrate=chIons[iz].RecombRate['rate']
#                ioneq[iz]=ionrate*ioneq[iz-1]/recrate
#            #
#            for iz in range(ionmax-1, -1, -1):
#                ionrate=chIons[iz].IonizRate['rate']
#                recrate=chIons[iz+1].RecombRate['rate']
#                ioneq[iz]=recrate*ioneq[iz+1]/ionrate
#            ionsum=ioneq.sum()
##            print ' ionsum = ', ionsum
#            ioneq=ioneq/ionsum
#            self.Ioneq=ioneq
#        #  ntemp >1
#        else:
#            ioneq=np.zeros((z+1,ntemp ), 'float32')
#            for it in range(ntemp):
#                factor=[]
#                for anIon in chIons:
#                    if type(anIon.IonizRate) != types.NoneType and type(anIon.RecombRate) != types.NoneType:
#                        rat=anIon.IonizRate['rate'][it]/anIon.RecombRate['rate'][it]
#                        factor.append(rat**2 + rat**(-2))
#                    else:
#                        factor.append(0.)
#                factor[0]=max(factor)
#                factor[-1]=max(factor)
#                ionmax=factor.index(min(factor))
#    #            print ' it, ionmax', it, ionmax
#                ioneq[ionmax, it]=1.
#                #
#                for iz in range(ionmax+1, z+1):
#                    ionrate=chIons[iz-1].IonizRate['rate'][it]
#                    recrate=chIons[iz].RecombRate['rate'][it]
#                    ioneq[iz, it]=ionrate*ioneq[iz-1, it]/recrate
#                #
#                for iz in range(ionmax-1, -1, -1):
#                    ionrate=chIons[iz].IonizRate['rate'][it]
#                    recrate=chIons[iz+1].RecombRate['rate'][it]
#                    ioneq[iz, it]=recrate*ioneq[iz+1, it]/ionrate
#                ionsum=ioneq[:, it].sum()
#    #            print ' ionsum = ', ionsum
#                ioneq[:, it]=ioneq[:, it]/ionsum
#            self.Ioneq=ioneq
##
#    def plot(self, stages=None, xr=None, yr=None, oplot=False, label=True, title=True,  bw=False):
#        '''Plots the ionization equilibria.
#
#        self.plot(xr=None, yr=None, oplot=False)
#        stages = sequence of ions to be plotted, neutral == 1, fully stripped == Z+1
#        xr = temperature range, yr = ion fraction range
#
#        for overplotting:
#        oplot="ioneqfilename" such as 'mazzotta'
#        or if oplot=True or oplot=1 and a widget will come up so that a file can be selected.'''
#        if bw:
#            linestyle=['k-','k--', 'k-.', 'k:']
#        else:
#            linestyle=['b-','r--', 'g-.', 'm:']
#        #
#        if type(stages) == types.NoneType:
#            stages=range(1, self.Z+2)
#        elif min(stages) < 1 or max(stages) > self.Z+1:
#            stages=range(1, self.Z+2)  #  spectroscopic notation
#        if type(xr) == types.NoneType:
#            xr=[self.Temperature.min(), self.Temperature.max()]
#        if type(yr) == types.NoneType:
#            yr=[0.01, 1.1]
#        xyr=list(xr)
#        xyr.extend(list(yr))
#        #
#        iz=stages[0]
#        pl.loglog(self.Temperature, self.Ioneq[iz-1])
#        if label:
#            idx=self.Ioneq[iz-1] == self.Ioneq[iz-1].max()
#            if idx.sum() > 1:
#                jdx=np.arange(len(idx))
#                idx=jdx[idx].max()
#            ann=const.Ionstage[iz-1]
#            pl.annotate(ann, [self.Temperature[idx], 0.7*self.Ioneq[iz-1, idx]], ha='center')
#        for iz in stages[1:]:
#            pl.plot(self.Temperature, self.Ioneq[iz-1], linestyle[0])
#            if label:
#                idx=self.Ioneq[iz-1] == self.Ioneq[iz-1].max()
#                if idx.sum() > 1:
#                    jdx=np.arange(len(idx))
#                    idx=jdx[idx].mean()
#                ann=const.Ionstage[iz-1]
#                pl.annotate(ann, [self.Temperature[idx], 0.7*self.Ioneq[iz-1, idx]], ha='center')
#        pl.xlabel('Temperature (K)')
#        pl.ylabel('Ion Fraction')
#        atitle='Chianti Ionization Equilibrium for '+El[self.Z-1].capitalize()
#        #
#        if oplot != False:
#            if type(oplot) == BooleanType:
#                result=self.ioneqRead(ioneqname='',default=False)
#                if result != False:
#                    atitle+='  & '+result['ioneqname'].replace('.ioneq', '')
#                    atitle+=' '+linestyle[0]
#                    for iz in ions:
#                        pl.plot(self.IoneqTemperature, self.IoneqAll[self.Z-1, iz-1],linestyle[0], linestyle[1])
#            elif type(oplot) == StringType:
#                atitle+='  & '+oplot+' '+linestyle[0]
#                atitle+=' '+linestyle[0]
#                result=self.ioneqRead(ioneqname=oplot,default=False)
#                if result != False:
#                    for iz in ions:
#                        pl.plot(self.IoneqTemperature, self.IoneqAll[self.Z-1, iz-1],linestyle[0], linestyle[1])
#            elif type(oplot) == ListType:
#                for iplot in range(len(oplot)):
#                    result=self.ioneqRead(ioneqname=oplot[iplot],default=False)
#                    if result != False:
#                        atitle+='  & '+oplot[iplot]+' '+linestyle[iplot%3]
#                        for iz in ions:
#                            pl.plot(self.IoneqTemperature, self.IoneqAll[self.Z-1, iz-1], linestyle[iplot%4])
#            else:
#                print ' oplot not understood ', oplot
#        if title:
#            pl.title(atitle)
#        pl.axis(xyr)
#    #
#    # -------------------------------------------------------------------------
#    #
