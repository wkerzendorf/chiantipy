#import chianti.core.ion as ion
#class ionWeb(ion):
#    """
#    a class that contains methods to be used for 'Chianti on the Web'
#    """
#    def gofntSelectLines(self,wvlRange=0, top=10,  saveFile=0):
#        """Provide a selection of lines for calculating the 'so-called' G(T) function.
#
#        Given as a function of both temperature and density.
#
#        Only the top( set by 'top') brightest lines are plotted."""
#        #
#        #self.emiss={"wvl":wvl,"emiss":em,"units":units,"plotLabels":plotLabels}
#        #
#        #
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
#        #
#        #
#        try:
#            ab=self.Abundance
#        except:
#            self.Abundance = util.abundanceRead()
#            ab=self.Abundance
#        #
#        fontsize=12
#        #
#        emiss=em["emiss"]
#        wvl=em["wvl"]
#        temperature=self.Temperature
#        density=self.Density
#        plotLabels=em["plotLabels"]
#        xLabel=plotLabels["xLabel"]
#        yLabel=plotLabels["yLabel"]
#        #
#        # find which lines are in the wavelength range if it is set
#        #
#        #
#        if type(wvlRange) != type(1):
#            igvl=util.between(wvl,wvlRange)
#        else:
#            igvl=range(len(wvl))
#        nlines=len(igvl)
#        if nlines ==0:
#            if chInteractive:
#                print ' no lines in selected interval'
#            else:
#                self.message = ' no lines in selected interval'
#            return
#        # find the top most intense lines
#        #
#        if (top > nlines) or (top == 0):
#            top=nlines
#        maxEmiss=np.zeros(nlines,'Float32')
#        for iline in range(nlines):
#            maxEmiss[iline]=emiss[igvl[iline]].max()
#        for iline in range(nlines):
#            if maxEmiss[iline]>=maxEmiss.max():
#                maxAll=emiss[igvl[iline]]
#                maxIndex = igvl[iline]
##        print ' maxIndex, maxAll = ', maxIndex,  maxAll
#        line=range(nlines)
#        igvlsort=np.take(igvl,np.argsort(maxEmiss))
#        topLines=igvlsort[-top:]
#        maxWvl='%5.3f' % wvl[topLines[-1]]
#        maxline=topLines[-1]
#        #
#        # need to make sure there are no negative values before plotting
#        good = np.where(emiss > 0.)
#        emissMin=emiss[good].min()
#        bad=np.where(emiss <= 0.)
#        emiss[bad]=emissMin
#        #
#        topLines=topLines[wvl[topLines].argsort()]
#        #
#        #
#        ntemp=self.Temperature.size
#        #
#        ndens=self.Density.size
#        #
#        ylabel = 'Emissivity relative to '+maxWvl
#        title = self.Spectroscopic
#        #
#        #  normally, ionWeb is only using in the non-interactive mode
#        if chInteractive:
#            pl.ion()
#        else:
#            pl.ioff()
#        #
#        #
#        if ndens==1 and ntemp==1:
#            print ' only a single temperature and density'
#            return
#        elif ndens == 1:
#            xlabel='Temperature (K)'
#            ngofnt = temperature.size
#            xvalues=temperature
#            outTemperature=temperature
#            outDensity=np.zeros(ntemp,'Float32')
#            outDensity.fill(density)
#            desc_str=' at Density = %10.2e' % density
#        elif ntemp == 1:
#            xvalues=density
#            ngofnt = density.size
#            outTemperature=np.zeros(ndens,'Float32')
#            outTemperature.fill(temperature)
#            outDensity=density
#            xlabel=r'$\rm{Electron Density (cm}^{-3}\rm{)}$'
#            desc_str=' at Temperature = %10.2e' % temperature
#        else:
#            outTemperature=temperature
#            outDensity=density
#            xlabel='Temperature (K)'
#            xvalues=temperature
#            ngofnt = temperature.size
#            desc_str=' for variable Density'
#            #
#        #
#        # put all actual plotting here
#        #
##        pl.ion()
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
#        if saveFile:
#            pl.savefig(saveFile)
#        else:
#            pl.draw()
#        #
##        print ' topLInes = ', wvl[topLines]
#        wvlChoices = []
#        for one in wvl[topLines]:
#            wvlChoices.append('%12.3f'%(one))
#        self.wvlChoices = wvlChoices
##        gline = gui.selectorDialog(wvlChoices,label='Select line(s)')
##        gline_idx=gline.selectedIndex
#        #
#        #
#        # -------------------------------------------------------------------------------------
#        #
#    def gofntShow(self,wvlRange=0,top=10,index=0, saveFile=0):
#        """Return a plot of the 'so-called' G(T) function fron the selected lines in index
#
#        Given as a function of both temperature and density.
#
#        Only the top( set by 'top') brightest lines are plotted."""
#        #
#        #self.emiss={"wvl":wvl,"emiss":em,"units":units,"plotLabels":plotLabels}
#        #
#        #
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
#        #
#        #
#        try:
#            ab=self.Abundance
#        except:
#            self.Abundance = util.abundanceRead()
#            ab=self.Abundance
#        #
#        fontsize=12
#        #
#        emiss=em["emiss"]
#        wvl=em["wvl"]
#        temperature=self.Temperature
#        density=self.Density
#        plotLabels=em["plotLabels"]
#        xLabel=plotLabels["xLabel"]
#        yLabel=plotLabels["yLabel"]
#        #
#        # find which lines are in the wavelength range if it is set
#        #
#        #
#        if type(wvlRange) != type(1):
#            igvl=util.between(wvl,wvlRange)
#        else:
#            igvl=range(len(wvl))
#        nlines=len(igvl)
#        if nlines ==0:
#            print ' no lines in selected interval'
#            self.Message = ' no lines in selected interval '
#            return
#        # find the top most intense lines
#        #
#        if top > nlines:
#            top=nlines
#        maxEmiss=np.zeros(nlines,'Float32')
#        for iline in range(nlines):
#            maxEmiss[iline]=emiss[igvl[iline]].max()
#        for iline in range(nlines):
#            if maxEmiss[iline]>=maxEmiss.max():
#                maxAll=emiss[igvl[iline]]
#                maxIndex = igvl[iline]
##        print ' maxIndex, maxAll = ', maxIndex,  maxAll
#        line=range(nlines)
#        igvlsort=np.take(igvl,np.argsort(maxEmiss))
#        topLines=igvlsort[-top:]
#        maxWvl='%5.3f' % wvl[topLines[-1]]
#        maxline=topLines[-1]
#        #
#        # need to make sure there are no negative values before plotting
#        good = np.where(emiss > 0.)
#        emissMin=emiss[good].min()
#        bad=np.where(emiss <= 0.)
#        emiss[bad]=emissMin
#        #
#        topLines=topLines[wvl[topLines].argsort()]
#        #
#        #
#        ntemp=self.Temperature.size
#        #
#        ndens=self.Density.size
#        #
#        ylabel = ' Gofnt '
#        title = self.Spectroscopic
#        #
#        #
#        if ndens==1 and ntemp==1:
#            print ' only a single temperature and density'
#            return
#        elif ndens == 1:
#            xlabel='Temperature (K)'
#            ngofnt = temperature.size
#            xvalues=temperature
#            outTemperature=temperature
#            outDensity=np.zeros(ntemp,'Float32')
#            outDensity.fill(density)
#            desc_str=' at Density = %10.2e' % density
#        elif ntemp == 1:
#            xvalues=density
#            ngofnt = density.size
#            outTemperature=np.zeros(ndens,'Float32')
#            outTemperature.fill(temperature)
#            outDensity=density
#            xlabel=r'$\rm{Electron Density (cm}^{-3}\rm{)}$'
#            desc_str=' at Temperature = %10.2e' % temperature
#        else:
#            outTemperature=temperature
#            outDensity=density
#            xlabel='Temperature (K)'
#            xvalues=temperature
#            ngofnt = temperature.size
#            desc_str=' for variable Density'
##            #
##        #
##        # put all actual plotting here
##        #
##        pl.ion()
##        pl.figure()
##        nxvalues=len(xvalues)
##        for iline in range(top):
##            tline=topLines[iline]
##            pl.loglog(xvalues,emiss[tline]/maxAll)
##            skip=2
##            start=divmod(iline,nxvalues)[1]
##            for ixvalue in range(start,nxvalues,nxvalues/skip):
##                pl.text(xvalues[ixvalue],emiss[tline,ixvalue]/maxAll[ixvalue],str(wvl[tline]))
##        pl.xlim(xvalues.min(),xvalues.max())
###       yl=pl.ylim()
###       pl.ylim(yl[0],1.2)
##        pl.xlabel(xlabel,fontsize=fontsize)
##        pl.ylabel(ylabel,fontsize=fontsize)
##        pl.title(title+desc_str,fontsize=fontsize)
##        pl.draw()
##        #
###        print ' topLInes = ', wvl[topLines]
##        wvlChoices = []
##        for one in wvl[topLines]:
##            wvlChoices.append('%12.3f'%(one))
##        gline = gui.selectorDialog(wvlChoices,label='Select line(s)')
##        gline_idx=gline.selectedIndex
##        #
#        gline_idx = index
#        # for now
#        ngofnt = 1
#        #
#        gAbund=self.Abundance
#        #
#        try:
#            thisIoneq=self.IoneqOne
#        except:
#            self.ioneqOne()
#        #        gioneq=np.where(thisIoneq > 0.)
#        #        y2=interpolate.splrep(np.log(self.IoneqAll['ioneqTemperature'][gioneq]),np.log(thisIoneq[gioneq]),s=0)  #allow smoothing,s=0)
#        #        gIoneq=interpolate.splev(np.log(temperature),y2)   #,der=0)
#        gIoneq=self.IoneqOne/density
#        #
#        #
#        #  ionWeb is normally only used in the non-interative mode
#        if chInteractive:
#            pl.ion()
#        else:
#            pl.ioff()
#        #
#        #
#        # plot the desired ratio
#        pl.figure()
#        g_line= topLines[gline_idx]#  [0]
#        #print ' g_line = ',g_line
#        #
#        gofnt=np.zeros(ngofnt,'float32')
#        if ngofnt > 1:
#            for aline in g_line:
#    #        for aline in gline_idx:
#                gofnt += gAbund*gIoneq*emiss[aline].squeeze()
#        else:
#            gofnt = gAbund*gIoneq*emiss[index].squeeze()
#
#        self.Gofnt={'temperature':outTemperature,'density':outDensity,'gofnt':gofnt}
#        #
#        pl.loglog(xvalues,gofnt)
#        pl.xlim(xvalues.min(),xvalues.max())
#        pl.xlabel(xlabel,fontsize=fontsize)
#        pl.ylabel('Gofnt',fontsize=fontsize)
#        pl.title(title+' '+str(wvl[g_line])+' '+desc_str, fontsize=fontsize)
#        if saveFile:
#            pl.savefig(saveFile)
#        else:
#            pl.show()
#        #pl.ioff()
#        #pl.show()
##        return
#    def intensityRatioSelectLines(self,wvlRange=0,top=10,  saveFile=0):
#        """Provide a selection of lines for calculating the 'so-called' G(T) function.
#
#        Given as a function of both temperature and density.
#
#        Only the top( set by 'top') brightest lines are plotted."""
#        #
#        #self.emiss={"wvl":wvl,"emiss":em,"units":units,"plotLabels":plotLabels}
#        #
#        #
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
#        #
#        #
#        try:
#            ab=self.Abundance
#        except:
#            self.Abundance = util.abundanceRead()
#            ab=self.Abundance
#        #
#        fontsize=12
#        #
#        emiss=em["emiss"]
#        wvl=em["wvl"]
#        temperature=self.Temperature
#        density=self.Density
#        plotLabels=em["plotLabels"]
#        xLabel=plotLabels["xLabel"]
#        yLabel=plotLabels["yLabel"]
#        #
#        # find which lines are in the wavelength range if it is set
#        #
#        #
#        if not isinstance(wvlRange, int):
#            igvl=util.between(wvl,wvlRange)
#        else:
#            igvl=range(len(wvl))
#        nlines=len(igvl)
#        if nlines ==0:
#            if chInteractive:
#                print ' no lines in selected interval'
#            else:
#                self.message = ' no lines in selected interval'
#            return
#        # find the top most intense lines
#        #
#        if (top > nlines) or (top == 0):
#            top=nlines
#        maxEmiss=np.zeros(nlines,'Float32')
#        for iline in range(nlines):
#            maxEmiss[iline]=emiss[igvl[iline]].max()
#        for iline in range(nlines):
#            if maxEmiss[iline]>=maxEmiss.max():
#                maxAll=emiss[igvl[iline]]
#                maxIndex = igvl[iline]
##        print ' maxIndex, maxAll = ', maxIndex,  maxAll
#        line=range(nlines)
#        igvlsort=np.take(igvl,np.argsort(maxEmiss))
#        topLines=igvlsort[-top:]
#        maxWvl='%5.3f' % wvl[topLines[-1]]
#        maxline=topLines[-1]
#        #
#        # need to make sure there are no negative values before plotting
#        good = np.where(emiss > 0.)
#        emissMin=emiss[good].min()
#        bad=np.where(emiss <= 0.)
#        emiss[bad]=emissMin
#        #
#        topLines=topLines[wvl[topLines].argsort()]
#        #
#        #
#        ntemp=self.Temperature.size
#        #
#        ndens=self.Density.size
#        #
#        ylabel = 'Emissivity relative to '+maxWvl
#        title = self.Spectroscopic
#        #
#        #
#        if ndens==1 and ntemp==1:
#            print ' only a single temperature and density'
#            return
#        elif ndens == 1:
#            xlabel='Temperature (K)'
#            ngofnt = temperature.size
#            xvalues=temperature
#            outTemperature=temperature
#            outDensity=np.zeros(ntemp,'Float32')
#            outDensity.fill(density)
#            desc_str=' at Density = %10.2e' % density
#        elif ntemp == 1:
#            xvalues=density
#            ngofnt = density.size
#            outTemperature=np.zeros(ndens,'Float32')
#            outTemperature.fill(temperature)
#            outDensity=density
#            xlabel=r'$\rm{Electron Density (cm}^{-3}\rm{)}$'
#            desc_str=' at Temperature = %10.2e' % temperature
#        else:
#            outTemperature=temperature
#            outDensity=density
#            xlabel='Temperature (K)'
#            xvalues=temperature
#            ngofnt = temperature.size
#            desc_str=' for variable Density'
#            #
#        #
#        # put all actual plotting here
#        #
##        pl.ion()
#        #  topLines are sorted by wavelength
#        ymax = np.max(1.2*emiss[topLines[0]]/maxAll)
#        ymin = ymax
#        pl.figure()
#        ax = pl.subplot(111)
#        nxvalues=len(xvalues)
#        for iline in range(top):
#            tline=topLines[iline]
#            pl.loglog(xvalues,emiss[tline]/maxAll)
#            if np.min(emiss[tline]/maxAll) < ymin:
#                ymin = np.min(emiss[tline]/maxAll)
#            if np.max(emiss[tline]/maxAll) > ymax:
#                ymax = np.max(emiss[tline]/maxAll)
#            skip=2
#            start=divmod(iline,nxvalues)[1]
#            for ixvalue in range(start,nxvalues,nxvalues/skip):
#                pl.text(xvalues[ixvalue],emiss[tline,ixvalue]/maxAll[ixvalue],str(wvl[tline]))
#        pl.xlim(xvalues.min(),xvalues.max())
##        print ' ymin, ymax = ', ymin, ymax
##        pl.ylim(ymin, ymax)
##       yl=pl.ylim()
##       pl.ylim(yl[0],1.2)
#        pl.xlabel(xlabel,fontsize=fontsize)
#        pl.ylabel(ylabel,fontsize=fontsize)
#        if ndens == ntemp and ntemp > 1:
#            pl.text(0.07, 0.5,title, horizontalalignment='left', verticalalignment='center', fontsize=fontsize,  transform = ax.transAxes)
#            #
#            ax2 = pl.twiny()
#            xlabelDen=r'Electron Density (cm$^{-3}$)'
#            pl.xlabel(xlabelDen, fontsize=fontsize)
#            pl.loglog(density,emiss[topLines[top-1]]/maxAll, visible=False)
#            ax2.xaxis.tick_top()
#            pl.ylim(ymin/1.2, 1.2*ymax)
#        else:
#            pl.ylim(ymin/1.2, 1.2*ymax)
#            pl.title(title+desc_str,fontsize=fontsize)
#        if saveFile:
#            pl.savefig(saveFile)
#        else:
#            pl.draw()
#        #
##        print ' topLInes = ', wvl[topLines]
#        wvlChoices = []
#        for one in wvl[topLines]:
#            wvlChoices.append('%12.3f'%(one))
#        self.wvlChoices = wvlChoices
#        self.topLines = topLines
#        #
#        #   -----------------------------------
#        #
#    def intensityRatioShow(self,numIdx, denIdx, saveFile=0):
#        """Plot the ratio of 2 lines or sums of lines.
#
#        Shown as a function of density and/or temperature."""
#        #
#        #        self.Emiss={"temperature":temperature,"density":density,"wvl":wvl,"emiss":em,
#        #        "plotLabels":plotLabels}
#        #
#        #
#        em = self.Emiss
#        #
##        doEmiss=False
##        try:
##            em=self.Emiss
##        except:
##            doEmiss=True
##        #
##        #
##        if doEmiss:
##            # new values of temperature or density
##            self.emiss()
##            em=self.Emiss
#        #
#        #
#        try:
#            ab=self.Abundance
#        except:
#            self.Abundance = util.abundanceRead()
#            ab=self.Abundance
#        emiss = em['emiss']
#        wvl = em["wvl"]
#        plotLabels=em["plotLabels"]
#        xLabel=plotLabels["xLabel"]
#        yLabel=plotLabels["yLabel"]
#        #
#        # find which lines are in the wavelength range if it is set
#        #
#        #
##        if not wvlRange:
##            igvl=range(len(wvl))
##        else:
##            igvl=util.between(wvl,wvlRange)
##        nlines=len(igvl)
#        #
##        print ' nlines = ',nlines
##        print ' iglv = ',igvl
##        igvl=np.take(igvl,wvl[igvl].argsort())
#        # find the top most intense lines
#        #
##        if (top > nlines) or (top == 0):
##            top=nlines
##        maxEmiss=np.zeros(nlines,'Float32')
##        for iline in range(nlines):
##            maxEmiss[iline]=emiss[igvl[iline]].max()
##        for iline in range(nlines):
##            if maxEmiss[iline]==maxEmiss.max():
##                maxAll=emiss[igvl[iline]]
##        line=range(nlines)
##        igvlsort=np.take(igvl,np.argsort(maxEmiss))
###        print 'igvlsort = ', igvlsort
##        topLines=igvlsort[-top:]
##        print ' topLines = ', topLines
##        topLines=topLines[wvl[topLines].argsort()]
#        topLines = self.topLines
#        maxWvl='%5.3f' % wvl[topLines[-1]]
#        maxline=topLines[-1]
#        #
#        #
#        #
#        # need to make sure there are no negative values before plotting
#        good = np.where(emiss > 0.)
#        emissMin=emiss[good].min()
#        bad=np.where(emiss <= 0.)
#        emiss[bad]=emissMin
#        #
#        #
#        ntemp=self.Temperature.size
#        #
#        ndens=self.Density.size
#        #
#        ylabel='Emissivity relative to '+maxWvl
#        title=self.Spectroscopic
#        #
#        #
#        if ndens==1 and ntemp==1:
#            print ' only a single temperature and density'
#            return
#        elif ndens == 1:
#            xlabel='Temperature (K)'
#            xvalues=self.Temperature
#            outTemperature=self.Temperature
#            outDensity=np.zeros(ntemp,'Float32')
#            outDensity.fill(self.Density)
#            desc_str=' at  Density = %10.2e (cm)$^{-3}$' % self.Density
#        elif ntemp == 1:
#            xvalues=self.Density
#            outTemperature=np.zeros(ndens,'Float32')
#            outTemperature.fill(self.Temperature)
#            outDensity=self.Density
#            xlabel=r'$\rm{Electron Density (cm)^{-3}}$'
#            desc_str=' at Temp = %10.2e (K)' % self.Temperature
#        else:
#            outTemperature=self.Temperature
#            outDensity=self.Density
#            xlabel='Temperature (K)'
#            xvalues=self.Temperature
#            desc_str=' for variable Density'
#            #
#        #
#        # put all actual plotting here
#        #
#        #
#        # num_idx and den_idx are tuples
#        #
#        if np.iterable(numIdx):
#            num_idx=numIdx
#        else:
#            num_idx = [numIdx]
#        if len(num_idx) == 0:
#            if chInteractive:
#                print ' no numerator lines were selected'
#            else:
#                self.Message = ' no numerator lines were selected'
#            return
#        #
#        if np.iterable(denIdx):
#            den_idx=denIdx
#        else:
#            den_idx = [denIdx]
#        #
#        if len(den_idx) == 0:
#            if chInteractive:
#                print ' no denominator lines were selected'
#            else:
#                self.Message = ' no denominator lines were selected'
#            return
#        #
##       print ' num_idx = ', num_idx
##       print ' toplines[num_idx] = ', topLines[num_idx]
##       num_line= topLines[num_idx]
##       den_line= topLines[den_idx]
##       #
##       print ' num_line = ', num_line
#        numEmiss=np.zeros(len(xvalues),'float32')
#        for aline in num_idx:
#            numEmiss+=emiss[topLines[aline]]
#        #
#        denEmiss=np.zeros(len(xvalues),'float32')
#        for aline in den_idx:
#            denEmiss+=emiss[topLines[aline]]
#        fontsize = 12
#        #
#        # plot the desired ratio
#        pl.figure()
#        ax = pl.subplot(111)
#        pl.loglog(xvalues,numEmiss/denEmiss)
#        pl.xlim(xvalues.min(),xvalues.max())
#        pl.xlabel(xlabel,fontsize=fontsize)
#        pl.ylabel('Ratio ('+self.Defaults['flux']+')',fontsize=fontsize)
#        desc = title + ':'
#        for aline in num_idx:
#            desc += ' ' + str(wvl[topLines[aline]])
#        desc +=' / '
#        for aline in den_idx:
#            desc += ' ' + str(wvl[topLines[aline]])
#        if ndens == ntemp and ntemp > 1:
#            pl.text(0.07, 0.5,desc, horizontalalignment='left', verticalalignment='center', fontsize=fontsize,  transform = ax.transAxes)
#            #
#            ax2 = pl.twiny()
#            xlabelDen=r'Electron Density (cm$^{-3}$)'
#            pl.xlabel(xlabelDen, fontsize=fontsize)
#            pl.loglog(outDensity,numEmiss/denEmiss, visible=False)
#            ax2.xaxis.tick_top()
#        else:
##            pl.ylim(ymin, ymax)
#            pl.title(desc,fontsize=fontsize)
##       desc=title+' '+str(wvl[num_line])+' / '+str(wvl[den_line])+' '+desc_str
##        pl.title(desc, fontsize=fontsize)
##       pl.title(title+' '+str(wvl[num_line])+' / '+str(wvl[den_line])+' '+desc_str,fontsize=fontsize)
##        pl.draw()
##        pl.ioff()
##        pl.show()
#        #
#        if saveFile:
#            pl.savefig(saveFile)
#        intensityRatioFileName=self.IonStr
#        for aline in num_idx:
#            intensityRatioFileName+= '_%3i'%(wvl[topLines[aline]])
#        intensityRatioFileName+='_2'
#        for aline in den_idx:
#            intensityRatioFileName+= '_%3i'%(wvl[topLines[aline]])
#        intensityRatioFileName+='.rat'
#        self.IntensityRatio={'ratio':numEmiss/denEmiss,'desc':desc,
#                'temperature':outTemperature,'density':outDensity,'filename':intensityRatioFileName}
#        #
#        # -------------------------------------------------------------------------------------
#        #
