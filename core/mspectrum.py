class mspectrum:
    ''' this is the multiprocessing version of spectrum
    Calculate the emission spectrum as a function of temperature and density.

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
    temperature/density.
    proc = the number of processors to use
    timeout - a small but non-zero value seems to be necessary
    '''
    def __init__(self, temperature, density, wavelength, filter=(chfilters.gaussianR, 1000.),  ionList = 0, minAbund=0., doContinuum=1, em = None,  proc=3,  verbose = 0,  timeout=0.1):
        t1 = datetime.now()
        masterlist = util.masterListRead()
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
        abundAll = self.AbundanceAll['abundance']
        nonzed = abundAll > 0.
        minAbundAll = abundAll[nonzed].min()
        if minAbund < minAbundAll:
            minAbund = minAbundAll
        ionInfo = util.masterListInfo()
        wavelength = np.asarray(wavelength)
        nWvl = wavelength.size
        self.Wavelength = wavelength
        wvlRange = [wavelength.min(), wavelength.max()]
        #
        proc = min([proc, mp.cpu_count()])
        #
        freeFree = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        freeBound = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        twoPhoton = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        lineSpectrum = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        #
        #  free-free multiprocessing setup
        ffWorkerQ = mp.Queue()
        ffDoneQ = mp.Queue()
        #
        #  free-bound multiprocessing setup
        #
        fbWorkerQ = mp.Queue()
        fbDoneQ = mp.Queue()
        #
        #  ion multiprocessing setup
        ionWorkerQ = mp.Queue()
        ionDoneQ = mp.Queue()
        #
        #
        self.Todo = []
        for iz in range(31):
            abundance = self.AbundanceAll['abundance'][iz-1]
            if abundance >= minAbund:
                if chInteractive and verbose:
                    print ' %5i %5s abundance = %10.2e '%(iz, const.El[iz-1],  abundance)
                #
                for ionstage in range(1, iz+2):
                    ionS = util.zion2name(iz, ionstage)
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
                        if chInteractive and verbose:
                            print ' setting up continuum calculation for :  ',  ionS
                        ffWorkerQ.put((ionS, temperature, wavelength))
                        fbWorkerQ.put((ionS, temperature, wavelength))
#                        fbInputs.append([ionS, temperature, wavelength])
                        #
                    if masterListTest and wvlTestMin and wvlTestMax and ioneqTest:
                        if chInteractive and verbose:
                            print ' setting up spectrum calculation for  :  ', ionS
                        ionWorkerQ.put((ionS, temperature, density, wavelength, filter))
                        self.Todo.append(ionS)
                    # get dielectronic lines
                    if masterListTestD and wvlTestMinD and wvlTestMaxD and ioneqTestD:
                        if chInteractive and verbose:
                            print ' setting up  spectrum calculation for  :  ', ionSd
#                        dielWorkerQ.put((ionSd, temperature, density, wavelength, filter))
                        ionWorkerQ.put((ionSd, temperature, density, wavelength, filter))
                        self.Todo.append(ionSd)
        #
        ffWorkerQSize = ffWorkerQ.qsize()
        fbWorkerQSize = fbWorkerQ.qsize()
        ionWorkerQSize = ionWorkerQ.qsize()
        if doContinuum:
            ffProcesses = []
            for i in range(proc):
                p = mp.Process(target=mputil.doFfQ, args=(ffWorkerQ, ffDoneQ))
                p.start()
                ffProcesses.append(p)
    #       timeout is not necessary
            for p in ffProcesses:
                if p.is_alive():
                    p.join(timeout=timeout)
#            for i in range(proc):
#                ffProcesses.append('STOP')
            #
            for iff in range(ffWorkerQSize):
                thisFreeFree =ffDoneQ.get()
                if nTempDen ==1:
                    freeFree += thisFreeFree['rate']
                else:
                    for iTempDen in range(nTempDen):
                        freeFree[iTempDen] += thisFreeFree['rate'][iTempDen]
            for p in ffProcesses:
                if not isinstance(p, str):
                    p.terminate()
        #
            fbProcesses = []
            for i in range(proc):
                p = mp.Process(target=mputil.doFbQ, args=(fbWorkerQ, fbDoneQ))
                p.start()
                fbProcesses.append(p)
    #       timeout is not necessary
            for p in fbProcesses:
                if p.is_alive():
                    p.join(timeout=timeout)
#            for i in range(proc):
#                fbProcesses.append('STOP')
            #
            for ifb in range(fbWorkerQSize):
                thisFreeBound = fbDoneQ.get()
                if thisFreeBound.has_key('rate'):
                    if nTempDen ==1:
                        freeBound += thisFreeBound['rate']
                    else:
                        for iTempDen in range(nTempDen):
                            freeBound[iTempDen] += thisFreeBound['rate'][iTempDen]
            for p in fbProcesses:
                if not isinstance(p, str):
                    p.terminate()
        #
        ionProcesses = []
        if ionWorkerQSize < proc:
            nproc = ionWorkerQSize
        for i in range(proc):
            p = mp.Process(target=mputil.doIonQ, args=(ionWorkerQ, ionDoneQ))
            p.start()
            ionProcesses.append(p)
#            ionWorkerQ.put('STOP')
#       timeout is not necessary
        for p in ionProcesses:
#            print' process is alive:  ', p.is_alive()
            if p.is_alive():
#                p.join()
                p.join(timeout=timeout)
#        for i in range(proc):
#            ionProcesses.append('STOP')
        self.Finished = []
        #
        for ijk in range(ionWorkerQSize):
            out = ionDoneQ.get()
            ions = out[0]
            self.Finished.append(ions)
            aspectrum = out[1]
            try:
                if nTempDen == 1:
                    lineSpectrum += aspectrum['intensity']
                else:
                    for iTempDen in range(nTempDen):
                        lineSpectrum[iTempDen] += aspectrum['intensity'][iTempDen]
               # check for two-photon emission
                if len(out) == 3:
                    tp = out[2]
                    if nTempDen == 1:
                        twoPhoton += tp['rate']
                    else:
                        for iTempDen in range(nTempDen):
                            twoPhoton[iTempDen] += tp['rate'][iTempDen]
            except:
                print '  error in ion pool'
        #
        for p in ionProcesses:
            if not isinstance(p, str):
                p.terminate()
        #
        #
        #
        self.FreeFree = {'wavelength':wavelength, 'intensity':freeFree.squeeze()}
        self.FreeBound = {'wavelength':wavelength, 'intensity':freeBound.squeeze()}
        self.LineSpectrum = {'wavelength':wavelength, 'intensity':lineSpectrum.squeeze()}
        self.TwoPhoton = {'wavelength':wavelength, 'intensity':twoPhoton.squeeze()}
        #
        total = freeFree + freeBound + lineSpectrum + twoPhoton
        t2 = datetime.now()
        dt=t2-t1
        if chInteractive and verbose:
            print ' elapsed seconds = ', dt.seconds
        if type(em) != types.NoneType:
            if nEm == 1:
                integrated = total*em
            else:
                integrated = np.zeros_like(wavelength)
                for iTempDen in range(nTempDen):
                    integrated += total[iTempDen]*em[iTempDen]
            self.Spectrum ={'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'em':em}
        else:
            self.Spectrum ={'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'worker':ionWorkerQ, 'done':ionDoneQ}
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
    #
    # -------------------------------------------------------------------------
    #
