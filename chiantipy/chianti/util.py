'''Utility functions, many for reading the CHIANTI database files.

Copyright 2009, 2010 Kenneth P. Dere

This software is distributed under the terms of the GNU General Public License
that is found in the LICENSE file


'''
import os
from types import *
from ConfigParser import *
import cPickle
from datetime import date
import numpy as np
from FortranFormat import *
import chianti.constants as const
#import chianti.util as util
#
def between(array,limits):
    '''returns an index array of elements of array where the values lie
    between the limits given as a 2 element list or tuple'''
    array=np.asarray(array)
    nlines=len(array)
    hi=np.where(array >= limits[0],range(1,nlines+1),0)
    lo=np.where(array <= limits[1],range(1,nlines+1),0)
    hilo=hi&lo
    out=[a -1  for a in hilo if a > 0]
    return out
    # -------------------------------------------------------------------------------------
    #
def ipRead(verbose=False):
    """ reads the ionization potential file, returns ip array in eV"""
    topdir=os.environ["XUVTOP"]
    ipname=os.path.join(topdir, 'ip','chianti.ip')
    ipfile=open(ipname)
    data=ipfile.readlines()
    ipfile.close()
    nip=0
    ndata=2
    maxz=0
    while ndata > 1:
        s1=data[nip]
        s2=s1.split()
        ndata=len(s2)
        nip=nip+1
        if int(s2[0]) > maxz:
            maxz=int(s2[0])
    if verbose:
        print ' maxz = ', maxz
    nip=nip-1
    ip=np.zeros((maxz, maxz), 'float32')
    for aline in data[0:nip]:
        s2=aline.split()
        iz=int(s2[0])
        ion=int(s2[1])
        ip[iz-1, ion-1]=float(s2[2])
    return ip*1.239841875e-4
    #
    # -------------------------------------------------------------------------------------
    #
def masterListInfo():
    """ returns information about ions in masterlist
    the reason for this file is to speed up multi-ion spectral calculations
    the information is stored in a pickled file 'masterlist_ions.pkl'
    if the file is not found, one will be created and the following information
    returned for each ion
    wmin, wmax :  the minimum and maximum wavelengths in the wgfa file
    tmin, tmax :  the minimum and maximum temperatures for which the ionization balance is nonzero"""
    dir=os.environ["XUVTOP"]
    infoPath = os.path.join(dir, 'masterlist')
    infoName=os.path.join(dir,'masterlist','masterlist_ions.pkl')
    masterName=os.path.join(dir,'masterlist','masterlist.ions')
    if os.path.isfile(infoName):
#       print ' file exists - ',  infoName
        pfile = open(infoName, 'r')
        masterListInfo = cPickle.load(pfile)
        pfile.close
    elif os.access(infoPath, os.W_OK):
        # the file does not exist but we have write access and will create it
        defaults = defaultsRead()
        print ' defaults = ', defaults
        ioneqName = defaults['ioneqfile']
        ioneq = ioneqRead(ioneqname = ioneqName)
        masterList = masterListRead()
        masterListInfo = {}
        haveZ = [0]*31
        haveStage = np.zeros((31, 31), 'Int32')
        haveDielectronic = np.zeros((31, 31), 'Int32')
        for one in masterList:
            ionInfo = convertName(one)
            z = ionInfo['Z']
            stage = ionInfo['Ion']
            haveZ[z] = 1
            dielectronic = ionInfo['Dielectronic']
            if dielectronic:
                haveDielectronic[z, stage] = 1
            else:
                haveStage[z, stage] = 1
            thisIoneq = ioneq['ioneqAll'][z- 1, stage - 1 + dielectronic]
            good = thisIoneq > 0.
            goodTemp = ioneq['ioneqTemperature'][good]
            tmin = goodTemp.min()
            tmax = goodTemp.max()
            vgood = thisIoneq == thisIoneq.max()
            vgoodTemp = ioneq['ioneqTemperature'][vgood][0]
            wgfa = wgfaRead(one)
            nZeros = wgfa['wvl'].count(0.)
            # two-photon transitions are denoted by a wavelength of zero (0.)
            while nZeros > 0:
                wgfa['wvl'].remove(0.)
                nZeros = wgfa['wvl'].count(0.)
            # unobserved lines are denoted with a negative wavelength
            wvl = np.abs(np.asarray(wgfa['wvl'], 'float64'))
            wmin = wvl.min()
            wmax = wvl.max()
            masterListInfo[one] = {'wmin':wmin, 'wmax':wmax, 'tmin':tmin, 'tmax':tmax, 'tIoneqMax':vgoodTemp}
        masterListInfo['haveZ'] = haveZ
        masterListInfo['haveStage'] = haveStage
        masterListInfo['haveDielectronic'] = haveDielectronic
        #  now do the bare ions from H thru Zn
        #  these are only involved in the continuum
        for iz in range(1, 31):
            ions = zion2name(iz, iz+1)
            thisIoneq = ioneq['ioneqAll'][iz-1, iz]
            good = thisIoneq > 0.
            goodTemp = ioneq['ioneqTemperature'][good]
            tmin = goodTemp.min()
            tmax = goodTemp.max()
            wmin=0.
            wmax = 1.e+30
            masterListInfo[ions] = {'wmin':wmin, 'wmax':wmax, 'tmin':tmin, 'tmax':tmax}
        pfile = open(infoName, 'w')
        cPickle.dump(masterListInfo, pfile)
        pfile.close
    else:
        # the file does not exist and we do NOT have write access to creat it
        # will just make an inefficient, useless version
        masterListInfo = {}
        for one in masterList:
            ionInfo = convertName(one)
            z = ionInfo['Z']
            stage = ionInfo['Ion']
            dielectronic = ionInfo['Dielectronic']
            wmin=0.
            wmax = 1.e+30
            masterListInfo[one] = {'wmin':wmin, 'wmax':wmax, 'tmin':1.e+4, 'tmax':1.e+9}
        #  now do the bare ions from H thru Zn
        #  these are only involved in the continuum
        for iz in range(1, 31):
            ions = zion2name(iz, iz+1)
            wmin=0.
            wmax = 1.e+30
            masterListInfo[ions] = {'wmin':wmin, 'wmax':wmax, 'tmin':1.e+4, 'tmax':1.e+9}
        pfile = open(infoName, 'w')
        cPickle.dump(masterListInfo, pfile)
        pfile.close
        masterListInfo = {'noInfo':'none'}
    return masterListInfo
    #
    # -------------------------------------------------------------------------------------
    #
def masterListRead():
    """ read a chianti masterlist file and return a list of files"""
    dir=os.environ["XUVTOP"]
    fname=os.path.join(dir,'masterlist','masterlist.ions')
    input=open(fname,'r')
    s1=input.readlines()
    dum=input.close()
    masterlist=[]
    for i in range(0,len(s1)):
        s1a=s1[i][:-1]
        s2=s1a.split(';')
        masterlist.append(s2[0].strip())
    return masterlist
    #
    # -------------------------------------------------------------------------------------
    #
def photoxRead(ions):
    """read chianti photoionization .photox files and return
        {"energy", "cross"} where energy is in Rydbergs and the
        cross section is in cm^2  """
    #
    zion=convertName(ions)
    if zion['Z'] < zion['Ion']:
        print ' this is a bare nucleus that has no ionization rate'
        return
    #
    fname=ion2filename(ions)
    paramname=fname+'.photox'
    input=open(paramname,'r')
    lines = input.readlines()
    input.close
    # get number of energies
#    neng = int(lines[0][0:6])
    dataEnd = 0
    lvl1 = []
    lvl2 = []
    energy = []
    cross = []
    icounter = 0
    while not dataEnd:
        lvl11 = int(lines[icounter][:8])
        lvl21 = int(lines[icounter][8:15])
        ener = lines[icounter][15:].split()
        energy1 = np.asarray(ener, 'float64')
        #
        icounter += 1
        irsl = int(lines[icounter][:8])
        ind0 = int(lines[icounter][8:15])
        if irsl != lvl11 or ind0 != lvl21:
            # this only happens if the file was written incorrectly
            print ' lvl1, lvl2 = ', lvl11, lvl21
            print ' irsl, indo = ', irsl,  ind0
            return
        crs = lines[icounter][15:].split()
        cross1 = np.asarray(crs, 'float64')
        lvl1.append(lvl11)
        lvl2.append(lvl21)
        energy.append(energy1)
        cross.append(cross1)
        icounter += 1
        dataEnd = lines[icounter].count('-1')
    ref = lines[icounter+1:-1]
    cross = np.asarray(cross, 'float64')
    energy = np.asarray(energy, 'float64')
    return {'lvl1':lvl1, 'lvl2':lvl2,'energy':energy, 'cross':cross,  'ref':ref}
    #
    # -------------------------------------------------------------------------------------
    #
def ionrecdatRead(filename):
    """ read chianti ionxdat, ionizdat, recombdat files and return
    {"ev":ev,"cross":cross,"crosserr":crosserr,"ref":ref}  not tested """
    #
    input=open(filename,'r')
    ionrec=input.readlines()
    dum=input.close()
    #
    # first get the number of data lines
    ndata=2
    iline=0
    while ndata > 1:
        s2=ionrec[iline].split()
        ndata=len(s2)
        iline=iline+1
    nline=iline-1
    #
    x=np.zeros(nline,'Float32')
    y=np.zeros(nline,'Float32')
    yerr=np.zeros(nline,'Float32')
#
    for iline in range(0,nline):
        ndata=len
        s2=ionrec[iline].split()
        ndata=len(s2)
        if ndata == 2:
            x[iline]=float(s2[0])
            y[iline]=float(s2[1])
            yerr[iline]=float(0.)
        else:
            x[iline]=float(s2[0])
            y[iline]=float(s2[1])
            yerr[iline]=float(s2[2])
    #
    ref=[]
    for iline in range(nline+1,len(ionrec)-1):
        s1a=ionrec[iline][:-1]
        ref.append(s1a.strip())


    ionrecdat={"x":x,"y":y,"yerr":yerr,"ref":ref}
    return ionrecdat
    #
    # --------------------------------------------------
    #
def z2element(z):
    """ convert Z to element string """
    if z-1 < len(const.El):
        thisel=const.El[z-1]
    else:
        thisel=''
    return thisel
    #
    # -------------------------------------------------------------------------------------
    #
def zion2name(z,ion, dielectronic=False):
    """ convert Z, ion to generic name  26, 13 -> fe_13 """
    if (z-1 < len(const.El)) and (ion <= z+1):
        thisone=const.El[z-1]+'_'+str(ion)
        if dielectronic:
            thisone+='d'
    else:
        thisone=''
    return thisone
    #
    # -------------------------------------------------------------------------------------
    #
def zion2filename(z,ion, dielectronic=False):
    """ convert Z to generic file name string """
    dir=os.environ["XUVTOP"]
    if (z-1 < len(const.El)) and (ion <= z+1):
        thisel=const.El[z-1]
    else:
        thisel=''
    if z-1 < len(const.El):
        thisone=const.El[z-1]+'_'+str(ion)
        if dielectronic:
            thisone+='d'
    else:
        thisone=''
    if thisel != '' :
        fname=os.path.join(dir,thisel,thisone,thisone)
    return fname
    #
    # -------------------------------------------------------------------------------------
    #
def zion2localFilename(z,ion, dielectronic=False):
    """ convert Z to generic file name string with current directory at top"""
    dir='.'
    if (z-1 < len(const.El)) and (ion <= z+1):
        thisel=const.El[z-1]
    else:
        thisel=''
    if z-1 < len(const.El):
        thisone=const.El[z-1]+'_'+str(ion)
        if dielectronic:
            thisone+='d'
    else:
        thisone=''
    if thisel != '' :
        fname=os.path.join(dir,thisel,thisone,thisone)
    return fname
    #
    # -------------------------------------------------------------------------------------
    #
def zion2spectroscopic(z,ion, dielectronic=False):
    """ convert Z and ion to spectroscopic notation string """
    if (z-1 < len(const.El)) and (ion <= z+1):
        spect=const.El[z-1].capitalize()+' '+const.Ionstage[ion-1]
        if dielectronic:
            spect+=' d'
    else:  spect = ''
    return spect
    #
    # -------------------------------------------------------------------------------------
    #
def convertName(name):
    """ convert ion name string to Z and Ion """
    s2=name.split('_')
    els=s2[0].strip()
    i1=const.El.index(els)+1
    ions=s2[1].strip()
    d=ions.find('d')
    if d >0 :
        dielectronic=True
        ions=ions.replace('d','')
    else: dielectronic=False
    return {'Z':int(i1),'Ion':int(ions),'Dielectronic':dielectronic, 'Element':els}
    #
    # -------------------------------------------------------------------------------------
    #
def defaultsRead(verbose=0):
    #
    #possibleDefaults = {'wavelength':['angstrom', 'kev', 'nm']}
    #symbolDefaults = {'wavelength':['A', 'keV', 'nm']}
    initDefaults={'abundfile': 'sun_photospheric','ioneqfile': 'chianti', 'wavelength': 'angstrom', 'flux': 'energy','gui':False}
    rcfile=os.path.join(os.environ['HOME'],'.chianti/chiantirc')
    if os.path.isfile(rcfile):
        print ' reading chiantirc file'
        config = RawConfigParser(initDefaults)
        config.read(rcfile)
        defaults = {}
        for anitem in config.items('chianti'):
            defaults[anitem[0]] = anitem[1]
        if defaults['gui'].lower() in ('t', 'y', 'yes', 'on', 'true', '1', 1, True):
            defaults['gui'] = True
        elif defaults['gui'].lower() in ('f', 'n', 'no', 'off', 'false', '0', 0, False):
            defaults['gui'] = False
    else:
        defaults = initDefaults
        if verbose:
            print ' chiantirc file (/HOME/.chianti/chiantirc) does not exist'
            print ' using the following defaults'
            for akey in defaults.keys():
                print akey, ' = ', defaults[akey]
    return defaults
    #
    # -------------------------------------------------------------------------------------
    #
def ion2filename(ions):
    """ convert ion string to generic file name string """
    dir=os.environ["XUVTOP"]
    zion=convertName(ions)
    el=z2element(zion['Z'])
    fname=os.path.join(dir,el,ions,ions)
    return fname
    #
    # -------------------------------------------------------------------------------------
    #
def el2z(els):
    """ from an the name of the element (1-2 letter) return Z"""
    z=const.El.index(els.lower())+1
    return z
    #
    # -------------------------------------------------------------------------------------
    #
def abundanceRead(abundancename=''):
    """ read an abundanc file and returns the abundance values relative to hydrogen"""
    abundance=np.zeros((50),'Float64')
    xuvtop=os.environ["XUVTOP"]
    if abundancename!='':
        # a specific abundance file name has been specified
        fname=os.path.join(xuvtop,'abundance',abundancename+'.abund')
    else:
        # the user will select an abundance file
        abundir=os.path.join(xuvtop,'abundance')
        abundlabel = 'ChiantiPy - Select an abundance file'
        fname = chpicker(abundir, filter='*.abund', label=abundlabel)
        if fname == None:
            print' no abundance file selected'
            return 0
        else:
            abundancefilename=os.path.basename(fname)
            abundancename,ext=os.path.splitext(abundancefilename)
#    else:
#        # the default abundance file will be used
#        abundancename=self.Defaults['abundfile']
#        fname=os.path.join(xuvtop,'abundance',abundancename+'.abund')
    input=open(fname,'r')
    s1=input.readlines()
    input.close()
    nlines=0
    idx=-1
    while idx <= 0:
        aline=s1[nlines][0:5]
        idx=aline.find('-1')
        nlines+=1
    nlines-=1
    for line in range(nlines):
        z,ab,element=s1[line].split()
        abundance[int(z)-1]=float(ab)
    gz=np.nonzero(abundance)
    abs=10.**(abundance[gz]-abundance[0])
    abundance.put(gz,abs)
    abundanceRef=s1[nlines+1:-2]
    return {'abundancename':abundancename,'abundance':abundance,'abundanceRef':abundanceRef}
    #
    # -------------------------------------------------------------------------------------
    #
def qrp(z,u):
    ''' qrp(Z,u)  u = E/IP
    calculate Qr-prime (equ. 2.12) of Fontes, Sampson and Zhang 1999'''
    #
    aa=1.13  # aa stands for A in equ 2.12
    #
    if z >= 16 :
        # use Fontes Z=20, N=1 parameters
        dd=3.70590
        c=-0.28394
        d=1.95270
        cc=0.20594
    else:
    # use Fontes Z=10, N=2 parameters
        dd=3.82652
        c=-0.80414
        d=2.32431
        cc=0.14424
    #
    if z > 20:
        cc+=((z-20.)/50.5)**1.11
    #
    bu=u <= 1.
    q=np.ma.array(u, 'float32', mask=bu, fill_value=0.)
    #
    #
    q=(aa*np.ma.log(u) + dd*(1.-1./u)**2 + cc*u*(1.-1./u)**4 + (c/u+d/u**2)*(1.-1/u))/u
    #
    q.set_fill_value(0.)  # I don't know why this is necessary
    return q  #  .set_fill_value(0.)
    #
    # -------------------------------------------------------------------------------------
    #

def elvlcRead(ions, filename = None, verbose=0,  useTh=0):
    """ read a chianti energy level file and returns
    {"lvl":lvl,"conf":conf,"term":term,"spin":spin,"l":l,"spd":spd,"j":j
    ,"mult":mult,"ecm":ecm,"eryd":eryd,"ecmth":ecmth,"erydth":erydth,"ref":ref,"pretty":pretty, 'ionS':ions}
    if a energy value for ecm or eryd is zero(=unknown), the theoretical values
    (ecmth and erydth) are inserted
    """
    #
    fstring='i3,i6,a15,i3,i3,a3,f4.1,i3,4f15.2'
    elvlcFormat=FortranFormat(fstring)
    #
    if type(filename) == NoneType:
        fname=ion2filename(ions)
        elvlname=fname+'.elvlc'
    else:
        elvlname = filename
        bname = os.path.basename(filename)
        ions = bname.split('.')[0]
    if not os.path.isfile(elvlname):
        print ' elvlc file does not exist:  ',elvlname
        return -1
    input=open(elvlname,'r')
    s1=input.readlines()
    input.close()
    nlvls=0
    ndata=2
    while ndata > 1:
        s1a=s1[nlvls][:-1]
        s2=s1a.split()
        ndata=len(s2)
        nlvls=nlvls+1
    nlvls-=1
    if verbose:
        print ' nlvls = ', nlvls
    lvl=[0]*nlvls
    conf=[0]*nlvls
    term=[0]*nlvls
    spin=[0]*nlvls
    l=[0]*nlvls
    spd=[0]*nlvls
    j=[0]*nlvls
    mult=[0]*nlvls
    ecm=[0]*nlvls
    eryd=[0]*nlvls
    ecmth=[0]*nlvls
    erydth=[0]*nlvls
    pretty=[0]*nlvls
    for i in range(0,nlvls):
        if verbose:
            print s1[i][0:115]
        inpt=FortranLine(s1[i][0:115],elvlcFormat)
        lvl[i]=inpt[0]
        conf[i]=inpt[1]
        term[i]=inpt[2]
        spin[i]=inpt[3]
        l[i]=inpt[4]
        spd[i]=inpt[5].strip()
        j[i]=inpt[6]
        mult[i]=inpt[7]
        ecm[i]=inpt[8]
        eryd[i]=inpt[9]
        ecmth[i]=inpt[10]
        erydth[i]=inpt[11]
        if ecm[i] == 0.:
            if useTh:
                ecm[i] = ecmth[i]
                eryd[i] = erydth[i]
        stuff = term[i].strip() + ' %1i%1s%3.1f'%( spin[i], spd[i], j[i])
        pretty[i] = stuff.strip()
    ref=[]
    for i in range(nlvls+1,len(s1)-1):
        s1a=s1[i][:-1]
        ref.append(s1a.strip())
#    self.const.Elvlc={"lvl":lvl,"conf":conf,"term":term,"spin":spin,"l":l,"spd":spd,"j":j
#            ,"mult":mult,"ecm":ecm,"eryd":eryd,"ecmth":ecmth,"erydth":erydth,"ref":ref}
    return {"lvl":lvl,"conf":conf,"term":term,"spin":spin,"l":l,"spd":spd,"j":j
            ,"mult":mult,"ecm":ecm,"eryd":eryd,"ecmth":ecmth,"erydth":erydth,"ref":ref,"pretty":pretty, 'ionS':ions}
    #
    # -------------------------------------------------------------------------------------
    #
def elvlcWrite(info, outfile=0, addLvl=0):
    ''' creates a .elvlc in the current directory
    info is a dictionary that must contain the following keys
    ionS, the Chianti style name of the ion such as c_4
    conf, an integer denoting the configuration - not too essential
    term, a string showing the configuration
    spin, an integer of the spin of the state in LS coupling
    l, an integer of the angular momentum quantum number
    spd, an string for the alphabetic symbol of the angular momemtum, S, P, D, etc
    j, a floating point number, the total angular momentum
    ecm, the observed energy in inverse cm, if unknown, the value is 0.
    eryd, the observed energy in Rydbergs, if unknown, the value is 0.
    ecmth, the calculated energy from the scattering calculation, in inverse cm
    erydth, the calculated energy from the scattering calculation in Rydbergs
    ref, the references in the literature to the data in the input info
    the output filename will be ionS+'.elvlc' unless outfile is specified
    addLvl is to add a constant value to the index of all levels
    '''
    gname = info['ionS']
    if outfile:
        elvlcName = outfile
    else:
        elvlcName = gname + '.elvlc'
    print ' elvlc file name = ', elvlcName
    out = open(elvlcName, 'w')
    for i,  conf in enumerate(info['conf']):
        mult = int(2.*info['j'][i]+1.)
        pstring = '%3i%6s%15s%3i%3i%2s%5.1f%3i%15.3f%15.6f%15.3f%15.6f \n'%(i+1+addLvl, conf, info['term'][i], info['spin'][i], info['l'][i], info['spd'][i], info['j'][i], mult, info['ecm'][i], info['eryd'][i], info['ecmth'][i], info['erydth'][i])
    #i3,a6,a15,2i3,a2,f5.1,i3,f15.3,f15.6,f15.3,f15.6
        out.write(pstring)
    out.write(' -1\n')
    out.write('%filename:  ' + elvlcName + '\n')
    info['ref'].append(' produced as a part of the George Mason University, University of Cambridge, University of Michigan \'CHIANTI\' atomic database for astrophysical spectroscopy consortium')
    today = date.today()
    info['ref'].append(' K. Dere (GMU) - ' + today.strftime('%Y %B %d'))
    for one in info['ref']:
        out.write(one+'\n')
    out.write(' -1\n')
    out.close()
    return
    #
    # -------------------------------------------------------------------------------------
    #
def wgfaRead(ions, filename = None):
    """ reads chianti wgfa file and returns
    {"lvl1":lvl1,"lvl2":lvl2,"wvl":wvl,"gf":gf,"avalue":avalue,"ref":ref}"""
    #
    if type(filename) == NoneType:
        fname=ion2filename(ions)
        wgfaname=fname+'.wgfa'
    else:
        wgfaname = filename
    input=open(wgfaname,'r')
    s1=input.readlines()
    dum=input.close()
    nwvl=0
    ndata=2
    while ndata > 1:
        s1a=s1[nwvl][:-1]
        s2=s1a.split()
        ndata=len(s2)
        nwvl=nwvl+1
    nwvl=nwvl-1
    lvl1=[0]*nwvl
    lvl2=[0]*nwvl
    wvl=[0.]*nwvl
    gf=[0.]*nwvl
    avalue=[0.]*nwvl
    #
    wgfaFormat='(2i5,f15.3,2e15.3)'
    for i in range(0,nwvl):
        inpt=FortranLine(s1[i],wgfaFormat)
        lvl1[i]=inpt[0]
        lvl2[i]=inpt[1]
        wvl[i]=inpt[2]
        gf[i]=inpt[3]
        avalue[i]=inpt[4]
    ref=[]
    for i in range(nwvl,len(s1)-1):
        s1a=s1[i][:-1]
        ref.append(s1a.strip())
    Wgfa={"lvl1":lvl1,"lvl2":lvl2,"wvl":wvl,"gf":gf,"avalue":avalue,"ref":ref, 'ionS':ions}
    return Wgfa
    #
    # --------------------------------------
    #
def wgfaWrite(info, outfile = 0, minBranch = 0.):
    ''' to write a wgfa file
    info is a dictionary the contains the following elements
    ionS, the Chianti style name of the ion such as c_4 for C IV
    lvl1 - the lower level, the ground level is 1
    lvl2 - the upper level
    wvl - the wavelength in Angstroms
    gf - the weighted oscillator strength
    avalue - the A value
    lower - descriptive text of the lower level (optional)
    upper - descriptive text of the upper level (optiona)
    ref - reference text, a list of strings
    minBranch:  the transition must have a branching ratio than the specified to be written to the file'''
    #
    gname = info['ionS']
    if outfile:
        wgfaname = outfile
    else:
        wgfaname = gname + '.wgfa'
    print ' wgfa file name = ', wgfaname
    if minBranch > 0.:
        info['ref'].append(' minimum branching ratio = %10.2e'%(minBranch))
    out = open(wgfaname, 'w')
    ntrans = len(info['lvl1'])
    nlvl = max(info['lvl2'])
    totalAvalue = np.zeros(nlvl, 'float64')
    if info.has_key('pretty1'):
        pformat = '%5i%5i%15.4f%15.3e%15.3e%30s - %30s'
    else:
        pformat = '%5i%5i%15.4f%15.3e%15.3e'
    for itrans, avalue in enumerate(info['avalue']):
        # for autoionization transitions, lvl1 can be less than zero
        if abs(info['lvl1'][itrans]) > 0 and info['lvl2'][itrans] > 0:
            totalAvalue[info['lvl2'][itrans] -1] += avalue

    for itrans, avalue in enumerate(info['avalue']):
        if avalue > 0.:
            branch = avalue/totalAvalue[info['lvl2'][itrans] -1]
        else:
            branch = 0.
        if branch > minBranch and abs(info['lvl1'][itrans]) > 0 and info['lvl2'][itrans] > 0:
            if info.has_key('pretty1'):
                pstring= pformat%(info['lvl1'][itrans], info['lvl2'][itrans], info['wvl'][itrans], info['gf'][itrans], avalue, info['pretty1'][itrans].rjust(30), info['pretty2'][itrans].ljust(30))
                out.write(pstring+'\n')
            else:
                pstring= pformat%(info['lvl1'][itrans], info['lvl2'][itrans], info['wvl'][itrans], info['gf'][itrans], avalue)
                out.write(pstring+'\n')
    out.write(' -1\n')
    out.write('%filename:  ' + wgfaname + '\n')
    for one in info['ref']:
        out.write(one+'\n')
    out.write(' -1\n')
    out.close()
    #
    # -------------------------------------------------------------------------------------
    #
def easplomRead(ions, extension='.splom'):
    """read chianti splom files and returns
    {"lvl1":lvl1,"lvl2":lvl2,"deryd":de,"gf":gf,"eryd":eout,"omega":omout}
    currently only works for 5 point spline fit files"""
    #
    #
    fname=ion2filename(ions)
    omname=fname+extension
    input=open(omname,'r')
    lines=input.readlines()
    input.close()
    format=FortranFormat('5i3,8e10.3')
    data=5
    iline=0
    lvl1=[]
    lvl2=[]
    ttype=[]
    gf=[]
    de=[]
    om=[]
    z=1
    while z > 0:
        omdat1=FortranLine(lines[iline],format)
        z=omdat1[0]
        if z > 0:
            l1=omdat1[2]
            l2=omdat1[3]
            ttype1=omdat1[4]
            gf1=omdat1[5]
            de1=omdat1[6]
            btf1=omdat1[7]
            om1=omdat1[8:]
            #
            lvl1.append(l1)
            lvl2.append(l2)
            ttype.append(ttype1)
            gf.append(gf1)
            de.append(de1)
            om.append(om1)
        iline=iline+1
    omout=np.asarray(om,'Float32')
    ref=lines[iline:-1]
#        omout=np.transpose(omout)
    if extension == '.omdat':
        Splom={"lvl1":lvl1,"lvl2":lvl2,'ttype':ttype,"gf":gf, "deryd":de,"omega":omout, 'ref':ref}
        return Splom
    elif  extension == '.easplom':
        Easplom={"lvl1":lvl1,"lvl2":lvl2,'ttype':ttype,"gf":gf, "deryd":de,"omega":omout, 'ref':ref}
        return Easplom

    return
    #
    # --------------------------------------------------
    #
def splomRead(ions, filename=None):
    """read chianti .splom files and return
    {"lvl1":lvl1,"lvl2":lvl2,"ttype":ttype,"gf":gf,"deryd":de,"f":f,"splom":splomout,"ref":hdr} not tested """
    #
    if type(filename) == NoneType:
        fname=ion2filename(ions)
        splomname=fname+'.splom'
    else:
        splomname = filename
    input=open(splomname,'r')
    #  need to read first line and see how many elements
    line1=input.readline()
    indices=line1[0:15]
    remainder=line1[16:]
    nom=remainder.split(' ')
    format=FortranFormat('5i3,'+str(len(nom))+'E10.2')
    #  go back to the beginning
    input.seek(0)
    lines=input.readlines()
    data=5
    iline=0
    lvl1=[]
    lvl2=[]
    ttype=[]
    gf=[]
    de=[]
    f=[]
    splom=[]
    ntrans=0
    while data > 1:
        splomdat=FortranLine(lines[iline],format)
        l1=splomdat[2]
        l2=splomdat[3]
        tt1=splomdat[4]
        gf1=splomdat[5]
        de1=splomdat[6]
        f1=splomdat[7]
        splom1=splomdat[8:]
        lvl1.append(int(l1))
        lvl2.append(int(l2))
        ttype.append(int(tt1))
        gf.append(float(gf1))
        de.append(float(de1))
        f.append(float(f1))
        splom.append(splom1)
        iline=iline+1
        data=len(lines[iline].split(' ',2))
    hdr=lines[iline+1:-1]
    de=np.asarray(de,'Float32')
    splomout=np.asarray(splom,'Float32')
    splomout=np.transpose(splomout)
    input.close()
    # note:  de is in Rydbergs
    splom={"lvl1":lvl1,"lvl2":lvl2,"ttype":ttype,"gf":gf,"deryd":de,"f":f
        ,"splom":splomout,"ref":hdr}
    return  splom
    #
    # --------------------------------------------------
    #
def splupsRead(ions, filename=None, prot=0, ci=0,  diel=0):
    """read a chianti splups file and return
    {"lvl1":lvl1,"lvl2":lvl2,"ttype":ttype,"gf":gf,"de":de,"cups":cups,"bsplups":bsplups,"ref":ref}
    if prot >0, then reads the psplups file
    if ci > 0, then reads cisplups file
    if diel > 0, then reads dielsplups file"""
    #
    if type(filename) == NoneType:
        fname=ion2filename(ions)
        if prot:
            splupsname=fname+'.psplups'
        elif ci:
            splupsname=fname+'.cisplups'
        elif diel:
            splupsname=fname+'.dielsplups'
        else:
            splupsname=fname+'.splups'
    else:
        splupsname = filename
    if not os.path.exists(splupsname):
        if prot:
            return None
        elif ci:
            return None
        else:
            return None
    # there is splups/psplups data
    else:
        input=open(splupsname,'r')
        s1=input.readlines()
        input.close()
        nsplups=0
        ndata=2
        while ndata > 1:
            s1a=s1[nsplups][:]
            s2=s1a.split()
            ndata=len(s2)
            nsplups=nsplups+1
        nsplups=nsplups-1
        lvl1=[0]*nsplups
        lvl2=[0]*nsplups
        ttype=[0]*nsplups
        gf = np.zeros(nsplups, 'float64')
        de= np.zeros(nsplups, 'float64')
        cups= np.zeros(nsplups, 'float64')
        nspl=[0]*nsplups
        splups=np.zeros((nsplups,9),'Float64')
        if prot:
            splupsFormat1='(3i3,8e10.3)'
            splupsFormat2='(3i3,12e10.3)'
        else:
#            splupsFormat1='(6x,3i3,8e10.3)'
            splupsFormat2='(6x,3i3,12e10.3)'
        #
        for i in range(0,nsplups):
            # checking for 5 or 9 point splines, need to check for 9 first!
#            try:
#                inpt=FortranLine(s1[i],splupsFormat2)
#                nspl[i] = 9
#                print ' 9 = ', inpt[6:]
#            except:
#                inpt=FortranLine(s1[i],splupsFormat1)
#                nspl[i] = 5
#                print ' 5 = ', inpt[6:]
#
            inpt=FortranLine(s1[i],splupsFormat2)
            lvl1[i]=inpt[0]
            lvl2[i]=inpt[1]
            ttype[i]=inpt[2]
            gf[i]=inpt[3]
            de[i]=inpt[4]
            cups[i]=inpt[5]
            spl=np.array(inpt[6:])
            good = spl > 0.
#            if good[5:].sum():
#                nspl[i] = 9
#            else:
#                nspl[i] = 5
            nspl[i] = good.sum()
            splups[i].put(range(nspl[i]), spl)
        #
        ref=[]
        for i in range(nsplups+1,len(s1)-1):
            s1a=s1[i][:-1]
            ref.append(s1a.strip())
        if prot:
#            self.Npsplups=nsplups
#            self.Psplups={"lvl1":lvl1,"lvl2":lvl2,"ttype":ttype,"gf":gf,"de":de,"cups":cups
#                ,"nspl":nspl,"splups":splups,"ref":ref}
            return {"lvl1":lvl1,"lvl2":lvl2,"ttype":ttype,"gf":gf,"de":de,"cups":cups
                ,"nspl":nspl,"splups":splups,"ref":ref}
        else:
#            self.Splups={"lvl1":lvl1,"lvl2":lvl2,"ttype":ttype,"gf":gf,"de":de,"cups":cups
#                ,"nspl":nspl,"splups":splups,"ref":ref}
            return {"lvl1":lvl1,"lvl2":lvl2,"ttype":ttype,"gf":gf,"de":de,"cups":cups
                ,"nspl":nspl,"splups":splups,"ref":ref}
    #
    # -------------------------------------------------------------------------------------
    #
def cireclvlRead(ions, type):
    ''' to read Chianti cilvl and reclvl files and return data
    must specify type as either cilvl or reclvl'''
    fname=ion2filename(ions)
    if type in ['cilvl', 'reclvl']:
        paramname=fname+'.' + type
    else:
        print('type must be one of "cilvl" or "reclvl" ')
        return {}
    if os.path.exists(paramname):
        input=open(paramname,'r')
        lines = input.readlines()
        input.close()
    else:
        return None
    #
    iline = 0
    idx = -1
    while idx <= 0:
        aline=lines[iline][0:5]
        idx=aline.find('-1')
        iline += 1
    ndata = iline - 1
    ntrans = ndata/2
    #
    nref = 0
    idx = -1
    while idx <= 0:
        aline=lines[iline][0:5]
        idx=aline.find('-1')
        iline += 1
        nref += 1
    nref -= 1
    #
    # need to find the maximum number of temperatures, not all lines are the same
    #
    ntemp = np.zeros(ntrans, 'int32')
    iline = 0
    for jline in range(0, ndata, 2):
        dummy = lines[jline].replace(os.linesep, '').split()
        ntemp[iline] = len(dummy[4:])
        iline += 1
    maxNtemp = ntemp.max()
#   print ' maxNtemp = ', maxNtemp
    temp = np.zeros((ntrans,maxNtemp), 'float64')
    iline = 0
    for jline in range(0, ndata, 2):
        recdat = lines[jline].replace(os.linesep, '').split()
        shortT = np.asarray(recdat[4:], 'float64')
        # the result of the next statement is to continue to replicate t
        t = np.resize(shortT, maxNtemp)
        temp[iline] = 10.**t
        iline += 1
    #
    lvl1 = np.zeros(ntrans, 'int64')
    lvl2 = np.zeros(ntrans, 'int64')
    ci = np.zeros((ntrans, maxNtemp), 'float64')
    #
    idat = 0
    for jline in range(1, ndata, 2):
        cidat = lines[jline].replace(os.linesep, '').split()
        shortCi = np.asarray(cidat[4:], 'float64')
        lvl1[idat] = int(cidat[2])
        lvl2[idat] = int(cidat[3])
        ci[idat] = np.resize(shortCi, maxNtemp)
        idat += 1
    return {'temperature':temp, 'ntemp':ntemp,'lvl1':lvl1, 'lvl2':lvl2, 'rate':ci,'ref':lines[ndata+1:-1], 'ionS':ions}
    #
def dilute(radius):
    ''' to calculate the dilution factor as a function distance from the center of a star in units of the stellar radius
    a radius of less than 1.0 (incorrect) results in a dilution factor of 0.'''
    if radius >= 1.:
        d = 0.5*(1. - np.sqrt(1. - 1./radius**2))
    else:
        d = 0.
    return d
    #
def diRead(ions):
    """read chianti direct ionization .params files and return
        {"info":info,"btf":btf,"ev1":ev1,"xsplom":xsplom,"ysplom":ysplom,"ref":hdr}
        info={"iz":iz,"ion":ion,"nspl":nspl,"neaev":neaev}"""
    #
    zion=convertName(ions)
    if zion['Z'] < zion['Ion']:
        print ' this is a bare nucleus that has no ionization rate'
        return
    #
    fname=ion2filename(ions)
    paramname=fname+'.diparams'
    input=open(paramname,'r')
    #  need to read first line and see how many elements
    line1=input.readline()
    indices=line1.split()
    iz=int(indices[0])
    ion=int(indices[1])
    nspl=indices[2]
    nfac=int(indices[3])
    neaev=int(indices[4])
    nspl=int(nspl)
    format=FortranFormat(str(nspl+1)+'E10.2')
    #
    ev1=np.zeros(nfac,'Float32')
    btf=np.zeros(nfac,'Float32')
    xsplom=np.zeros([nfac, nspl],'Float32')
    ysplom=np.zeros([nfac, nspl],'Float32')
    #
    for ifac in range(nfac):
        line=input.readline()
        paramdat=FortranLine(line,format)
        btf[ifac]=paramdat[0]
        xsplom[ifac]=paramdat[1:]
        line=input.readline()
        paramdat=FortranLine(line,format)
        ev1[ifac]=paramdat[0]
        ysplom[ifac]=paramdat[1:]
    if neaev:
        line=input.readline()
        eacoef=line.split()
#            print ' eaev = ', type(eacoef), eacoef
        eaev=[float(avalue) for avalue in eacoef]
#            print ' eaev = ', type(eaev), eaev
#            print ' eaev = ', type(eaev), eaev
#            if len(eaev) == 1:
#                eaev=float(eaev[0])
#                eaev=np.asarray(eaev, 'float32')
#            else:
#                eaev=np.asarray(eaev, 'float32')
    else:
        eaev=0.
    hdr=input.readlines()
    input.close()
    info={"iz":iz,"ion":ion,"nspl":nspl,"neaev":neaev, 'nfac':nfac}
    DiParams={"info":info,"btf":btf,"ev1":ev1,"xsplom":xsplom,"ysplom":ysplom, 'eaev':eaev,"ref":hdr}
    return DiParams
    #
    # -------------------------------------------------------------------------------------
    #
def eaRead(ions):
    '''read a chianti excitation-autoionization file  return ionization rate
    derived from splupsRead
    {"lvl1":lvl1,"lvl2":lvl2,"ttype":ttype,"gf":gf,"de":de,"cups":cups,"bsplups":bsplups,"ref":ref}'''
    fname=ion2filename(ions)
    splupsname=fname+'.easplups'
    if not os.path.exists(splupsname):
        print ' could not find file:  ', splupsname
        self.Splups={"lvl1":-1}
        return {"lvl1":-1}
    # there is splups/psplups data
    else:
        input=open(splupsname,'r')
        s1=input.readlines()
        dum=input.close()
        nsplups=0
        ndata=2
        while ndata > 1:
            s1a=s1[nsplups][:]
            s2=s1a.split()
            ndata=len(s2)
            nsplups=nsplups+1
        nsplups=nsplups-1
        lvl1=[0]*nsplups
        lvl2=[0]*nsplups
        ttype=[0]*nsplups
        gf=[0.]*nsplups
        de=[0.]*nsplups
        cups=[0.]*nsplups
        nspl=[0]*nsplups
        splups=np.zeros((nsplups,9),'Float32')
        splupsFormat1='(6x,3i3,8e10.0)'
        splupsFormat2='(6x,3i3,12e10.0)'
        #
        for i in range(0,nsplups):
            try:
                inpt=FortranLine(s1[i],splupsFormat1)
            except:
                inpt=FortranLine(s1[i],splupsFormat2)
            lvl1[i]=inpt[0]
            lvl2[i]=inpt[1]
            ttype[i]=inpt[2]
            gf[i]=inpt[3]
            de[i]=inpt[4]
            cups[i]=inpt[5]
            if len(inpt)  > 13:
                nspl[i]=9
                splups[i].put(range(9),inpt[6:])
            else:
                nspl[i]=5
                splups[i].put(range(5),inpt[6:])
        #
        ref=[]
        for i in range(nsplups+1,len(s1)-1):
            s1a=s1[i][:-1]
            ref.append(s1a.strip())
#        self.EaParams={"lvl1":lvl1,"lvl2":lvl2,"ttype":ttype,"gf":gf,"de":de,"cups":cups
#                ,"nspl":nspl,"splups":splups,"ref":ref}
        return {"lvl1":lvl1,"lvl2":lvl2,"ttype":ttype,"gf":gf,"de":de,"cups":cups
                ,"nspl":nspl,"splups":splups,"ref":ref}
    #
    # -------------------------------------------------------------------------------------
    #
def rrRead(ions):
    """read chianti radiative recombination .rrparams files and return
        {'rrtype','params','ref'}"""
    #
    #
    fname=ion2filename(ions)
    paramname=fname+'.rrparams'
    if os.path.isfile(paramname):
        input=open(paramname,'r')
        #  need to read first line and see how many elements
        lines=input.readlines()
        input.close()
        rrtype=int(lines[0])
        ref=lines[3:-2]
        #
        if rrtype == 1:
            # a Badnell type
            fmt=FortranFormat('3i5,e12.4,f10.5,2e12.4')
            params=FortranLine(lines[1],fmt)
            RrParams={'rrtype':rrtype, 'params':params, 'ref':ref}
        elif rrtype == 2:
            # a Badnell type
            fmt=FortranFormat('3i5,e12.4,f10.5,2e11.4,f10.5,e12.4')
            params=FortranLine(lines[1],fmt)
            RrParams={'rrtype':rrtype, 'params':params, 'ref':ref}
        elif rrtype == 3:
            # a Shull type
            fmt=FortranFormat('2i5,2e12.4')
            params=FortranLine(lines[1],fmt)
            RrParams={'rrtype':rrtype, 'params':params, 'ref':ref}
        else:
            RrParams=None
            print ' for ion %5s unknown RR type = %5i' %(ions, rrtype)
        return RrParams
    else:
        return {'rrtype':-1}

    #
    # -------------------------------------------------------------------------------------
    #
def drRead(ions):
    """read chianti dielectronic recombination .drparams files and return
        {'rrtype','params','ref'}"""
    #
    #
    fname=ion2filename(ions)
    paramname=fname+'.drparams'
    if os.path.isfile(paramname):
        input=open(paramname,'r')
        #  need to read first line and see how many elements
        lines=input.readlines()
        input.close()
        drtype=int(lines[0])
        ref=lines[4:-1]
        #
        if drtype == 1:
            # a Badnell type
            fmt=FortranFormat('2i5,8e12.4')
            eparams=np.asarray(FortranLine(lines[1],fmt)[2:], 'float64')
            cparams=np.asarray(FortranLine(lines[2],fmt)[2:], 'float64')
            DrParams={'drtype':drtype, 'eparams':eparams,'cparams':cparams,  'ref':ref}
        elif drtype == 2:
            # shull type
            fmt=FortranFormat('2i5,4e12.4')
            params=np.asarray(FortranLine(lines[1],fmt)[2:], 'float64')
            DrParams={'drtype':drtype, 'params':params, 'ref':ref}
        else:
            DrParams = None
            print ' for ion %5s unknown DR type = %5i' %(ions, drtype)
    else:
        DrParams=None
    return DrParams
    #
    # -------------------------------------------------------------------------------------
    #
def ioneqRead(ioneqname=''):
    """ reads an ioneq file and stores temperatures and ionization
    equilibrium values in self.IoneqTemperature and self.Ioneq and returns
    a dictionary containing these value and the reference to the literature"""
    dir=os.environ["XUVTOP"]
    if ioneqname == '':
        # the user will select an ioneq file
        ioneqdir=os.path.join(dir,'ioneq')
        fname=chpicker(ioneqdir,filter='*.ioneq',label = 'Select an Ionization Equilibrium file')
        if fname ==None:
            print' no ioneq file selected'
            return False
        else:
            ioneqfilename=os.path.basename(fname)
            ioneqname,ext=os.path.splitext(ioneqfilename)
    else:
        fname=os.path.join(dir,'ioneq',ioneqname+'.ioneq')
        if not os.path.isfile(fname):
            print ' file does not exist:  ', fname
            path=os.path.join(dir,'ioneq')
            filelist=os.listdir(path)
            ioneqlist=[]
            for aname in fnmatch.filter(filelist,'*.ioneq'):
                ioneqlist.append(aname)
            print ' the following files do exist:'
            for one in ioneqlist:
                print ' - ', one.replace('.ioneq', '')
            return False
    #
    input=open(fname,'r')
    s1=input.readlines()
    input.close()
    ntemp,nele=s1[0].replace('  ',' ').strip().split(' ')
    nTemperature=int(ntemp)
    nElement=int(nele)
    tformat=FortranFormat(str(nTemperature)+'f6.2')
    ioneqTemperature=FortranLine(s1[1],tformat)
    ioneqTemperature=np.asarray(ioneqTemperature[:],'Float32')
    ioneqTemperature=10.**ioneqTemperature
    nlines=0
    idx=-1
    while idx <= 0:
        aline=s1[nlines][0:5]
        idx=aline.find('-1')
        nlines+=1
    nlines-=1
    #
    ioneqformat=FortranFormat('2i3,'+str(nTemperature)+'e10.2')
    #
    ioneqAll=np.zeros((nElement,nElement+1,nTemperature),'Float32')
    for iline in range(2,nlines):
        out=FortranLine(s1[iline],ioneqformat)
        iz=out[0]
        ion=out[1]
        ioneqAll[iz-1,ion-1].put(range(nTemperature),np.asarray(out[2:],'Float32'))
    ioneqRef=s1[nlines+1:-2]
    del s1
    return {'ioneqname':ioneqname,'ioneqAll':ioneqAll,'ioneqTemperature':ioneqTemperature,'ioneqRef':ioneqRef}
    #
    # -------------------------------------------------------------------------------------
    #
def gffRead():
    '''to read the free-free gaunt factors of Sutherland, 1998, MNRAS, 300, 321.
    this function reads the file and reverses the values of g2 and u'''
    xuvtop = os.environ['XUVTOP']
    fileName = os.path.join(xuvtop, 'continuum','gffgu.dat' )
    input = open(fileName)
    lines = input.readlines()
    input.close()
    #
    #  the 1d stuff below is to make it easy to use interp2d
    ngamma=41
    nu=81
    nvalues = ngamma*nu
    g2 = np.zeros(ngamma, 'float64')
    g21d = np.zeros(nvalues, 'float64')
    u = np.zeros(nu, 'float64')
    u1d = np.zeros(nvalues, 'float64')
    gff = np.zeros((ngamma, nu), 'float64')
    gff1d = np.zeros(nvalues, 'float64')
    #
    iline = 5
    ivalue = 0
    for ig2 in range(ngamma):
        for iu in range(nu):
            values = lines[iline].split()
            u[iu] = float(values[1])
            u1d[ivalue] = float(values[1])
            g2[ig2] = float(values[0])
            g21d[ivalue] = float(values[0])
            gff[ig2, iu] = float(values[2])
            gff1d[ivalue] = float(values[2])
            iline+=1
            ivalue += 1
    #
    return {'g2':g2, 'g21d':g21d,  'u':u, 'u1d':u1d,  'gff':gff,  'gff1d':gff1d}
    #
    # ----------------------------------------------------------------------------------------
    #
def gffintRead():
    '''to read the integrated free-free gaunt factors of Sutherland, 1998, MNRAS, 300, 321.'''
    xuvtop = os.environ['XUVTOP']
    fileName = os.path.join(xuvtop, 'continuum','gffint.dat' )
    input = open(fileName)
    lines = input.readlines()
    input.close()
    #
    ngamma=41
    g2 = np.zeros(ngamma, 'float64')
    gffint = np.zeros(ngamma, 'float64')
    s1 = np.zeros(ngamma, 'float64')
    s2 = np.zeros(ngamma, 'float64')
    s3 = np.zeros(ngamma, 'float64')
    #
    ivalue = 0
    start = 4
    for iline in range(start,start+ngamma):
        values = lines[iline].split()
        g2[ivalue] = float(values[0])
        gffint[ivalue] = float(values[1])
        s1[ivalue] = float(values[2])
        s2[ivalue] = float(values[3])
        s3[ivalue] = float(values[4])
        ivalue += 1
    #
    return {'g2':g2, 'gffint':gffint, 's1':s1, 's2':s2, 's3':s3}
    #
    # ----------------------------------------------------------------------------------------
    #
def itohRead():
    '''to read in the free-free gaunt factors of Itoh et al. (ApJS 128, 125, 2000)'''
    xuvtop = os.environ['XUVTOP']
    itohName = os.path.join(xuvtop, 'continuum', 'itoh.dat')
    input = open(itohName)
    lines = input.readlines()
    input.close()
    gff = np.zeros((30, 121), 'float64')
    for iline in range(30):
        gff[iline]= np.asarray(lines[iline].split(), 'float64')
    return {'itohCoef':gff}
    #
    #
    # ----------------------------------------------------------------------------------------
    #
def klgfbRead():
    '''to read CHIANTI files file containing the free-bound gaunt factors for n=1-6 from Karzas and Latter, 1961, ApJSS, 6, 167
    returns {pe, klgfb}, the photon energy and the free-bound gaunt factors'''
    xuvtop = os.environ['XUVTOP']
    fname = os.path.join(xuvtop, 'continuum', 'klgfb.dat')
    input = open(fname)
    lines = input.readlines()
    input.close()
    #
    ngfb = int(lines[0].split()[0])
    nume = int(lines[0].split()[1])

    gfb = np.zeros((ngfb, ngfb, nume), 'float64')
    nlines = len(lines)
#        print 'nlines, nume, ngfb = ', nlines,  nume, ngfb
    pe = np.asarray(lines[1].split(), 'float64')
    for iline in range(2, nlines):
        data = lines[iline].split()
        n = int(data[0])
        l = int(data[1])
        gfb[n-1, l] = np.array(data[2:], 'float64')
    return {'pe':pe, 'klgfb':gfb}
    #
    # ----------------------------------------------------------------------------------------
    #
def fblvlRead(filename, verbose=False):
    """ read a chianti energy level file and returns
    {"lvl":lvl,"conf":conf,"term":term,"spin":spin,"l":l,"spd":spd,"j":j
    ,"mult":mult,"ecm":ecm,"eryd":eryd,"ref":ref}"""
#        #  ,format='(i5,a20,2i5,a3,i5,2f20.3)'
    fstring='i5,a20,2i5,a3,i5,2f20.3'
    elvlcFormat=FortranFormat(fstring)
    #
    if os.path.exists(filename):
        input=open(filename,'r')
        s1=input.readlines()
        input.close()
        nlvls=0
        ndata=2
        while ndata > 1:
            s1a=s1[nlvls][:-1]
            s2=s1a.split()
            ndata=len(s2)
            nlvls=nlvls+1
        nlvls-=1
        if verbose:
            print ' nlvls = ', nlvls
        lvl=[0]*nlvls
        conf=[0]*nlvls
        pqn=[0]*nlvls
        l=[0]*nlvls
        spd=[0]*nlvls
        mult=[0]*nlvls
        ecm=[0]*nlvls
        ecmth=[0]*nlvls
        for i in range(0,nlvls):
            if verbose:
                print s1[i]
            inpt=FortranLine(s1[i],elvlcFormat)
            lvl[i]=inpt[0]
            conf[i]=inpt[1].strip()
            pqn[i]=inpt[2]
            l[i]=inpt[3]
            spd[i]=inpt[4].strip()
            mult[i]=inpt[5]
            if inpt[6] == 0.:
                ecm[i]=inpt[7]
            else:
                ecm[i]=inpt[6]
                ecmth[i]=inpt[7]
        ref=[]
        for i in range(nlvls+1,len(s1)-1):
            s1a=s1[i][:-1]
            ref.append(s1a.strip())
        return {"lvl":lvl,"conf":conf,'pqn':pqn,"l":l,"spd":spd,"mult":mult,
            "ecm":ecm,'ecmth':ecmth, 'ref':ref}
    else:
        return {'errorMessage':' fblvl file does not exist'}
    #
    # -----------------------------------------------------------------
    #
def vernerRead():
    '''Reads the Verner & Yakovlev (A&AS 109, 125, 1995) photoionization cross-section data'''
    xuvtop = os.environ['XUVTOP']
    fname = os.path.join(xuvtop, 'continuum', 'verner_short.txt')
    input = open(fname)
    lines = input.readlines()
    input.close()
    #
    nlines=465
    maxZ = 30+1
    maxNel = 30 +1# also equal max(stage)
    #
    #z = np.array(nlines,'int32')
    #nel = np.array(nlines,'int32')
    pqn = np.zeros((maxZ,maxNel),'int32')
    l = np.zeros((maxZ,maxNel),'int32')
    eth = np.zeros((maxZ,maxNel),'float64')
    e0 = np.zeros((maxZ,maxNel),'float64')
    sig0 = np.zeros((maxZ,maxNel),'float64')
    ya = np.zeros((maxZ,maxNel),'float64')
    p = np.zeros((maxZ,maxNel),'float64')
    yw = np.zeros((maxZ,maxNel),'float64')
    #
    fstring='i2,i3,i2,i2,6f11.3'
    vernerFormat=FortranFormat(fstring)
    #
    for iline in range(nlines):
        out=FortranLine(lines[iline],vernerFormat)
        z = out[0]
        nel = out[1]
        stage = z - nel + 1
        pqn[z,stage] = out[2]
        l[z,stage] = out[3]
        eth[z,stage] = out[4]
        e0[z,stage] = out[5]
        sig0[z,stage] = out[6]
        ya[z,stage] = out[7]
        p[z,stage] = out[8]
        yw[z,stage] = out[9]
    #
    return {'pqn':pqn, 'l':l, 'eth':eth, 'e0':e0, 'sig0':sig0, 'ya':ya, 'p':p, 'yw':yw}
    #
    #-----------------------------------------------------------
    #
def twophotonHRead():
    ''' to read the two-photon A values and distribution function for the H seq'''
    xuvtop = os.environ['XUVTOP']
    fName = os.path.join(xuvtop, 'continuum', 'hseq_2photon.dat')
    dFile = open(fName, 'r')
    a = dFile.readline()
    y0 = np.asarray(a.split())
    a = dFile.readline()
    z0 = np.asarray(a.split())
    nz = 30
    avalue = np.zeros(nz, 'float64')
    asum = np.zeros(nz, 'float64')
    psi0 = np.zeros((nz, 17), 'float64')
    for iz in range(nz):
        a=dFile.readline().split()
        avalue[iz] = float(a[1])
        asum[iz] = float(a[2])
        psi = np.asarray(a[3:])
        psi0[iz] = psi
    dFile.close()
    return {'y0':y0, 'z0':z0, 'avalue':avalue, 'asum':asum, 'psi0':psi0.reshape(30, 17)}
    #
    #-----------------------------------------------------------
    #
def twophotonHeRead():
    ''' to read the two-photon A values and distribution function for the He seq'''
    xuvtop = os.environ['XUVTOP']
    fName = os.path.join(xuvtop, 'continuum', 'heseq_2photon.dat')
    dFile = open(fName, 'r')
    a = dFile.readline()
    y0 = np.asarray(a.split())
    nz = 30
    avalue = np.zeros(nz, 'float64')
    psi0 = np.zeros((nz, 41), 'float64')
    for iz in range(1, nz):
        a=dFile.readline().split()
        avalue[iz] = float(a[1])
        psi = np.asarray(a[2:])
        psi0[iz] = psi
    dFile.close()
    return {'y0':y0, 'avalue':avalue, 'psi0':psi0.reshape(30, 41)}
    #
    #-----------------------------------------------------------
    #
def versionRead():
    """ read the version number of the CHIANTI database"""
    xuvtop = os.environ['XUVTOP']
    vFileName = os.path.join(xuvtop, 'VERSION')
    vFile = open(vFileName)
    versionStr = vFile.readline()
    vFile.close()
    return versionStr.strip()
    #
    #-----------------------------------------------------------
    #
def scale_bti(evin,crossin,f,ev1):
    """
    apply BT ionization scaling to (energy, cross-section)
    returns [bte,btx]
    """
    u=evin/ev1
    bte=1.-np.log(f)/np.log(u-1.+f)
    btx=u*crossin*(ev1**2)/(np.log(u)+1.)
    return [bte,btx]
    #
    #-----------------------------------------------------------
    #
def descale_bti(bte,btx,f,ev1):
    """
    descale BT ionization scaling
    returns [energy,cross-section]
    """
    u=1.-f+np.exp(np.log(f)/(1.-bte))
    energy=u*ev1
    cross=(np.log(u)+1.)*btx/(u*ev1**2)
    return [energy,cross]
    #
    #-----------------------------------------------------------
    #
def descale_bt(bte,btomega,f,ev1):
    """
    descale BT excitation scaling
    returns [energy,collision strength]
    """
    u=1.-f+np.exp(np.log(f)/(1.-bte))
    energy=u*ev1
    omega=(np.log(u)-1.+np.exp(1.))*btomega
    return [energy,omega]
    #
    #-----------------------------------------------------------
    #
def scale_bt(evin,omega,f,ev1):
    """
    apply BT excitation scaling to (energy, collision strength)
    returns [bte,btomega]
    """
    u=evin/ev1
    bte=1.-np.log(f)/np.log(u-1.+f)
    btomega=omega/(np.log(u)-1.+np.exp(1.))
    return [bte,btomega]
