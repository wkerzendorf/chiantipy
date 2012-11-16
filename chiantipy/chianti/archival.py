from FortranFormat import *
import chianti.constants as const
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
        return {'status':0}
    status = 1
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
        term[i]=inpt[2].strip()
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
            ,"mult":mult,"ecm":ecm,"eryd":eryd,"ecmth":ecmth,"erydth":erydth,"ref":ref,"pretty":pretty, 'ionS':ions, 'status':status}
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
        thisTerm = info['term'][i].ljust(14)
        pstring = '%3i%6s%15s%3i%3i%2s%5.1f%3i%15.3f%15.6f%15.3f%15.6f \n'%(i+1+addLvl, conf, thisTerm, info['spin'][i], info['l'][i], info['spd'][i], info['j'][i], mult, info['ecm'][i], info['eryd'][i], info['ecmth'][i], info['erydth'][i])
    #i3,a6,a15,2i3,a2,f5.1,i3,f15.3,f15.6,f15.3,f15.6
        out.write(pstring)
    out.write(' -1\n')
    out.write('%filename:  ' + elvlcName + '\n')
    for one in info['ref']:
        out.write(one+'\n')
    out.write(' -1\n')
    out.close()
    return
