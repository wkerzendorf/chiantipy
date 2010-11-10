import chianti
import chianti.core as ch
import numpy as np
def splups2rate(ionS, temperature):
    ''' this is primarily to get rates from the dielectronic splups files'''
    ion = ch.ion(ionS, temperature)
    ion.upsilonDescale()
    # since this is a dielectronic ion, the ionization potential is already subtracted
    # from the excitation energy
    tform = '%3i%3i%3i%4i'
    if isinstance(temperature, float):
        nTemp = 1
    else:
        nTemp = len(temperature)
    tempForm = tform+nTemp*'%10.3f'
    rateForm = tform+nTemp*'%10.3e'
    for isplups in range(ion.Nsplups):
        l1 = ion.Splups['lvl1'][isplups]
        l2 = ion.Splups['lvl2'][isplups]
        print tempForm%(ion.Z, ion.Ion, l1, l2, np.log10(temperature))
        print rateForm%(ion.Z, ion.Ion, l1, l2, ion.Upsilon['exRate'][isplups])
