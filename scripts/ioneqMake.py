# a python script to create ioneq files
#
from datetime import date
import numpy as np
import chianti
import chianti.core as ch
#
chRef = ['produced as a part of the George Mason University, University of Cambridge, University of Michigan \'CHIANTI\' atomic database for astrophysical spectroscopy consortium']
today = date.today()
chRef.append(' K. Dere (GMU) - ' + today.strftime('%Y %B %d'))
#
nt = 107
nz = 30
tlog = 3.7 + 0.05*np.arange(nt)
t=10.**tlog
print ' tmin, tmax = ', t.min(), t.max()
#
out = open('chianti_new.ioneq', 'w')
pstring = '%7i %7i  \n'%(nt, nz)
print ' pstring = ', pstring
out.write(pstring)
#
pstring = ''
for it in range(nt):
    pstring += '%6.2f'%(tlog[it])
pstring += '   \n'
out.write(pstring)
#
for iz in range(nz):
    ion = ch.ioneq(iz+1, t)
#    print ' ioneq.shape = ', ion.Ioneq.shape
    for stage in range(iz+2):
        pstring = '%3i%3i'%(iz+1, stage+1)
        for it in range(nt):
            pstring += '%10.3e'%(ion.Ioneq[stage, it])
        pstring += '  \n'
        out.write(pstring)
#
out.write('-1   \n')
out.write('%filename = chianti.ioneq \n')
out.write('%ionization equilibrium:  CHIANTI \n')
for one in chRef:
    out.write(one + '\n')
out.write('-1   \n')
out.close()
#
