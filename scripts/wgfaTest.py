'''
a python script to time the two different wgfaRead functions
'''
import chianti.util as util
import futil
from datetime import datetime


filename = util.zion2filename(26, 13)+'.wgfa'

print ' filename = ', filename

t1 = datetime.now()

wgfax = futil.wgfaReadx(filename)

t2 = datetime.now()
dt1=t2-t1

wgfa = util.wgfaRead('', filename)

t3 = datetime.now()
dt2 = t3 - t2
print ' elapsed seconds = ', dt1.microseconds,  dt2.microseconds
