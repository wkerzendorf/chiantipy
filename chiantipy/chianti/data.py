'''
module for collecting various top-level Chianti data
for the keyword arguments below
temperature = temperature in Kelvin
eDensity is the electron density per cubic cm
hDensity is the hydrogen density per cubic cm
pDensity is the proton density per cubic cm
radTemperature is the radiation temperature of central source
rStar is the distance of the plasma from the source in units of the sources radius
distance is the distance from the central source
'''
import os
import util
#import constants
#import filters
#import mputil
#
if os.environ.has_key('CHIANTIPY_INTERACTIVE'):
    chInteractive = int(os.environ['CHIANTIPY_INTERACTIVE'])
else:
    chInteractive = 1

#if chInteractive:
#    import pylab as pl
#else:
#    import matplotlib
#    matplotlib.use('Agg')
#    import matplotlib.pyplot as pl
###
xuvtop = os.environ['XUVTOP']
#chInteractive=1
Defaults = util.defaultsRead()
Ip = util.ipRead()
MasterList = util.masterListRead()
AbundanceAll = util.abundanceRead(abundancename = Defaults['abundfile'])
IoneqAll = util.ioneqRead(ioneqname = Defaults['ioneqfile'])
import version
__version__ = version.__version__
__version_info__ = version.__version_info__
keywordArgs = ['temperature','eDensity','hDensity', 'pDensity','radTemperature', 'rStar', 'distance']
