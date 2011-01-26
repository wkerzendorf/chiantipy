'''the ChiantiPy - CHIANTI Python package
calculates various aspects of emission line and continua from the
CHIANTI atomic database for astrophysical spectroscopy'''
import os
import util
import constants
import filters
import mputil
#
try:
    chInteractive = int(os.environ['CHIANTIPY_INTERACTIVE'])
except:
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
Defaults = util.defaultsRead(verbose = chInteractive)
Ip = util.ipRead()
MasterList = util.masterListRead()
AbundanceAll = util.abundanceRead(abundancename = Defaults['abundfile'])
IoneqAll = util.ioneqRead(ioneqname = Defaults['ioneqfile'])
import version
__version__ = version.__version__
__version_info__ = version.__version_info__
#import core
