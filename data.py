'''module for collecting various top-level Chianti data'''
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
if chInteractive:
    IoneqAll = util.ioneqRead(ioneqname = Defaults['ioneqfile'])
import version
__version__ = version.__version__
__version_info__ = version.__version_info__
