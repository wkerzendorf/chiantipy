'''the ChiantiPy - CHIANTI Python package
calculates various aspects of emission line and continua from the
CHIANTI atomic database for astrophysical spectroscopy'''
import os
import constants
import filters
import mputil
#
#try:
#    chInteractive = int(os.environ['CHIANTIPY_INTERACTIVE'])
#except:
#    chInteractive = 1

#if chInteractive:
#    import pylab as pl
#else:
#    import matplotlib
#    matplotlib.use('Agg')
#    import matplotlib.pyplot as pl
###
#xuvtop = os.environ['XUVTOP']
##chInteractive=1
#Defaults = util.defaultsRead(verbose = chInteractive)
#Ip = util.ipRead()
#MasterList = util.masterListRead()
#AbundanceAll = util.abundanceRead(abundancename = Defaults['abundfile'])
#IoneqAll = util.ioneqRead(ioneqname = Defaults['ioneqfile'])
#import version
#__version__ = version.__version__
#__version_info__ = version.__version_info__
#import core
import pylab as pl
if pl.rcParams['backend'].lower() == 'qt4agg':
    import gui_qt.gui as gui
elif pl.rcParams['backend'].lower() == 'wxagg':
    import gui_wx.gui as gui
elif pl.rcParams['backend'].lower() == 'gtkagg':
    import gui_cl.gui as gui
elif pl.rcParams['backend'].lower() == 'agg':
    import gui_cl.gui as gui
elif pl.rcParams['backend'].lower() == 'agg':
    import gui_cl.gui as gui
elif pl.rcParams['backend'].lower() == 'macosx':
    import gui_cl.gui as gui
else:
    print ' - Warning - '
    print ' - in order to use the various gui dialogs, the matlpotlib/pylab backend needs'
    print ' - to be either Qt4Agg or WXAgg - '
    print ' - in order to use the command line dialogs, the matlpotlib/pylab backend needs'
    print ' - to be GTKAgg or MacOSX - '
    print ' - current backend is ',pl.rcParams['backend']
    print ' - the full functionality of the chianti.core.ion class may not be available'
    print ' - it would probably be better to set your matplotlib backend to either'
    print ' - Qt4Agg, WXAgg, GTKAgg, or MacOSX'
    print ' - using the command line dialogs for now but there could be problems -'
    import gui_cl.gui as gui
#
# placed here because util needs gui
import util

