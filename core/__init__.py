''' chianti.core - contains the main classes for ChiantiPy users.

Copyright 2009, 2010 Kenneth P. Dere

This software is distributed under the terms of the GNU General Public License
that is found in the LICENSE file


'''
#import chianti
#import os
#import types
#from datetime import datetime
##from ConfigParser import *
#import numpy as np
#from scipy import interpolate
##
## the following is necessary to make chiantipy non interactive for the web
#try:
#    chInteractive = int(os.environ['CHIANTIPY_INTERACTIVE'])
#except:
#    chInteractive = 1
##
#try:
##    from Queue.Queue import *
#    import multiprocessing as mp
#    from chianti import mputil
#    import chianti.mputil as mputil
##    from multiprocessing import Pool as mp.Pool
##    from multiprocessing import Process,  Queue
#except:
#    if chInteractive:
#        print ' your version of Python does not support multiprocessing \n you will not be able to use mspectrum'
##
#if chInteractive:
#    import pylab as pl
#else:
#    import matplotlib
#    matplotlib.use('Agg')
#    import matplotlib.pyplot as pl
#try:
#    from matplotlib.delaunay.triangulate import Triangulation
#except:
#    from scikits.delaunay.triangulate import Triangulation
##import chianti
##from chianti import util
##import chianti.util as util
##import chianti.constants as const
##import chianti.filters as chfilters
##
#xuvtop=os.environ['XUVTOP']
##
#import chianti.util as util
#ip = util.ipRead()
#MasterList = util.masterListRead()
#defaults = util.defaultsRead(verbose = chInteractive)
#print ' core __init_ defaults = ', defaults
#
#from chianti import *
from Spectrum import spectrum
from Mspectrum import mspectrum
from Continuum import continuum
from RadLoss import radLoss
#from ion import ioneq
from Ion import ion
from Ion import ionWeb
from Ion import ioneq
from Ion import photoioneq
##from ionweb import ionWeb
#
#if chInteractive:
#    import matplotlib as pl
#else:
#    import pylab as pl
##
###
##
##if defaults['gui'] == 'qt':
##   import chianti.gui_qt.gui as gui
##elif defaults['gui'] == 'wx':
##    import chianti.gui_wx.gui as gui
##else:
##    print ' unknown gui - ',defaults['gui']
##    print ' the full functionality of the chiant.core.ion class will not be available'
#    #
#if pl.rcParams['backend'].lower() == 'qt4agg':
#    import chianti.gui_qt.gui as gui
#elif pl.rcParams['backend'].lower() == 'wxagg':
#    import chianti.gui_wx.gui as gui
#elif pl.rcParams['backend'].lower() == 'gtkagg':
#    import chianti.gui_cl.gui as gui
#elif pl.rcParams['backend'].lower() == 'agg':
#    import chianti.gui_cl.gui as gui
#elif pl.rcParams['backend'].lower() == 'agg':
#    import chianti.gui_cl.gui as gui
#elif pl.rcParams['backend'].lower() == 'macosx ':
#    import chianti.gui_cl.gui as gui
#else:
#    print ' - Warning - '
#    print ' - in order to use the various gui dialogs, the matlpotlib/pylab backend needs'
#    print ' - to be either Qt4Agg or WXAgg - '
#    print ' - in order to use the command line dialogs, the matlpotlib/pylab backend needs'
#    print ' - to be GTKAgg or MacOSX - '
#    print ' - current backend is ',pl.rcParams['backend']
#    print ' - the full functionality of the chianti.core.ion class may not be available'
#    print ' - it would probably be better to set your matplotlib backend to either'
#    print ' - Qt4Agg, WXAgg, GTKAgg, or MacOSX'
#    print ' - using the command line dialogs for now but there could be problems -'
#    import chianti.gui_cl.gui as gui
#    #
