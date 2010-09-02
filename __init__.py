'''the CHIANTI Python package'''
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
##
xuvtop = os.environ['XUVTOP']
chInteractive=1
Defaults = util.defaultsRead(verbose = chInteractive)
Ip = util.ipRead()
MasterList = util.masterListRead()
AbundanceAll = util.abundanceRead(abundancename = Defaults['abundfile'])
IoneqAll = util.ioneqRead(ioneqname = Defaults['ioneqfile'])
import core
__version__ = '0.3'
