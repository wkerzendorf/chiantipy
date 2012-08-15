'''
utilities for testing, improving CHIANTI atomic data files
'''
import chianti.util as util
#
def uniqueConfigurations(elvlc):
    '''
    to get a list of the unique configurations
    '''
    unique = []
    for aterm in elvlc['term']:
        if aterm not in unique:
             unique.append(aterm)
    #
    return unique
    #
    #  --------------------------------------------------------
    #
def energiesCompare(elvlc):
    '''
    to compare observed versus theoretical for each configuration
    '''
    return
