from distutils.core import setup
setup(name = 'ChiantiPy',
    version = '0.1.1',
    description = 'a Python interface to the CHIANTI atomic database for astrophysical spectroscopy',
    long_description = open('README').read(),
    author = 'Ken Dere',
    author_email = 'kdere@gmu.edu',
    url = 'http://chiantipy.sourceforge.net',
    download_url = 'http://sourceforge.net/projects/chiantipy',
    package_dir = {'chianti':''},
    packages = ['chianti','chianti.gui_qt','chianti.gui_wx','chianti.gui_cl'],
    classifiers = [
    'Development Status :: 3 - Alpha',
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'Intended Audience :: End Users/Desktop',
    'License :: OSI Approved :: GNU General Public License (GPL)',
    'Operating System :: POSIX :: Linux',
    'Programming Language :: Python',
    'Topic :: Scientific/Engineering :: Astronomy',
    'Topic :: Scientific/Engineering :: Physics'
    ]
    )
