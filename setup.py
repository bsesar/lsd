#!/usr/bin/env python

import os, os.path, sys

try:
        import numpy
except ImportError:
        raise Exception("LSD requires numpy")

numpy_include = numpy.get_include()

# ------ no changes below! If you need to change, it's a bug! -------
from distutils.core import setup, Extension
from sys import platform

import numpy
inc = [numpy_include]

longdesc = """Large Survey Database"""

args = { 
	'name'			: "lsd",
	'version'		: "0.3.1", # When you change this, modify __version__ in __init__.py
	'description'	 	: "Large Survey Database Python Module",
	'long_description'	: longdesc,
	'license'		: "GPLv2",
	'author'		: "Mario Juric",
	'author_email'		: "mjuric@cfa.harvard.edu",
	'maintainer'		: "Mario Juric",
	'maintainer_email'	: "mjuric@cfa.harvard.edu",
	'url'			: "http://mwscience.net/lsd",
	'download_url'		: "http://mwscience.net/lsd/download",
	'classifiers'		: [
					'Development Status :: 3 - Alpha',
					'Intended Audience :: Science/Research', 
					'Intended Audience :: Developers',
					'License :: OSI Approved :: GNU General Public License (GPL)', 
					'Programming Language :: C++', 
					'Programming Language :: Python :: 2', 
					'Programming Language :: Python :: 2.7',
					'Operating System :: POSIX :: Linux',
					'Topic :: Database',
					'Topic :: Scientific/Engineering :: Astronomy'
	],
	'scripts'	: ['src/lsd-footprint', 'src/lsd-import-sdss', 'src/lsd-make-object-catalog',
	 			'src/lsd-import-dvo', 'src/lsd-import-smf',
	 			'src/lsd-query', 'src/lsd-xmatch', 'src/mr-peer'],
	'packages'	: ['lsd', 'mr'],
	'package_dir'	: {'': 'src'},
	'ext_modules'	: [Extension('lsd.native', ['src/native/main.cpp'], include_dirs=inc)],
	'data_files'    : [('share/lsd/examples', ['src/examples/latitude_histogram.py', 'src/examples/count_rows.py'])]
}

setup(**args)
