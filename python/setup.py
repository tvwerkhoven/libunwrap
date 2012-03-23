#!/usr/bin/env python
# encoding: utf-8
"""
setup.py -- setup file for the Python part of libunwrap

Created by Tim van Werkhoven (werkhoven@strw.leidenuniv.nl) on 2012-03-20
Copyright (c) 2012 Tim van Werkhoven. All rights reserved.
"""
import sys

# Try importing to see if we have NumPy available (we need this)
try:
	import numpy
	from numpy.distutils.core import setup, Extension
	from numpy.distutils.misc_util import Configuration
except:
	print "Could not load NumPy (numpy.distutils.{core,misc_util}), required by this package. Aborting."
	sys.exit(1)

# Setup extension module for C stuff
extlibs = []
extlibs.append(Extension('unwrap_c',
			define_macros =
				[('MAJOR_VERSION', '0'),
				('MINOR_VERSION', '1')],
			include_dirs = [numpy.get_include(), "../c/unwrap"],
			libraries = ["m"],
			library_dirs = [],
			extra_compile_args=["-O3", "-ffast-math", "-Wall", "-DHAVE_DEBUGPRINT"],
			extra_link_args=None,
			sources = ['unwrap-c/libunwrap-py.cc', 'unwrap-c/libunwrap-py.h',
			'../c/unwrap/libunwrap.cc', '../c/unwrap/libunwrap.h']))

# Setup
setup(name = 'unwrap',
	version = '0.1',
	description = 'Phase unwrapping with various algorithms',
	author = 'Tim van Werkhoven',
	author_email = 'werkhoven@strw.leidenuniv.nl',
	url = '',
	license = "GPL",
	packages = ['unwrap'],
	# This is for the C module
	ext_package = 'unwrap',
	ext_modules = extlibs)
