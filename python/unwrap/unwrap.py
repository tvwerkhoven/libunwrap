#!/usr/bin/env python
# encoding: utf-8
"""
This is unwrap.py, a Python wrapper for phase unwrapping routines.
"""

## @file unwrap.py
# @brief Python wrapper for unwrap.c library
# @author Tim van Werkhoven (werkhoven@strw.leidenuniv.nl)
# @date 20120319
#
# Created by Tim van Werkhoven on 2012-03-19.
# Copyright (c) 2012 Tim van Werkhoven (werkhoven@strw.leidenuniv.nl)
#
# This file is licensed under the Creative Commons Attribution-Share Alike
# license versions 3.0 or higher, see
# http://creativecommons.org/licenses/by-sa/3.0/


## @brief Pure Python wrapper for C-Python library
#
# This is a Python wrapper for the Python-C library. It mostly checks and
# sanitizes the parameters before passing it on to the Python-C library.
#
# Available methods include:
# - Itoh
# - Floodfill
def unwrap(phase, method="itoh")
	print "unwrap(...)"

	# Check method validity
	val_meths = ['itoh', 'flood']
	if (method.lower() not in val_meths):
		raise ValueError("Invalid method, must be one of " + str(val_meths)).

	# Check phase dimension

	# Check phase datatype

	# Wrap to C-library

	# Return results

	return
