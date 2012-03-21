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

# Importe C routines
import unwrap_c

import numpy as N
import unittest

## @brief Unwrap phase, wrapper for C-function
#
# This is a Python wrapper for the Python-C library. It mostly checks and
# sanitizes the parameters before passing it on to the Python-C library.
#
# Available methods include:
# - Itoh
# - Floodfill
def unwrap(phase, method="itoh", *args):
	# Check method validity
	val_meths = ['itoh', 'flood']
	if (method.lower() not in val_meths):
		raise ValueError("Invalid method, must be one of " + str(val_meths))

	# Check sanity of phase
	__check_sanity(phase)

	# Check method and call function
	if (method == 'flood'):
		return floodfill(phase, *args)
	elif (method == 'itoh'):
		return phase
	else:
		return phase

## @brief Floodfill phase unwrapping
#
def floodfill(phase, startx=0, starty=0):
	# Check sanity of phase
	__check_sanity(phase)

	return unwrap_c.floodfill(phase, startx, starty)

## @brief Check if phase seems ok
#
def __check_sanity(phase, ndim=2, dtlist=[N.float]):
	if (phase.ndim != ndim):
		raise ValueError("Phase should be %d-dimensional." % (ndim))
	if (phase.dtype not in dtlist):
		raise ValueError("Phase should be on of %s." % str(dtlist))

### Test routines ###########################################################

class TestSanityCheck(unittest.TestCase):
	def setUp(self):
		"""Generate fake phase"""
		self.sz = (257, 509)
		self.val_meths = ['itoh', 'flood']
		grid = N.indices(self.sz) / N.r_[self.sz].reshape(-1,1,1)
		# Random phase
		self.phase = N.sin(grid[0] * N.pi * 2)*3. + N.cos(grid[1] * N.pi * 8)*4.
		# Wrapped phase
		self.phase_wr = (self.phase % (2*N.pi)) - N.pi

	# Shallow function tests
	def test0b_sanity(self):
		"""Test __check_sanity()"""
		self.assertRaisesRegexp(ValueError, 'Phase should be 2-dimensional.',
                        unwrap, N.arange(100.0))
		self.assertRaisesRegexp(ValueError, 'Phase should be on of.*',
                        unwrap, N.arange(100).reshape(10,10))

	def test0c_unwrap(self):
		"""Test unwrap methods check"""
		self.assertRaisesRegexp(ValueError, 'Invalid method.*',
                        unwrap, N.arange(100.0), 'invalid')

	def test1a(self):
		"""Testing unwrapped unwrapping (dummy run)"""
		for m in self.val_meths:
			test_uw = unwrap(self.phase, method=m)
			self.assertTrue(N.allclose(test_uw, self.phase))


# This must be the final part of the file, code after this won't be executed
if __name__ == "__main__":
	unittest.main()
