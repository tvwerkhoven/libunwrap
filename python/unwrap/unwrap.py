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

## @brief Quality-guided floodfill phase unwrapping
#
def floodfill(phase, quality):
	# Check sanity of phase & quality
	__check_sanity(phase)
	__check_sanity(quality)

	return unwrap_c.floodfill(phase, quality)

## @brief Check if a map seems ok
#
def __check_sanity(data, ndim=2, dtlist=[N.float]):
	if (data.ndim != ndim):
		raise ValueError("Data should be %d-dimensional." % (ndim))
	if (data.dtype not in dtlist):
		raise ValueError("Data should be on of %s." % str(dtlist))

### Test routines ###########################################################

class TestSanityCheck(unittest.TestCase):
	def setUp(self):
		"""Generate fake phase"""
		self.sz = (257, 257)
		self.verb = 2
		grid = N.indices(self.sz, dtype=N.float) / N.r_[self.sz].reshape(-1,1,1)
		# Random phase
		self.phase = N.sin(grid[0] * N.pi * 2)*3. + N.cos(grid[1] * N.pi * 8)*4.
		# Wrapped phase
		self.phase_wr = (self.phase % (2*N.pi)) - N.pi
		# Quality map
		self.qualmap = 1.0*grid.sum(0)

	# Shallow function tests
	def test0a_inspect(self):
		"""Inspect data"""
# 		plt.figure()
# 		plt.imshow(self.phase)
# 		plt.figure()
# 		plt.imshow(self.phase_wr)
# 		raw_input()
		pass

	def test0b_sanity(self):
		"""Test __check_sanity()"""
		self.assertRaisesRegexp(ValueError, '.*should be 2-dimensional.',
                        floodfill, N.arange(100.0), N.arange(100.0))
		self.assertRaisesRegexp(ValueError, '.*should be on of.*',
                        floodfill, N.arange(100).reshape(10,10), N.arange(100).reshape(10,10))

	def test1a_flood_dummy(self):
		"""Testing unwrapped unwrapping (dummy run)"""
		test_uw = floodfill(self.phase, self.qualmap)
		self.assertAlmostEqual(test_uw.sum(), self.phase.sum())

	def test2a_flood_qual(self):
		"""Testing quality guided floodfill unwrapping"""
		test_uw = floodfill(self.phase_wr, self.qualmap)
		test_uw -= test_uw.mean()
		if (self.verb > 1):
			plt.figure()
			plt.title("Original phase")
			plt.imshow(self.phase)
			plt.colorbar()

			plt.figure()
			plt.title("Wrapped phase")
			plt.imshow(self.phase_wr)
			plt.colorbar()

			plt.figure()
			plt.title("Recovered phase")
			plt.imshow(test_uw)
			plt.colorbar()
			raw_input()

		self.assertAlmostEqual((test_uw - self.phase).sum(), 0)


# This must be the final part of the file, code after this won't be executed
if __name__ == "__main__":
	import pylab as plt
	unittest.main()
