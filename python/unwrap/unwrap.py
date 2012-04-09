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
def flood_quality(phase, quality):
	# Check sanity of phase & quality
	__check_sanity(phase)
	__check_sanity(quality)

	return unwrap_c.flood_quality(phase, quality)

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
		self.sz = (33, 57)
		self.verb = 2
		grid = N.indices(self.sz, dtype=N.float) / N.r_[self.sz].reshape(-1,1,1)
		# Random phase
		self.phase = N.sin(grid[0] * N.pi)*3. + N.cos(grid[1] * N.pi * 2)*4.
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
                        flood_quality, N.arange(100.0), N.arange(100.0))
		self.assertRaisesRegexp(ValueError, '.*should be on of.*',
                        flood_quality, N.arange(100).reshape(10,10), N.arange(100).reshape(10,10))

	def test1a_flood_dummy(self):
		"""Testing unwrapped unwrapping (dummy run)"""
		test_uw = flood_quality(self.phase, self.qualmap)
		self.assertAlmostEqual(test_uw.sum(), self.phase.sum())

	def test2a_flood_qual(self):
		"""Testing quality guided floodfill unwrapping"""
		test_uw = flood_quality(self.phase_wr, self.qualmap)
		if (self.verb > 1):
			plt.figure()
			plt.title("Original phase")
			plt.imshow(self.phase)
			plt.colorbar()

			plt.figure()
			plt.title("Quality map")
			plt.imshow(self.qualmap)
			plt.colorbar()
			print N.argmax(self.qualmap), N.max(self.qualmap), N.argwhere(self.qualmap == self.qualmap.max())

			plt.figure()
			plt.title("Wrapped phase")
			plt.imshow(self.phase_wr)
			plt.colorbar()

			plt.figure()
			plt.title("Recovered phase")
			plt.imshow(test_uw)
			plt.colorbar()

			raw_input()

		diff = test_uw - self.phase
		diff -= diff.mean()
		self.assertAlmostEqual(diff.sum(), 0)

class TestUnwrap(unittest.TestCase):
	def setUp(self):
		"""Generate phase and wrapped phase of shape (257, 509)"""
		sz = (257, 509)
		grid = N.indices(sz) *2.5 / N.r_[sz].reshape(-1,1,1)
		self.phase = ((grid*2.0)**2.0).sum(0)
		self.phase_wr = self.phase % (2*N.pi)
		# Flat image to test if unwrap() does not introduce errors
		self.flat = N.ones(sz)

	# Shallow data tests
	def test0a_wrapping(self):
		"""Phase amp should be > 2 pi, wrapped phase amp should be <= 2 pi"""
		self.assertTrue(self.phase.ptp() > 2*N.pi)
		self.assertTrue(self.phase_wr.ptp() <= 2*N.pi)

	def test0b_wrapping(self):
		"""Wrapped phase should be unequal to original phase"""
		self.assertFalse(N.allclose(self.phase_wr, self.phase))
		self.assertTrue(self.phase.ptp() > 2*N.pi)
		self.assertTrue(self.phase_wr.ptp() <= 2*N.pi)

	# Shallow function test
	def test1a_ret_shape_type(self):
		"""Return shape and type should be sane"""
		test_unwr = flood_quality(self.phase_wr, self.flat)
		self.assertEqual(test_unwr.shape, self.phase.shape)
		self.assertEqual(test_unwr.dtype, self.phase.dtype)

	# Deep function test
	def test2a_unwrap(self):
		"""Unwrapping flat data should do nothing"""
		test_unwr = flood_quality(self.flat, self.flat)
		self.assertTrue(N.allclose(test_unwr, self.flat))

	def test2b_unwrap(self):
		"""Unwrapping original phase should do nothing"""
		test_unwr = flood_quality(self.phase, self.flat)
		self.assertTrue(N.allclose(test_unwr, self.phase))

	def test2c_unwrab(self):
		"""Unwrapped phase should be same as original phase"""
		test_unwr = flood_quality(self.phase_wr, self.flat)
		self.assertTrue(N.allclose(test_unwr, self.phase))


# This must be the final part of the file, code after this won't be executed
if __name__ == "__main__":
	import pylab as plt
	unittest.main()
