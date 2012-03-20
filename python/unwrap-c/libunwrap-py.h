/**
@file libshifts-c.h
@brief Header file for libshifts-c.c
@author Tim van Werkhoven (tim@astro.su.se)
@date 20090507

Created by Tim van Werkhoven on 2009-05-07.
Copyright (c) 2009 Tim van Werkhoven (tim@astro.su.se)

This file is licensed under the Creative Commons Attribution-Share Alike
license versions 3.0 or higher, see
http://creativecommons.org/licenses/by-sa/3.0/
*/

#ifndef __LIBSHIFTS_C_H__  // __LIBSHIFTS_C_H__
#define __LIBSHIFTS_C_H__

//
// Headers
//

#include <Python.h>				// For python extension
#include <numpy/arrayobject.h> 	// For numpy
//#include <Numeric/arrayobject.h> 	// For numpy (old)
//#include <numarray/arrayobject.h> 	// For numpy (old)
#include <stdint.h>
#include <sys/time.h>			// For timestamps
#include <time.h>					// For timestamps
#include <math.h>					// For pow()

//
// Defines
//

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

// Defines for comparison algorithms
#define COMPARE_XCORR 0				// Direct cross correlation
#define COMPARE_SQDIFF 1			// Square difference
#define COMPARE_ABSDIFFSQ 2		// Absolute difference squared
#define COMPARE_FFT 3					// Fourier method

// Defines for extremum finding algorithms
#define EXTREMUM_MAXVAL 0		 	// maximum value (no interpolation)
#define EXTREMUM_2D9PTSQ 1		// 2d 9 point parabola interpolation
#define EXTREMUM_2D5PTSQ 2		// 2d 5 point parabola interpolation

// Defines for reference usage
#define REF_BESTRMS 0         // Use subimages with best RMS as reference, 
								 							// 'refopt' should be an integer indicating how 
								 							// many references should be used.
#define REF_STATIC 1					// Use static reference subapertures, pass a 
															// list to the 'refopt' parameter to specify
															// which subaps should be used.

#define NTHREADS 4            // Number of threads to work with

//#define DEBUG               // Enable debug output

//
// Types
//

typedef float float32_t;  // 'Standard' 32 bit float type
typedef double float64_t; // 'Standard' 64 bit float type

// Thread data type
struct thread_data32 {
	float32_t *img;         // Big image and stride
	int32_t stride;
	int32_t *mask;          // Correlation mask and stride
	int32_t maskstride;
  float32_t *shifts;      // Store shifts here
	int32_t (*sapos)[2];	  // Subaperture positions
  int32_t nsa;
	int32_t (*sfpos)[2];	  // Subfield positions
  int32_t nsf;
	float32_t *ref;         // Reference subaperture (already normalized)
  int32_t refsa;
	int32_t *sasize;  	    // Subaperture size
	int32_t *sfsize;	      // Subfield size
	int32_t *shran;         // Shift range to test
	int32_t dosa[2];        // Subapertures this thread should process
  int compmeth;           // Comparison method to use
  int extmeth;            // Interpolation method to use
};

//
// Prototypes
//

// Python-accessible functions
static PyObject * libshifts_calcshifts(PyObject *self, PyObject *args);

// Helper routines
int _findrefidx_float32(float32_t *image, int32_t stride, int32_t sapos[][2], int npos, int32_t sasize[2], int refmode, int refopt, int32_t **list, int32_t *nref);

// Main 'glueing' routine
int _calcshifts_float32(float32_t *image, int32_t stride, int32_t *mask, int32_t maskstride, int32_t sapos[][2], int nsa, int32_t sasize[2], int32_t sfpos[][2], int nsf, int32_t sfsize[2], int32_t shran[2], int compmeth, int extmeth, int32_t *reflist, int32_t nref, float32_t **shifts);

void *_procsubaps_float32(void* args);

// Square difference image comparison
int _sqdiff(float32_t *img, int32_t imgsize[2], int32_t imstride, float32_t *ref, int32_t refsize[2], int32_t refstride, int32_t *mask, int32_t maskstride, float32_t *diffmap, int32_t pos[2], int32_t range[2], int bigref);

// Absolute difference squared image comparison
int _absdiffsq(float32_t *img, int32_t imgsize[2], int32_t imstride, float32_t *ref, int32_t refsize[2], int32_t refstride, int32_t *mask, int32_t maskstride, float32_t *diffmap, int32_t pos[2], int32_t range[2], int bigref);

// Direct cross-correlation image comparison
int _crosscorr(float32_t *img, int32_t imgsize[2], int32_t imstride, float32_t *ref, int32_t refsize[2], int32_t refstride, int32_t *mask, int32_t maskstride, float32_t *diffmap, int32_t pos[2], int32_t range[2], int bigref);

// 9-point quadratic interpolation
int _9pquadint(float32_t *diffmap, int32_t diffsize[2], float32_t shvec[2], int32_t shran[2]);

// 5-point quadratic interpolation
int _5pquadint(float32_t *diffmap, int32_t diffsize[2], float32_t shvec[2], int32_t shran[2]);

// One-liner help functions
int _comp_dbls(const double *a, const double *b) {
  if ((*b - *a) > 0) return 1;
  if ((*b - *a) < 0) return -1;
  return 0;
  //return (*b - *a) > 0 ? 1 : ( (*b - *a) < 0) ? -1 : 0 ;
}

#endif // __LIBSHIFTS_C_H__
