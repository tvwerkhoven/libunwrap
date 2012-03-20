/*
 libunwrap-py.c -- Python bindings for libunwrap C library.
 Copyright (C) 2012 Tim van Werkhoven <werkhoven@strw.leidenuniv.nl>
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*!
 @file libunwrap-py.c
 @brief Python bindings for libunwrap C library.
 @author Tim van Werkhoven (werkhoven@strw.leidenuniv.nl)
 @date 20120320 
*/

//
// Headers
//

#include <Python.h>				// For python extension
#include <numpy/arrayobject.h> 	// For numpy
//#include <Numeric/arrayobject.h> 	// For numpy
//#include <numarray/arrayobject.h> 	// For numpy
#include <sys/time.h>			// For timestamps
#include <time.h>					// For timestamps
#include <math.h>					// For pow()
#include "libunwrap-py.h" 	// For this file

//
// Methods table for this module
//

static PyMethodDef LibunwrapMethods[] = {
	{"floodfill",  libunwrap_floodfill, METH_VARARGS, "Flood-fill unwrap."},
	{NULL, NULL, 0, NULL}        /* Sentinel */
};

//
// Init module methods
//
PyMODINIT_FUNC init_libshifts(void) {
	(void) Py_InitModule("_libunwrap", LibunwrapMethods);
	// Init numpy usage
	import_array();
}

//
// Main python function
//
static PyObject *libshifts_calcshifts(PyObject *self, PyObject *args) {
	// Python function arguments
	PyArrayObject *image, *saccdpos, *saccdsize, *sfccdpos, *sfccdsize, *shrange, *corrmask;
	// Optional arguments with default settings
	int compmeth = COMPARE_SQDIFF, extmeth = EXTREMUM_2D9PTSQ;
	int refmode = REF_BESTRMS, refopt=1;
	// Generic variables
	int i, j, ret;
	// Return variables
	float32_t *shifts;
	int32_t *reflist, nref;
	PyArrayObject *retshifts, *retreflist;
	
	//
	// Parse arguments from Python function
	//
	if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!O!|iiii", 
		&PyArray_Type, &image, 								// Image data
		&PyArray_Type, &saccdpos, 						// Subaperture positions
		&PyArray_Type, &saccdsize, 						// SA size
		&PyArray_Type, &sfccdpos, 						// Subfield positions
		&PyArray_Type, &sfccdsize, 						// SF size
		&PyArray_Type, &shrange, 							// Shift range
		&PyArray_Type, &corrmask, 						// Correlation mask
		&compmeth,														// Compmethod to use
		&extmeth, 														// Interpolation method to use
		&refmode,															// Subap reference mode
		&refopt)) {
		PyErr_SetString(PyExc_SyntaxError, "In calcshifts: failed to parse arguments.");
		return NULL; 
		}
	
#ifdef DEBUG
	printf("libshifts_calcshifts(): img: 0x%p, sapos: 0x%p, sasize: 0x%p, sfpos: 0x%p, sfsize: 0x%p, shran: 0x%p.", image, saccdpos, saccdsize, sfccdpos, sfccdsize, shrange);
#endif
	
	//
	// Verify arguments (most should be int32)
	//
	
	if (PyArray_TYPE((PyObject *) saccdpos) != NPY_INT32 ||
	 	PyArray_TYPE((PyObject *) saccdsize) != NPY_INT32 ||
	 	PyArray_TYPE((PyObject *) sfccdpos) != NPY_INT32 || 
		PyArray_TYPE((PyObject *) sfccdsize) != NPY_INT32 ||
		PyArray_TYPE((PyObject *) corrmask) != NPY_INT32 ||
		PyArray_TYPE((PyObject *) shrange) != NPY_INT32) {
#ifdef DEBUG
		printf("...: coordinates, ranges or masks are not int32.\n");
#endif
		PyErr_SetString(PyExc_ValueError, "In calcshifts: coordinates, ranges and masks should be int32.");
		return NULL;
		}
	
	// Refopt should be somewhere between 0 and 100
	if (refopt < 0 || refopt > 100) {
#ifdef DEBUG
		printf("...: refopt value invalid.\n");
#endif
		PyErr_SetString(PyExc_ValueError, "In calcshifts: refopt value invalid.");
		return NULL;	
	}
	
	//
	// Parse arguments
	//
	
	// Get number of subapertures and subfields
	int nsa = (int) PyArray_DIM((PyObject*) saccdpos, 0);
	int nsf = (int) PyArray_DIM((PyObject*) sfccdpos, 0);
#ifdef DEBUG
		printf("...: nsa: %d, nsf: %d.\n", nsa, nsf);
#endif	
	// Convert options
	int32_t sapos[nsa][2];
	for (i=0; i<nsa; i++) {
		sapos[i][0] = *((uint32_t *)PyArray_GETPTR2((PyObject *) saccdpos, i, 0));
		sapos[i][1] = *((uint32_t *)PyArray_GETPTR2((PyObject *) saccdpos, i, 1));
#ifdef DEBUG
		printf("...: sa %d: %d,%d.\n", i, sapos[i][0], sapos[i][1]);
#endif
	}
	int32_t sfpos[nsf][2];
	for (i=0; i<nsf; i++) {
		sfpos[i][0] = *((uint32_t *)PyArray_GETPTR2((PyObject *) sfccdpos, i, 0));
		sfpos[i][1] = *((uint32_t *)PyArray_GETPTR2((PyObject *) sfccdpos, i, 1));
#ifdef DEBUG
		printf("...: sf %d: %d,%d.\n", i, sfpos[i][0], sfpos[i][1]);
#endif
	}
	int32_t sasize[2];
	sasize[0] = *((uint32_t *) PyArray_GETPTR1((PyObject *) saccdsize, 0));
	sasize[1] = *((uint32_t *) PyArray_GETPTR1((PyObject *) saccdsize, 1));
#ifdef DEBUG
	printf("...: sasize: %d,%d.\n", sasize[0], sasize[1]);
#endif
	int32_t sfsize[2];
	sfsize[0] = *((uint32_t *) PyArray_GETPTR1((PyObject *) sfccdsize, 0));
	sfsize[1] = *((uint32_t *) PyArray_GETPTR1((PyObject *) sfccdsize, 1));
#ifdef DEBUG
	printf("...: sfsize: %d,%d.\n", sfsize[0], sfsize[1]);
#endif
	
	int32_t shran[2];
	shran[0] = *((uint32_t *) PyArray_GETPTR1((PyObject *) shrange, 0));
	shran[1] = *((uint32_t *) PyArray_GETPTR1((PyObject *) shrange, 1));
#ifdef DEBUG
	printf("...: shran: %d,%d.\n", shran[0], shran[1]);
#endif
	
	//
	// Process data using various subroutines, depending on datatype
	//
	
	switch (PyArray_TYPE((PyObject *) image)) {
		// Image is of float32 type, use appropriate routines
		case (NPY_FLOAT32): {
#ifdef DEBUG
			printf("...: Found type NPY_FLOAT32\n");
#endif
			// Convenience pointers
			float32_t *im32 = (float32_t *) PyArray_DATA((PyObject *) image);
			int32_t imstride = PyArray_STRIDES((PyObject *) image)[0] / PyArray_ITEMSIZE((PyObject *) image);
			int32_t *mask = (int32_t *) PyArray_DATA((PyObject *) corrmask);
			int32_t maskstride = PyArray_STRIDES((PyObject *) corrmask)[0] / PyArray_ITEMSIZE((PyObject *) corrmask);
			// First get the reference subapertures
#ifdef DEBUG
			printf("...: Calling _findrefidx_float32().\n");
#endif
			ret = _findrefidx_float32(im32, imstride, sapos, nsa, sasize, refmode, refopt, &reflist, &nref);

#ifdef DEBUG
			printf("...: Calling _calcshifts_float32().\n");
#endif
			ret = _calcshifts_float32(im32, imstride, mask, maskstride, sapos, nsa, sasize, sfpos, nsf, sfsize, shran, compmeth, extmeth, reflist, nref, &shifts);
			break;
		}
		case (NPY_FLOAT64): {
			PyErr_SetString(PyExc_NotImplementedError, "In calcshifts: NPY_FLOAT64 not supported (yet).");
			return NULL;
			break;
		}
		default: {
#ifdef DEBUG
			printf("...: unsupported type.\n");
#endif
			PyErr_SetString(PyExc_NotImplementedError, "In calcshifts: datatype not supported.");
			return NULL;
		}
	}
	
	//
	// Reformat results to NumPy arrays
	//
		
	// Build numpy vector from 'reflist', a 1-d array
	npy_intp rr_dims[] = {nref};
	retreflist = (PyArrayObject*) PyArray_SimpleNewFromData(1, rr_dims,
		NPY_INT32, (void *) reflist);
	PyArray_FLAGS(retreflist) |= NPY_OWNDATA;
	
	// Build numpy vector from 'shifts', a 4-d array
	npy_intp s_dims[] = {nref, nsa, nsf, 2};
	retshifts = (PyArrayObject*) PyArray_SimpleNewFromData(4, s_dims,
		NPY_FLOAT32, (void *) shifts);
	PyArray_FLAGS(retshifts) |= NPY_OWNDATA;
	
	if (!PyArray_CHKFLAGS(retreflist, NPY_OWNDATA) ||
	 	!PyArray_CHKFLAGS(retshifts, NPY_OWNDATA)) {
		PyErr_SetString(PyExc_RuntimeError, "In calcshifts: unable to own 'reflist' or 'retshifts' data, aborting");
		free(reflist);
		free(shifts);
		return NULL;
		}
	
	return Py_BuildValue("{s:N,s:N}", "shifts", retshifts, "refapts", retreflist);
}

//
// Wrapper functions
//

int _calcshifts_float32(float32_t *image, int32_t stride, int32_t *mask, int32_t maskstride, int32_t sapos[][2], int nsa, int32_t sasize[2], int32_t sfpos[][2], int nsf, int32_t sfsize[2], int32_t shran[2], int compmeth, int extmeth, int32_t *reflist, int32_t nref, float32_t **shifts) {
#ifdef DEBUG
	printf("_calcshifts_float32(): im: 0x%p\n", image);
	printf("...: stride: %d, nsa: %d, nsf: %d, shran: %d,%d, nref: %d\n", stride, nsa, nsf, shran[0], shran[1], nref);
#endif
	int refsa, i, j, ret=0;
	float32_t pix, mean;
	int32_t *refpos;
	float32_t ref[sasize[0]*sasize[1]];				// Reference subap
#ifdef DEBUG
	printf("...: allocating memory (%dx%dx%dx%d).\n", nref, nsa, nsf, 2);
#endif
 	// allocate memory for shifts
	(*shifts) = (float32_t *) malloc(nref * nsa * nsf * 2 * sizeof(float32_t));

	//
	// Init different threads
	//
	pthread_t threads[NTHREADS];
	struct thread_data32 thr_dat[NTHREADS];
	int thr, thrret=0;
	void *thrstat;
#ifdef DEBUG
		printf ("Initializing thread data\n");
#endif
	for (thr=0; thr<NTHREADS; thr++) {
		thr_dat[thr].img = image;
		thr_dat[thr].stride = stride;
		thr_dat[thr].mask = mask;
		thr_dat[thr].maskstride = maskstride;
		thr_dat[thr].shifts = (*shifts);
		thr_dat[thr].sapos = sapos;
		thr_dat[thr].nsa = nsa;
		thr_dat[thr].sfpos = sfpos;
		thr_dat[thr].nsf = nsf;
		thr_dat[thr].ref = ref;
		thr_dat[thr].sasize = sasize;
		thr_dat[thr].sfsize = sfsize;
		thr_dat[thr].shran = shran;
		thr_dat[thr].dosa[0] = thr*nsa/NTHREADS;
		thr_dat[thr].dosa[1] = (thr+1)*nsa/NTHREADS;
		thr_dat[thr].compmeth = compmeth;
		thr_dat[thr].extmeth = extmeth;
	}
	
	//
	// Loop over reference subapertures
	//
	for (refsa=0; refsa<nref; refsa++) {
		refpos = sapos[reflist[refsa]];
#ifdef DEBUG
		printf("...: parsing ref %d at (%d,%d).\n", reflist[refsa], refpos[0], refpos[1]);
#endif
		// Cut out reference subap, calculate mean
		mean = 0;
		for (j=0; j<sasize[1]; j++) {
			for (i=0; i<sasize[0]; i++) {
				pix = image[(refpos[1] + j) * stride + refpos[0] + i];
				ref[sasize[0] * j + i] = pix;
				mean += pix;
			}
		}
		mean = mean/(sasize[1]*sasize[0]);
#ifdef DEBUG
		printf("...: ref mean was: %g.\n", mean);
#endif
		// Divide reference subap by mean
		for (j=0; j<sasize[1]; j++)
			for (i=0; i<sasize[0]; i++)
				ref[sasize[0] * j + i] /= mean;

#ifdef DEBUG
		printf("...: starting different threads.\n");
#endif
		for (thr=0; thr<NTHREADS; thr++) {
			thr_dat[thr].refsa = refsa;
			pthread_create(&threads[thr], 0, _procsubaps_float32, (void *)(thr_dat+thr));
		}
		
#ifdef DEBUG
		printf("...: joining threads.\n");
#endif
		for (thr=0; thr<NTHREADS; thr++)
			thrret += pthread_join(threads[thr], &thrstat);

	}
	
	return ret;
}

//
// Worker thread function
//

void *_procsubaps_float32(void* args) {
	// Re-cast argument to right type
	struct thread_data32 *dat = (struct thread_data32 *) args;
	int sa, sf, i, j, ret=0;
	int nsa = dat->nsa, nsf = dat->nsf;
	float32_t pix, mean;
	int _refoff = dat->refsa * (nsa*nsf*2);
	float32_t _subimg[dat->sasize[0]*dat->sasize[1]];		// Subap to test
	float32_t *_subfield;											// Pointer to subfield
	int32_t diffsize[] = {(dat->shran[0]*2+1), (dat->shran[1]*2+1)};
#ifdef DEBUG
	printf("...: diffmap: %dx%d.\n", diffsize[0], diffsize[1]);
#endif
	float32_t *diffmap = (float32_t *) calloc(diffsize[0] * diffsize[1], \
		sizeof(float32_t));
	float32_t shvec[2];

#ifdef DEBUG
	printf ("Thread working on %d--%d\n", dat->dosa[0], dat->dosa[1]);
#endif
	//
	// Loop over subapertures
	//
	
	for (sa=dat->dosa[0]; sa<dat->dosa[1]; sa++) {
#ifdef DEBUG
		printf("...: parsing sa %d at (%d,%d).. ", sa, dat->sapos[sa][0], dat->sapos[sa][1]);
#endif
		// Cut out subaperture, calculate mean
		mean = 0;
		for (j=0; j<dat->sasize[1]; j++) {
			for (i=0; i<dat->sasize[0]; i++) {
				pix = dat->img[(dat->sapos[sa][1] + j) * dat->stride + dat->sapos[sa][0] + i];
				_subimg[dat->sasize[0] * j + i] = pix;
				mean += pix;
			}
		}
		mean = mean/(dat->sasize[1]*dat->sasize[0]);
#ifdef DEBUG
		printf("mean: %g... ", mean);
#endif
		// Divide subap by mean
		for (j=0; j<dat->sasize[1]; j++)
			for (i=0; i<dat->sasize[0]; i++)
				_subimg[dat->sasize[0] * j + i] /= mean;

		//
		// Loop over subfields
		//
		for (sf=0; sf<dat->nsf; sf++) {
#ifdef DEBUG
			printf("sf %d... ", sf);
#endif
			// Cut out subfield
			_subfield = _subimg + (dat->sfpos[sf][1] * dat->sasize[0]) + dat->sfpos[sf][0];
			mean = 0.0;
			for (j=0; j<dat->sfsize[1]; j++)
				for (i=0; i<dat->sfsize[0]; i++)
					mean += _subfield[dat->sasize[0] * j + i];
#ifdef DEBUG					
			printf("sf mean: %g... ", mean/(dat->sfsize[1] * dat->sfsize[0]));
#endif
	
			// Calculate correlation map
			switch (dat->compmeth) {
				case (COMPARE_ABSDIFFSQ): {
					ret = _absdiffsq(_subfield, dat->sfsize, dat->sasize[0], dat->ref, dat->sasize, dat->sasize[0], dat->mask, dat->maskstride, diffmap, dat->sfpos[sf], dat->shran, 1);
					break;
				}
				case (COMPARE_SQDIFF): {
					ret = _sqdiff(_subfield, dat->sfsize, dat->sasize[0], dat->ref, dat->sasize, dat->sasize[0], dat->mask, dat->maskstride, diffmap, dat->sfpos[sf], dat->shran, 1);
					break;
				}
				case (COMPARE_XCORR): {
					ret = _crosscorr(_subfield, dat->sfsize, dat->sasize[0], dat->ref, dat->sasize, dat->sasize[0], dat->mask, dat->maskstride, diffmap, dat->sfpos[sf], dat->shran, 1);
					break;
				}
				default: {
					return NULL;
					break;
				}
			}
			if (ret == 2)
				printf("_procsubaps_float32(): err in COMPARE @ sa: %d, sf :%d\n", \
					sa, sf);
			
			// Find subpixel maximum
			switch (dat->extmeth) {
				case (EXTREMUM_2D9PTSQ): {
					ret = _9pquadint(diffmap, diffsize, shvec, dat->shran);
					break;
				}
				case (EXTREMUM_2D5PTSQ): {
					ret = _5pquadint(diffmap, diffsize, shvec, dat->shran);					
					break;
				}
				default: {
					return NULL;
					break;
				}
			}
			if (ret == 2)
				printf("_procsubaps_float32(): err in EXTREMUM @ sa: %d, sf :%d\n", \
				 	sa, sf);

			dat->shifts[_refoff + sa * (nsf * 2) + sf * (2) + 0] = (float32_t) shvec[0]-dat->shran[0];
			dat->shifts[_refoff + sa * (nsf * 2) + sf * (2) + 1] = (float32_t) shvec[1]-dat->shran[1];
#ifdef DEBUG
			printf("sh: (%.3g, %.3g) ", shvec[0]-dat->shran[0], shvec[1]-dat->shran[1]);
#endif
		}
#ifdef DEBUG
		printf("\n");
#endif
	}
	return NULL;
}

//
// Find reference subapertures
//

int _findrefidx_float32(float32_t *image, int32_t stride, int32_t sapos[][2], int npos, int32_t sasize[2], int refmode, int refopt, int32_t **list, int32_t *nref) {
#ifdef DEBUG
	printf("_findrefidx_float32() im: 0x%p, pos: 0x%p, size: 0x%p.\n", image, sapos, sasize);
	printf("_findrefidx_float32() refmode: %d, refopt: %d.\n", refmode, refopt);
#endif
	// Find a reference subaperture in 'image'
	int sa, ssa, i, j;
	
	if (refmode == REF_STATIC) {
		// Only return a static reference subaperture
		*list = (int32_t*) malloc(1 * sizeof(int32_t));
		(*list)[0] = refopt;
		*nref = 1;
		return 0;
	} // refmode == REF_STATIC
	else if (refmode == REF_BESTRMS){
		if (refopt < 1) refopt = 1;
		*list = malloc(refopt * sizeof(int32_t));
		*nref = refopt;
		double rmslist[npos], rmslists[npos];

		// Calculate RMS values
		for (sa=0; sa < npos; sa++) {
#ifdef DEBUG
			printf("...: checking subaperture %d...", sa);
#endif
			rmslist[sa] = 0;
			double mean=0;
			// Calculate mean
			// ONLY USE THE CENTRAL 50% OF THE SUBAP!
			for (j=sasize[1]/4; j<3*sasize[1]/4; j++)
				for (i=sasize[0]/4; i<3*sasize[0]/4; i++)
					mean += image[(sapos[sa][1] + j) * stride + sapos[sa][0] + i];
			mean /= ((sasize[0]/2) * (sasize[1]/2));
			for (j=sasize[1]/4; j<3*sasize[1]/4; j++)
				for (i=sasize[0]/4; i<3*sasize[0]/4; i++)
					rmslist[sa] += pow(image[(sapos[sa][1] + j) * stride + sapos[sa][0] + i] - mean, 2.0);

			rmslist[sa] = 100.0*pow(rmslist[sa]/((sasize[0]/2) * (sasize[1]/2)), \
			 	0.5)/mean;
			rmslists[sa] = rmslist[sa];
#ifdef DEBUG
			printf(" rms is: %g\n", rmslist[sa]);
#endif
		}

		// Get first refopt best values
#ifdef DEBUG
		printf("...: sorting RMS values\n");
#endif
		qsort((void*) rmslists, npos, sizeof(double), _comp_dbls);

		for (sa=0; sa < npos; sa++) {
#ifdef DEBUG
			printf("...: sa %d rms: %g.\n", sa, rmslists[sa]);
#endif
		}

		for (ssa=0; ssa < refopt; ssa++) {
#ifdef DEBUG
			printf("...: Searching #%d rms: %g... ", ssa, rmslists[ssa]);
#endif
			for (sa=0; sa < npos; sa++) {
				if (rmslist[sa] == rmslists[ssa]) {
#ifdef DEBUG
					printf("found at sa %d, %g == %g\n", sa, rmslists[ssa], rmslist[sa]);
#endif
					(*list)[ssa] = sa;
					break;
				}
			}
		}
		return 0;
	} // refmode == REF_BESTRMS
	// We only arrive here when things go wrong
	return -1;
}

//
// Interpolation functions
//

int _9pquadint(float32_t *diffmap, int32_t diffsize[2], float32_t shvec[2], int32_t shran[2]) {
	// Find maximum
	float32_t max = diffmap[0], pix;
	int32_t maxidx[] = {0,0};
	int i, j;
	for (j=0; j<diffsize[1]; j++) {
		for (i=0; i<diffsize[0]; i++) {
			pix = diffmap[j * diffsize[0] + i];
			if (pix > max) {
				max = pix;
				maxidx[0] = i;
				maxidx[1] = j;
			}
		}
	}
#ifdef DEBUG
	printf("max @: %g (%d,%d) ", max, maxidx[0], maxidx[1]);
#endif
	if (maxidx[0] == 0 || maxidx[0] == diffsize[0]-1 ||
		maxidx[1] == 0 || maxidx[1] == diffsize[1]-1) {
			// Out of bound, interpolation failed
			shvec[0] = maxidx[0];
			shvec[1] = maxidx[1];
			return 1;
	}
	// Now interpolate around the maximum	
	float32_t a2 = 0.5 * (diffmap[(maxidx[1]) * diffsize[0] + maxidx[0] + 1] - \
		diffmap[(maxidx[1]) * diffsize[0] + maxidx[0]-1]);
	float32_t a3 = 0.5 * diffmap[(maxidx[1]) * diffsize[0] + maxidx[0] + 1] - \
		diffmap[(maxidx[1]) * diffsize[0] + maxidx[0]] + \
		0.5 * diffmap[(maxidx[1]) * diffsize[0] + maxidx[0]-1];
	float32_t a4 = 0.5 * (diffmap[(maxidx[1] + 1) * diffsize[0] + maxidx[0]] -
		diffmap[(maxidx[1] - 1) * diffsize[0] + maxidx[0]]);
	float32_t a5 = 0.5 * diffmap[(maxidx[1]+1) * diffsize[0] + maxidx[0]] -
		diffmap[(maxidx[1]) * diffsize[0] + maxidx[0]] + \
		0.5 * diffmap[(maxidx[1]-1) * diffsize[0] + maxidx[0]];
	float32_t a6 = 0.25 * (diffmap[(maxidx[1]+1) * diffsize[0] + maxidx[0]+1] -\
		diffmap[(maxidx[1]-1) * diffsize[0] + maxidx[0]+1] - \
		diffmap[(maxidx[1]+1) * diffsize[0] + maxidx[0]-1] + \
		diffmap[(maxidx[1]-1) * diffsize[0] + maxidx[0]-1]);
	
	shvec[0] = maxidx[0] + (2.0*a2*a5-a4*a6)/(a6*a6-4.0*a3*a5);
	shvec[1] = maxidx[1] + (2.0*a3*a4-a2*a6)/(a6*a6-4.0*a3*a5);
#ifdef DEBUG
	if (!(shvec[0] >= 0 && shvec[0] <= diffsize[0]) || !(shvec[1] >= 0 && shvec[1] <= diffsize[1])) {
		printf("NAN: %g,%g (%d,%d). %g - %g - %g - %g - %g\n", 
			shvec[0], shvec[1], maxidx[0], maxidx[1], a2, a3, a4, a5, a6);
		printf("a4: 0.5 * (%g - %g)\n", \
			diffmap[(maxidx[1] + 1) * diffsize[0] + maxidx[0]], \
			diffmap[(maxidx[1] - 1) * diffsize[0] + maxidx[0]]);
		printf("a5: 0.5 * %g - %g + 0.5 * %g\n", \
			diffmap[(maxidx[1] + 1) * diffsize[0] + maxidx[0]], \
			diffmap[(maxidx[1]) * diffsize[0] + maxidx[0]], \
			diffmap[(maxidx[1]-1) * diffsize[0] + maxidx[0]]);
		printf("a6: 0.25 * (%g - %g - %g + %g)\n", \
			diffmap[(maxidx[1]+1) * diffsize[0] + maxidx[0]+1], \
			diffmap[(maxidx[1]-1) * diffsize[0] + maxidx[0]+1], \
			diffmap[(maxidx[1]+1) * diffsize[0] + maxidx[0]-1], \
			diffmap[(maxidx[1]-1) * diffsize[0] + maxidx[0]-1]);
		printf("diff: %dx%d\n", diffsize[0], diffsize[1]);
		for (i=0; i<diffsize[0]*diffsize[1]; i++)
			printf("%g ", diffmap[i]);
		printf("\n");
		return 2;
	}
#endif
	return 0;
}

int _5pquadint(float32_t *diffmap, int32_t diffsize[2], float32_t shvec[2], int32_t shran[2]) {
	// Find maximum
	float32_t max = diffmap[0], pix;
	int32_t maxidx[] = {0,0};
	int i, j;
	for (j=0; j<diffsize[1]; j++) {
		for (i=0; i<diffsize[0]; i++) {
			pix = diffmap[j * diffsize[0] + i];
			if (pix > max) {
				max = pix;
				maxidx[0] = i;
				maxidx[1] = j;
			}
		}
	}
#ifdef DEBUG
	printf("max: %g (%d,%d) ", max, maxidx[0], maxidx[1]);
#endif
	if (maxidx[0] == 0 || maxidx[0] == diffsize[0]-1 ||
		maxidx[1] == 0 || maxidx[1] == diffsize[1]-1) {
			// Out of bound, interpolation failed
			shvec[0] = maxidx[0];
			shvec[1] = maxidx[1];
			return 1;
	}
	// Now interpolate around the maximum	
	// a2 = f(x)
	float32_t a2 = 0.5 * (diffmap[(maxidx[1]) * diffsize[0] + maxidx[0] + 1] - \
		diffmap[(maxidx[1]) * diffsize[0] + maxidx[0]-1]);
	// a3 = f(x)
	float32_t a3 = 0.5 * diffmap[(maxidx[1]) * diffsize[0] + maxidx[0] + 1] - \
		diffmap[(maxidx[1]) * diffsize[0] + maxidx[0]] + \
		0.5 * diffmap[(maxidx[1]) * diffsize[0] + maxidx[0]-1];
	// a4 = f(y) -- 0.5 * map[x, y+1] - 0.5 * map[x, y-1]
	float32_t a4 = 0.5 * (diffmap[(maxidx[1] + 1) * diffsize[0] + maxidx[0]] -
		diffmap[(maxidx[1] - 1) * diffsize[0] + maxidx[0]]);
	// a5 = f(y) -- 0.5 * map[x, y+1] - map[x, y] + 0.5 * map[x, y-1]
	float32_t a5 = 0.5 * diffmap[(maxidx[1]+1) * diffsize[0] + maxidx[0]] -
		diffmap[(maxidx[1]) * diffsize[0] + maxidx[0]] + \
		0.5 * diffmap[(maxidx[1]-1) * diffsize[0] + maxidx[0]];
	
	shvec[0] = maxidx[0] + (2.0*a2)/(-4.0*a3);
	shvec[1] = maxidx[1] + (2.0*a4)/(-4.0*a5);
#ifdef DEBUG
	if (!(shvec[0] >= 0 && shvec[0] <= diffsize[0]) || !(shvec[1] >= 0 && shvec[1] <= diffsize[1])) {
		printf("NAN: %g,%g (%d,%d). %g - %g - %g - %g\n", 
			shvec[0], shvec[1], maxidx[0], maxidx[1], a2, a3, a4, a5);
		printf("a4: 0.5 * (%g - %g)\n", \
			diffmap[(maxidx[1] + 1) * diffsize[0] + maxidx[0]], \
			diffmap[(maxidx[1] - 1) * diffsize[0] + maxidx[0]]);
		printf("a5: 0.5 * %g - %g + 0.5 * %g\n", \
			diffmap[(maxidx[1] + 1) * diffsize[0] + maxidx[0]], \
			diffmap[(maxidx[1]) * diffsize[0] + maxidx[0]], \
			diffmap[(maxidx[1]-1) * diffsize[0] + maxidx[0]]);
		printf("a6: 0.25 * (%g - %g - %g + %g)\n", \
			diffmap[(maxidx[1]+1) * diffsize[0] + maxidx[0]+1], \
			diffmap[(maxidx[1]-1) * diffsize[0] + maxidx[0]+1], \
			diffmap[(maxidx[1]+1) * diffsize[0] + maxidx[0]-1], \
			diffmap[(maxidx[1]-1) * diffsize[0] + maxidx[0]-1]);
		printf("diff: %dx%d\n", diffsize[0], diffsize[1]);
		for (i=0; i<diffsize[0]*diffsize[1]; i++)
			printf("%g ", diffmap[i]);
		printf("\n");
		return 2;
	}
#endif
	return 0;
}

//
// Image comparison functions
//

int _sqdiff(float32_t *img, int32_t imgsize[2], int32_t imstride, float32_t *ref, int32_t refsize[2], int32_t refstride, int32_t *mask, int32_t maskstride, float32_t *diffmap, int32_t pos[2], int32_t range[2], int bigref) {
	double tmpsum, diff;
	// Loop ranges
	int sh0min = -range[0] + pos[0];
	int sh0max =  range[0] + pos[0];
	int sh1min = -range[1] + pos[1];
	int sh1max =  range[1] + pos[1];
	
	// If bigref is 1, we expect that the reference is large enough to naively 
	// shift img around.
	int sh0, sh1, i, j;
	if (bigref == 1) {
		// Loop over all shifts to be tested
		for (sh0=sh0min; sh0 <= sh0max; sh0++) {
			for (sh1=sh1min; sh1 <= sh1max; sh1++) {
				// Loop over all pixels within the img and refimg, and compute
				// the cross correlation between the two.
				tmpsum = 0.0;
				for (j=0; j<imgsize[1]; j++) {
					for (i=0; i<imgsize[0]; i++) {
						// Check the mask if we need to correlate this pixel
//						if (mask[j*maskstride + i] == 1) {
							// First get the difference...
							diff = img[j*imstride + i] - ref[(j+sh1)*refstride + i+sh0];
							// ...then square this
						tmpsum += mask[j*maskstride + i]*diff*diff;
//						}
					}
				}
				// Store the current correlation value in the map. Use
				// negative value to ensure that we get a maximum for best
				// match in diffmap (this allows to use a more general 
				// maximum-finding method, instead of splitting between maxima
				// and minima)
				diffmap[(sh1-sh1min) * (range[0]*2+1) + (sh0-sh0min)] = -tmpsum;
			}
		}
	}
	// If bigref is zero, we use clipping of ref and img to only use the 
	// intersection of the two datasets for comparison, using normalisation to 
	// make the results consistent. This does increase the noise in the shift 
	// measurement though.
	else {
		for (sh0=sh0min; sh0 <= sh0max; sh0++) {
			for (sh1=sh1min; sh1 <= sh1max; sh1++) {
				tmpsum = 0.0;
				// If the shift to compare is negative, we must make sure 
				// i+sh0 in 'ref' will not be negative, so we start i at -sh0.
				// If the shift is positive, we must make sure that i+sh0 in 
				// 'ref' will not go out of bound, so we stop i at Nimg[0] - 
				// sh0.
				// N<array>[<index>] gives the <index>th size of <array>, i.e. 
				// Nimg[1] gives the second dimension of array 'img'
				for (j=0 - min(sh0, 0); j<imgsize[1] - max(sh0, 0); j++){
					for (i=0 - min(sh1, 0); i<imgsize[0] - max(sh1, 0); i++) {
						// Check the mask if we need to correlate this pixel
//						if (mask[j*maskstride + i] == 1) {
							// First get the difference...
							diff = img[j*imstride + i] - ref[(j+sh1)*refstride + i+sh0];
							//diff = img(i,j) - ref(i+sh1,j+sh0);
							// ...then square this
							tmpsum += mask[j*maskstride + i]*diff*diff;
//						}
					}
				}
				// Scale the value found by dividing it by the number of 
				// pixels we compared.
				diffmap[(sh1-sh1min) * (range[0]*2+1) +(sh0-sh0min)] = -tmpsum /
					((imgsize[0]-abs(sh0)) * (imgsize[1]-abs(sh1)));
			}
		}
	}
	
	return 0;
}

int _absdiffsq(float32_t *img, int32_t imgsize[2], int32_t imstride, float32_t *ref, int32_t refsize[2], int32_t refstride, int32_t *mask, int32_t maskstride, float32_t *diffmap, int32_t pos[2], int32_t range[2], int bigref) {
	double tmpsum;
	// Loop ranges
	int sh0min = -range[0] + pos[0];
	int sh0max =  range[0] + pos[0];
	int sh1min = -range[1] + pos[1];
	int sh1max =  range[1] + pos[1];
	
	// If bigref is 1, we expect that the reference is large enough to naively 
	// shift img around.
	int sh0, sh1, i, j;
	double max = 0;
	int maxidx[] = {-1,-1};
	if (bigref == 1) {
		// Loop over all shifts to be tested
		for (sh1=sh1min; sh1 <= sh1max; sh1++) {
			for (sh0=sh0min; sh0 <= sh0max; sh0++) {
				tmpsum = 0.0;
				for (j=0; j<imgsize[1]; j++) {
					for (i=0; i<imgsize[0]; i++) {
//						if (mask[j*maskstride + i] == 1) {
							tmpsum += mask[j*maskstride + i] *\
								fabs(img[j*imstride + i] - ref[(j+sh1)*refstride + i+sh0]);
//						}
					}
				}
				// if (max == 0)
				// 	max = -(tmpsum*tmpsum);
				// if (-(tmpsum*tmpsum) > max) {
				// 	max = -(tmpsum*tmpsum);
				// 	maxidx[0] = sh0-sh0min;
				// 	maxidx[1] = sh1-sh1min;
				// }
				diffmap[(sh1-sh1min) * (range[0]*2+1) + (sh0-sh0min)] = \
				 	-(tmpsum*tmpsum);
			}
		}
	}
	// If bigref is zero, we use clipping of ref and img to only use the 
	// intersection of the two datasets for comparison, using normalisation to 
	// make the results consistent. This does increase the noise in the shift 
	// measurement though.
	else {
		for (sh0=sh0min; sh0 <= sh0max; sh0++) {
			for (sh1=sh1min; sh1 <= sh1max; sh1++) {
				tmpsum = 0.0;
				for (j=0 - min(sh0, 0); j<imgsize[1] - max(sh0, 0); j++){
					for (i=0 - min(sh1, 0); i<imgsize[0] - max(sh1, 0); i++) {
//						if (mask[j*maskstride + i] == 1) {
							tmpsum += mask[j*maskstride + i] *\
								fabsf(img[j*imstride + i] - ref[(j+sh1)*refstride + i+sh0]);
//						}
					}
				}
				// Scale the value found by dividing it by the number of 
				// pixels we compared.
				diffmap[(sh1-sh1min) * (range[0]*2+1) +(sh0-sh0min)] = \
				 	-(tmpsum*tmpsum) / ((imgsize[0]-abs(sh0)) * (imgsize[1]-abs(sh1)));
			}
		}
	}
	
	return 0;
}

int _crosscorr(float32_t *img, int32_t imgsize[2], int32_t imstride, float32_t *ref, int32_t refsize[2], int32_t refstride, int32_t *mask, int32_t maskstride, float32_t *diffmap, int32_t pos[2], int32_t range[2], int bigref) {
	double tmpsum;
	// Loop ranges
	int sh0min = -range[0] + pos[0];
	int sh0max =  range[0] + pos[0];
	int sh1min = -range[1] + pos[1];
	int sh1max =  range[1] + pos[1];
	
	// If bigref is 1, we expect that the reference is large enough to naively 
	// shift img around.
	int sh0, sh1, i, j;
	if (bigref == 1) {
		// Loop over all shifts to be tested
		for (sh0=sh0min; sh0 <= sh0max; sh0++) {
			for (sh1=sh1min; sh1 <= sh1max; sh1++) {
				tmpsum = 0.0;
				for (j=0; j<imgsize[1]; j++) {
					for (i=0; i<imgsize[0]; i++) {
//						if (mask[j*maskstride + i] == 1) {
							tmpsum += mask[j*maskstride + i] * \
								(img[j*imstride + i] * ref[(j+sh1)*refstride + i+sh0]) * \
								(img[j*imstride + i] * ref[(j+sh1)*refstride + i+sh0]);
//						}
					}
				}
				diffmap[(sh1-sh1min) * (range[0]*2+1) + (sh0-sh0min)] = tmpsum;
			}
		}
	}
	// If bigref is zero, we use clipping of ref and img to only use the 
	// intersection of the two datasets for comparison, using normalisation to 
	// make the results consistent. This does increase the noise in the shift 
	// measurement though.
	else {
		for (sh0=sh0min; sh0 <= sh0max; sh0++) {
			for (sh1=sh1min; sh1 <= sh1max; sh1++) {
				tmpsum = 0.0;
				for (j=0 - min(sh0, 0); j<imgsize[1] - max(sh0, 0); j++){
					for (i=0 - min(sh1, 0); i<imgsize[0] - max(sh1, 0); i++) {
//						if (mask[j*maskstride + i] == 1) {
							tmpsum += mask[j*maskstride + i] * \
								(img[j*imstride + i] * ref[(j+sh1)*refstride + i+sh0]) * \
								(img[j*imstride + i] * ref[(j+sh1)*refstride + i+sh0]);
//						}
					}
				}
				// Scale the value found by dividing it by the number of 
				// pixels we compared.
				diffmap[(sh1-sh1min) * (range[0]*2+1) +(sh0-sh0min)] = \
					tmpsum / ((imgsize[0]-abs(sh0)) * (imgsize[1]-abs(sh1)));
			}
		}
	}
	
	return 0;
}
