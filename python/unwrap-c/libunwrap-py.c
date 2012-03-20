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
// #include <sys/time.h>			// For timestamps
//#include <time.h>					// For timestamps
//#include <math.h>					// For pow()

#include "debugprint.h"
#include "libunwrap-py.h"

//
// Methods table for this module
//

static PyMethodDef LibunwrapMethods[] = {
	{"helloworld",  libunwrap_helloworld, METH_VARARGS, "Hello World routine."},
	{"floodfill",  libunwrap_floodfill, METH_VARARGS, "Flood-fill unwrap."},
	{NULL, NULL, 0, NULL}        /* Sentinel */
};

//
// Init module methods
//
PyMODINIT_FUNC initunwrap_c(void) {
	(void) Py_InitModule("unwrap_c", LibunwrapMethods);
	// Init numpy usage
	import_array();
}

//
// Main python functions
//
static PyObject * libunwrap_helloworld(PyObject *self, PyObject *args) {
  fprintf(stderr, "Hello world!\n");
  Py_INCREF(Py_None);
  return Py_None; 
}

static PyObject *libunwrap_floodfill(PyObject *self, PyObject *args) {
  PyArrayObject* phase;     // Wrappted input phase
	int startx, starty;       // Starting coordinates
  PyArrayObject *ret_phase; // Unwrapped output phase
  
  //
	// Parse arguments from Python function
	//
	if (!PyArg_ParseTuple(args, "O!ii", 
                        &PyArray_Type, &phase,    // Wrapped phase
                        &startx,									// Start x-coordinate
                        &starty                   // Start y-coordinate
                        )) {
		PyErr_SetString(PyExc_SyntaxError, "In floodfill: failed to parse arguments.");
		return NULL; 
  }

  DEBUGPRINT("phase: 0x%p, start: %d,%d\n", phase, startx, starty);
 
  int nd = PyArray_NDIM(phase);
  if (nd != 2) {
		PyErr_SetString(PyExc_RuntimeError, "In floodfill: can only work with 2-d phase.");
		return NULL;
  }
  int ph_w = (int) PyArray_DIM((PyObject*) phase, 0);
	int ph_h = (int) PyArray_DIM((PyObject*) phase, 1);

  DEBUGPRINT("#dim: %d, dims: %d, %d, size: %d\n", nd, ph_w, ph_h, ph_w*ph_h);

  Py_INCREF(Py_None);
  return Py_None; 
//  // Build numpy vector from C array 'phase_uw', a 2-d array
//	npy_intp s_dims[] = {512, 256};
//	ret_phase = (PyArrayObject*) PyArray_SimpleNewFromData(4, s_dims, NPY_FLOAT32, (void *) phase_uw);
//	PyArray_FLAGS(ret_phase) |= NPY_OWNDATA;
//	
//	if (!PyArray_CHKFLAGS(ret_phase, NPY_OWNDATA)) {
//		PyErr_SetString(PyExc_RuntimeError, "In floodfill: unable to set dat ownership for 'phase_uw'.");
//		free(ret_phase);
//		return NULL;
//  }
//	
//	return Py_BuildValue("{s:N,s:N}", "phase", ret_phase, "refapts", retreflist);  
}