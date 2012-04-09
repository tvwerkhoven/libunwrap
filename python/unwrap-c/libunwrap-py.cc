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

#include <stdlib.h>

#include <Python.h>				// For python extension
#include <numpy/arrayobject.h> 	// For numpy
//#include <Numeric/arrayobject.h> 	// For numpy
//#include <numarray/arrayobject.h> 	// For numpy

#include "debugprint.h"
#include "libunwrap-py.h"

#include "libunwrap.h"

//
// Methods table for this module
//

static PyMethodDef LibunwrapMethods[] = {
	{"helloworld",  libunwrap_helloworld, METH_VARARGS, "Hello World routine."},
	{"flood_quality",  libunwrap_flood_quality, METH_VARARGS, "Quality guided flood-fill unwrap."},
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

static PyObject *libunwrap_flood_quality(PyObject *self, PyObject *args) {
  PyArrayObject* phase;     // Wrapped input phase
  PyArrayObject* qual;      // Quality map
  
  PyObject *ph_unwrap_obj;
  //
	// Parse arguments from Python function.
	//
	if (!PyArg_ParseTuple(args, "O!O!", 
                        &PyArray_Type, &phase,    // Wrapped phase
                        &PyArray_Type, &qual      // Quality map
                        )) {
		PyErr_SetString(PyExc_SyntaxError, "floodfill: failed to parse args.");
		return NULL; 
  }

  DEBUGPRINT("phase: 0x%p, qual: 0x%p\n", phase, qual);

  // Inspect array dimensions
  int ph_0 = (int) PyArray_DIM((PyObject*) phase, 0);
  int ph_1 = (int) PyArray_DIM((PyObject*) phase, 1);

  
  DEBUGPRINT("#dim: %d, dims: %d, %d, size: %d\n", nd, ph_0, ph_1, ph_0*ph_1);

  switch (PyArray_TYPE((PyObject *) phase)) {
		case (NPY_FLOAT64): {
      DEBUGPRINT("%s\n", "NPY_FLOAT64");

      // Copy input phase to ensure contiguous memory
			ph_unwrap_obj = PyArray_FROM_OTF((PyObject *) phase, NPY_DOUBLE, NPY_ENSURECOPY | NPY_OUT_ARRAY);
      
      // Make sure the quality map is well-behaved
      PyObject *qual64_obj = PyArray_FROM_OTF((PyObject *) qual, NPY_DOUBLE, NPY_IN_ARRAY);
      
      // Get pointers to data
      double *ph_uw = (double *) PyArray_DATA(ph_unwrap_obj);
      double *qual64 = (double *) PyArray_DATA(qual64_obj);
      
      // Got data, call floodfill now. We reverse the dimensions (why?)
      unwrap_flood_quality(ph_uw, qual64, ph_1, ph_0);
      
      // We don't need a quality map reference anymore
      Py_DECREF(qual64_obj);
			break;
		}
		default: {
      DEBUGPRINT("%s\n", "unsupported type");
			PyErr_SetString(PyExc_NotImplementedError, "In floodfill: datatype not supported.");
			return NULL;
		}
	}
  
  return Py_BuildValue("N", ph_unwrap_obj);
}

