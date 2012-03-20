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

#ifndef HAVE_LIBUNWRAP_PY_H
#define HAVE_LIBUNWRAP_PY_H

//
// Headers
//

#include <Python.h>				// For python extension
#include <numpy/arrayobject.h> 	// For numpy

static PyObject *libunwrap_helloworld(PyObject *self, PyObject *args);
static PyObject *libunwrap_floodfill(PyObject *self, PyObject *args);

#endif // HAVE_LIBUNWRAP_PY_H
