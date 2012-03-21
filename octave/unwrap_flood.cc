/*
 unwrap_flood.c -- Octave binding for libunwrap
 Copyright (C) 2012 Visa Korkiakoski <korkiakoski@strw.leidenuniv.nl>
 
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

#include "string.h"
#include "oct.h"

extern "C"
{
#include "libunwrap.h"
}  /* end extern "C" */




// void
// mexFunction (int nlhs, mxArray *plhs[], int nrhs, 
//              const mxArray *prhs[])


/*!
 @brief Octave interface for the floodfill unwrapping routine. 
 
 OBSOLETE!

 @author Visa Korkiakoski
*/
DEFUN_DLD(unwrap_flood, args, nargout, 
	  "unwrapped = unwrap_flood(phase, pup, po1, po2, usedonemask, borderlimitfrac)")
{
  const double  *phsm;
  int            po1, po2;

  int nargin = args.length();
  if (nargin != 6) {
    print_usage();
    return octave_value_list();
  }

  // Get the phase dimension
  int ndims = args(0).dims().length();
  if (ndims != 2) {
    error("Phase needs to be a matrix.");
    return octave_value_list();
  }
  int dim1 = (args(0).dims())(0);
  int dim2 = (args(0).dims())(1);
  if (ndims != 2 || dim1!= dim2) {
    error("Phase needs to be a square matrices.");
    return octave_value_list();
  }
  g_phdim = dim1;

  ndims = args(1).dims().length();
  if (ndims != 2) {
    error("Pupil needs to be a matrix.");
    return octave_value_list();
  }
  dim1 = (args(1).dims())(0);
  dim2 = (args(1).dims())(1);
  if (dim1!= dim2) {
    error("Pupil needs to be a square matrices.");
    return octave_value_list();
  }
  if (g_phdim != dim1) {
    error("Pupil needs same dimension.");
    return octave_value_list();
  }


  dim1 = (args(2).dims())(0);
  dim2 = (args(2).dims())(1);
  if (dim1 != 1 || dim2 != 1) {
    error("Position 1 needs to be scalar.");
    return octave_value_list();
  }
  dim1 = (args(3).dims())(0);
  dim2 = (args(3).dims())(1);
  if (dim1 != 1 || dim2 != 1) {
    error("Position 2 needs to be scalar.");
    return octave_value_list();
  }

  dim1 = (args(4).dims())(0);
  dim2 = (args(4).dims())(1);
  if (dim1 != 1 || dim2 != 1) {
    error("Usedonemask needs to be scalar.");
    return octave_value_list();
  }

  dim1 = (args(5).dims())(0);
  dim2 = (args(5).dims())(1);
  if (dim1 != 1 || dim2 != 1) {
    error("borderlimitfrac needs to be a scalar.");
    return octave_value_list();
  }


  po1 = (int)(args(2).array_value())(0);
  po2 = (int)(args(3).array_value())(0);
  g_usedonemask = (int)(args(4).array_value())(0);
  g_borderlimitfrac = (args(5).array_value())(0);

  phsm  = args(0).array_value().data();
  g_pup = args(1).array_value().data();
  
  Matrix mxWrapped = Matrix(g_phdim, g_phdim);
  g_wrapped = mxWrapped.fortran_vec();
  memcpy(g_wrapped, phsm, sizeof(double)*g_phdim*g_phdim);
  g_donemask = (double *)calloc(g_phdim*g_phdim, sizeof(double));

#if(0)
  octave_stdout << "po1: " << po1 << "  po2: "<< po2<< "\n";
  octave_stdout << "phsm: " << phsm <<"\n";
  octave_stdout << "g_pup: " << g_pup <<"\n";
  octave_stdout << "g_wrapped: " << g_wrapped <<"\n";
  octave_stdout << "phdim: " << g_phdim <<"\n";
  octave_stdout << "donemask: " << g_donemask <<"\n";
#endif

  // Start the recursive unwrapping
  unwrap_flood(po1, po2, 0);

  return octave_value(mxWrapped);
}
