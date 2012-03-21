/*
 unwrapflood.c -- Octave binding for libunwrap
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





/*!
 @brief Octave interface for the quality based unwrapping function. 
 @author Visa Korkiakoski
*/
DEFUN_DLD(unwrap_qua, args, nargout, 
	  "unwrapped = unwrap_oct(phase, quality)")
{
  const double  *phsm;
  const double  *qua;
  int      phdim;
  double  *wrapped;

  int nargin = args.length();
  if (nargin != 2) {
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
    error("Phase needs to be a square matrix.");
    return octave_value_list();
  }
  phdim = dim1;

  ndims = args(1).dims().length();
  if (ndims != 2) {
    error("Quality needs to be a matrix.");
    return octave_value_list();
  }
  dim1 = (args(1).dims())(0);
  dim2 = (args(1).dims())(1);
  if (dim1!= dim2) {
    error("Quality needs to be a square matrix.");
    return octave_value_list();
  }
  if (phdim != dim1) {
    error("Pupil needs same dimension.");
    return octave_value_list();
  }

  phsm = args(0).array_value().data();
  qua  = args(1).array_value().data();
  
  Matrix mxWrapped = Matrix(phdim, phdim);
  wrapped = mxWrapped.fortran_vec();
  memcpy(wrapped, phsm, sizeof(double)*phdim*phdim);

#if(0)
  octave_stdout << "phsm: " << phsm <<"\n";
  octave_stdout << "qua: " << qua <<"\n";
  octave_stdout << "wrapped: " << wrapped <<"\n";
#endif

  // Start the recursive unwrapping
  unwrap_quality(wrapped, qua, phdim);

  return octave_value(mxWrapped);
} // unwrap_qua
