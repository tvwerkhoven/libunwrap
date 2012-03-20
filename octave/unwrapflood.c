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

#include <string.h>
#include <math.h>
#include "mex.h"
#define BORDERLIMITFRAC_DEFAULT 0.9

// #define M_PI  3.141592653589793238462643383279502884197169399



double *g_wrapped  = NULL;
double *g_pup      = NULL;
double *g_donemask = NULL;
double  g_borderlimitfrac = 0;
long    g_count = 0;
long    g_phdim = 0;
int     g_usedonemask = 0;

//
// Usage: unwrapflood(phsm, pup, po1, po2)
//
void
mexFunction (int nlhs, mxArray *plhs[], int nrhs, 
	     const mxArray *prhs[])
{
  mxArray *mxWrapped;
  double  *phsm;
  int      po1, po2;

  if (nrhs != 6) {
    mexErrMsgTxt("Usage: unwrapflood phase pup po1 po2 usedonemask borderlimitfrac.");
  }


  // Get the phase dimension
  int ndims = (int)(mxGetNumberOfDimensions(prhs[0]));
  if (ndims != 2) {
    mexErrMsgTxt("Phase needs to be a matrix.");
  }
  int dim1 = (int)(mxGetDimensions(prhs[0]))[0];
  int dim2 = (int)(mxGetDimensions(prhs[0]))[1];
  if (ndims != 2 || dim1!= dim2) {
    mexErrMsgTxt("Phase needs to be a square matrices.");
  }
  g_phdim = dim1;

  ndims = (int)(mxGetNumberOfDimensions(prhs[1]));
  if (ndims != 2) {
    mexErrMsgTxt("Pupil needs to be a matrix.");
  }
  dim1 = (int)(mxGetDimensions(prhs[1]))[0];
  dim2 = (int)(mxGetDimensions(prhs[1]))[1];
  if (dim1!= dim2) {
    mexErrMsgTxt("Pupil needs to be a square matrices.");
  }
  if (g_phdim != dim1) {
    mexErrMsgTxt("Pupil needs same dimension.");
  }


  dim1 = (int)(mxGetDimensions(prhs[2]))[0];
  dim2 = (int)(mxGetDimensions(prhs[2]))[1];
  if (dim1 != 1 || dim2 != 1) {
    mexErrMsgTxt("Position 1 needs to be scalar.");
  }
  dim1 = (int)(mxGetDimensions(prhs[3]))[0];
  dim2 = (int)(mxGetDimensions(prhs[3]))[1];
  if (dim1 != 1 || dim2 != 1) {
    mexErrMsgTxt("Position 2 needs to be scalar.");
  }

  dim1 = (int)(mxGetDimensions(prhs[4]))[0];
  dim2 = (int)(mxGetDimensions(prhs[4]))[1];
  if (dim1 != 1 || dim2 != 1) {
    mexErrMsgTxt("Usedonemask needs to be scalar.");
  }

  dim1 = (int)(mxGetDimensions(prhs[5]))[0];
  dim2 = (int)(mxGetDimensions(prhs[5]))[1];
  if (dim1 != 1 || dim2 != 1) {
    mexErrMsgTxt("borderlimitfrac needs to be a scalar.");
  }


  po1 = (int)(mxGetPr(prhs[2]))[0];
  po2 = (int)(mxGetPr(prhs[3]))[0];
  g_usedonemask = (int)(mxGetPr(prhs[4]))[0];
  g_borderlimitfrac = (double)(mxGetPr(prhs[5]))[0];

  // mexPrintf("pos: %d, %d, g_phdim: %ld, usedonemask: %d\n", po1, po2, g_phdim, g_usedonemask);


  phsm  = mxGetPr(prhs[0]);
  g_pup = mxGetPr(prhs[1]);
  
  mxWrapped = mxCreateDoubleMatrix(g_phdim, g_phdim, mxREAL);
  g_wrapped = mxGetPr(mxWrapped);
  memcpy(g_wrapped, phsm, sizeof(double)*g_phdim*g_phdim);

  g_donemask = calloc(g_phdim*g_phdim, sizeof(double));

  // Start the recursive unwrapping
  unwrapflood(po1, po2, 0);

  plhs[0] = mxWrapped;
} // mexFunction




//
// Try to make a floodfill unwrap.
//
void unwrapflood(int po1, int po2, int itco)
{

  //int i1;
  //for (i1=0; i1<128;i1++)
  //g_wrapped[i1 + 60*128]= 0;

  int    i1, i2, itcolim;
  int    curpo1, curpo2;
  double thestep;

  // How many recursions we permit?
  itco = itco+1;
  if (g_usedonemask == 0)
    itcolim = 300;
  else
    itcolim = 300; // g_phdim*g_phdim;
  if (itco > itcolim)
    return;
    
  // Handle neighbors
  for (i1=-1; i1<=1; i1++) {
    curpo1=po1-i1;

    if (curpo1 >= 0 && curpo1 < g_phdim) {

      for (i2=-1; i2<=1; i2++) {
	curpo2=po2-i2;
    
	if (curpo2 >= 0 && curpo2 < g_phdim &&
	    g_pup[     curpo1 + curpo2*g_phdim] > 1e-9 &&
	    (g_donemask[curpo1 + curpo2*g_phdim] != 1 || g_usedonemask==0)) {

	  // && pup(curpo1,curpo2)~=0 && donemask(curpo1,curpo2)~=1
	
	  thestep = g_wrapped[po1+po2*g_phdim] - g_wrapped[curpo1+curpo2*g_phdim];

	  // mexPrintf("thestep: %e\n", thestep);
	  if (fabs(thestep) >= g_borderlimitfrac*2*M_PI) {

	    g_wrapped[curpo1+curpo2*g_phdim] = 
	      g_wrapped[curpo1+curpo2*g_phdim] + 2*M_PI*round(thestep/(2*M_PI));
	    g_donemask[curpo1+curpo2*g_phdim]=1;

	    g_count++;

	    unwrapflood(curpo1, curpo2, itco+1);
	  }	
	}
      }
    }
  }
} // unwrapflood
