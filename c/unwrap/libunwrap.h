/*
 libunwrap.h -- Unwrap library header file
 Copyright (C) 2012 Visa Korkiakoski <korkiakoski@strw.leidenuniv.nl>
   & Tim van Werkhoven <werkhoven@strw.leidenuniv.nl>
 
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

#ifndef HAVE_LINBUNWRAP_H
#define HAVE_LINBUNWRAP_H

// Global variables -- they are used only by the obsolete flood
// filling routine.
extern const double *g_pup;
extern double *g_wrapped;
extern double *g_donemask;
extern double  g_borderlimitfrac;
extern long    g_count;
extern long    g_phdim;
extern int     g_usedonemask;

/*!
  Data structure used internally by the quality guided unwrapping.
 */
typedef struct {
  int  *doneMask;     //!< Mask counting the pixels that are already unwrapped 
  int   phdim;        //!< Dimension of the phase
  int   listsz;       //!< Number of elements in the list
  int   unwcount;     //!< Number of unwrapped points

  int  *borderListPrevs; //!< Array contining indexes to the previous nodes in the list linking the border pixels. Has same dimension as the phase.

  int  *borderListNexts; //!< Array contining indexes to the next nodes in the list linking the border pixels. Has same dimension as the phase

  int   borderListFirst; //!< Index to the first node in the list linking the border pixels. */
} unwrapqdata_t;


void unwrap_flood(int po1, int po2, int itco);

void unwrap_flood_quality(double *ph, const double *quality, int phdim);

// These functions are used internally 
int findmax(const double *arr, int len);
double valid_neighs_getmean(int po1, int po2, double *ph, int *doneMask, int phdim);
void floodborder_add(unwrapqdata_t *uwd, int po1, int po2);
void floodborder_remove(unwrapqdata_t *uwd, const double *quality, int po1, int po2);
int floodborder_findmaxneighbor(unwrapqdata_t *uwd, const double *quality, 
				int *maxpo1, int *maxpo2);

#endif // HAVE_LINBUNWRAP_H
