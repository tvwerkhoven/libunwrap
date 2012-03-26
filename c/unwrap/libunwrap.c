/*
 libunwrap.cc -- Phase unwrapping library
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

#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "libunwrap.h"

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399
#endif


// These are used by unwrap_flood()
#define BORDERLIMITFRAC_DEFAULT 0.9
const double *g_pup      = NULL;
double *g_wrapped  = NULL;
double *g_donemask = NULL;
double  g_borderlimitfrac = 0;
long    g_count = 0;
long    g_phdim = 0;
int     g_usedonemask = 0;

/*!
 @brief Try to make a floodfill unwrap.
 @author Visa Korkiakoski

 This function is obsolete, please use unwrap_quality() instead.
*/
void unwrap_flood(int po1, int po2, int itco) 
{  
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
            
            unwrap_flood(curpo1, curpo2, itco+1);
          }	
        }
      }
    }
  }
} // unwrap_flood(int po1, int po2, int itco) 


/*!
 @brief Quality guided unwrapping
 @author Visa Korkiakoski

 Unwraps the phase in the order of quality array. The highest quality
 pixels will be unwrapped first. To find out which pixel is unwrapped
 (based on detecting a (n x 2PI) phase shift comparing to the values
 of the neighboring pixels), this algorithm goes through all the
 neighbors of already unwrapped pixels and finds the one with highest
 quality.

 TODO: add references to literature & publications.

 @param ph [in, out] phase to be unwrapped -- this array will be modified

 @param quality [in] quality array same size as phase. 0 means no phase information available, higher values indicate higher reliability

 @param phdim [in] dimension of ph and quality. They are assumed to be square matrices.
*/
void unwrap_flood_quality(double *ph, const double *quality, int phdim) 
{
  unwrapqdata_t  uwd = {0};
  int  i, i1, i2, maxi;
  int neighbo1, neighbo2;
  double meaval, thestep;

  uwd.phdim    = phdim;
  uwd.listsz   = 0;
  uwd.unwcount = 0;
  uwd.doneMask = calloc(phdim*phdim, sizeof(int));
  uwd.borderListPrevs = calloc(phdim*phdim, sizeof(int));
  uwd.borderListNexts = calloc(phdim*phdim, sizeof(int));
  uwd.borderListFirst = -1; // No border
  for (i=0; i<phdim*phdim; i++) {
    uwd.borderListPrevs[i] = -1; // point to nil
    uwd.borderListNexts[i] = -1; // point to nil
  }
  
  // First starting point: start with the pixel with the highest quality
  maxi = findmax(quality, phdim*phdim);
  uwd.doneMask[maxi] = 1;
  uwd.unwcount++;
  neighbo1 = maxi % phdim;
  neighbo2 = maxi / phdim;
  // This point (neighbo1, neighbo2) is now unwrapped, add to the border list
  floodborder_add(&uwd, neighbo1, neighbo2);
  
  // Find the neighbor of the border of the unwrapped phase with the highest quality
  while (floodborder_findmaxneighbor(&uwd, quality, &neighbo1, &neighbo2) >= 0) {
    // Pixel under study: neighbo1, neighbo2
    
    // If the phase of the current testpixel differs more than pi from the unwrapped pixels, we need to unwrap the current pixel
    meaval = valid_neighs_getmean(neighbo1, neighbo2, ph, uwd.doneMask, phdim);
    thestep = meaval - ph[neighbo1+neighbo2*phdim];
    if (fabs(thestep) > M_PI) {
      ph[neighbo1 + phdim*neighbo2] += 2*M_PI*round(thestep/(2*M_PI));
    }
    // Pixel is now unwrapped, mark as completed
    uwd.doneMask[neighbo1 + phdim*neighbo2] = 1;
    uwd.unwcount++;
    
    // The new neighboring point needs to be added to the flooding border
    floodborder_add(&uwd, neighbo1, neighbo2);

    // When enough pixels are unwrapped, some unwrapped pixels are not at the 
    // border anymore, remove these from the list to prevent needless computations
    for (i1=neighbo1-1; i1<=neighbo1+1; i1++) {
      if (i1>=0 && i1 < phdim) {
        for (i2=neighbo2-1; i2<=neighbo2+1; i2++) {
          if (i2>=0 && i2 < phdim) {
            // Removes a point from the flooding border, if all neighbors are 
            // already unwrapped.
            floodborder_remove(&uwd, quality, i1, i2);
          }
        }
      }
    }    

  } // flooding border has valid neighbors

  // Free the resources
  free(uwd.doneMask);
  free(uwd.borderListPrevs);
  free(uwd.borderListNexts);

} // unwrap_flood_quality(double *ph, const double *quality, int phdim) 



/*!
 @brief Finds the index of maximum value in an array
 @author Visa Korkiakoski
 */
int findmax(const double *arr, int len)
{
  double maxval=arr[0];
  int i, maxi=0;
  for (i=1; i<len; i++) {
    if (arr[i] > maxval) {
      maxval = arr[i];
      maxi = i;
    }
  }
  return maxi;
} // findmax



/*!
 @brief Finds the mean value of neighbour points for which doneMask equals 1
 @author Visa Korkiakoski
  
 @todo Document parameters
 
 @param [in] po1 Coordinate to check (x)
 @param [in] po2 Coordinate to check (y)
 @param [in] *ph 2D phase information
 @param [in] *doneMask 
 @param [in] phdim Size of phase 
 */
double valid_neighs_getmean(int po1, int po2, double *ph, int *doneMask, int phdim)
{
  int    i1, i2, inds=0;
  double meaval=0;

  // Loop over all neighbours of point (po1, po2). For all neighbours with 
  // doneMask equal to 1, sum the values.
  for (i1=po1-1; i1<=po1+1; i1++) {
    if (i1>=0 && i1 < phdim) {
      for (i2=po2-1; i2<=po2+1; i2++) {
        if (i2>=0 && i2 < phdim) {
          if (doneMask[i1 + i2*phdim] == 1) {
            meaval += ph[i1 + i2*phdim];
            inds ++;
          }
        }
      }
    }
  }

  // No valid neighbours found?
  if (inds==0) {
    fprintf(stderr, "valid_neighs_getmean() INTERNAL ERROR.\n");
    return -1;
  }
  return ((double)meaval)/((double)inds);
}

/*!
  @brief  Adds a point to the floodborder.
  @author Visa Korkiakoski

  @param uwd [in, out] unwrapping data structure
  @param po1 [in] point to add, 1st coordinate
  @param po2 [in] point to add, 2nd coordinate
 */
void floodborder_add(unwrapqdata_t *uwd, int po1, int po2)
{
  int newIndex = po1+po2*uwd->phdim;

  if (uwd->borderListFirst >= 0) {
    uwd->borderListPrevs[ uwd->borderListFirst ] = newIndex;
  }

  uwd->borderListNexts[newIndex] = uwd->borderListFirst;
  uwd->borderListPrevs[newIndex] = -1; // this is the first
  
  // The added item is first at the list
  uwd->borderListFirst = newIndex;
  uwd->listsz++;
} // floodborder_add



/*!
  @brief  Removes a point from the floodborder, if necessary.
  @author Visa Korkiakoski

  @param uwd [in, out] unwrapping data structure
  @param quality [in] quality array (the points where quality is 0, are considered done)
  @param po1 [in] point to add, 1st coordinate
  @param po2 [in] point to add, 2nd coordinate
 */
void floodborder_remove(unwrapqdata_t *uwd, const double *quality, int po1, int po2)
{
  int i1, i2;
  int unfinished=0;
  int remIndex = po1 + po2*uwd->phdim;

  // Is point (po1, po2) included in the borderList?
  if (uwd->borderListNexts[remIndex] < 0 &&  uwd->borderListPrevs[remIndex] < 0)
    return; // no, it isn't  (or it is just the sole remaining)


  // Point is removed, if all the neighboring points are done, or if
  // they contain no phase information.
  for (i1=po1-1; i1<=po1+1; i1++) {
    if (i1>=0 && i1 < uwd->phdim) {
      for (i2=po2-1; i2<=po2+1; i2++) {
        if (i2>=0 && i2 < uwd->phdim) {
          if (uwd->doneMask[i1 + i2*uwd->phdim] == 0 &&
              quality[i1 + i2*uwd->phdim] > 0) {
            unfinished = 1;
            break;
          }
        }
      }
    }
  }
  if (unfinished == 0) {
    int prev = uwd->borderListPrevs[remIndex];
    int next = uwd->borderListNexts[remIndex];

    if (next >= 0)
      uwd->borderListPrevs[next] = prev;

    if (prev >= 0)
      uwd->borderListNexts[prev] = next;
    else
      uwd->borderListFirst = next;

    // Mark this as non-used
    uwd->borderListPrevs[remIndex] = -1;
    uwd->borderListNexts[remIndex] = -1;
    uwd->listsz--;
  } // unfinished ==  0, point to remove
} // floodborder_remove



/*!
  @brief Finds the maximum quality neighboring points of the floodborder that are not marked by doneMask.
  @author Visa Korkiakoski

  @param uwd [in, out] unwrapping data structure
  @param quality [in] quality array, should be positive (the points where quality is 0, are considered done)
  @param maxpo1 [out] the point where quality is highest, 1st coordinate
  @param maxpo2 [out] the point where quality is highest, 2nd coordinate
  @return the index of (maxpo1, maxpo2), or -1 if no valid maximum is found
 */
int floodborder_findmaxneighbor(unwrapqdata_t *uwd, const double *quality, 
				int *maxpo1, int *maxpo2)
{
  int     curIndex = uwd->borderListFirst;
  int     po1, po2, i1, i2, maxi=-1, itcount=0;
  double  maxQuality = -1;

  // Go through all the neighbors of the borderList
  while (curIndex >= 0) {
    itcount++;
    if (itcount > uwd->listsz) {
      fprintf(stderr, "floodborder_findmaxneighbor() fail: itcount > uwd->listsz, lists broken\n");
      return -1;
    }
    
    // Compute the coordinates corresponding to the index
    po1 = curIndex % uwd->phdim;
    po2 = curIndex / uwd->phdim;


    // Enumerate neighbors
    for (i1=po1-1; i1<=po1+1; i1++) {
      if (i1>=0 && i1 < uwd->phdim) {
        for (i2=po2-1; i2<=po2+1; i2++) {
          if (i2>=0 && i2 < uwd->phdim) {
            if (uwd->doneMask[i1 + i2*uwd->phdim] == 0 &&
                quality[i1 + i2*uwd->phdim] > maxQuality) {
              maxQuality = quality[i1 + i2*uwd->phdim];
              *maxpo1 = i1;
              *maxpo2 = i2;
              maxi = i1 + i2*uwd->phdim;
            }	  
          }
        }
      }
    }
    
    // Take the next node
    curIndex = uwd->borderListNexts[curIndex];
  } // enumerate borderList nodes
  
  return maxi;
} // floodborder_findmaxneighbor
