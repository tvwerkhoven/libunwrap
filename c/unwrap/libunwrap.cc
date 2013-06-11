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
#include <queue>
#include <algorithm>

#include "libunwrap.h"
#include "debugprint.h"

using namespace std;

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399
#endif

// Unbelievable!!
#ifdef _WIN32
inline double round( double d )
{
  return floor( d + 0.5 );
}
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
 @param phx [in] width of 'ph' and 'quality'
 @param dim2 [in] height of 'ph' and 'quality'
*/
void unwrap_flood_quality_slow(double *ph, const double *quality, size_t dim1, size_t dim2) 
{
DEBUGPRINT("ph: 0x%p, quality: 0x%p, dim1: %ld, dim2: %ld\n", ph, quality, dim1, dim2);
  unwrapqdata_t uwd;
  size_t  i, i1, i2;
  ssize_t maxi;
  size_t  neighx, neighy;
  double  meaval, thestep;

  uwd.dim1    = dim1;
  uwd.dim2    = dim2;
  uwd.nel    = dim1 * dim2;
  uwd.listsz   = 0;
  uwd.unwcount = 0;
  uwd.doneMask = (int *) calloc(dim1*dim2, sizeof( *(uwd.doneMask) ));
  uwd.borderListPrevs = (ssize_t *) calloc(dim1*dim2, sizeof( *(uwd.borderListPrevs) ));
  uwd.borderListNexts = (ssize_t *) calloc(dim1*dim2, sizeof( *(uwd.borderListPrevs) ));
  uwd.borderListFirst = -1; // No border
  for (i=0; i<dim1*dim2; i++) {
    uwd.borderListPrevs[i] = -1; // point to nil
    uwd.borderListNexts[i] = -1; // point to nil
  }
  
  // First starting point
  maxi = findmax(quality, dim1 * dim2);
  uwd.doneMask[maxi] = 1;
  uwd.unwcount++;
  neighx = maxi % dim1;
  neighy = maxi / dim2;
DEBUGPRINT("found maximum at %ld, or (%ld, %ld)\n", maxi, neighx, neighy);
  
  // This point (neighbo1, neighbo2) is now unwrapped, add to the border list
  floodborder_add(&uwd, neighx, neighy);
  
  while (floodborder_findmaxneighbor(&uwd, quality, &neighx, &neighy) >= 0) {
    
    // Is it necessary to wrap the neighboring point?
    meaval = valid_neighs_getmean(neighx, neighy, ph, uwd.doneMask, dim1, dim2);
    thestep = meaval - ph[neighx + neighy * dim1];

    if (fabs(thestep) > M_PI) {
      ph[neighx + neighy*dim1] += 2.0 * M_PI * round(thestep/(2*M_PI));
    }
    
    // Pixel is now unwrapped, mark as completed
    uwd.doneMask[neighx + neighy*dim1] = 1;
    uwd.unwcount++;
    
    // The new neighboring point needs to be added to the flooding border
    floodborder_add(&uwd, neighx, neighy);

    // Check if some points need to be removed from the flooding border
    for (i1=max(size_t(1), neighx)-1; i1<min(neighx+2, dim1); i1++) {
      for (i2=max(size_t(1), neighy)-1; i2<min(neighy+2, dim2); i2++) {
        // Removes a point from the flooding border, if the point
        // has no potential neighbors
        floodborder_remove(&uwd, quality, i1, i2);
      }
    }    

  } // flooding border has valid neighbors

  // Free the resources
  free(uwd.doneMask);
  free(uwd.borderListPrevs);
  free(uwd.borderListNexts);

} // unwrap_flood_quality(double *ph, const double *quality, int phdim) 



class QualityComparison
{
  const double *const quality;
public:
  QualityComparison(const double *const qualitypo) : quality(qualitypo) {;}

  bool operator() (const int& lhs, const int&rhs) const
  {
    return (quality[lhs] < quality[rhs]);
  }
};




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
 @param dim1 [in] width of 'ph' and 'quality'
 @param dim2 [in] height of 'ph' and 'quality'
*/
void unwrap_flood_quality(double *ph, const double *quality, ssize_t dim1, ssize_t dim2) 
{
DEBUGPRINT("ph: 0x%p, quality: 0x%p, dim1: %ld, dim2: %ld\n", ph, quality, dim1, dim2);
  ssize_t     maxi;
  int         ind1, ind2, ind3, ind4, po1, po2;
  double      meaval, thestep;
  int         unwcount=0;
  int        *doneMask   = (int *)calloc(dim1*dim2, sizeof(int));
  int        *borderMask = (int *)calloc(dim1*dim2, sizeof(int));

  // priority_queue< int, vector<int>, QualityComparison> Adjoint;
  typedef priority_queue< int, vector<int>, QualityComparison> mypq_type;
  mypq_type Adjoint(quality);

  // typedef priority_queue<int,vector<int>,mycomparison> 
  // mypq_type fifth (mycomparison());
  // mypq_type sixth (mycomparison(true));



  // Adjoint.reserve(dim1*dim2);

  // Indexes to the phase map are computed this way:
  //   ind = po1 + dim1*po2;
  // Coordinates from index are computed this way:
  //   po1 =       ind % dim1
  //   po2 = floor(ind / dim1)

  // First starting point
  maxi = findmax(quality, dim1 * dim2);
  doneMask[maxi] = 1;
  // Add the neighbors to adjoint
  po1 = maxi % dim1;
  po2 = maxi / dim1;

  ind1 = (po1-1) + (po2)*dim1;
  ind2 = (po1+1) + (po2)*dim1;
  ind3 = (po1) + (po2-1)*dim1;
  ind4 = (po1) + (po2+1)*dim1;
  if (po1-1 >= 0)    {Adjoint.push(ind1);  borderMask[ind1]=1; }
  if (po1+1 < dim1)  {Adjoint.push(ind2);  borderMask[ind2]=1; }
  if (po2-1 >= 0)    {Adjoint.push(ind3);  borderMask[ind3]=1; }
  if (po2+1 < dim2)  {Adjoint.push(ind4);  borderMask[ind4]=1; }
    
  
  while (Adjoint.empty() == false) {

    // Get the phase map index with highest quality
    maxi = Adjoint.top();
    po1 = maxi % dim1;
    po2 = maxi / dim1;
    
    // Remove the added point from the border
    Adjoint.pop();
    borderMask[maxi] = 0;

    // Is it necessary to wrap it?
    meaval = valid_neighs_getmean(po1, po2, ph, doneMask, dim1, dim2);
    thestep = meaval - ph[po1 + po2*dim1];
    if (fabs(thestep) > M_PI) {
      ph[po1 + po2*dim1] += 2.0 * M_PI * round(thestep/(2*M_PI));
    }

    // printf("adding (%d,%d), meaval=%f, step=%f\n", po1, po2, meaval, thestep);
    
    // Pixel is now unwrapped, mark as completed
    doneMask[po1 + po2*dim1] = 1;
    unwcount++;
    
    // New neighbors need to be added to the adjoint
    ind1 = (po1-1) + (po2)*dim1;
    ind2 = (po1+1) + (po2)*dim1;
    ind3 = (po1) + (po2-1)*dim1;
    ind4 = (po1) + (po2+1)*dim1;
    if (po1-1 >= 0   && borderMask[ind1]==0 && doneMask[ind1]==0 && quality[ind1]>0) {
      Adjoint.push(ind1);  borderMask[ind1]=1; 
    }
    if (po1+1 < dim1 && borderMask[ind2]==0 && doneMask[ind2]==0 && quality[ind2]>0) {
      Adjoint.push(ind2);  borderMask[ind2]=1; 
    }
    if (po2-1 >= 0   && borderMask[ind3]==0 && doneMask[ind3]==0 && quality[ind3]>0) {
      Adjoint.push(ind3);  borderMask[ind3]=1; 
    }
    if (po2+1 < dim2 && borderMask[ind4]==0 && doneMask[ind4]==0 && quality[ind4]>0) {
      Adjoint.push(ind4);  borderMask[ind4]=1; 
    }

  } // flooding border has valid neighbors


  free(doneMask);
  free(borderMask);
} // unwrap_flood_quality_heap(double *ph, const double *quality, int phdim) 




/*!
 @brief Finds the index of maximum value in an array
 @author Visa Korkiakoski
 */
size_t findmax(const double *arr, size_t len)
{
  double maxval=arr[0];
  size_t i, maxi=0;
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
double valid_neighs_getmean(size_t po1, size_t po2, const double * const ph, const int * const doneMask, size_t dim1, size_t dim2)
{
  size_t inds = 0;
  double meaval = 0.0;

  // Loop over all neighbours of point (po1, po2). For all neighbours with 
  // doneMask equal to 1, sum the values.
  for (size_t i1 = max(size_t(1), po1)-1; i1<min(po1+2, dim1); i1++) {
    for (size_t i2 = max(size_t(1), po2)-1; i2<min(po2+2, dim2); i2++) {
      if (doneMask[i1 + i2*dim1] == 1) {
        meaval += ph[i1 + i2*dim1];
        inds++;
      }
    }
  }

  // No valid neighbours found?
  if (inds == 0) {
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
void floodborder_add(unwrapqdata_t *uwd, size_t pox, size_t poy)
{
  // Convert to one-dimensional index
  size_t newIndex = pox + poy * uwd->dim1;

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
void floodborder_remove(unwrapqdata_t *uwd, const double *quality, size_t pox, size_t poy)
{
  int unfinished = 0;
  
  // Convert to one-dimensional index
  size_t remIndex = pox + poy * uwd->dim1;

  // Is point (pox, poy) included in the borderList?
  if (uwd->borderListNexts[remIndex] < 0 &&  uwd->borderListPrevs[remIndex] < 0)
    return; // no, it isn't  (or it is just the sole remaining)

  // Point is removed, if all the neighboring points are done, or if
  // they contain no phase information.
  for (size_t i1 = max(size_t(1), pox)-1; i1<min(pox+2, uwd->dim1); i1++) {
    for (size_t i2 = max(size_t(1), poy)-1; i2<min(poy+2, uwd->dim2); i2++) {
      if (uwd->doneMask[i1 + i2*uwd->dim1] == 0 &&
          quality[i1 + i2*uwd->dim1] > 0) {
        unfinished = 1;
        break;
      }
    }
  }
  if (unfinished == 0) {
    ssize_t prev = uwd->borderListPrevs[remIndex];
    ssize_t next = uwd->borderListNexts[remIndex];

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
ssize_t floodborder_findmaxneighbor(unwrapqdata_t *uwd, const double *quality, 
				size_t *maxpo1, size_t *maxpo2)
{
  ssize_t curIndex = uwd->borderListFirst;
  size_t  po1, po2, itcount=0;;
  ssize_t maxi = -1;
  double  maxQuality = -1;

  // Go through all the neighbors of the borderList
  while (curIndex >= 0) {
    itcount++;
    if (itcount > uwd->listsz) {
      fprintf(stderr, "floodborder_findmaxneighbor() fail: itcount > uwd->listsz, lists broken\n");
      return -1;
    }
    
    // Compute the coordinates corresponding to the index
    po1 = curIndex % uwd->dim1;
    po2 = curIndex / uwd->dim1;

    // Enumerate neighbors
    for (size_t i1 = max(size_t(1), po1)-1; i1<min(po1+2, uwd->dim1); i1++) {
      for (size_t i2 = max(size_t(1), po2)-1; i2<min(po2+2, uwd->dim2); i2++) {
        if (uwd->doneMask[i1 + i2*uwd->dim1] == 0 &&
            quality[i1 + i2*uwd->dim1] > maxQuality) {
          maxQuality = quality[i1 + i2*uwd->dim1];
          *maxpo1 = i1;
          *maxpo2 = i2;
          maxi = i1 + i2*uwd->dim1;
        }
      }
    }
    
    // Take the next node
    curIndex = uwd->borderListNexts[curIndex];
  } // enumerate borderList nodes
  
  return maxi;
} // floodborder_findmaxneighbor
