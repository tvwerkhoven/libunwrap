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
#define BORDERLIMITFRAC_DEFAULT 0.9

// #define M_PI  3.141592653589793238462643383279502884197169399

double *g_wrapped  = NULL;
double *g_pup      = NULL;
double *g_donemask = NULL;
double  g_borderlimitfrac = 0;
long    g_count = 0;
long    g_phdim = 0;
int     g_usedonemask = 0;

/*!
 @brief Try to make a floodfill unwrap.
 @author Visa Korkiakoski
*/
void unwrapflood(int po1, int po2, int itco) {
  
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
