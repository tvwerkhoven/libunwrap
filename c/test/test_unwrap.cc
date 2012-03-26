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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "libunwrap.h"

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399
#endif

/*!
 @brief  Tests the unwrapping library.

 To test the unwrap library, either takes files from the command line to 
 unwrap, or generates random data on the fly to test. Note that in the latter 
 case the results will not make sense, and are merely for robustness testing.
 
 The supplied files are assumed to be a binary blob of 64 bit double data, 
 representing a square array. Their dimension is then inferred from the
 sizes the take at disk.

 The test input is unwrapped, and stored to disk.
 
 Run without arguments for help.
 */
int main(int argc, char *argv[])
{
  char   *fname, *qualname, *resname;
  FILE   *stream;
  long   fsize=0, readed=0;
  double *wrapped, *quality;
  size_t phdim=0;
  bool gendata = false;


  // 2 arguments: generate data ourselves
  if (argc == 2) {
    gendata = true;
    resname  = argv[1];
    printf("Running %s, auto-generating data, storing to %s.\n", argv[0], resname);
  }
  // 4 arguments: we get input data
  else if (argc == 4) {
    fname    = argv[1];
    qualname = argv[2];
    resname  = argv[3];
    printf("Running %s, phase: %s, quality: %s, storing to %s.\n", argv[0], fname, qualname, resname);
  }
  else {
    printf("Usage:\n%s <phname> <qualityname> <resname>.\nOR\n", argv[0]);
    printf("%s <resname>.\n", argv[0]);
    return -1;
  }
  
  if (gendata) {
    phdim = 257;

    printf("Generating test data (phdim %ld)...\nOR\n", phdim);
    
    // Allocate memory for quality and phase
    wrapped  = (double *) malloc(phdim * phdim * sizeof(*wrapped));
    quality  = (double *) malloc(phdim * phdim * sizeof(*quality));
    
    // Randomly generate phase and quality
//    for (size_t i0=0; i0<phdim; i0++) {
//      for (size_t i1=0; i1<phdim; i1++) {
//        wrapped[i0 + i1*phdim] = drand48();
//        quality[i0 + i1*phdim] = drand48();
//      }
//    }
    // Smooth test data (like in Python)
    for (size_t i0=0; i0<phdim; i0++) {
      for (size_t i1=0; i1<phdim; i1++) {
        wrapped[i0 + i1*phdim] = sin(M_PI*2.0*i0/phdim)*3.0 + cos(M_PI*8.0*i0/phdim)*4.0;
        quality[i0 + i1*phdim] = i0 + i1;
      }
    }
  }
  else {
    // Read the phase to be unwrapped
    printf("Reading file %s...\n", fname);
    
    stream = fopen(fname, "rb");
    if (stream == NULL) {
      printf("Could not open file %s.\n", fname);
      return -1;
    }
    
    fseek(stream, 0, SEEK_END); // seek to end of file
    fsize = ftell(stream); // get current file pointer
    fseek(stream, 0, SEEK_SET); // seek back to beginning of file
    
    // Check if the file size makes sense
    if (fabs(sqrt(fsize/8.) - round(sqrt(fsize/8.))) > 1e-9) {
      printf("The phase array needs to be a square matrix in double format.\n");
      fclose(stream);
      return -1;
    }
    
    printf("File size %ld, dimension %g\n", fsize, sqrt(fsize/8.));
    phdim = (size_t) sqrt(fsize/8.);
    
    // wrapped and quality are of type double, but data is read as char
    wrapped  = (double *) calloc(fsize, sizeof(char));
    quality  = (double *) calloc(fsize, sizeof(char));
    
    // Read in data
    fread(wrapped, sizeof(char), fsize, stream);
    fclose(stream);
    
    stream = fopen(qualname, "rb");  
    if (stream == NULL) {
      printf("File %s not found.\n", qualname);
      return -1;
    }
    readed = fread(quality, sizeof(char), fsize, stream);
    if (readed != fsize) {
      printf("Error reading %s (%ld/%ld)\n", qualname, readed, fsize);
    }
    fclose(stream);
  }
  
  printf("Unwrapping...");

  // Do the recursive unwrap
  unwrap_flood_quality(wrapped, quality, phdim, phdim);
  
  // Save result
  stream = fopen(resname, "wb");
  fwrite(wrapped, sizeof(*wrapped), phdim*phdim, stream);
  fclose(stream);

  free(wrapped);
  free(quality);

  return 0;
} // main
