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


/*!
  @brief  Tests the unwrapping library.

  As arguments takes the names of wrapped phase, quality array and a
  result name. The arrays are assumed to be square, and they are
  assumed to be double type. Their dimension is then inferred from the
  sizes the take at disk.

  Then, the arrays are read, unwrapped, and then saved to disk.
 */
int main(int argc, char *argv[])
{
  char   *fname, *qualname, *resname;
  FILE   *stream;
  long    fsize, readed;
  double *wrapped, *quality;
  int     phdim;

  if (argc != 4) {
    printf("Usage: test_unwrap <phname> <qualityname> <resname>.\n");
    return -1;
  }
    
  fname   = argv[1];
  qualname = argv[2];
  resname = argv[3];

  // Read the phase to be unwrapped
  printf("Reading file %s...\n", fname);

  stream = fopen(fname, "rb");
  if (stream == NULL) {
    printf("File %s not found.\n", fname);
    return -1;
  }
  
  fseek(stream, 0, SEEK_END); // seek to end of file
  fsize = ftell(stream); // get current file pointer
  fseek(stream, 0, SEEK_SET); // seek back to beginning of file

  printf("File size %ld, dimension %f\n", fsize, sqrt(fsize/8.));
  phdim = sqrt(fsize/8.);

  // Check if the file size makes sense
  if (fabs(sqrt(fsize/8.) - round(sqrt(fsize/8.))) > 1e-9) {
    printf("The phase array needs to be a square matrix in double format.\n");
    fclose(stream);
    return -1;
  }

  wrapped  = calloc(fsize, sizeof(char));
  quality  = calloc(fsize, sizeof(char));

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


  printf("Unwrapping...");

  // Do the recursive unwrap
  unwrap_quality(wrapped, quality, phdim);
  
  // Save result
  stream = fopen(resname, "wb");
  fwrite(wrapped, sizeof(char), fsize, stream);
  fclose(stream);


  free(wrapped);
  free(quality);

  return 0;
} // main
