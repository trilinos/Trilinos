/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "az_aztec.h"

void write_vec(const char *filename, int n_equations,double *x)
 
/*  write ASCII data file:
*/

{
  FILE *out_file ;
  int i;

 
  if ( (out_file = fopen( filename, "w")) == NULL ) {
    fprintf(stderr,"Error: Cannot open file: %s\n",filename);
    return;
  }
 
  for (i=0; i< n_equations; i++)
  fprintf(out_file, "%20.15e\n", x[i]) ;

  fclose(out_file);

/* end write_vec */
}
