/*@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#include <stdlib.h>
#include <stdio.h>

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
