/*
// @HEADER
// ***********************************************************************
// 
//                Pliris: Parallel Dense Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#include "mpi.h"
#include "stdio.h"
#include "defines.h"
#include "x_factor.h"
#include "x_solve.h"
#include "permute.h"
#include "distribute.h"
#include "clean_code.h"


main (int argc, char **argv) {
  int nproc_perrow,my_rhs,nrows,ncols;
  int first_row,first_col,proc_row,proc_col;
  int nrhs_all;
  double matrix[66];
  double rhs[8];
  int permute[8];
  double time_factor;
  int i,n;

  MPI_Init(&argc,&argv);
  n = 8;
  nproc_perrow = 1;
  nrhs_all = 1;

  distmat_(&nproc_perrow,&n,&nrhs_all,&nrows,&ncols,
           &first_row,&first_col,&my_rhs,&proc_row,&proc_col);

  i = 0;
  matrix[i++] = 1;
  matrix[i++] = -0.0432915;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = -0.0774953;
  matrix[i++] = 0;
  matrix[i++] = -0.228296;
  matrix[i++] = -0.285212;
  matrix[i++] = -0.015859;
  matrix[i++] = 1;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = -0.142507;
  matrix[i++] = -0.0836315;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 1;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 1;
  matrix[i++] = -0.00273798;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = -0.00246323;
  matrix[i++] = -0.00246323;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = -0.00273798;
  matrix[i++] = 1;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 1;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = -0.0836315;
  matrix[i++] = -0.142507;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 1;
  matrix[i++] = -0.015859;
  matrix[i++] = -0.285212;
  matrix[i++] = -0.228296;
  matrix[i++] = 0;
  matrix[i++] = -0.0774953;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = -0.0432915;
  matrix[i++] = 1;

  dfactor_(matrix,&n,&nproc_perrow,permute,&time_factor);
    

  i = 0;
  rhs[i++] = 0.288703;
  rhs[i++] = 0.268113;
  rhs[i++] = 0.553465;
  rhs[i++] = 0.536645;
  rhs[i++] = 0.536645;
  rhs[i++] = 0.553465;
  rhs[i++] = 0.268113;
  rhs[i++] = 0.288703;

  printf("Rhs:\n");
  for (i = 0; i < 8; i++) printf("  %d %g\n",i,rhs[i]);

  dpermute_(matrix, permute);

  dsolve_(matrix,permute,rhs,&nrhs_all);

  printf("Soln:\n");
  for (i = 0; i < 8; i++) printf("  %d %g\n",i,rhs[i]);

  i = 0;
  rhs[i++] = 0.288703;
  rhs[i++] = 0.268113;
  rhs[i++] = 0.553465;
  rhs[i++] = 0.536645;
  rhs[i++] = 0.536645;
  rhs[i++] = 0.553465;
  rhs[i++] = 0.268113;
  rhs[i++] = 0.288703;

  dsolve_(matrix,permute,rhs,&nrhs_all);

  printf("Soln: #2 \n");
  for (i = 0; i < 8; i++) printf("  %d %g\n",i,rhs[i]);

  cleancode_ ( );

  distmat_(&nproc_perrow,&n,&nrhs_all,&nrows,&ncols,
           &first_row,&first_col,&my_rhs,&proc_row,&proc_col);
  i = 0;
  matrix[i++] = 1;
  matrix[i++] = -0.0432915;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = -0.0774953;
  matrix[i++] = 0;
  matrix[i++] = -0.228296;
  matrix[i++] = -0.285212;
  matrix[i++] = -0.015859;
  matrix[i++] = 1;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = -0.142507;
  matrix[i++] = -0.0836315;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 1;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 1;
  matrix[i++] = -0.00273798;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = -0.00246323;
  matrix[i++] = -0.00246323;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = -0.00273798;
  matrix[i++] = 1;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 1;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = -0.0836315;
  matrix[i++] = -0.142507;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = 1;
  matrix[i++] = -0.015859;
  matrix[i++] = -0.285212;
  matrix[i++] = -0.228296;
  matrix[i++] = 0;
  matrix[i++] = -0.0774953;
  matrix[i++] = 0;
  matrix[i++] = 0;
  matrix[i++] = -0.0432915;
  matrix[i++] = 1;
  matrix[i++] = 0;
  matrix[i++] = 0;

  dfactor_(matrix,&n,&nproc_perrow,permute,&time_factor);
    

  i = 0;
  rhs[i++] = 0.288703;
  rhs[i++] = 0.268113;
  rhs[i++] = 0.553465;
  rhs[i++] = 0.536645;
  rhs[i++] = 0.536645;
  rhs[i++] = 0.553465;
  rhs[i++] = 0.268113;
  rhs[i++] = 0.288703;

  printf("Rhs:\n");
  for (i = 0; i < 8; i++) printf("  %d %g\n",i,rhs[i]);

  dpermute_(matrix, permute);

  dsolve_(matrix,permute,rhs,&nrhs_all);

  printf("Soln: problem redone\n");
  for (i = 0; i < 8; i++) printf("  %d %g\n",i,rhs[i]);

  MPI_Finalize ( ) ;

  exit(0);

}
