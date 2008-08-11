// @HEADER
// ***********************************************************************
// 
//                 TriUtils: Trilinos Utilities Package
//                 Copyright (2001) Sandia Corporation
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

#include "Trilinos_Util.h"

void Trilinos_Util_read_coo(char *data_file, int MyPID,
	      int *N_global, int *n_nonzeros,
	      double **val, int **bindx,
	      double **x, double **b, double **xexact)
#undef DEBUG 
     /*  read ASCII data file:
	 line 1: N_global, number of entries (%d,%d)
	 line 2-...: i,j,real (%d, %d, %f)
     */

{
  FILE *data ;


  int i, n_entries, N_columns;
  int ii, jj ;
  int kk = 0;
  int max_ii = 0, max_jj = 0;
  double value;
  int *pntr, *indx1, *pntr1;
  double *val1;

  int MAXBLOCKSIZE = 1;
  
  if(MyPID == 0) 
    { 

 
      data = fopen(data_file,"r") ;
 
      fscanf(data, "%d %d %d", N_global, &N_columns, &n_entries) ;
      if (N_columns != *N_global)
	perror("Matrix dimensions must be the same");
      printf("Reading from file: %s\n",data_file);
      printf("Number of equations = %d\n",*N_global);
      printf("Number of entries   = %d\n",n_entries);
    

      *bindx = (int   *) calloc(n_entries+1,sizeof(int)) ;
      *val = (double *) calloc(n_entries+1,sizeof(double)) ;

      pntr1 = (int   *) calloc(n_entries+1,sizeof(int)) ;
      indx1 = (int   *) calloc(n_entries+1,sizeof(int)) ;
      val1 = (double *) calloc(n_entries+1,sizeof(double)) ;
  
      pntr = (int   *) calloc(n_entries+1,sizeof(int)) ;

      if ((pntr) == NULL) 
	perror("Error: Not enough space to create matrix");

      while(!feof(data))
	{
	  fscanf(data, "%d %d %lf", &ii, &jj, &value) ;
	  max_ii = Trilinos_Util_max(max_ii,ii);
	  max_jj = Trilinos_Util_max(max_jj,jj);
#ifdef DEBUG
	  printf("Entry %d, %d = %lf.\n",ii,jj,value);
#endif
	  (*bindx)[kk] = ii-1;
	  pntr[kk] = jj-1;
	  (*val)[kk] = value;
	  kk++;
	}
      *n_nonzeros = kk-1;
      *N_global = max_ii;
      if (max_ii != max_jj) perror("Error: Number of rows and columns not equal");

      printf("Number of nonzeros = %d\n",*n_nonzeros);

      /* Convert real part in the following way:
	 - Convert COO to CSR
	 - CSR to CSC
	 - CSC to CSR (columns are now in ascending order)
	 - CSR to MSR
      */
      Trilinos_Util_coocsr(*N_global, *n_nonzeros, *val, *bindx, pntr, val1, indx1, pntr1);

      Trilinos_Util_csrcsc(*N_global, *N_global, 0, 0, val1, indx1, pntr1,*val, *bindx, pntr);

      Trilinos_Util_csrcsc(*N_global, *N_global, 0, 0, *val,*bindx, pntr, val1, indx1, pntr1);

      Trilinos_Util_csrmsr(*N_global, val1, indx1, pntr1, *val,*bindx, *val,*bindx);


      *b = (double   *) calloc((*N_global)*MAXBLOCKSIZE,sizeof(double)) ;
      *x = (double   *) calloc((*N_global)*MAXBLOCKSIZE,sizeof(double)) ;

      if ((*x) == NULL) 
	perror("Error: Not enough space to create matrix");

  
      /* Set RHS to a random vector, initial guess to zero */
      for (i=0;i<*N_global;i++)
	{
	  (*b)[i] = ((double) rand())/ ((double) RAND_MAX);
	  (*x)[i] = 0.0;
	}
  }
   
      /* Release unneeded space */

  free((void *) pntr);
  free((void *) val1);
  free((void *) indx1);
  free((void *) pntr1);
  
  
  /* end read_coo */
}
