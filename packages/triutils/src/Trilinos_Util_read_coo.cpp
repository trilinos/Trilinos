// @HEADER
// ***********************************************************************
// 
//                 TriUtils: Trilinos Utilities Package
//                 Copyright (2011) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
 
      int count = fscanf(data, "%d %d %d", N_global, &N_columns, &n_entries) ;
      assert(count > 0);
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
	  int count = fscanf(data, "%d %d %lf", &ii, &jj, &value) ;
    assert(count > 0);
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
