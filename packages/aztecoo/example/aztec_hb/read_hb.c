/*
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
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
//@HEADER
*/

#include <stdlib.h>
#include <stdio.h>
#include "az_aztec.h"
#include "iohb.h"
#include "prototypes.h"

void read_hb(char *data_file, int *proc_config,
	      int *N_global, int *n_nonzeros, 
	      double **val, int **bindx,
	      double **x, double **b, double **xexact)
#undef DEBUG 
     /*  read ASCII data file:
	 line 1: N_global, number of entries (%d,%d)
	 line 2-...: i,j,real (%d, %d, %f)
     */

{
  FILE *in_file ;
  char Title[73], Key[9], Rhstype[4];
  char Type[4];
  char Ptrfmt[17], Indfmt[17], Valfmt[21], Rhsfmt[21];
  int Ptrcrd, Indcrd, Valcrd, Rhscrd;

  int i, n_entries, N_columns, Nrhs;
  int ii, jj ;
  int kk = 0;
  int isym;
  int ione = 1;
  double value, res;
  double *cnt;
  int *pntr, *indx1, *pntr1;
  double *val1;

  int MAXBLOCKSIZE = 1;
  
  if(proc_config[AZ_node] == 0) 
    { 
      in_file = fopen( data_file, "r");
      if (in_file == NULL)
	{
	  printf("Error: Cannot open file: %s\n",data_file);
	  exit(1);
	}

      /* Get information about the array stored in the file specified in the  */
      /* argument list:                                                       */

      printf("Reading matrix info from %s...\n",data_file);
      
      in_file = fopen( data_file, "r");
      if (in_file == NULL)
	{
	  printf("Error: Cannot open file: %s\n",data_file);
	  exit(1);
	}

      readHB_header(in_file, Title, Key, Type, N_global, &N_columns, 
		    &n_entries, &Nrhs,
		    Ptrfmt, Indfmt, Valfmt, Rhsfmt,
		    &Ptrcrd, &Indcrd, &Valcrd, &Rhscrd, Rhstype);
      fclose(in_file);

      if (Nrhs < 0 ) Nrhs = 0;

      printf("***************************************************************\n");
      printf("Matrix in file %s is %d x %d, \n",data_file, *N_global, N_columns);
      printf("with %d nonzeros with type %3s;\n", n_entries, Type);
      printf("***************************************************************\n");
      printf("Title: %72s\n",Title);
      printf("***************************************************************\n");
      /*Nrhs = 0; */
      printf("%d right-hand-side(s) available.\n",Nrhs);

      if (Type[0] != 'R') perror("Can only handle real valued matrices");
      isym = 0;
      if (Type[1] == 'S') 
	{
	  printf("Converting symmetric matrix to nonsymmetric storage\n");
	  n_entries = 2*n_entries - N_columns;
	  isym = 1;
	}
      if (Type[2] != 'A') perror("Can only handle assembled matrices");
      if (N_columns != *N_global) perror("Matrix dimensions must be the same");
      *n_nonzeros = n_entries;
      
      /* Read the matrix information, generating the associated storage arrays  */
      printf("Reading the matrix from %s...\n",data_file);

      /* Allocate space.  Note that we add extra storage in case of zero
	 diagonals.  This is necessary for conversion to MSR format. */

      pntr   = (int    *) calloc(N_columns+1,sizeof(int)) ;
      *bindx = (int    *) calloc(n_entries+N_columns+1,sizeof(int)) ;
      *val   = (double *) calloc(n_entries+N_columns+1,sizeof(double)) ;
      
      readHB_mat_double(data_file, pntr, *bindx, *val);

      /* If a rhs is specified in the file, read one, 
	 generating the associate storage */
      if (Nrhs > 0 && Rhstype[2] =='X')
	{
	  printf("Reading right-hand-side vector(s) from %s...\n",data_file);
	  *b = (double *) calloc(N_columns,sizeof(double));
	  readHB_aux_double(data_file, 'F', (*b));
	  printf("Reading exact solution  vector(s) from %s...\n",data_file);
	  *xexact = (double *) calloc(N_columns,sizeof(double));
	      readHB_aux_double(data_file, 'X', (*xexact));

	}
      else
	{
	  
	  /* Set Xexact to a random vector */

	  printf("Setting  random exact solution  vector\n");
	  *xexact = (double *) calloc(N_columns,sizeof(double));
	  
	  for (i=0;i<*N_global;i++)	 (*xexact)[i] = drand48();
	  
	  /* Compute b to match xexact */
	  
	  *b = (double   *) calloc(N_columns*MAXBLOCKSIZE,sizeof(double)) ;
	  if ((*b) == NULL) perror("Error: Not enough space to create rhs");
 
	  for (i = 0; i <= *N_global; i++) pntr[i]--;
	  for (i = 0; i <= n_entries; i++) (*bindx)[i]--;

      scscmv (isym, N_columns, N_columns, (*val), (*bindx), pntr, (*xexact), (*b));

	  for (i = 0; i <= *N_global; i++) pntr[i]++;
	  for (i = 0; i <= n_entries; i++) (*bindx)[i]++;
	}

      /* Compute residual using CSC format */

      for (i = 0; i <= *N_global; i++) pntr[i]--;
      for (i = 0; i <= n_entries; i++) (*bindx)[i]--;
      res = scscres(isym, *N_global, *N_global, (*val), (*bindx), pntr, 
		    (*xexact), (*b));
      printf(
	      "The residual using CSC format and exact solution is %12.4g\n",
	      res);
      for (i = 0; i <= *N_global; i++) pntr[i]++;
      for (i = 0; i <= n_entries; i++) (*bindx)[i]++;

      
      /* Set initial guess to zero */
      
      *x = (double   *) calloc((*N_global)*MAXBLOCKSIZE,sizeof(double)) ;
      
      if ((*x) == NULL) 
	perror("Error: Not enough space to create guess");
      
      
      /* Set RHS to a random vector, initial guess to zero */
      for (i=0;i<*N_global;i++) (*x)[i] = 0.0;
      
      
      /* Allocate temporary space */
      
      pntr1 = (int   *) calloc(N_columns+1,sizeof(int)) ;
      indx1 = (int   *) calloc(n_entries+N_columns+1,sizeof(int)) ;
      val1 = (double *) calloc(n_entries+N_columns+1,sizeof(double)) ;
      
      
      /* Convert in the following way:
	 - CSC to CSR 
	 - CSR to MSR
      */
      csrcsc_(N_global,&ione,&ione,
	      *val,*bindx,pntr,
	      val1,indx1,pntr1);
      
      if (Type[1] == 'S') 
	{
	  int *indu;
	  int ierr;
	  indu = indx1+n_entries; /* Use end of bindx for workspace */
	  ssrcsr_(&N_columns, val1, indx1, pntr1, &n_entries,
				  val1, indx1, pntr1, indu, &ierr);
	  if (ierr !=0 ) 
	    {
	    printf(" Error in converting from symmetric form\n  IERR = %d\n",ierr);
	    abort();
	    }
	}

      csrmsr_(N_global,val1,indx1,pntr1,
	      *val,*bindx,
	      *val,*bindx);
      
      /* Recompute number of nonzeros in case there were zero diagonals */
      
      *n_nonzeros = (*bindx)[*N_global] - 2; /* Still in Fortran mode so -2 */
      
      /* Finally, convert bindx vectors to zero base */
      
      for (i=0;i<*n_nonzeros+1;i++) (*bindx)[i] -= 1;
      

      printf("The residual using MSR format and exact solution is %12.4g\n",
	      smsrres (*N_global, *N_global, (*val), (*bindx), 
		       (*xexact), (*xexact), (*b)));

      /* Release unneeded space */

      free((void *) val1);
      free((void *) indx1);
      free((void *) pntr1);
      free((void *) pntr);
    }
  
  /* end read_hb */
}
