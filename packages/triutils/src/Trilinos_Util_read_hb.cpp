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
#include "iohb.h"

void Trilinos_Util_read_hb(char *data_file, int MyPID,
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
  char Type[4] = "XXX";
  char Ptrfmt[17], Indfmt[17], Valfmt[21], Rhsfmt[21];
  int Ptrcrd, Indcrd, Valcrd, Rhscrd;

  int i, n_entries=0, N_columns=0, Nrhs=0;
  int isym;
  double res;
  int *pntr, *indx1, *pntr1;
  double *val1;

  int MAXBLOCKSIZE = 1;
  
  if(MyPID == 0)  { 
      
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

      printf("%s", "***************************************************************\n");
      printf("Matrix in file %s is %d x %d, \n",data_file, *N_global, N_columns);
      printf("with %d nonzeros with type %3s;\n", n_entries, Type);
      printf("%s", "***************************************************************\n");
      printf("Title: %72s\n",Title);
      printf("%s", "***************************************************************\n");
      /*Nrhs = 0; */
      printf("%d right-hand-side(s) available.\n",Nrhs);

      if (Type[0] != 'R') perror("Can only handle real valued matrices");
      isym = 0;
      if (Type[1] == 'S') 
	{
	  printf("%s", "Converting symmetric matrix to nonsymmetric storage\n");
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

      for (i = 0; i <= *N_global; i++) pntr[i]--;
      for (i = 0; i <= n_entries; i++) (*bindx)[i]--;

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

	  printf("%s", "Setting  random exact solution  vector\n");
	  *xexact = (double *) calloc(N_columns,sizeof(double));
	  
	  for (i=0;i<*N_global;i++)	 (*xexact)[i] =
                                       ((double)rand())/((double) RAND_MAX);
	  
	  /* Compute b to match xexact */
	  
	  *b = (double   *) calloc(N_columns*MAXBLOCKSIZE,sizeof(double)) ;
	  if ((*b) == NULL) perror("Error: Not enough space to create rhs");
 

      Trilinos_Util_scscmv (isym, N_columns, N_columns, (*val), (*bindx), pntr, (*xexact), (*b));
	}

      /* Compute residual using CSC format */

      res = Trilinos_Util_scscres(isym, *N_global, *N_global, (*val), (*bindx), pntr, 
		    (*xexact), (*b));
      printf(
	      "The residual using CSC format and exact solution is %12.4g\n",
	      res);

      
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
      Trilinos_Util_csrcsc(*N_global, *N_global, 0, 0, *val, *bindx, pntr, val1, indx1, pntr1);
      
      if (Type[1] == 'S') 
	{
	  int *indu, *iwk;
	  int ierr;
	  indu = new int[N_columns];
	  iwk = new int[N_columns+1];
	  ierr = Trilinos_Util_ssrcsr(3, 1, N_columns, val1, indx1, pntr1, n_entries,
	  		  val1, indx1, pntr1, indu, iwk);
	  delete [] indu;
	  delete [] iwk;
	  if (ierr !=0 ) 
	    {
	    printf(" Error in converting from symmetric form\n  IERR = %d\n",ierr);
	    abort();
	    }
	}

      Trilinos_Util_csrmsr(*N_global, val1, indx1, pntr1, *val, *bindx, *val, *bindx);
      
      /* Recompute number of nonzeros in case there were zero diagonals */
      
      *n_nonzeros = (*bindx)[*N_global] - 1; /* Still in Fortran mode so -2 */
            

      printf("The residual using MSR format and exact solution is %12.4g\n",
	      Trilinos_Util_smsrres (*N_global, *N_global, (*val), (*bindx), 
		       (*xexact), (*xexact), (*b)));

      /* Release unneeded space */

      free((void *) val1);
      free((void *) indx1);
      free((void *) pntr1);
      free((void *) pntr);
    }
  
  /* end read_hb */
}
//
//----------------------------------------------------------------------
// Alternate version of Trilinos_Util_read_hb.cpp
// This version only reads in the matrix from an HB
// data set. It ignores any right-hand side, initial guess,
// and exact solution information that may be in the file.
//---------------------------------------------------------------------- 
//

void Trilinos_Util_read_hb(char *data_file, int MyPID,
	      int *N_global, int *n_nonzeros, 
	      double **val, int **bindx)
#undef DEBUG 
     /*  read ASCII data file:
	 line 1: N_global, number of entries (%d,%d)
	 line 2-...: i,j,real (%d, %d, %f)
     */

{
  FILE *in_file ;
  char Title[73], Key[9], Rhstype[4];
  char Type[4] = "XXX";
  char Ptrfmt[17], Indfmt[17], Valfmt[21], Rhsfmt[21];
  int Ptrcrd, Indcrd, Valcrd, Rhscrd;

  int i, n_entries=0, N_columns=0, Nrhs=0;
  int isym;
  int *pntr, *indx1, *pntr1;
  double *val1;

  
  if(MyPID == 0)  { 
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

      printf("%s", "***************************************************************\n");
      printf("Matrix in file %s is %d x %d, \n",data_file, *N_global, N_columns);
      printf("with %d nonzeros with type %3s;\n", n_entries, Type);
      printf("%s", "***************************************************************\n");
      printf("Title: %72s\n",Title);
      printf("%s", "***************************************************************\n");
      /*Nrhs = 0; */
      //printf("%d right-hand-side(s) available.\n",Nrhs);

      if (Type[0] != 'R') perror("Can only handle real valued matrices");
      isym = 0;
      if (Type[1] == 'S') 
	{
	  printf("%s", "Converting symmetric matrix to nonsymmetric storage\n");
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

      for (i = 0; i <= *N_global; i++) pntr[i]--;
      for (i = 0; i <= n_entries; i++) (*bindx)[i]--;
   
      /* Allocate temporary space */
      
      pntr1 = (int   *) calloc(N_columns+1,sizeof(int)) ;
      indx1 = (int   *) calloc(n_entries+N_columns+1,sizeof(int)) ;
      val1 = (double *) calloc(n_entries+N_columns+1,sizeof(double)) ;
      
      
      /* Convert in the following way:
	 - CSC to CSR 
	 - CSR to MSR
      */
      Trilinos_Util_csrcsc(*N_global, *N_global, 0, 0, *val, *bindx, pntr, val1, indx1, pntr1);
      
      if (Type[1] == 'S') 
	{
	  int *indu, *iwk;
	  int ierr;
	  indu = new int[N_columns];
	  iwk = new int[N_columns+1];
	  ierr = Trilinos_Util_ssrcsr(3, 1, N_columns, val1, indx1, pntr1, n_entries,
	  		  val1, indx1, pntr1, indu, iwk);
	  delete [] indu;
	  delete [] iwk;
	  if (ierr !=0 ) 
	    {
	    printf(" Error in converting from symmetric form\n  IERR = %d\n",ierr);
	    abort();
	    }
	}

      Trilinos_Util_csrmsr(*N_global, val1, indx1, pntr1, *val, *bindx, *val, *bindx);
      
      /* Recompute number of nonzeros in case there were zero diagonals */
      
      *n_nonzeros = (*bindx)[*N_global] - 1; /* Still in Fortran mode so -2 */
            
      /* Release unneeded space */

      free((void *) val1);
      free((void *) indx1);
      free((void *) pntr1);
      free((void *) pntr);
    }
  
  /* end read_hb */
}
