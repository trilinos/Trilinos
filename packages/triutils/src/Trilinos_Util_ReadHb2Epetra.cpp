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
#include "Epetra_ConfigDefs.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include <vector>
#include <algorithm>

void Trilinos_Util_ReadHb2Epetra_internal(char *data_file,
				 const Epetra_Comm  &comm, 
				 Epetra_Map *& map, 
				 Epetra_CrsMatrix *& A, 
				 Epetra_Vector *& x, 
				 Epetra_Vector *& b,
				 Epetra_Vector *&xexact,
				 bool FakeLongLong) {
  FILE *in_file ;
  int numGlobalEquations=0, N_columns=0, n_entries=0, Nrhs=0;
  char Title[73], Key[9], Rhstype[4];
  char Type[4] = "XXX";
  char Ptrfmt[17], Indfmt[17], Valfmt[21], Rhsfmt[21];
  int Ptrcrd, Indcrd, Valcrd, Rhscrd;
  
  for(int ii=0; ii<73; ++ii) Title[ii] = '\0';

  int * bindx, * pntr, * indx1, * pntr1;
  double * val, * val1, * hbx, * hbxexact, * hbb;

  hbx = 0; hbb = 0; hbxexact = 0; hbb = 0;

  if(comm.MyPID() == 0)  { 
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
      
      readHB_header(in_file, Title, Key, Type, &numGlobalEquations, &N_columns, 
		    &n_entries, &Nrhs,
		    Ptrfmt, Indfmt, Valfmt, Rhsfmt,
		    &Ptrcrd, &Indcrd, &Valcrd, &Rhscrd, Rhstype);
      fclose(in_file);

      if (Nrhs < 0 ) Nrhs = 0;

      printf("%s", "***************************************************************\n");
      printf("Matrix in file %s is %d x %d, \n",data_file, numGlobalEquations, N_columns);
      printf("with %d nonzeros with type %3s;\n", n_entries, Type);
      printf("%s", "***************************************************************\n");
      printf("Title: %72s\n",Title);
      printf("%s", "***************************************************************\n");
      /*Nrhs = 0; */
      printf("%d right-hand-side(s) available.\n",Nrhs);

      if (Type[0] != 'R') perror("Can only handle real valued matrices");
      int isym = 0;
      if (Type[1] == 'S') 
	{
	  printf("%s", "Converting symmetric matrix to nonsymmetric storage\n");
	  n_entries = 2*n_entries - N_columns;
	  isym = 1;
	}
      if (Type[2] != 'A') perror("Can only handle assembled matrices");
      if (N_columns != numGlobalEquations) perror("Matrix dimensions must be the same");
      
      /* Read the matrix information, generating the associated storage arrays  */
      printf("Reading the matrix from %s...\n",data_file);

      /* Allocate space.  Note that we add extra storage in case of zero
	 diagonals.  This is necessary for conversion to MSR format. */

      pntr   = (int    *) calloc(N_columns+1,sizeof(int)) ;
      bindx = (int    *) calloc(n_entries+N_columns+1,sizeof(int)) ;
      val   = (double *) calloc(n_entries+N_columns+1,sizeof(double)) ;
      
      readHB_mat_double(data_file, pntr, bindx, val);

      /* Translate integer arrays to zero base */
      for (int i = 0; i <= numGlobalEquations; i++) pntr[i]--;
      {for (int i = 0; i <= n_entries; i++) bindx[i]--;}

      /* If a rhs is specified in the file, read one, 
	 generating the associate storage */
      if (Nrhs > 0 && Rhstype[2] =='X')
	{
	  printf("Reading right-hand-side vector(s) from %s...\n",data_file);
	  hbb = (double *) calloc(N_columns,sizeof(double));
	  readHB_aux_double(data_file, 'F', hbb);
	  printf("Reading exact solution  vector(s) from %s...\n",data_file);
	  hbxexact = (double *) calloc(N_columns,sizeof(double));
	      readHB_aux_double(data_file, 'X', hbxexact);

	}
      else
	{
	  
	  /* Set Xexact to a random vector */

	  printf("%s", "Setting  random exact solution  vector\n");
	  hbxexact = (double *) calloc(N_columns,sizeof(double));
	  
	  for (int i=0;i<numGlobalEquations;i++)	 hbxexact[i] = 
                                              ((double)
                                               rand())/((double) RAND_MAX);
	  
	  /* Compute b to match xexact */
	  
	 hbb = (double   *) calloc(N_columns,sizeof(double)) ;
	  if (hbb == NULL) perror("Error: Not enough space to create rhs");
 

      Trilinos_Util_scscmv (isym, N_columns, N_columns, val, bindx, pntr, hbxexact, hbb);
	}

      /* Compute residual using CSC format */

      double res = Trilinos_Util_scscres(isym, numGlobalEquations, numGlobalEquations, val, bindx, pntr, 
		    hbxexact, hbb);
      printf(
	      "The residual using CSC format and exact solution is %12.4g\n",
	      res);

      
      /* Set initial guess to zero */
      
      hbx = (double   *) calloc(numGlobalEquations,sizeof(double)) ;
      
      if (hbx == NULL) 
	perror("Error: Not enough space to create guess");
      
      
      /* Set RHS to a random vector, initial guess to zero */
      {for (int i=0;i<numGlobalEquations;i++) hbx[i] = 0.0;}
      
      
      /* Allocate temporary space */
      
      pntr1 = (int   *) calloc(N_columns+1,sizeof(int)) ;
      indx1 = (int   *) calloc(n_entries+N_columns+1,sizeof(int)) ;
      val1 = (double *) calloc(n_entries+N_columns+1,sizeof(double)) ;
      
      
      /* Convert in the following way:
	 - CSC to CSR 
	 - CSR to MSR
      */
      Trilinos_Util_csrcsc(numGlobalEquations, numGlobalEquations, 0, 0, val, bindx, pntr, val1, indx1, pntr1);
      
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
  }
  comm.Broadcast(&numGlobalEquations, 1, 0);
  int nlocal = 0;
  if (comm.MyPID()==0) nlocal = numGlobalEquations;

  // Note: We avoid converting all functions called by this function (and functions
  // called by them) to support genuine long longs.  Instead we (optionally) fake
  // long long in this function by having all int based data converted to long long
  // so that all the tests/code that uses this function at least run.  Even though this
  // code will not be able to handle "large" maps, it will at least work with
  // long long. -- jhurani@txcorp.com

  if(FakeLongLong) {
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    map = new Epetra_Map((long long) numGlobalEquations, nlocal, (long long) 0, comm); // Create map with all elements on PE 0
    A = new Epetra_CrsMatrix(Copy, *map, 0); // Construct matrix on PE 0

    std::vector<long long> LL_indices;
    if (comm.MyPID()==0) {
      for (int i=0; i<numGlobalEquations; i++) {
        if(pntr1[i+1]-pntr1[i] > 0) {
          LL_indices.resize(pntr1[i+1]-pntr1[i]);
		  const int* begin = indx1+pntr1[i];
		  const int* end = begin + LL_indices.size();
		  std::copy(begin, end, LL_indices.begin());
		  A->InsertGlobalValues((long long) i, pntr1[i+1]-pntr1[i], val1+pntr1[i], &LL_indices.front());
        }
      }
    }
#else
    throw "Trilinos_Util_ReadHb2Epetra_internal: long long requested but no long long maps";
#endif
  }
  else
  {
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    map = new Epetra_Map(numGlobalEquations, nlocal, 0, comm); // Create map with all elements on PE 0
    A = new Epetra_CrsMatrix(Copy, *map, 0); // Construct matrix on PE 0

    if (comm.MyPID()==0)
      for (int i=0; i<numGlobalEquations; i++)
        A->InsertGlobalValues(i, pntr1[i+1]-pntr1[i], val1+pntr1[i], indx1+pntr1[i]);
#else
    throw "Trilinos_Util_ReadHb2Epetra_internal: int requested but no int maps";
#endif
  }

  A->FillComplete();
  
  x = new Epetra_Vector(Copy, *map, hbx);
  b = new Epetra_Vector(Copy, *map, hbb);
  xexact = new Epetra_Vector(Copy, *map, hbxexact);

  Epetra_Vector bcomp(*map);

  A->Multiply(false, *xexact, bcomp);
  double residual;
  bcomp.Norm2(&residual);
  if (comm.MyPID()==0) cout << "Norm of computed b = " << residual << endl;
  b->Norm2(&residual);
  if (comm.MyPID()==0) cout << "Norm of given b    = " << residual << endl;
  bcomp.Update(-1.0, *b, 1.0);
  bcomp.Norm2(&residual);
  if (comm.MyPID()==0) cout << "Norm of difference between computed b and given b for xexact = " << residual << endl;
  
  /* Release unneeded space */
  if (comm.MyPID()==0) {
    if (hbb!=0) free((void *) hbb);
    if (hbx!=0) free((void *) hbx);
    if (hbxexact!=0) free((void *) hbxexact);
    free((void *) val);
    free((void *) bindx);
    free((void *) val1);
    free((void *) indx1);
    free((void *) pntr1);
    free((void *) pntr);
  }
  return;
}

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES

void Trilinos_Util_ReadHb2Epetra(char *data_file,
				 const Epetra_Comm  &comm, 
				 Epetra_Map *& map, 
				 Epetra_CrsMatrix *& A, 
				 Epetra_Vector *& x, 
				 Epetra_Vector *& b,
				 Epetra_Vector *&xexact) {
  Trilinos_Util_ReadHb2Epetra_internal(data_file, comm, map, A, x, b, xexact, false);
}

#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES

void Trilinos_Util_ReadHb2Epetra64(char *data_file,
				 const Epetra_Comm  &comm, 
				 Epetra_Map *& map, 
				 Epetra_CrsMatrix *& A, 
				 Epetra_Vector *& x, 
				 Epetra_Vector *& b,
				 Epetra_Vector *&xexact) {
  Trilinos_Util_ReadHb2Epetra_internal(data_file, comm, map, A, x, b, xexact, true);
}

#endif
