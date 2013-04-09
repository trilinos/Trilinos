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

template<typename int_type>
void Trilinos_Util_ReadHpc2Epetra_internal(char *data_file,
				 const Epetra_Comm  &comm, 
				 Epetra_Map *& map, 
				 Epetra_CrsMatrix *& A, 
				 Epetra_Vector *& x, 
				 Epetra_Vector *& b,
				 Epetra_Vector *&xexact) {


  FILE *in_file ;

  int l;
  int * lp = &l;
  double v;
  double * vp = &v;
#ifdef DEBUG
  bool debug = true;
#else
  bool debug = false;
#endif

  int size = comm.NumProc();
  int rank = comm.MyPID();
  printf("Reading matrix info from %s...\n",data_file);
  
  in_file = fopen( data_file, "r");
  if (in_file == NULL)
    {
      printf("Error: Cannot open file: %s\n",data_file);
      exit(1);
    }
  int_type numGlobalEquations, total_nnz;
  int cnt;
  if(sizeof(int) == sizeof(int_type)) {
    cnt = fscanf(in_file,"%d",&numGlobalEquations);
    assert(cnt > 0);
    cnt = fscanf(in_file,"%d",&total_nnz);
    assert(cnt > 0);
  }
  else if(sizeof(long long) == sizeof(int_type)) {
    cnt = fscanf(in_file,"%lld",&numGlobalEquations);
    assert(cnt > 0);
    cnt = fscanf(in_file,"%lld",&total_nnz);
    assert(cnt > 0);
  }
  else
    assert(false);

  map = new Epetra_Map(numGlobalEquations, (int_type) 0, comm); // Create map with uniform distribution
  
  A = new Epetra_CrsMatrix(Copy, *map, 0); // Construct matrix

  x = new Epetra_Vector(*map);
  b = new Epetra_Vector(*map);
  xexact = new Epetra_Vector(*map);
  int numMyEquations = map->NumMyPoints();

  // Allocate arrays that are of length numMyEquations

  // Find max nnz per row for this processor

  int max_nnz = 0;

  for (int i=0; i<numGlobalEquations; i++) {
      cnt = fscanf(in_file, "%d",lp); /* row #, nnz in row */
      assert(cnt > 0);
      if (map->MyGID(i)) max_nnz = EPETRA_MAX(max_nnz,l);
    }


  // Allocate arrays that are of length local_nnz
  double * list_of_vals = new double[max_nnz];
  int_type *list_of_inds = new int_type   [max_nnz];

  {for (int_type i=0; i<numGlobalEquations; i++)
    {
      int cur_nnz;
      cnt = fscanf(in_file, "%d",&cur_nnz);
      assert(cnt > 0);
      if (map->MyGID(i)) // See if nnz for row should be added
	{
	  if (debug) cout << "Process "<<rank
			  <<" of "<<size<<" getting row "<<i<<endl;
	  int nnz_kept = 0;
	  for (int j=0; j<cur_nnz; j++) 
	    {
	      cnt = fscanf(in_file, "%lf %d",vp,lp);
        assert(cnt > 0);
	      if (v!=0.0) {
		list_of_vals[nnz_kept] = v;
		list_of_inds[nnz_kept] = l;
		nnz_kept++;
	      }
	    }
	  A->InsertGlobalValues(i, nnz_kept, list_of_vals, list_of_inds);
	}
      else
	for (int j=0; j<cur_nnz; j++) {
    cnt = fscanf(in_file, "%lf %d",vp,lp); // otherwise read and discard
    assert(cnt > 0);
  }
    }}

  double xt, bt, xxt;
  {for (int_type i=0; i<numGlobalEquations; i++) 
    {
      if (map->MyGID(i)) // See if entry should be added
	{
	  if (debug) cout << "Process "<<rank<<" of "
                       <<size<<" getting RHS "<<i<<endl;
	  cnt = fscanf(in_file, "%lf %lf %lf",&xt, &bt, &xxt);
    assert(cnt > 0);
	  int cur_local_row = map->LID(i);
	  (*x)[cur_local_row] = xt;
	  (*b)[cur_local_row] = bt;
	  (*xexact)[cur_local_row] = xxt;
	}
      else
      {
	  cnt = fscanf(in_file, "%lf %lf %lf",vp, vp, vp); // or thrown away
    assert(cnt > 0);
    }
    }}

  fclose(in_file);

  
  if (debug)
    cout << "Process "<<rank<<" of "<<size<<" has "<<numMyEquations
	 << " rows. Min global row "<< map->MinMyGID64()
	 <<" Max global row "<< map->MaxMyGID64() <<endl
	 <<" and "<<A->NumMyNonzeros()<<" nonzeros."<<endl;

  A->FillComplete();
  

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
  
  delete [] list_of_vals;
  delete []list_of_inds;

  return;
}


#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES

void Trilinos_Util_ReadHpc2Epetra(char *data_file,
				 const Epetra_Comm  &comm, 
				 Epetra_Map *& map, 
				 Epetra_CrsMatrix *& A, 
				 Epetra_Vector *& x, 
				 Epetra_Vector *& b,
				 Epetra_Vector *&xexact) {
  Trilinos_Util_ReadHpc2Epetra_internal<int>(data_file, comm, map, A, x, b, xexact);
}

#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES

void Trilinos_Util_ReadHpc2Epetra64(char *data_file,
				 const Epetra_Comm  &comm, 
				 Epetra_Map *& map, 
				 Epetra_CrsMatrix *& A, 
				 Epetra_Vector *& x, 
				 Epetra_Vector *& b,
				 Epetra_Vector *&xexact) {
  Trilinos_Util_ReadHpc2Epetra_internal<long long>(data_file, comm, map, A, x, b, xexact);
}

#endif
