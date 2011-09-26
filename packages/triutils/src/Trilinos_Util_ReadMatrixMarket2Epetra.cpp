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

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "Trilinos_Util_CountMatrixMarket.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

int Trilinos_Util_ReadMatrixMarket2Epetra( char *data_file,
				      const Epetra_Comm  &comm, 
				      Epetra_Map *& map, 
				      Epetra_CrsMatrix *& A, 
				      Epetra_Vector *& x, 
				      Epetra_Vector *& b,
				      Epetra_Vector *&xexact ) {
  FILE *in_file ;
  int N_rows, nnz ; 

  bool diag = true;   // If set to false, this prevents any entries on the diagonal 
                      
  const int BUFSIZE = 800 ; 
  char buffer[BUFSIZE] ; 
  std::vector<int> non_zeros;   // Number of non-zeros in each row
  Trilinos_Util_CountMatrixMarket( data_file, non_zeros, N_rows, nnz, comm ) ; 

  std::vector<int> ptrs(N_rows+1) ; // Pointers into inds and vals for the start of each row
  std::vector<int> inds(nnz);     //  Column Indices
  std::vector<double> vals(nnz);  //  Matrix values
  std::vector<int> iptrs ;        //  Current pointers into inds and vals for each row

  if(comm.MyPID() == 0)  { 
    //  ptrs, inds and vals together constitute a compressed row storage of the matrix.


    in_file = fopen( data_file, "r");
    assert (in_file != NULL) ;  // Checked in Trilinos_Util_CountMatrixMarket() 

    ptrs[0] = 0 ; 
    for( int i=0; i< N_rows; i++ ) { 
      ptrs[i+1] = ptrs[i] + non_zeros[i]; 
    }

    iptrs = ptrs ; //  Current pointers into inds and vals for each row

    fgets( buffer, BUFSIZE, in_file ) ;  // Pick symmetry info off of this string 
    bool symmetric = false ; 
    string headerline1 = buffer;
    if ( headerline1.find("symmetric") != string::npos) symmetric = true; 
    fgets( buffer, BUFSIZE, in_file ) ;

    while ( fgets( buffer, BUFSIZE, in_file ) ) { 
      int i, j; 
      double val ; 
      i = -13 ;   // Check for blank lines 
      sscanf( buffer, "%d %d %lg", &i, &j, &val ) ; 
      assert( i != -13) ; 
      if ( diag || i != j ) { 
	//	if ( i == j && i == 1 ) val *= 1.0001 ;
	int iptr = iptrs[i-1] ; 
	iptrs[i-1]++ ;
	vals[iptr] = val ; 
	inds[iptr] = j-1 ; 
	//
	//  If this is a symmetric matrix, we need to enter the entry 
	//  for the other triangular half
	//
	if (symmetric && i != j ) {
	  iptr = iptrs[j-1] ; 
	  iptrs[j-1]++;
	  vals[iptr] = val ; 
	  inds[iptr] = i-1 ; 
	}
      }
    }
    fclose(in_file);


    if ( diag ) { 
      for (int i=0; i<N_rows; i++)
	assert( iptrs[i] == ptrs[i+1] ) ; 
    }
  }


  int nlocal = 0;
  if (comm.MyPID()==0) nlocal = N_rows;
  map = new Epetra_Map(N_rows, nlocal, 0, comm); // Create map with all elements on PE 0
  
  A = new Epetra_CrsMatrix(Copy, *map, 0); // Construct matrix on PE 0
  if (comm.MyPID()==0)
    for (int i=0; i<N_rows; i++) {
      A->InsertGlobalValues(i, iptrs[i]-ptrs[i], &vals[ptrs[i]], &inds[ptrs[i]]);
    }
  A->FillComplete();

  Epetra_Vector diagA(*map);

  A->ExtractDiagonalCopy( diagA ) ; 

  assert( diag || A->NoDiagonal() ) ;

  std::vector<double> xyz(N_rows);
  
  x = new Epetra_Vector(Copy, *map, &xyz[0]);
  b = new Epetra_Vector(Copy, *map, &xyz[0]);
  xexact = new Epetra_Vector(Copy, *map, &xyz[0]);

  EPETRA_CHK_ERR( x->PutScalar( 0.0 ) );
  EPETRA_CHK_ERR( xexact->Random( ) ) ; 
  EPETRA_CHK_ERR( A->Multiply( false, *xexact, *b ) ); 

  assert( map->SameAs(A->RowMap()) ) ;
 
  return 0;
}
