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
#include <limits>
#include "Trilinos_Util_CountMatrixMarket.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

template<typename int_type>
int Trilinos_Util_ReadMatrixMarket2Epetra_internal( char *data_file,
				      const Epetra_Comm  &comm, 
				      Epetra_Map *& map, 
				      Epetra_CrsMatrix *& A, 
				      Epetra_Vector *& x, 
				      Epetra_Vector *& b,
				      Epetra_Vector *&xexact ) {
  FILE *in_file ;
  int_type N_rows, nnz ; 

  bool diag = true;   // If set to false, this prevents any entries on the diagonal 
                      
  const int BUFSIZE = 800 ; 
  char buffer[BUFSIZE] ; 
  std::vector<int> non_zeros;   // Number of non-zeros in each row
  Trilinos_Util_CountMatrixMarket( data_file, non_zeros, N_rows, nnz, comm ) ; 

  std::vector<int_type> ptrs(N_rows+1) ; // Pointers into inds and vals for the start of each row
  std::vector<int_type> inds(nnz);     //  Column Indices
  std::vector<double> vals(nnz);  //  Matrix values
  std::vector<int_type> iptrs ;        //  Current pointers into inds and vals for each row

  if(comm.MyPID() == 0)  { 
    //  ptrs, inds and vals together constitute a compressed row storage of the matrix.


    in_file = fopen( data_file, "r");
    assert (in_file != NULL) ;  // Checked in Trilinos_Util_CountMatrixMarket() 

    ptrs[0] = 0 ; 
    for( int_type i=0; i< N_rows; i++ ) { 
      ptrs[i+1] = ptrs[i] + non_zeros[i]; 
    }

    iptrs = ptrs ; //  Current pointers into inds and vals for each row

    // Pick symmetry info off of this string
    if (fgets( buffer, BUFSIZE, in_file) == NULL)
        assert(false);
    bool symmetric = false ; 
    string headerline1 = buffer;
    if ( headerline1.find("symmetric") != string::npos) symmetric = true; 
    if (fgets( buffer, BUFSIZE, in_file) == NULL)
        assert(false);

    while ( fgets( buffer, BUFSIZE, in_file ) ) { 
      int_type i, j; 
      double val ; 
      i = -13 ;   // Check for blank lines 
      if(sizeof(int) == sizeof(int_type))
        sscanf( buffer, "%d %d %lg", &i, &j, &val ) ; 
      else if(sizeof(long long) == sizeof(int_type))
        sscanf( buffer, "%lld %lld %lg", &i, &j, &val ) ; 
      else
        assert(false);
      assert( i != -13) ; 
      if ( diag || i != j ) { 
	//	if ( i == j && i == 1 ) val *= 1.0001 ;
	int_type iptr = iptrs[i-1] ; 
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
      for (int_type i=0; i<N_rows; i++)
	assert( iptrs[i] == ptrs[i+1] ) ; 
    }
  }


  int nlocal = 0;
  if (comm.MyPID()==0) {
    if(N_rows > std::numeric_limits<int>::max())
      throw "Triutils: Trilinos_Util_ReadMatrixMarket2Epetra_internal: N_rows > std::numeric_limits<int>::max()";
    nlocal = static_cast<int>(N_rows);
  }
  map = new Epetra_Map(N_rows, nlocal, (int_type) 0, comm); // Create map with all elements on PE 0
  
  A = new Epetra_CrsMatrix(Copy, *map, 0); // Construct matrix on PE 0
  if (comm.MyPID()==0)
    for (int_type i=0; i<N_rows; i++) {
      if(iptrs[i]-ptrs[i] > std::numeric_limits<int>::max())
        throw "Triutils: Trilinos_Util_ReadMatrixMarket2Epetra_internal: iptrs[i]-ptrs[i] > std::numeric_limits<int>::max()";
      A->InsertGlobalValues(i, (int) (iptrs[i]-ptrs[i]), &vals[ptrs[i]], &inds[ptrs[i]]);
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

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES

int Trilinos_Util_ReadMatrixMarket2Epetra( char *data_file,
				      const Epetra_Comm  &comm, 
				      Epetra_Map *& map, 
				      Epetra_CrsMatrix *& A, 
				      Epetra_Vector *& x, 
				      Epetra_Vector *& b,
				      Epetra_Vector *&xexact ) {
  return Trilinos_Util_ReadMatrixMarket2Epetra_internal<int>(data_file, comm, map, A, x, b, xexact);
}

#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES

int Trilinos_Util_ReadMatrixMarket2Epetra64( char *data_file,
				      const Epetra_Comm  &comm, 
				      Epetra_Map *& map, 
				      Epetra_CrsMatrix *& A, 
				      Epetra_Vector *& x, 
				      Epetra_Vector *& b,
				      Epetra_Vector *&xexact ) {
  return Trilinos_Util_ReadMatrixMarket2Epetra_internal<long long>(data_file, comm, map, A, x, b, xexact);
}

#endif
