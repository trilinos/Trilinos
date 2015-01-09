// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

//
//  CrsMatrixTranspose( Epetra_CrsMatrix *In,  Epetra_CrsMatrix *Out ) 
//  fills the matrix "Out" with the transpose of the matrix "In".
//
//  Notes:
//    Only works on pseudo-distributed matrices, where all rows are stored
//    on process 0.
//    Testing has been limited to square matrices
//    This implementation requires an extra copy of the matrix as an 
//    intermediate. 
//
//  One of the bizarre aspects of this code is that it uses the 
//  results of calls to ExtractMyRowView() to populate the input 
//  to calls to InsertGlobalValues().  This only works because
//  the code is only designed for matrices that all stored entirely on 
//  process 0.
//
//
#include "CrsMatrixTranspose.h"
#include <vector>
#include "Epetra_Comm.h"

int CrsMatrixTranspose( Epetra_CrsMatrix *In,  Epetra_CrsMatrix *Out ) { 

  int ierr = 0;
   
  int iam = In->Comm().MyPID() ;

  long long numentries = In->NumGlobalNonzeros64();
  int NumRowEntries = 0;
  double *RowValues = 0;
  int *ColIndices = 0;

  long long numrows = In->NumGlobalRows64();
  long long numcols = In->NumGlobalCols64();

  std::vector <int> Ap( numcols+1 );       // Column i is stored in Aval(Ap[i]..Ap[i+1]-1)
  std::vector <int> nextAp( numcols+1 );   // Where to store next value in Column i
#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
  std::vector <long long> Ai( EPETRA_MAX( numcols, numentries) ) ; //  Row indices
#else
  std::vector <int> Ai( EPETRA_MAX( numcols, numentries) ) ; //  Row indices
#endif
  std::vector <double> Aval( EPETRA_MAX( numcols, numentries) ) ; 

  if ( iam == 0 ) { 

    assert( In->NumMyRows() == In->NumGlobalRows64() ) ; 
    //
    //  Count the number of entries in each column
    //
    std::vector <int>RowsPerCol( numcols ) ; 
    for ( int i = 0 ; i < numcols ; i++ ) RowsPerCol[i] = 0 ; 
    for ( int MyRow = 0; MyRow <numrows; MyRow++ ) {
      ierr = In->ExtractMyRowView( MyRow, NumRowEntries, RowValues, ColIndices );
      assert( ierr == 0 ) ;
      for ( int j = 0; j < NumRowEntries; j++ ) { 
	RowsPerCol[ ColIndices[j] ] ++ ; 
      }
    }
    //
    //  Set Ap and nextAp based on RowsPerCol
    //
    Ap[0] = 0 ; 
    for ( int i = 0 ; i < numcols ; i++ ) {
      Ap[i+1]= Ap[i] + RowsPerCol[i] ; 
      nextAp[i] = Ap[i];
    }
    //
    //  Populate Ai and Aval 
    //
    for ( int MyRow = 0; MyRow <numrows; MyRow++ ) {
      ierr = In->ExtractMyRowView( MyRow, NumRowEntries, RowValues, ColIndices );
      assert( ierr == 0 ) ;
      for ( int j = 0; j < NumRowEntries; j++ ) { 
	Ai[ nextAp[ ColIndices[j] ] ] = MyRow ; 
	Aval[ nextAp[ ColIndices[j] ] ] = RowValues[j] ; 
	nextAp[ ColIndices[j] ] ++ ; 
      }
    }

    //
    //  Insert values into Out 
    //
    for ( int MyRow = 0; MyRow <numrows; MyRow++ ) {
      int NumInCol = Ap[MyRow+1] -  Ap[MyRow] ;
#if !defined(EPETRA_NO_32BIT_GLOBAL_INDICES) || !defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
      Out->InsertGlobalValues( MyRow, NumInCol, &Aval[Ap[MyRow]], 
			   &Ai[Ap[MyRow]] );
#endif
      assert( Out->IndicesAreGlobal() ) ; 
    }
  } else {
    assert( In->NumMyRows() == 0 ) ; 
  }

  ierr = Out->FillComplete();
  assert( ierr==0 ) ;
  return 0 ; 
}
