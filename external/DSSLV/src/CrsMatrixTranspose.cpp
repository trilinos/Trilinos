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

   
  int iam = In->Comm().MyPID() ;

  int numentries = In->NumGlobalNonzeros();
  int NumRowEntries;
  double *RowValues;
  int *ColIndices;

  int numrows = In->NumGlobalRows();
  int numcols = In->NumGlobalCols();

  vector <int> Ap( numcols+1 );       // Column i is stored in Aval(Ap[i]..Ap[i+1]-1)
  vector <int> nextAp( numcols+1 );   // Where to store next value in Column i
  vector <int> Ai( EPETRA_MAX( numcols, numentries) ) ; //  Row indices
  vector <double> Aval( EPETRA_MAX( numcols, numentries) ) ; 

  if ( iam == 0 ) { 

    assert( In->NumMyRows() == In->NumGlobalRows() ) ; 
    //
    //  Count the number of entries in each column
    //
    vector <int>RowsPerCol( numcols ) ; 
    for ( int i = 0 ; i < numcols ; i++ ) RowsPerCol[i] = 0 ; 
    for ( int MyRow = 0; MyRow <numrows; MyRow++ ) {
      assert( In->ExtractMyRowView( MyRow, NumRowEntries, RowValues, ColIndices ) == 0 ) ;
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
      assert( In->ExtractMyRowView( MyRow, NumRowEntries, RowValues, ColIndices ) == 0 ) ;
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
      Out->InsertGlobalValues( MyRow, NumInCol, &Aval[Ap[MyRow]], 
			   &Ai[Ap[MyRow]] );
      assert( Out->IndicesAreGlobal() ) ; 
    }
  } else {
    assert( In->NumMyRows() == 0 ) ; 
  }


  assert( Out->TransformToLocal()==0 ) ;
  return 0 ; 
}
