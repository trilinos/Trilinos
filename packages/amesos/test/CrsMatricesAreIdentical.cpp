#include "CrsMatricesAreIdentical.h"
#include "Epetra_Comm.h"

int CrsMatricesAreIdentical( Epetra_CrsMatrix *A, Epetra_CrsMatrix *B ) {

  int numrows = A->NumMyRows();
  int iam = A->Comm().MyPID() ; 

  assert( A->NumGlobalNonzeros() == B->NumGlobalNonzeros() ) ; 
  assert( A->NumGlobalRows() == B->NumGlobalRows() ) ; 
  assert( A->NumGlobalCols() == B->NumGlobalCols() ) ; 
  if( ! ( A->NumGlobalDiagonals() == B->NumGlobalDiagonals() ) ) { 
    cout << "  A->NumGlobalDiagonals() = " <<  A->NumGlobalDiagonals() << endl ; 
    cout << "  B->NumGlobalDiagonals() = " <<  B->NumGlobalDiagonals() << endl ; 
  }
  assert( A->NumGlobalDiagonals() == B->NumGlobalDiagonals() ) ; 
  assert( A->NumMyNonzeros() == B->NumMyNonzeros() ) ; 
  assert( A->NumMyRows() == B->NumMyRows() ) ; 
  if( ! ( A->NumMyCols() == B->NumMyCols() ) ) { 
    cout << " iam = " << iam << "  A->NumMyRows() = " <<  A->NumMyRows() << endl ; 
    cout << " iam = " << iam << "  A->NumMyCols() = " <<  A->NumMyCols() << endl ; 
    cout <<  " iam = " << iam << "  B->NumMyCols() = " <<  B->NumMyCols() << endl ; 
  }
  //  assert( A->NumMyCols() == B->NumMyCols() ) ; 
  assert( A->NumMyDiagonals() == B->NumMyDiagonals() ) ; 




  if( A->MaxNumEntries() != B->MaxNumEntries() ) {
    cout << "  A->MaxNumEntries() = " <<  A->MaxNumEntries() 
	 << "  B->MaxNumEntries() = " <<  B->MaxNumEntries() << endl;
  }
  assert( A->GlobalMaxNumEntries() == B->GlobalMaxNumEntries() ) ; 

  assert( A->IndexBase() == B->IndexBase() ) ; 
  int *Aindices, *Bindices;
  double *Avalues, *Bvalues;
  const int MaxVals = 1000;
  double Aval[MaxVals], Bval[MaxVals];
  int Aind[MaxVals], Bind[MaxVals];
  int Anumentries, Bnumentries;
  for ( int i = 0; i < numrows; i++ ) { 
    //    cout << " i = " << i << endl ; 
    // cout << " CrsMatrices AreIdentical::  A->NumGlobalEntries = " 
    //	 << A->NumGlobalEntries(i) << endl ; 
    // cout << " CrsMatrices AreIdentical::  B->NumGlobalEntries = " 
    // << B->NumGlobalEntries(i) << endl ; 
    Anumentries =  A->NumGlobalEntries(i) ;
    Bnumentries =  B->NumGlobalEntries(i) ;
    // cout << " i = " << i << endl ; 
    // cout << " CrsMatrices AreIdentical::  A->NumGlobalEntries = " 
    // << A->NumGlobalEntries(i) << endl ; 
    // cout << " CrsMatrices AreIdentical::  B->NumGlobalEntries = " 
    // << B->NumGlobalEntries(i) << endl ; 
    // cout << " ABOVE row number = " << i << " A num entries = " << Anumentries << " iam = " << iam << endl ; 
    // cout << " ABOVE row number = " << i << " B num entries = " << Bnumentries << " iam = " << iam << endl ; 
    //    assert( A->NumAllocatedGlobalEntries(i) == B->NumAllocatedGlobalEntries(i) ) ; 
    //
    //  I am suspicious that this will not work in parallel because
    //  my guess is that i is a local row number, not a global row number
    //
    A->ExtractMyRowView( i, Anumentries, Avalues, Aindices ) ;
    // cout << " row number = " << i << " A num entries = " << Anumentries << " iam = " << iam << endl ; 
    B->ExtractMyRowView( i, Bnumentries, Bvalues, Bindices ) ;
    // cout << " row number = " << i << " B num entries = " << Bnumentries << " iam = " << iam << endl ; 
    //    cout << " Avals = " << Avalues[0] << " " <<  Avalues[1] << " " <<  Avalues[2] << " " << endl ; 
    //    cout << " Bvals = " << Bvalues[0] << " " <<  Bvalues[1] << " " <<  Bvalues[2] << " " << Bvalues[3] <<  endl ; 
    //    cout << " Aindices = " << Aindices[0] << " " <<  Aindices[1] << " " <<  Aindices[2] << " " << endl ; 
    //    cout << " Bindices = " << Bindices[0] << " " <<  Bindices[1] << " " <<  Bindices[2] << " " << Bindices[3] <<  endl ; 

    A->ExtractMyRowCopy( i, MaxVals, Anumentries, (double *) &Aval, (int *)&Aind ) ;
    //    cout << " row number = " << i << " num entries = " << Anumentries << " iam = " << iam << endl ; 
    B->ExtractMyRowCopy( i, MaxVals, Bnumentries, (double *)&Bval, (int *)&Bind ) ;
    //    cout << " row number = " << i << " num entries = " << Bnumentries << " iam = " << iam << endl ; 
    //    cout << " Avals = " << Aval[0] << " " <<  Aval[1] << " " <<  Aval[2] << " " << endl ; 
    //    cout << " Bvals = " << Bval[0] << " " <<  Bval[1] << " " <<  Bval[2] << " " << Bval[3] <<  endl ; 
    //    cout << " Aind = " << Aind[0] << " " <<  Aind[1] << " " <<  Aind[2] << " " << endl ; 
    //    cout << " Bind = " << Bind[0] << " " <<  Bind[1] << " " <<  Bind[2] << " " << Bind[3] <<  endl ; 
    assert( Anumentries == Bnumentries ) ; 
    for ( int j = 0 ; j < Anumentries ; j++ ) { 
      assert( Avalues[j] == Bvalues[j] ) ; 
      assert( Aindices[j] == Bindices[j] ) ; 
      //      cout << " iam = " << iam << " Avalues[" << i << "," << Aindices[j] << "]=" << Avalues[j] << endl ; 
    }
    // cout << " i = " << i << endl ; 
    // cout << " CrsMatrices AreIdentical::  A->NumGlobalEntries = " 
    //	 << A->NumGlobalEntries(i) << endl ; 
    // cout << " CrsMatrices AreIdentical::  B->NumGlobalEntries = " 
    //	 << B->NumGlobalEntries(i) << endl ; 
    assert( A->NumGlobalEntries(i) == B->NumGlobalEntries(i) ) ; 
    assert( A->NumMyEntries(i) == B->NumMyEntries(i) ) ; 

  }

  return 1;
}
