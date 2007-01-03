#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_Vector.h"
#include "Teuchos_RefCountPtr.hpp"
#include <vector>
#include <assert.h>

using namespace Teuchos;
int SmallRowPermute( int in ) { return in + 2 ; } 
int BigRowPermute( int in ) { return 5*in + 1 ; } 
int NoPermute( int in ) { return in ; } 
int SmallColPermute( int in ) { return in + 4 ; } 
  
//
//  Diagonal:  0=no change, 1=eliminate entry
//             from the map for the largest row element in process 0
//             2=add diagonal entries to the matrix, with a zero value 
//             (assume row map contains all diagonal entries). 
//
//  ReindexRowMap:  
//    0=no change, 1= add 2 (still contiguous), 2=non-contiguous
//  
//  ReindexColMap
//    0=same as RowMap, 1=add 4 - Different From RowMap, but contiguous) 
//
//  RangeMap:
//    0=no change, 1=serial map, 2=bizarre distribution, 3=replicated map
//
//  DomainMap:
//    0=no change, 1=serial map, 2=bizarre distribution, 3=replicated map
//
RefCountPtr<Epetra_CrsMatrix> NewMatNewMap(Epetra_CrsMatrix& In, 
					   int Diagonal, 
					   int ReindexRowMap,
					   int ReindexColMap,
					   int RangeMapType,
					   int DomainMapType
					   )
{

  //
  //  If we are making no change, return the original matrix (which has a linear map) 
  //
#if 0
  cout << __FILE__ << "::" << __LINE__ << " " 
       << Diagonal << " " 
       << ReindexRowMap << " " 
       << ReindexColMap << " " 
       << RangeMapType << " " 
       << DomainMapType << " " << endl ; 
#endif

  if ( Diagonal + ReindexRowMap + ReindexColMap + RangeMapType + DomainMapType == 0 ) {
    RefCountPtr<Epetra_CrsMatrix> ReturnOrig = rcp( &In, false );
    return ReturnOrig ;
  }

  //
  //  Diagonal==2 is used for a different purpose - 
  //    Making sure that the diagonal of the matrix is non-empty.
  //  Note:  The diagonal must exist in In.RowMap().
  //
  if ( Diagonal == 2 ) { 
    assert( ReindexRowMap==0 && ReindexColMap == 0 ) ; 
  }

  int (*RowPermute)(int in) = 0;
  int (*ColPermute)(int in) = 0;

  assert( Diagonal >= 0  && Diagonal <= 2 ); 
  assert( ReindexRowMap>=0 && ReindexRowMap<=2 );
  assert( ReindexColMap>=0 && ReindexColMap<=1 );
  assert( RangeMapType>=0 && RangeMapType<=3 );
  assert( DomainMapType>=0 && DomainMapType<=3 );

  Epetra_Map DomainMap = In.DomainMap();
  Epetra_Map RangeMap = In.RangeMap();
  Epetra_Map ColMap = In.ColMap();
  Epetra_Map RowMap = In.RowMap();
  int NumMyRowElements = RowMap.NumMyElements();
  int NumMyColElements = ColMap.NumMyElements();
  int NumMyRangeElements = RangeMap.NumMyElements();
  int NumMyDomainElements = DomainMap.NumMyElements();

  int NumGlobalRowElements = RowMap.NumGlobalElements();
  int NumGlobalColElements = ColMap.NumGlobalElements();
  int NumGlobalRangeElements = RangeMap.NumGlobalElements();
  int NumGlobalDomainElements = DomainMap.NumGlobalElements();
  assert( NumGlobalRangeElements == NumGlobalDomainElements ) ; 

  vector<int> MyGlobalRowElements( NumMyRowElements ) ; 
  vector<int> NumEntriesPerRow( NumMyRowElements ) ; 
  vector<int> MyPermutedGlobalRowElements( NumMyRowElements ) ; 
  vector<int> MyGlobalColElements( NumMyColElements ) ; 
  vector<int> MyPermutedGlobalColElements( NumMyColElements ) ; // Used to create the column map
  vector<int> MyPermutedGlobalColElementTable( NumMyColElements ) ; // To convert local indices to global
  vector<int> MyGlobalRangeElements( NumMyRangeElements ) ; 
  vector<int> MyPermutedGlobalRangeElements( NumMyRangeElements ) ; 
  vector<int> MyGlobalDomainElements( NumMyDomainElements ) ; 
  vector<int> MyPermutedGlobalDomainElements( NumMyDomainElements ) ; 
  RowMap.MyGlobalElements(&MyGlobalRowElements[0]);
  ColMap.MyGlobalElements(&MyGlobalColElements[0]);
  RangeMap.MyGlobalElements(&MyGlobalRangeElements[0]);
  DomainMap.MyGlobalElements(&MyGlobalDomainElements[0]);

  switch( ReindexRowMap ) {
  case 0:
    RowPermute = &NoPermute ;
    break; 
  case 1:
    RowPermute = &SmallRowPermute ;
    break; 
  case 2:
    RowPermute = BigRowPermute ;
    break; 
  }
  switch( ReindexColMap ) {
  case 0:
    ColPermute = RowPermute ;
    break; 
  case 1:
    ColPermute = &SmallColPermute ;
    break; 
  }

  //
  //  Create Serial Range and Domain Maps based on the permuted indexing
  //
  int nlocal = 0;
  if (In.Comm().MyPID()==0) nlocal = NumGlobalRangeElements;
  vector<int> AllIDs( NumGlobalRangeElements ) ; 
  for ( int i = 0; i < NumGlobalRangeElements ; i++ ) AllIDs[i] = (*RowPermute)( i ) ; 
  Epetra_Map SerialRangeMap( -1, nlocal, &AllIDs[0], 0, In.Comm()); 
  vector<int> AllIDBs( NumGlobalRangeElements ) ; 
  for ( int i = 0; i < NumGlobalRangeElements ; i++ ) AllIDBs[i] = (*ColPermute)( i ) ; 
  Epetra_Map SerialDomainMap( -1, nlocal, &AllIDBs[0], 0, In.Comm()); 

  //
  //  Create Bizarre Range and Domain Maps based on the permuted indexing
  //  These are nearly serial, having all but one element on process 0
  //  The goal here is to make sure that we can use Domain and Range maps 
  //  that are neither serial, nor distributed in the normal manner.
  //
  vector<int> AllIDCs( NumGlobalRangeElements ) ; 
  for ( int i = 0; i < NumGlobalRangeElements ; i++ ) AllIDCs[i] = (*ColPermute)( i ) ; 
  if ( In.Comm().NumProc() > 1 ) { 
    if (In.Comm().MyPID()==0) nlocal = NumGlobalRangeElements-1;
    if (In.Comm().MyPID()==1) {
      nlocal = 1;
      AllIDCs[0] = (*ColPermute)( NumGlobalRangeElements - 1 );
    }
  } 
  int iam = In.Comm().MyPID();
  Epetra_Map BizarreDomainMap( -1, nlocal, &AllIDCs[0], 0, In.Comm()); 

  vector<int> AllIDDs( NumGlobalRangeElements ) ; 
  for ( int i = 0; i < NumGlobalRangeElements ; i++ ) AllIDDs[i] = (*RowPermute)( i ) ; 
  if ( In.Comm().NumProc() > 1 ) { 
    if (In.Comm().MyPID()==0) nlocal = NumGlobalRangeElements-1;
    if (In.Comm().MyPID()==1) {
      nlocal = 1;
      AllIDDs[0] = (*RowPermute)( NumGlobalRangeElements -1 ) ;
    }
  } 
  Epetra_Map BizarreRangeMap( -1, nlocal, &AllIDDs[0], 0, In.Comm()); 


  //
  //  Compute the column map 
  //
  //  If Diagonal==1, remove the column corresponding to the last row owned 
  //  by process 0.  Removing this column from a tridiagonal matrix, leaves
  //  a disconnected, but non-singular matrix.  
  //
  int NumMyColElementsOut = 0 ; 
  int NumGlobalColElementsOut ; 
  if ( Diagonal == 1 ) 
    NumGlobalColElementsOut = NumGlobalColElements-1; 
  else
    NumGlobalColElementsOut = NumGlobalColElements; 
  if ( Diagonal == 1 && iam==0 ) { 
    for ( int i=0; i < NumMyColElements  ; i++ ) {
      if ( MyGlobalColElements[i] != MyGlobalRowElements[NumMyRowElements-1] ) {
	MyPermutedGlobalColElements[NumMyColElementsOut++] = 
	  (*ColPermute)( MyGlobalColElements[i] ) ; 
      }
    }
    assert( NumMyColElementsOut == NumMyColElements-1 );
  } else {
    for ( int i=0; i < NumMyColElements  ; i++ )  
      MyPermutedGlobalColElements[i] = 
	(*ColPermute)( MyGlobalColElements[i] ) ; 
    NumMyColElementsOut = NumMyColElements ; 
    if ( Diagonal == 2 ) {
      //  For each row, make sure that the column map has this row in it, 
      //    if it doesn't, add it to the column map.  
      //  Note:  MyPermutedGlobalColElements == MyGlobalColElements when 
      //  Diagonal==2 because  ( Diagonal == 2 ) implies:
      //     ReindexRowMap==0 && ReindexColMap == 0  - see assert above
      for ( int i=0; i < NumMyRowElements  ; i++ ) {
	bool MissingDiagonal = true; 
	for ( int j=0; j < NumMyColElements; j++ ) { 
	  if ( MyGlobalRowElements[i] == MyGlobalColElements[j] ) {
	    MissingDiagonal = false; 
	  }
	}
	if ( MissingDiagonal ) {
	  MyPermutedGlobalColElements.resize(NumMyColElements+1);
	  MyPermutedGlobalColElements[NumMyColElementsOut] = MyGlobalRowElements[i];
	  NumMyColElementsOut++;
	}
      }
      In.Comm().SumAll(&NumMyColElementsOut,&NumGlobalColElementsOut,1); 
    }
  }

  //
  //  These tables are used both as the permutation tables and to create the maps.
  //
  for ( int i=0; i < NumMyColElements  ; i++ ) 
    MyPermutedGlobalColElementTable[i] = 
      (*ColPermute)( MyGlobalColElements[i] ) ; 
  for ( int i=0; i < NumMyRowElements  ; i++ ) 
    MyPermutedGlobalRowElements[i] = 
      (*RowPermute)( MyGlobalRowElements[i] ) ; 
  for ( int i=0; i < NumMyRangeElements  ; i++ ) 
    MyPermutedGlobalRangeElements[i] = 
      (*RowPermute)( MyGlobalRangeElements[i] ) ; 
  for ( int i=0; i < NumMyDomainElements  ; i++ ) 
    MyPermutedGlobalDomainElements[i] = 
      (*ColPermute)( MyGlobalDomainElements[i] ) ; 

  RefCountPtr<Epetra_Map> PermutedRowMap = 
    rcp( new Epetra_Map( NumGlobalRowElements, NumMyRowElements, 
			 &MyPermutedGlobalRowElements[0], 0, In.Comm() ) ); 
									
  RefCountPtr<Epetra_Map> PermutedColMap = 
    rcp( new Epetra_Map( NumGlobalColElementsOut, NumMyColElementsOut, 
			 &MyPermutedGlobalColElements[0], 0, In.Comm() ) ); 
									
  RefCountPtr<Epetra_Map> PermutedRangeMap = 
    rcp( new Epetra_Map( NumGlobalRangeElements, NumMyRangeElements, 
			 &MyPermutedGlobalRangeElements[0], 0, In.Comm() ) ); 
									
  RefCountPtr<Epetra_Map> PermutedDomainMap = 
    rcp( new Epetra_Map( NumGlobalDomainElements, NumMyDomainElements, 
			 &MyPermutedGlobalDomainElements[0], 0, In.Comm() ) ); 
									
  //
  //  These vectors are filled and then passed to InsertGlobalValues 
  //
  vector<int> ThisRowIndices( In.MaxNumEntries() );
  vector<double> ThisRowValues( In.MaxNumEntries() );
  vector<int> PermutedGlobalColIndices( In.MaxNumEntries() );

  //cout << __FILE__ << "::" <<__LINE__ << endl ; 
  RefCountPtr<Epetra_CrsMatrix> Out = 
    rcp( new Epetra_CrsMatrix( Copy, *PermutedRowMap, *PermutedColMap, 0 ) );

  for (int i=0; i<NumMyRowElements; i++)
    {

      int NumIndicesThisRow;
      assert( In.ExtractMyRowCopy( i, 
				   In.MaxNumEntries(),
				   NumIndicesThisRow,
				   &ThisRowValues[0],
				   &ThisRowIndices[0] ) == 0 ) ;
      for (int j = 0 ; j < NumIndicesThisRow ; j++ )
	{
	  PermutedGlobalColIndices[j] = MyPermutedGlobalColElementTable[ ThisRowIndices[j] ]  ;
	}
      bool MissingDiagonal = false; 
      if ( Diagonal==2 ) { 
	//
	assert( MyGlobalRowElements[i] == MyPermutedGlobalRowElements[i] );
	MissingDiagonal = true; 
	for( int j =0 ; j < NumIndicesThisRow ; j++ ) {
	  if ( PermutedGlobalColIndices[j] == MyPermutedGlobalRowElements[i] ) {
	    MissingDiagonal = false ; 
	  }
	}
#if 0
	cout  << __FILE__ << "::" << __LINE__ 
	      << " i = " << i 
	      << " MyPermutedGlobalRowElements[i]  = " << MyPermutedGlobalRowElements[i] 
	      <<   " MissingDiagonal = " << MissingDiagonal << endl ; 
#endif

      }
      if ( MissingDiagonal ) { 
	ThisRowValues.resize(NumIndicesThisRow+1) ; 
	ThisRowValues[NumIndicesThisRow] = 0.0;
	PermutedGlobalColIndices.resize(NumIndicesThisRow+1);
	PermutedGlobalColIndices[NumIndicesThisRow] = MyPermutedGlobalRowElements[i] ;
	
#if 0
	cout  << __FILE__ << "::" << __LINE__ 
	      << " i = " << i 
	      << "NumIndicesThisRow = " << NumIndicesThisRow 
	      << "ThisRowValues[NumIndicesThisRow = " << ThisRowValues[NumIndicesThisRow] 
	      << " PermutedGlobalColIndices[NumIndcesThisRow] = " << PermutedGlobalColIndices[NumIndicesThisRow] 
	      << endl ; 
#endif

	NumIndicesThisRow++  ;

      } 
      assert( Out->InsertGlobalValues( MyPermutedGlobalRowElements[i], 
				       NumIndicesThisRow,
				       &ThisRowValues[0],
				       &PermutedGlobalColIndices[0] ) >= 0 ); 
    }

  //

  Epetra_LocalMap ReplicatedMap( NumGlobalRangeElements, 0, In.Comm() );

  RefCountPtr<Epetra_Map> OutRangeMap ;
  RefCountPtr<Epetra_Map> OutDomainMap ;
  
  switch( RangeMapType ) {
  case 0:
    OutRangeMap = PermutedRangeMap ;
    break;
  case 1:
    OutRangeMap = rcp(&SerialRangeMap, false); 
    break;
  case 2:
    OutRangeMap = rcp(&BizarreRangeMap, false); 
    break;
  case 3:
    OutRangeMap = rcp(&ReplicatedMap, false); 
    break;
  }
  //  switch( DomainMapType ) {
  switch( DomainMapType ) {
  case 0:
    OutDomainMap = PermutedDomainMap ;
    break;
  case 1:
    OutDomainMap = rcp(&SerialDomainMap, false); 
    break;
  case 2:
    OutDomainMap = rcp(&BizarreDomainMap, false); 
    break;
  case 3:
    OutDomainMap = rcp(&ReplicatedMap, false); 
    break;
  }
#if 0
  assert(Out->FillComplete( *PermutedDomainMap, *PermutedRangeMap )==0);
#else
  assert(Out->FillComplete( *OutDomainMap, *OutRangeMap )==0);
#endif

#if 0
  cout << __FILE__ << "::" << __LINE__ << endl ;
  Out->Print( cout ) ; 
#endif

  return Out;
}



