
#include <EpetraExt_SolverMap_CrsMatrix.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>

#include <vector>

namespace EpetraExt {

CrsMatrix_SolverMap::
~CrsMatrix_SolverMap()
{
  if( newObj_ && newObj_ != origObj_ ) delete newObj_;
  if( NewGraph_ ) delete NewGraph_;
  if( NewColMap_ ) delete NewColMap_;
}

CrsMatrix_SolverMap::NewTypeRef
CrsMatrix_SolverMap::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  assert( !orig.IndicesAreGlobal() );

  //test if matrix has missing local columns in its col map
  const Epetra_Map & RowMap = orig.RowMap();
  const Epetra_Map & ColMap = orig.ColMap();
  const Epetra_Comm & Comm = RowMap.Comm();
  int NumMyRows = RowMap.NumMyElements();
  int Match = 0;
  for( int i = 0; i < NumMyRows; ++i )
    if( RowMap.GID(i) != ColMap.GID(i) )
    {
      Match = 1;
      break;
    }

  int MatchAll = 0;
  Comm.SumAll( &Match, &MatchAll, 1 );

  if( !MatchAll )
  {
    newObj_ = origObj_;
  }
  else
  {
//    cout << "RowMap\n";
//    cout << RowMap;
//    cout << "ColMap\n";
//    cout << ColMap;
    //create ColMap with all local rows included
    vector<int> Cols(NumMyRows);
    for( int i = 0; i < NumMyRows; ++i )
      Cols[i] = RowMap.GID(i);

    int NumMyCols = ColMap.NumMyElements();
    for( int i = 0; i < NumMyCols; ++i )
      if( !RowMap.MyGID( ColMap.GID(i) ) ) Cols.push_back( ColMap.GID(i) );
    
    int NewNumMyCols = Cols.size();
    int NewNumGlobalCols;
    Comm.SumAll( &NewNumMyCols, &NewNumGlobalCols, 1 );
    NewColMap_ = new Epetra_Map( NewNumGlobalCols, NewNumMyCols,  &Cols[0], RowMap.IndexBase(), Comm );

    cout << RowMap;
    Comm.Barrier();
    cout << ColMap;
    Comm.Barrier();
    cout << *NewColMap_;
    Comm.Barrier();

    //New Graph
    vector<int> NumIndicesPerRow( NumMyRows );
    for( int i = 0; i < NumMyRows; ++i )
      NumIndicesPerRow[i] = orig.NumMyEntries(i);
    NewGraph_ = new Epetra_CrsGraph( Copy, RowMap, *NewColMap_, &NumIndicesPerRow[0] );

    int MaxNumEntries = orig.MaxNumEntries();
    int NumEntries;
    vector<int> Indices( MaxNumEntries );
    for( int i = 0; i < NumMyRows; ++i )
    {
      int RowGID = RowMap.GID(i);
      orig.Graph().ExtractGlobalRowCopy( RowGID, MaxNumEntries, NumEntries, &Indices[0] );
      NewGraph_->InsertGlobalIndices( RowGID, NumEntries, &Indices[0] );
    }
    NewGraph_->TransformToLocal();

//    cout << *NewGraph_;

    //intial construction of matrix 
    Epetra_CrsMatrix * NewMatrix = new Epetra_CrsMatrix( View, *NewGraph_ );

    //insert views of row values
    int * myIndices;
    double * myValues;
    int indicesCnt;
    int numMyRows = NewMatrix->NumMyRows();
    for( int i = 0; i < numMyRows; ++i )
    {
      orig.ExtractMyRowView( i, indicesCnt, myValues, myIndices );
      NewGraph_->ExtractMyRowView( i, indicesCnt, myIndices );

      NewMatrix->InsertMyValues( i, indicesCnt, myValues, myIndices );
    }

    NewMatrix->TransformToLocal();

    newObj_ = NewMatrix;
  }

  return *newObj_;
}

} // namespace EpetraExt

