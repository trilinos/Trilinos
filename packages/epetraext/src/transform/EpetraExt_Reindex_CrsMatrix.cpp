
#include <EpetraExt_Reindex_CrsMatrix.h>

#include <vector>

#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_IntVector.h>

#include <Epetra_Export.h>
#include <Epetra_Import.h>

namespace EpetraExt {

CrsMatrix_Reindex::
~CrsMatrix_Reindex()
{
  if( newObj_ ) delete newObj_;
  if( NewColMap_ ) delete NewColMap_;
}

CrsMatrix_Reindex::NewTypeRef
CrsMatrix_Reindex::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  //test map, must have same number of local and global elements as original row map
  Epetra_Map & OldRowMap = const_cast<Epetra_Map&>(orig.RowMap());
  Epetra_Map & OldColMap = const_cast<Epetra_Map&>(orig.ColMap());
  int NumMyElements = OldRowMap.NumMyElements();
  int NumGlobalElements = OldRowMap.NumGlobalElements();
  assert( OldRowMap.NumMyElements() == NewRowMap_.NumMyElements() );

  //Construct new Column Map
  Epetra_IntVector Cols( OldRowMap );
  Epetra_IntVector NewCols( OldColMap );
  Epetra_Import Importer( OldColMap, OldRowMap );

  for( int i = 0; i < NumMyElements; ++i )
    Cols[i] = NewRowMap_.GID(i);

  NewCols.Import( Cols, Importer, Insert );

  vector<int*> NewColIndices(1);
  NewCols.ExtractView( &NewColIndices[0] );

  int NumMyColElements = OldColMap.NumMyElements();
  int NumGlobalColElements = OldColMap.NumGlobalElements();

  NewColMap_ = new Epetra_Map( NumGlobalColElements, NumMyColElements, NewColIndices[0], 0, OldColMap.Comm() );

  //intial construction of matrix 
  Epetra_CrsMatrix * NewMatrix = new Epetra_CrsMatrix( View, NewRowMap_, *NewColMap_, 0 );

  //insert views of row values
  int * myIndices;
  double * myValues;
  int indicesCnt;
  int numMyRows = NewMatrix->NumMyRows();
  for( int i = 0; i < numMyRows; ++i )
  {
    orig.ExtractMyRowView( i, indicesCnt, myValues, myIndices );
    NewMatrix->InsertMyValues( i, indicesCnt, myValues, myIndices );
  }

  NewMatrix->TransformToLocal();

  newObj_ = NewMatrix;

  return *NewMatrix;
}

} // namespace EpetraExt

