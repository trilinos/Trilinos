
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>

#include <EDT_CrsMatrix_View.h>

namespace EpetraExt {

CrsMatrix_View::
~CrsMatrix_View()
{
  if( newObj_ ) delete newObj_;
}

CrsMatrix_View::NewTypeRef
CrsMatrix_View::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  if( orig.IndicesAreGlobal() ) cout << "EDT_CrsMatrix_View: Indices must be LOCAL!\n";
  assert( !orig.IndicesAreGlobal() );

  //test graph, new graph must be contiguous subset of old

  //intial construction of matrix 
  Epetra_CrsMatrix * newMatrix( new Epetra_CrsMatrix( View, NewGraph_ ) );

  //insert views of row values
  int * myIndices;
  double * myValues;
  int indicesCnt;
  int numMyRows = newMatrix->NumMyRows();
  for( int i = 0; i < numMyRows; ++i )
  {
    orig.ExtractMyRowView( i, indicesCnt, myValues, myIndices );

    int newIndicesCnt = indicesCnt;
    bool done = false;
    for( int j = 0; j < indicesCnt; ++j )
      if( !done && NewGraph_.GCID( myIndices[j] ) == -1 )
      {
        newIndicesCnt = j;
        done = true;
      }

    newMatrix->InsertMyValues( i, newIndicesCnt, myValues, myIndices );
  }

  newMatrix->TransformToLocal();

  newObj_ = newMatrix;

  return *newMatrix;
}

} // namespace EpetraExt

