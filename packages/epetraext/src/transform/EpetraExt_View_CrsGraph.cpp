
#include <EpetraExt_View_CrsGraph.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_BlockMap.h>

#include <vector>

using std::vector;

namespace EpetraExt {

CrsGraph_View::
~CrsGraph_View()
{
  if( newObj_ ) delete newObj_;
}

CrsGraph_View::NewTypeRef
CrsGraph_View::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  //Error, must be local indices
  assert( !orig.IndicesAreGlobal() );

  //test maps, new map must be left subset of old
  const Epetra_BlockMap & oRowMap = orig.RowMap();
  const Epetra_BlockMap & oColMap = orig.ColMap();

  int oNumRows = oRowMap.NumMyElements();
  int oNumCols = oRowMap.NumMyElements();
  int nNumRows = NewRowMap_->NumMyElements();
  int nNumCols = 0;
  if( NewColMap_ ) nNumCols = NewColMap_->NumMyElements();

  bool matched = true;
  for( int i = 0; i < nNumRows; ++i )
    matched = matched && ( oRowMap.GID(i) == NewRowMap_->GID(i) );
  if( nNumCols )
    for( int i = 0; i < nNumCols; ++i )
      matched = matched && ( oColMap.GID(i) == NewColMap_->GID(i) );

  if( !matched ) cout << "EDT_CrsGraph_View: Bad Row or Col Mapping\n";
  assert( matched );

  //intial construction of graph
  vector<int> numIndices( nNumRows );
  vector<int*> indices( nNumRows );
  for( int i = 0; i < nNumRows; ++i )
  {
    orig.ExtractMyRowView( i, numIndices[i], indices[i] );
    int j = 0;
    if( nNumCols )
    {
      while( j < numIndices[i] && NewColMap_->GID(indices[i][j]) != -1 ) ++j;
      numIndices[i] = j;
    }
  }

  Epetra_CrsGraph * newGraph( new Epetra_CrsGraph( View,
                                                   *NewRowMap_,
                                                   *NewColMap_,
                                                   &numIndices[0] ) );

  //insert views of row indices
  for( int i = 0; i < nNumRows; ++i )
    newGraph->InsertMyIndices( i, numIndices[i], indices[i] );

  newGraph->TransformToLocal();

  newObj_ = newGraph;

  return *newGraph;
}

} // namespace EpetraExt

