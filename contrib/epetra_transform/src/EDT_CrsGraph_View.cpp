
#include <vector>

#include <EDT_CrsGraph_View.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_BlockMap.h>

namespace Epetra_Transform {

std::auto_ptr<Epetra_CrsGraph> CrsGraph_View::operator()( const Epetra_CrsGraph & original )
{
  //Error, must be local indices
  if( original.IndicesAreGlobal() ) return std::auto_ptr<Epetra_CrsGraph>(0);

  //test maps, new map must be left subset of old

  //intial construction of graph
  int numMyRows = NewRowMap_.NumMyElements();
  int numMyCols = -1;
  if( NewDomainMap_ ) numMyCols = NewDomainMap_->NumMyElements();

  vector<int> numIndices( numMyRows );
  vector<int*> indices( numMyRows );
  for( int i = 0; i < numMyRows; ++i )
  {
    original.ExtractMyRowView( i, numIndices[i], indices[i] );
    int j = 0;
    if( numMyCols != -1 )
    {
      while( j < numIndices[i] && NewDomainMap_->GID(indices[i][j]) != -1 ) ++j;
      numIndices[i] = j;
    }
  }
  std::auto_ptr<Epetra_CrsGraph> newGraph( new Epetra_CrsGraph( View,
                                                                NewRowMap_,
                                                                &numIndices[0] ) );

  //insert views of row indices
  for( int i = 0; i < numMyRows; ++i )
    newGraph->InsertMyIndices( i, numIndices[i], indices[i] );

  if( NewDomainMap_ ) newGraph->TransformToLocal( const_cast<Epetra_BlockMap*>(NewDomainMap_),
                                                  const_cast<Epetra_BlockMap*>(&NewRowMap_) );
  else                newGraph->TransformToLocal();

  return newGraph;
}

} //namespace Epetra_Transform
