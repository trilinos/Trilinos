
#include <vector>

#include <EDT_CrsGraph_View.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_BlockMap.h>

namespace Epetra_Transform {

std::auto_ptr<Epetra_CrsGraph> CrsGraph_View::operator()( const Epetra_CrsGraph & original )
{
  if( original.IndicesAreGlobal() ) 
    //Error, must be local indices
    return std::auto_ptr<Epetra_CrsGraph>(0);

  //test maps, new map must be left subset of old

  //intial construction of graph
  int numMyRows = NewMap_.NumMyElements();
  vector<int> numIndices( numMyRows );
  vector<int*> myIndices( numMyRows );
  for( int i = 0; i < numMyRows; ++i )
  {
    original.ExtractMyRowView( i, numIndices[i], myIndices[i] );
    int j = 0;
    while( j < numMyRows && myIndices[i][j] < numMyRows ) ++j;
    numIndices[i] = j;
  }
  std::auto_ptr<Epetra_CrsGraph> newGraph( new Epetra_CrsGraph( View,
                                                    NewMap_,
                                                    &numIndices[0] ) );
  newGraph->TransformToLocal();

  //insert views of row indices
  for( int i = 0; i < numMyRows; ++i )
    newGraph->InsertMyIndices( i, numIndices[i], myIndices[i] );

  return newGraph;
}

} //namespace Epetra_Transform
