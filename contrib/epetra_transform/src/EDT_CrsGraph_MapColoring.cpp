
#include <EDT_CrsGraph_MapColoring.h>

#include <EDT_CrsGraph_Transpose.h>

#include <vector>
#include <set>
#include <map>

#include <Epetra_CrsGraph.h>
#include <Epetra_MapColoring.h>
#include <Epetra_Map.h>

namespace Epetra_Transform {

std::auto_ptr<Epetra_MapColoring> CrsGraph_MapColoring::operator()( const Epetra_CrsGraph & original )
{
  int err;

  const Epetra_BlockMap & RowMap = original.RowMap();
  int nRows = RowMap.NumMyElements();

  int NumIndices;
  int * Indices;

  //Generate Local Distance-1 Adjacency Graph
  //(Transpose of original excluding non-local column indices)
  std::auto_ptr<Epetra_CrsGraph> Adj1 = CrsGraph_Transpose( true )( original );
  if( verbose_ ) cout << "Adjacency 1 Graph!\n" << *Adj1;

  int Delta = Adj1->MaxNumIndices();
  cout << endl << "Delta: " << Delta << endl;

  //Generation of Local Distance-2 Adjacency Graph
  Epetra_CrsGraph Adj2( Copy, RowMap, RowMap, 0 );
  int NumAdj1Indices;
  int * Adj1Indices;
  for( int i = 0; i < nRows; ++i )
  {
    assert( Adj1->ExtractMyRowView( i, NumAdj1Indices, Adj1Indices ) == 0 );

    for( int j = 0; j < NumAdj1Indices; ++j )
    {
      assert( original.ExtractMyRowView( Adj1Indices[j], NumIndices, Indices ) == 0 );
      int NumLocalIndices = 0;
      for( int k = 0; k < NumIndices; ++k )
        if( Indices[k] < nRows ) NumLocalIndices++; 
      assert( Adj2.InsertMyIndices( i, NumLocalIndices, Indices ) >= 0 );
    }
  }
  assert( Adj2.TransformToLocal() == 0 );
  if( verbose_ ) cout << "Adjacency 2 Graph!\n" << Adj2;

  vector<int> rowOrder( nRows );
  //Simple reordering
  {
    multimap<int,int> adjMap;
    for( int i = 0; i < nRows; ++i )
      adjMap.insert( pair<int,int>( Adj2.NumMyIndices(i), i ) );
    multimap<int,int>::iterator iter = adjMap.begin();
    multimap<int,int>::iterator end = adjMap.end();
    for( int i = 1; iter != end; ++iter, ++i )
      rowOrder[ nRows - i ] = iter->second;
  }

  Epetra_MapColoring * ColorMap = new Epetra_MapColoring( RowMap );

  //Application of Greedy Algorithm to generate Color Map
  int Size = Delta * Delta + 1;
  set<int> allowedColors;
  for( int col = 0; col < nRows; ++col )
  {
    for( int i = 0; i < Size; ++i ) allowedColors.insert( i+1 ); 

    Adj2.ExtractMyRowView( rowOrder[col], NumIndices, Indices );

    for( int i = 0; i < NumIndices; ++i )
      if( (*ColorMap)[ Indices[i] ] > 0 ) allowedColors.erase( (*ColorMap)[ Indices[i] ] );

    (*ColorMap)[ rowOrder[col] ] = *(allowedColors.begin());
  }

  if( verbose_ ) cout << "ColorMap!\n" << *ColorMap;

  return std::auto_ptr<Epetra_MapColoring> ( ColorMap );
}

}
