
#include <EDT_CrsGraph_MapColoring.h>

#include <EDT_CrsGraph_Transpose.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_MapColoring.h>
#include <Epetra_Map.h>

#include <vector>
#include <set>
#include <map>

using std::vector;
using std::set;
using std::map;

namespace EpetraExt {
namespace Transform {

NewTypePtr CrsGraph_MapColoring::operator()( OriginalTypeRef original )
{
  int err;

  const Epetra_BlockMap & RowMap = original.RowMap();
  int nRows = RowMap.NumMyElements();
  const Epetra_BlockMap & ColMap = original.ColMap();
  int nCols = ColMap.NumMyElements();

  if( verbose_ ) cout << "RowMap:\n" << RowMap;
  if( verbose_ ) cout << "ColMap:\n" << ColMap;

  Epetra_CrsGraph * base = &( const_cast<Epetra_CrsGraph&>(original) );

  int NumIndices;
  int * Indices;

  // For parallel applications, add in boundaries to coloring
  bool distributedGraph = RowMap.DistributedGlobal();
  if( distributedGraph )
  {
    base = new Epetra_CrsGraph( Copy, ColMap, ColMap, 0 );

    for( int i = 0; i < nRows; ++i )
    {
      assert( original.ExtractMyRowView( i, NumIndices, Indices ) == 0 );
      assert( base->InsertMyIndices( i, NumIndices, Indices ) >= 0 );

      for( int j = 0; j < NumIndices; ++j )
        if( Indices[j] >= nRows )
          assert( base->InsertMyIndices( Indices[j], 1, &i ) >= 0 );
    } 

    base->TransformToLocal();
  }

  if( verbose_ ) cout << "Base Graph:\n" << *base << endl;

  //Generate Local Distance-1 Adjacency Graph
  //(Transpose of original excluding non-local column indices)
  Epetra_CrsGraph* Adj1 = CrsGraph_Transpose( true )( *base );
  if( verbose_ ) cout << "Adjacency 1 Graph!\n" << *Adj1;

  int Delta = Adj1->MaxNumIndices();
  cout << endl << "Delta: " << Delta << endl;

  //Generation of Local Distance-2 Adjacency Graph
  Epetra_CrsGraph Adj2( Copy, ColMap, ColMap, 0 );
  int NumAdj1Indices;
  int * Adj1Indices;
  for( int i = 0; i < nCols; ++i )
  {
    assert( Adj1->ExtractMyRowView( i, NumAdj1Indices, Adj1Indices ) == 0 );

    for( int j = 0; j < NumAdj1Indices; ++j )
    {
      assert( base->ExtractMyRowView( Adj1Indices[j], NumIndices, Indices ) == 0 );
      int NumLocalIndices = 0;
      for( int k = 0; k < NumIndices; ++k )
        if( Indices[k] < nCols ) NumLocalIndices++; 
      assert( Adj2.InsertMyIndices( i, NumLocalIndices, Indices ) >= 0 );
    }
  }
  assert( Adj2.TransformToLocal() == 0 );
  if( verbose_ ) cout << "Adjacency 2 Graph!\n" << Adj2;

  vector<int> rowOrder( nCols );
  //Simple reordering
  {
    multimap<int,int> adjMap;
    for( int i = 0; i < nCols; ++i )
      adjMap.insert( pair<const int,int>( Adj2.NumMyIndices(i), i ) );
    multimap<int,int>::iterator iter = adjMap.begin();
    multimap<int,int>::iterator end = adjMap.end();
    for( int i = 1; iter != end; ++iter, ++i )
      rowOrder[ nCols - i ] = iter->second;
  }

  Epetra_MapColoring * ColorMap = new Epetra_MapColoring( ColMap );

  //Application of Greedy Algorithm to generate Color Map
  int Size = Delta * Delta + 1;
  set<int> allowedColors;
  for( int col = 0; col < nCols; ++col )
  {
    for( int i = 0; i < Size; ++i ) allowedColors.insert( i+1 ); 

    Adj2.ExtractMyRowView( rowOrder[col], NumIndices, Indices );

    for( int i = 0; i < NumIndices; ++i )
      if( (*ColorMap)[ Indices[i] ] > 0 ) allowedColors.erase( (*ColorMap)[ Indices[i] ] );

    (*ColorMap)[ rowOrder[col] ] = *(allowedColors.begin());
  }

  if( verbose_ ) cout << "ColorMap!\n" << *ColorMap;

  if( distributedGraph ) delete base;

  return NewTypePtr( ColorMap );
}

}
