
#include <EDT_CrsGraph_MapColoringIndex.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_MapColoring.h>
#include <Epetra_IntVector.h>
#include <Epetra_Map.h>

#include <map>

using std::map;

EpetraExt::CrsGraph_MapColoringIndex::NewTypePtr EpetraExt::CrsGraph_MapColoringIndex::operator()( EpetraExt::CrsGraph_MapColoringIndex::OriginalTypeRef original )
{
  int err;

  const Epetra_BlockMap & RowMap = original.RowMap();
  int nRows = RowMap.NumMyElements();

  int NumColors = ColorMap_.NumColors();
  int * ListOfColors = ColorMap_.ListOfColors();

  map<int,int> MapOfColors;
  for( int i = 0; i < NumColors; ++i ) MapOfColors[ ListOfColors[i] ] = i;

  //initial setup of stl vector of IntVectors for indexing
  vector<int> dummy( nRows, -1 );
  NewTypePtr IndexVec( new vector<Epetra_IntVector>( NumColors, Epetra_IntVector( Copy, RowMap, &dummy[0] ) ) );

  int MaxNumIndices = original.MaxNumIndices();
  int NumIndices;
  vector<int> Indices( MaxNumIndices );

  for( int i = 0; i < nRows; ++i )
  {
    original.ExtractGlobalRowCopy( original.GRID(i), MaxNumIndices, NumIndices, &Indices[0] );

    for( int j = 0; j < NumIndices; ++j )
     (*IndexVec)[ MapOfColors[ColorMap_(Indices[j])] ][i] = Indices[j];
  }

  return IndexVec;
}

