
#include <EDT_CrsGraph_Transpose.h>

#include <vector>

#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>

namespace Epetra_Transform {

std::auto_ptr<Epetra_CrsGraph> CrsGraph_Transpose::operator()( const Epetra_CrsGraph & original )
{
  int err;

  int nCols = original.NumMyCols();
  int nRows = original.NumMyRows();

  int numIndices;
  int * Indices;

  vector<int> TransNumNZ( nCols, 0 );
  for( int i = 0; i < nRows; ++i )
  {
    original.ExtractMyRowView( i, numIndices, Indices );
    for( int j = 0; j < numIndices; ++j )
      ++TransNumNZ[ Indices[j] ];
  }

  int ** TransIndices = new int*[nCols];
  for( int i = 0; i < nCols; ++i )
    if( TransNumNZ[i] )
    {
      TransIndices[i] = new int[ TransNumNZ[i] ];
      TransNumNZ[i] = 0;
    }
    else
      TransIndices[i] = 0;

  for( int i = 0; i < nRows; ++i )
  {
    original.ExtractMyRowView( i, numIndices, Indices );
    for( int j = 0; j < numIndices; ++j )
    {
      int TransRow = Indices[j];
      int TransLoc = TransNumNZ[TransRow]++;
      TransIndices[TransRow][TransLoc] = i;
    }
  }

  Epetra_CrsGraph SharedTransGraph( View, original.ImportMap(), 0 );
  for( int i = 0; i < nCols; ++i )
    SharedTransGraph.InsertMyIndices( i, TransNumNZ[i], TransIndices[i] );
  SharedTransGraph.TransformToLocal();

  Epetra_CrsGraph * TransposeGraph = new Epetra_CrsGraph( Copy, original.RowMap(), 0 );
  Epetra_Export Exporter( original.ImportMap(), original.RowMap() ); 
  TransposeGraph->Export( SharedTransGraph, Exporter, Add );
  TransposeGraph->TransformToLocal();

  for( int i = 0; i < nCols; ++i ) if( TransIndices[i] ) delete [] TransIndices[i];
  if( TransIndices ) delete [] TransIndices;

  return std::auto_ptr<Epetra_CrsGraph>( TransposeGraph );
}

}
