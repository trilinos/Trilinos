
#include <EDT_CrsGraph_Transpose.h>

#include <vector>

#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>

namespace Epetra_Transform {

std::auto_ptr<Epetra_CrsGraph> CrsGraph_Transpose::operator()( const Epetra_CrsGraph & original )
{
  int err;

  int nRows = original.NumMyRows();
  int nCols = original.NumMyCols();

  const Epetra_BlockMap & RowMap = original.RowMap();

  int numIndices;
  int * Indices;

  Epetra_CrsGraph * TransposeGraph = 0;

  if( !ignoreNonLocalCols_ && original.DistributedGlobal() )
  {
    vector<int> TransNumNZ( nCols, 0 );
    for( int i = 0; i < nRows; ++i )
    {
      original.ExtractMyRowView( i, numIndices, Indices );
      for( int j = 0; j < numIndices; ++j ) ++TransNumNZ[ Indices[j] ];
    }

    vector< vector<int> > TransIndices( nCols );
    for( int i = 0; i < nCols; ++i )
      if( TransNumNZ[i] )
      {
        TransIndices[i].resize( TransNumNZ[i] );
        TransNumNZ[i] = 0;
      }

    for( int i = 0; i < nRows; ++i )
    {
      original.ExtractMyRowView( i, numIndices, Indices );
      for( int j = 0; j < numIndices; ++j )
        TransIndices[ Indices[j] ][ TransNumNZ[ Indices[j] ]++ ] = i;
    }

    Epetra_CrsGraph SharedTransGraph( View, original.ImportMap(), &TransNumNZ[0] );
    for( int i = 0; i < nCols; ++i )
      if( TransNumNZ[i] ) SharedTransGraph.InsertMyIndices( i, TransNumNZ[i], &TransIndices[i][0] );
    SharedTransGraph.TransformToLocal();

    TransposeGraph = new Epetra_CrsGraph( Copy, RowMap, 0 );
    Epetra_Export Exporter( original.ImportMap(), RowMap ); 
    TransposeGraph->Export( SharedTransGraph, Exporter, Add );
    TransposeGraph->TransformToLocal();
  }
  else
  {
    vector<int> TransNumNZ( nRows, 0 );
    for( int i = 0; i < nRows; ++i )
    {
      original.ExtractMyRowView( i, numIndices, Indices );
      for( int j = 0; j < numIndices; ++j )
        if( Indices[j] < nRows ) ++TransNumNZ[ Indices[j] ];
    }

    vector< vector<int> > TransIndices( nRows );
    for( int i = 0; i < nRows; ++i )
      if( TransNumNZ[i] )
      {
        TransIndices[i].resize( TransNumNZ[i] );
        TransNumNZ[i] = 0;
      }

    for( int i = 0; i < nRows; ++i )
    {
      original.ExtractMyRowView( i, numIndices, Indices );
      for( int j = 0; j < numIndices; ++j )
        if( Indices[j] < nRows ) TransIndices[ Indices[j] ][ TransNumNZ[ Indices[j] ]++ ] = i;
    }

    TransposeGraph = new Epetra_CrsGraph( Copy, RowMap, &TransNumNZ[0] );

    for( int i = 0; i < nRows; ++i )
      if( TransNumNZ[i] ) TransposeGraph->InsertMyIndices( i, TransNumNZ[i], &TransIndices[i][0] );

    TransposeGraph->TransformToLocal();
  }

  return std::auto_ptr<Epetra_CrsGraph>( TransposeGraph );
}

}
