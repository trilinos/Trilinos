
#include <EDT_CrsGraph_Transpose.h>

#include <Epetra_Export.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>

#include <vector>

using std::vector;

namespace EpetraExt {

CrsGraph_Transpose::
~CrsGraph_Transpose()
{
  delete newObj_;
}

CrsGraph_Transpose::NewTypeRef
CrsGraph_Transpose::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  int err;

  int nRows = orig.NumMyRows();
  int nCols = orig.NumMyCols();

  const Epetra_BlockMap & RowMap = orig.RowMap();

  int numIndices;
  int * Indices;

  Epetra_CrsGraph * TransposeGraph = 0;

  if( !ignoreNonLocalCols_ && orig.DistributedGlobal() )
  {
    vector<int> TransNumNZ( nCols, 0 );
    for( int i = 0; i < nRows; ++i )
    {
      orig.ExtractMyRowView( i, numIndices, Indices );
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
      orig.ExtractMyRowView( i, numIndices, Indices );
      for( int j = 0; j < numIndices; ++j )
        TransIndices[ Indices[j] ][ TransNumNZ[ Indices[j] ]++ ] = i;
    }

    Epetra_CrsGraph SharedTransGraph( View, orig.ImportMap(), RowMap, &TransNumNZ[0] );
    for( int i = 0; i < nCols; ++i )
      if( TransNumNZ[i] ) SharedTransGraph.InsertMyIndices( i, TransNumNZ[i], &TransIndices[i][0] );
    SharedTransGraph.TransformToLocal();

    TransposeGraph = new Epetra_CrsGraph( Copy, RowMap, 0 );
    Epetra_Export Exporter( orig.ImportMap(), RowMap ); 
    TransposeGraph->Export( SharedTransGraph, Exporter, Add );
    TransposeGraph->TransformToLocal();
  }
  else
  {
    vector<int> TransNumNZ( nRows, 0 );
    for( int i = 0; i < nRows; ++i )
    {
      orig.ExtractMyRowView( i, numIndices, Indices );
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
      orig.ExtractMyRowView( i, numIndices, Indices );
      for( int j = 0; j < numIndices; ++j )
        if( Indices[j] < nRows ) TransIndices[ Indices[j] ][ TransNumNZ[ Indices[j] ]++ ] = i;
    }

    TransposeGraph = new Epetra_CrsGraph( Copy, RowMap, RowMap, &TransNumNZ[0] );

    for( int i = 0; i < nRows; ++i )
      if( TransNumNZ[i] ) TransposeGraph->InsertMyIndices( i, TransNumNZ[i], &TransIndices[i][0] );

    TransposeGraph->TransformToLocal();
  }

  newObj_ = TransposeGraph;

  return *TransposeGraph;
}

} // namespace EpetraExt

