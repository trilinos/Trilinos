//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include <EpetraExt_Transpose_CrsGraph.h>

#include <Epetra_Export.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>

#include <vector>

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

  int nRows = orig.NumMyRows();
  int nCols = orig.NumMyCols();

  const Epetra_BlockMap & RowMap = orig.RowMap();

  int numIndices;
  int * Indices;

  Epetra_CrsGraph * TransposeGraph = 0;

  if( !ignoreNonLocalCols_ && orig.DistributedGlobal() )
  {
    std::vector<int> TransNumNZ( nCols, 0 );
    for( int i = 0; i < nRows; ++i )
    {
      orig.ExtractMyRowView( i, numIndices, Indices );
      for( int j = 0; j < numIndices; ++j ) ++TransNumNZ[ Indices[j] ];
    }

    std::vector< std::vector<int> > TransIndices( nCols );
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
    SharedTransGraph.FillComplete();

    TransposeGraph = new Epetra_CrsGraph( Copy, RowMap, 0 );
    Epetra_Export Exporter( orig.ImportMap(), RowMap ); 
    TransposeGraph->Export( SharedTransGraph, Exporter, Add );
    TransposeGraph->FillComplete();
  }
  else
  {
    std::vector<int> TransNumNZ( nRows, 0 );
    for( int i = 0; i < nRows; ++i )
    {
      orig.ExtractMyRowView( i, numIndices, Indices );
      for( int j = 0; j < numIndices; ++j )
        if( Indices[j] < nRows ) ++TransNumNZ[ Indices[j] ];
    }

    std::vector< std::vector<int> > TransIndices( nRows );
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

    TransposeGraph->FillComplete();
  }

  newObj_ = TransposeGraph;

  return *TransposeGraph;
}

} // namespace EpetraExt

