//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef HAVE_NOX_EPETRAEXT

// Use local version of needed coloring file
#include <nox_EDT_CrsGraph_Transpose.H>

#include <vector>

#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>

namespace Epetra_Transform {

Epetra_CrsGraph* CrsGraph_Transpose::operator()( const Epetra_CrsGraph & original )
{
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

    Epetra_CrsGraph SharedTransGraph( View, original.ImportMap(), RowMap, &TransNumNZ[0] );
    for( int i = 0; i < nCols; ++i )
      if( TransNumNZ[i] ) SharedTransGraph.InsertMyIndices( i, TransNumNZ[i], &TransIndices[i][0] );
    SharedTransGraph.FillComplete();

    TransposeGraph = new Epetra_CrsGraph( Copy, RowMap, 0 );
    Epetra_Export Exporter( original.ImportMap(), RowMap ); 
    TransposeGraph->Export( SharedTransGraph, Exporter, Add );
    TransposeGraph->FillComplete();
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

    TransposeGraph = new Epetra_CrsGraph( Copy, RowMap, RowMap, &TransNumNZ[0] );

    for( int i = 0; i < nRows; ++i )
      if( TransNumNZ[i] ) TransposeGraph->InsertMyIndices( i, TransNumNZ[i], &TransIndices[i][0] );

    TransposeGraph->FillComplete();
  }

  return TransposeGraph;
}

}
#endif
