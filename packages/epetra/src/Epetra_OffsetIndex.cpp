//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include "Epetra_OffsetIndex.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Distributor.h"
#include "Epetra_Comm.h"

//==============================================================================
// Epetra_OffsetIndex constructor from Importer
Epetra_OffsetIndex::Epetra_OffsetIndex( const Epetra_CrsGraph & SourceGraph,
                                        const Epetra_CrsGraph & TargetGraph,
                                        Epetra_Import & Importer )
  : Epetra_Object("Epetra::OffsetIndex"),
    NumSame(0),
    SameOffsets(0),
    NumPermute(0),
    PermuteOffsets(0),
    NumExport(0),
    NumRemote(0),
    RemoteOffsets(0)
    DataOwned(true)
{
  NumSame = Importer.NumSameIDs();
  int * SameLIDs = Importer.SameLIDs();
  NumPermute = Importer.NumPermuteIDs();
  int * PermuteLIDs = Importer.PermuteLIDs();

  NumExport = Importer.NumExportIDs();
  int * ExportLIDs = Importer.ExportLIDs();

  NumRemote = Importer.NumRemoteIDs();
  int * RemoteLIDs = Importer.RemoteLIDs();

  GenerateLocalOffsets_( SourceGraph, TargetGraph,
                         SameLIDs, PermuteLIDs );

  GenerateRemoteOffsets_( SourceGraph, TargetGraph,
                          ExportLIDs, RemoteLIDs,
                          Importer.Distributor() );
}

//==============================================================================
// Epetra_OffsetIndex constructor from Exporter
Epetra_OffsetIndex::Epetra_OffsetIndex( const Epetra_CrsGraph & SourceGraph,
                                        const Epetra_CrsGraph & TargetGraph,
                                        Epetra_Export & Exporter )
  : Epetra_Object("Epetra::OffsetIndex"),
    NumSame(0),
    SameOffsets(0),
    NumPermute(0),
    PermuteOffsets(0),
    NumRemote(0),
    RemoteOffsets(0)
    DataOwned(true)
{
  NumSame = Exporter.NumSameIDs();
  int * SameLIDs = Exporter.SameLIDs();
  NumPermute = Exporter.NumPermuteIDs();
  int * PermuteLIDs = Exporter.PermuteLIDs();

  NumExport = Exporter.NumExportIDs();
  int * ExportLIDs = Exporter.ExportLIDs();

  NumRemote = Exporter.NumRemoteIDs();
  int * RemoteLIDs = Exporter.RemoteLIDs();

  GenerateLocalOffsets_( SourceGraph, TargetGraph,
                         SameLIDs, PermuteLIDs );

  GenerateRemoteOffsets_( SourceGraph, TargetGraph,
                          ExportLIDs, RemoteLIDs,
                          Exporter.Distributor() );
}

//==============================================================================
// Epetra_OffsetIndex copy constructor 
Epetra_OffsetIndex::Epetra_OffsetIndex(const Epetra_OffsetIndex& Indexor)
  : Epetra_Object(Indexor),
    NumSame(Indexor.NumSame),
    SameOffsets(Indexor.SameOffsets),
    NumPermute(Indexor.NumPermute),
    PermuteOffsets(Indexor.PermuteOffsets),
    NumRemote(Indexor.NumRemote),
    RemoteOffsets(Indexor.RemoteOffsets),
    DataOwned(false)
{
}

//==============================================================================
Epetra_OffsetIndex::GenerateLocalOffsets_( const Epetra_CrsGraph & SourceGraph,
                                           const Epetra_CrsGraph & TargetGraph,
                                           const int * SameLIDs,
                                           const int * PermuteLIDs )
{
  int GlobalMaxNumSourceIndices = SourceGraph.GlobalMaxNumIndices();
  int GlobalMaxNumTargetIndices = TargetGraph.GlobalMaxNumIndices();

  int NumSourceIndices;
  int SourceIndices[GlobalMaxNumSourceIndices];

  //setup Same Offsets
  SameOffsets = new int*[NumSame];
  for( int i = 0; i < NumSame; ++i ) SameOffsets[i] = 0;

  for( int i = 0; i < NumSame; ++i ) {
    int GID = SourceGraph.GRID(SameLIDs[i]);
    SourceGraph.ExtractGlobalRowCopy( GID,
                                      GlobalMaxNumSourceIndices,
                                      NumSourceIndices,
                                      SourceIndices );

    if( NumSourceIndices > 0 ) SameOffsets[i] = new int[NumSourceIndices];

    int Loc = 0;
    int Start = 0;
    for( int j = 0; j < NumSourceIndices; ++j ) {
      Start = Loc;
      if( TargetGraph->FindGlobalIndexLoc(GID,SourceIndices[j],Start,Loc) )
        SameOffsets[i][j] = Loc;
      else
        SameOffsets[i][j] = -1;
    }
  }

  //do also for permuted ids
  PermuteOffsets = new int*[NumPermute];
  for( int i = 0; i < NumPermute; ++i ) PermuteOffsets[i] = 0;

  for( int i = 0; i < NumPermute; ++i ) {
    int GID = SourceGraph.GRID(PermuteLIDs[i]);
    SourceGraph.ExtractGlobalRowCopy( GID,
                                      GlobalMaxNumSourceIndices,
                                      NumSourceIndices,
                                      SourceIndices );

    if( NumSourceIndices > 0 ) PermuteOffsets[i] = new int[NumSourceIndices];

    int Loc = 0;
    int Start = 0;
    for( int j = 0; j < NumSourceIndices; ++j ) {
      Start = Loc;
      if( TargetGraph->FindGlobalIndexLoc(GID,SourceIndices[j],Start,Loc) )
        PermuteOffsets[i][j] = Loc;
      else
        PermuteOffsets[i][j] = -1;
    }
  }
}

//==============================================================================
Epetra_OffsetIndex::GenerateRemoteOffsets_( const Epetra_CrsGraph & SourceGraph,
                                            const Epetra_CrsGraph & TargetGraph,
                                            const int * ExportLIDs,
                                            const int * RemoteLIDs,
                                            Epetra_Distributor & Distor )
{
  int GlobalMaxNumIndices = SourceGraph.GlobalMaxNumIndices()

  int NumIndices;
  int Indices[GlobalMaxNumIndices];

  //Pack Source Rows
  int * Sizes = 0;
  if( NumExport > 0 ) Sizes = new int[NumExport];
  int TotalSize = 0;
  for( int i = 0; i < NumExport; ++i ) {
    Sizes[i] = SourceGraph.NumMyIndices(ExportLIDs[i]) + 1;
    TotalSize += Sizes[i];
  }

  int * SourceArray = new int[TotalSize+1];
  int Loc = 0;
  for( int i = 0; i < NumExport; ++i ) {
    int GID = SourceGraph.GRID(ExportLIDs[i]);
    SourceArray[Loc] = Sizes[i]-1;
    SourceGraph.ExtractGlobalRowCopy( GID,
                                      GlobalMaxNumIndices,
                                      NumIndices,
                                      SourceArray[Loc+1] );
    Loc += Sizes[i];
  }

  //Push to Target
  int * RecvArray = 0;
  Distor.Do( SourceArray, sizeof(int), Sizes, RecvArray );

  //Construct RemoteOffsets
  if( NumRemote > 0 ) RemoteOffsets = new int*[NumRemote];
  for( int i = 0; i < NumRemote; ++i ) RemoteOffsets[i] = 0;

  Loc = 0;
  for( int i = 0; i < NumRemote; ++i ) {
    NumIndices = RecvArray[Loc];
    ++Loc;
    int GID = TargetGraph.GRID(RemoteLIDs[i]);
    int FLoc = 0;
    int Start = 0;
    for( int j = 0; j < NumIndices; ++j ) {
      Start = FLoc;
      if( TargetGraph->FindGlobalIndexLoc(GID,RecvArray[Loc],Start,FLoc) )
        RemoteOffsets[i][j] = FLoc;
      else
        RemoteOffsets[i][j] = -1;
      ++Loc;
    }
  }

  if( Sizes ) delete [] Sizes;
  if( SourceArray ) delete [] SourceArray;
  if( RecvArray ) delete [] RecvArray;
}

//==============================================================================
// Epetra_Import destructor 
Epetra_Import::~Epetra_Import() {

  if( DataOwned )
  {
    for( int i = 0; i < NumSame; ++i )
      if( SameOffsets[i] ) delete [] SameOffsets[i];
    for( int i = 0; i < NumPermute; ++i )
      if( PermuteOffsets[i] ) delete [] PermuteOffsets[i];
    for( int i = 0; i < NumRemote; ++i )
      if( RemoteOffsets[i] ) delete [] RemoteOffsets[i];
  }
}

//=============================================================================
void Epetra_OffsetIndex::Print(ostream & os) const
{
  os << "Number of Same IDs = " << NumSame << endl;
  os << "Number of Permute IDs = " << NumPermute << endl;
  os << "Number of Remote IDs = " << NumRemote << endl;
  
  return;
}

