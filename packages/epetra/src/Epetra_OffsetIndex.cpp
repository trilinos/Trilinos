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
    NumSame_(0),
    SameOffsets_(0),
    NumPermute_(0),
    PermuteOffsets_(0),
    NumExport_(0),
    NumRemote_(0),
    RemoteOffsets_(0),
    DataOwned_(true)
{
  NumSame_ = Importer.NumSameIDs();

  NumPermute_ = Importer.NumPermuteIDs();
  int * PermuteLIDs = Importer.PermuteToLIDs();

  NumExport_ = Importer.NumExportIDs();
  int * ExportLIDs = Importer.ExportLIDs();

  NumRemote_ = Importer.NumRemoteIDs();
  int * RemoteLIDs = Importer.RemoteLIDs();

  GenerateLocalOffsets_( SourceGraph, TargetGraph,
                         PermuteLIDs );

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
    NumSame_(0),
    SameOffsets_(0),
    NumPermute_(0),
    PermuteOffsets_(0),
    NumExport_(0),
    NumRemote_(0),
    RemoteOffsets_(0),
    DataOwned_(true)
{
  NumSame_ = Exporter.NumSameIDs();

  NumPermute_ = Exporter.NumPermuteIDs();
  int * PermuteLIDs = Exporter.PermuteToLIDs();

  NumExport_ = Exporter.NumExportIDs();
  int * ExportLIDs = Exporter.ExportLIDs();

  NumRemote_ = Exporter.NumRemoteIDs();
  int * RemoteLIDs = Exporter.RemoteLIDs();

  GenerateLocalOffsets_( SourceGraph, TargetGraph,
                         PermuteLIDs );

  GenerateRemoteOffsets_( SourceGraph, TargetGraph,
                          ExportLIDs, RemoteLIDs,
                          Exporter.Distributor() );
}

//==============================================================================
// Epetra_OffsetIndex copy constructor 
Epetra_OffsetIndex::Epetra_OffsetIndex(const Epetra_OffsetIndex& Indexor)
  : Epetra_Object(Indexor),
    NumSame_(Indexor.NumSame_),
    SameOffsets_(Indexor.SameOffsets_),
    NumPermute_(Indexor.NumPermute_),
    PermuteOffsets_(Indexor.PermuteOffsets_),
    NumRemote_(Indexor.NumRemote_),
    RemoteOffsets_(Indexor.RemoteOffsets_),
    DataOwned_(false)
{
}

//==============================================================================
// Epetra_OffsetIndex destructor
Epetra_OffsetIndex::~Epetra_OffsetIndex()
{
  if( DataOwned_ )
  {
    for( int i = 0; i < NumSame_; ++i )
      if( SameOffsets_[i] ) delete [] SameOffsets_[i];
    delete [] SameOffsets_;
    for( int i = 0; i < NumPermute_; ++i )
      if( PermuteOffsets_[i] ) delete [] PermuteOffsets_[i];
    delete [] PermuteOffsets_;
    for( int i = 0; i < NumRemote_; ++i )
      if( RemoteOffsets_[i] ) delete [] RemoteOffsets_[i];
    delete [] RemoteOffsets_;
  }
}

//==============================================================================
void Epetra_OffsetIndex::GenerateLocalOffsets_( const Epetra_CrsGraph & SourceGraph,
                                                const Epetra_CrsGraph & TargetGraph,
                                                const int * PermuteLIDs )
{
  const int GlobalMaxNumSourceIndices = SourceGraph.GlobalMaxNumIndices();
  const int GlobalMaxNumTargetIndices = TargetGraph.GlobalMaxNumIndices();

  int NumSourceIndices;
  int * SourceIndices = 0;
  if( GlobalMaxNumSourceIndices>0 ) SourceIndices = new int[GlobalMaxNumSourceIndices];

  //setup Same Offsets
  SameOffsets_ = new int*[NumSame_];
  for( int i = 0; i < NumSame_; ++i ) SameOffsets_[i] = 0;

  for( int i = 0; i < NumSame_; ++i ) {
    int GID = SourceGraph.GRID(i);
    SourceGraph.ExtractGlobalRowCopy( GID,
                                      GlobalMaxNumSourceIndices,
                                      NumSourceIndices,
                                      SourceIndices );

    if( NumSourceIndices > 0 ) SameOffsets_[i] = new int[NumSourceIndices];

    int Loc = 0;
    int Start = 0;
    for( int j = 0; j < NumSourceIndices; ++j ) {
      Start = Loc;
      if( TargetGraph.FindGlobalIndexLoc(GID,SourceIndices[j],Start,Loc) )
        SameOffsets_[i][j] = Loc;
      else
        SameOffsets_[i][j] = -1;
    }
  }

  //do also for permuted ids
  PermuteOffsets_ = new int*[NumPermute_];
  for( int i = 0; i < NumPermute_; ++i ) PermuteOffsets_[i] = 0;

  for( int i = 0; i < NumPermute_; ++i ) {
    int GID = SourceGraph.GRID(PermuteLIDs[i]);
    SourceGraph.ExtractGlobalRowCopy( GID,
                                      GlobalMaxNumSourceIndices,
                                      NumSourceIndices,
                                      SourceIndices );

    if( NumSourceIndices > 0 ) PermuteOffsets_[i] = new int[NumSourceIndices];

    int Loc = 0;
    int Start = 0;
    for( int j = 0; j < NumSourceIndices; ++j ) {
      Start = Loc;
      if( TargetGraph.FindGlobalIndexLoc(PermuteLIDs[i],SourceIndices[j],Start,Loc) )
        PermuteOffsets_[i][j] = Loc;
      else
        PermuteOffsets_[i][j] = -1;
    }
  }

  if( GlobalMaxNumSourceIndices>0 ) delete [] SourceIndices;
}

//==============================================================================
void Epetra_OffsetIndex::GenerateRemoteOffsets_( const Epetra_CrsGraph & SourceGraph,
                                                 const Epetra_CrsGraph & TargetGraph,
                                                 const int * ExportLIDs,
                                                 const int * RemoteLIDs,
                                                 Epetra_Distributor & Distor )
{
  int numProcs = SourceGraph.RowMap().Comm().NumProc();
  if (numProcs < 2) {
    return;
  }

  const int GlobalMaxNumIndices = SourceGraph.GlobalMaxNumIndices();

  int NumIndices;
  int * Indices = 0;
  if( GlobalMaxNumIndices>0 ) Indices = new int[GlobalMaxNumIndices];

  //Pack Source Rows
  int * Sizes = 0;
  if( NumExport_ > 0 ) Sizes = new int[NumExport_];
  int TotalSize = 0;
  for( int i = 0; i < NumExport_; ++i ) {
    Sizes[i] = SourceGraph.NumMyIndices(ExportLIDs[i]) + 1;
    TotalSize += Sizes[i];
  }

  int * SourceArray = new int[TotalSize+1];
  int Loc = 0;
  for( int i = 0; i < NumExport_; ++i ) {
    int GID = SourceGraph.GRID(ExportLIDs[i]);
    SourceArray[Loc] = Sizes[i]-1;
    SourceGraph.ExtractGlobalRowCopy( GID,
                                      GlobalMaxNumIndices,
                                      NumIndices,
                                      &(SourceArray[Loc+1]) );
    Loc += Sizes[i];
  }

  //Push to Target
  char * cRecvArray = 0;
  int * RecvArray = 0;
  int RecvArraySize = 0;
  Distor.Do( reinterpret_cast<char *>(SourceArray),
             sizeof(int),
             Sizes,
             RecvArraySize,
             cRecvArray );
  RecvArray = reinterpret_cast<int*>(cRecvArray);

  //Construct RemoteOffsets
  if( NumRemote_ > 0 ) RemoteOffsets_ = new int*[NumRemote_];
  for( int i = 0; i < NumRemote_; ++i ) RemoteOffsets_[i] = 0;

  Loc = 0;
  for( int i = 0; i < NumRemote_; ++i ) {
    NumIndices = RecvArray[Loc];
    RemoteOffsets_[i] = new int[NumIndices];
    ++Loc;
    int GID = TargetGraph.GRID(RemoteLIDs[i]);
    int FLoc = 0;
    int Start = 0;
    for( int j = 0; j < NumIndices; ++j ) {
      Start = FLoc;
      if( TargetGraph.FindGlobalIndexLoc(RemoteLIDs[i],RecvArray[Loc],Start,FLoc) )
        RemoteOffsets_[i][j] = FLoc;
      else
        RemoteOffsets_[i][j] = -1;
      ++Loc;
    }
  }

  if( GlobalMaxNumIndices>0 ) delete [] Indices;
  if( Sizes ) delete [] Sizes;
  if( SourceArray ) delete [] SourceArray;
  if( RecvArraySize ) delete [] cRecvArray;
}

//=============================================================================
void Epetra_OffsetIndex::Print(ostream & os) const
{
  os << "Number of Same IDs = " << NumSame_ << endl;
  os << "Number of Permute IDs = " << NumPermute_ << endl;
  os << "Number of Remote IDs = " << NumRemote_ << endl;
  
  return;
}

