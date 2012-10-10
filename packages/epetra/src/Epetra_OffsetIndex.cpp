//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
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
// ************************************************************************
//@HEADER

#include "Epetra_ConfigDefs.h"
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

  if(!SourceGraph.RowMap().GlobalIndicesTypeMatch(TargetGraph.RowMap()))
     throw ReportError("Epetra_OffsetIndex::Epetra_OffsetIndex: SourceGraph and TargetGraph global indices type mismatch", -1);
  if(SourceGraph.RowMap().GlobalIndicesInt()) {
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
     GenerateLocalOffsets_<int>( SourceGraph, TargetGraph,
                            PermuteLIDs );

     GenerateRemoteOffsets_<int>( SourceGraph, TargetGraph,
                             ExportLIDs, RemoteLIDs,
                             Importer.Distributor() );
#else
    throw ReportError("Epetra_OffsetIndex::Epetra_OffsetIndex: ERROR, GlobalIndicesInt but no API for it.",-1);
#endif
  }
  else if(SourceGraph.RowMap().GlobalIndicesLongLong()) {
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
     GenerateLocalOffsets_<long long>( SourceGraph, TargetGraph,
                            PermuteLIDs );

     GenerateRemoteOffsets_<long long>( SourceGraph, TargetGraph,
                             ExportLIDs, RemoteLIDs,
                             Importer.Distributor() );
#else
    throw ReportError("Epetra_OffsetIndex::Epetra_OffsetIndex: ERROR, GlobalIndicesLongLong but no API for it.",-1);
#endif
  }
  else
     throw ReportError("Epetra_OffsetIndex::Epetra_OffsetIndex: SourceGraph global indices type unknown", -1);
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

  if(!SourceGraph.RowMap().GlobalIndicesTypeMatch(TargetGraph.RowMap()))
     throw ReportError("Epetra_OffsetIndex::Epetra_OffsetIndex: SourceGraph and TargetGraph global indices type mismatch", -1);
  if(SourceGraph.RowMap().GlobalIndicesInt()) {
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
     GenerateLocalOffsets_<int>( SourceGraph, TargetGraph,
                            PermuteLIDs );

     GenerateRemoteOffsets_<int>( SourceGraph, TargetGraph,
                             ExportLIDs, RemoteLIDs,
                             Exporter.Distributor() );
#else
    throw ReportError("Epetra_OffsetIndex::Epetra_OffsetIndex: ERROR, GlobalIndicesInt but no API for it.",-1);
#endif
  }
  else if(SourceGraph.RowMap().GlobalIndicesLongLong()) {
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
     GenerateLocalOffsets_<long long>( SourceGraph, TargetGraph,
                            PermuteLIDs );

     GenerateRemoteOffsets_<long long>( SourceGraph, TargetGraph,
                             ExportLIDs, RemoteLIDs,
                             Exporter.Distributor() );
#else
    throw ReportError("Epetra_OffsetIndex::Epetra_OffsetIndex: ERROR, GlobalIndicesLongLong but no API for it.",-1);
#endif
  }
  else
     throw ReportError("Epetra_OffsetIndex::Epetra_OffsetIndex: SourceGraph global indices type unknown", -1);
}

//==============================================================================
// Epetra_OffsetIndex copy constructor 
Epetra_OffsetIndex::Epetra_OffsetIndex(const Epetra_OffsetIndex& Indexor)
  : Epetra_Object(Indexor),
    NumSame_(Indexor.NumSame_),
    SameOffsets_(Indexor.SameOffsets_),
    NumPermute_(Indexor.NumPermute_),
    PermuteOffsets_(Indexor.PermuteOffsets_),
    NumExport_(0),
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
template<typename int_type>
void Epetra_OffsetIndex::GenerateLocalOffsets_( const Epetra_CrsGraph & SourceGraph,
                                                const Epetra_CrsGraph & TargetGraph,
                                                const int * PermuteLIDs )
{
  const int GlobalMaxNumSourceIndices = SourceGraph.GlobalMaxNumIndices();

  int NumSourceIndices;
  int_type * SourceIndices = 0;
  if( GlobalMaxNumSourceIndices>0 ) SourceIndices = new int_type[GlobalMaxNumSourceIndices];

  //setup Same Offsets
  SameOffsets_ = new int*[NumSame_];
  for( int i = 0; i < NumSame_; ++i ) SameOffsets_[i] = 0;

  for( int i = 0; i < NumSame_; ++i ) {
    int_type GID = (int_type) SourceGraph.GRID64(i);
    SourceGraph.ExtractGlobalRowCopy( GID,
                                      GlobalMaxNumSourceIndices,
                                      NumSourceIndices,
                                      SourceIndices );

    if( NumSourceIndices > 0 ) SameOffsets_[i] = new int[NumSourceIndices];

    int Loc = 0;
    int Start = 0;
    for( int j = 0; j < NumSourceIndices; ++j ) {
      Start = Loc;
      if( TargetGraph.FindGlobalIndexLoc(i,SourceIndices[j],Start,Loc) )
        SameOffsets_[i][j] = Loc;
      else
        SameOffsets_[i][j] = -1;
    }
  }

  //do also for permuted ids
  PermuteOffsets_ = new int*[NumPermute_];
  for( int i = 0; i < NumPermute_; ++i ) PermuteOffsets_[i] = 0;

  for( int i = 0; i < NumPermute_; ++i ) {
    int_type GID = (int_type) SourceGraph.GRID64(PermuteLIDs[i]);
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
template<typename int_type>
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
  /* "Indices" appears to be unused -- jhurani@txcorp.com
  int * Indices = 0;
  if( GlobalMaxNumIndices>0 ) Indices = new int[GlobalMaxNumIndices];
  */

  //Pack Source Rows
  int * Sizes = 0;
  if( NumExport_ > 0 ) Sizes = new int[NumExport_];
  int TotalSize = 0;
  for( int i = 0; i < NumExport_; ++i ) {
    Sizes[i] = SourceGraph.NumMyIndices(ExportLIDs[i]) + 1;
    TotalSize += Sizes[i];
  }

  int_type * SourceArray = new int_type[TotalSize+1];
  int Loc = 0;
  for( int i = 0; i < NumExport_; ++i ) {
    int_type GID = (int_type) SourceGraph.GRID64(ExportLIDs[i]);
    SourceArray[Loc] = Sizes[i]-1;
    SourceGraph.ExtractGlobalRowCopy( GID,
                                      GlobalMaxNumIndices,
                                      NumIndices,
                                      &(SourceArray[Loc+1]) );
    Loc += Sizes[i];
  }

  //Push to Target
  char * cRecvArray = 0;
  int_type * RecvArray = 0;
  int RecvArraySize = 0;

  Distor.Do( reinterpret_cast<char *>(SourceArray),
             (int)sizeof(int_type),
             Sizes,
             RecvArraySize,
             cRecvArray );
  RecvArray = reinterpret_cast<int_type*>(cRecvArray);

  //Construct RemoteOffsets
  if( NumRemote_ > 0 ) RemoteOffsets_ = new int*[NumRemote_];
  for( int i = 0; i < NumRemote_; ++i ) RemoteOffsets_[i] = 0;

  Loc = 0;
  for( int i = 0; i < NumRemote_; ++i ) {
    NumIndices = (int) RecvArray[Loc];
    RemoteOffsets_[i] = new int[NumIndices];
    ++Loc;
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

  /* "Indices" appears to be unused -- jhurani@txcorp.com
  if( GlobalMaxNumIndices>0 ) delete [] Indices;
  */

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

