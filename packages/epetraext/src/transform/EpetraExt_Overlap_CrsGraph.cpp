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

#include <EpetraExt_Overlap_CrsGraph.h>

#include <Epetra_Import.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>

namespace EpetraExt {

CrsGraph_Overlap::
~CrsGraph_Overlap()
{
  if( newObj_ ) delete newObj_;

  if( OverlapMap_ ) delete OverlapMap_;
}

CrsGraph_Overlap::NewTypeRef
CrsGraph_Overlap::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  //check that this is a distributed graph and overlap level is not zero
  if( orig.DomainMap().DistributedGlobal() && levelOverlap_ )
  {
    Epetra_CrsGraph * OverlapGraph = new Epetra_CrsGraph( orig );
    OverlapMap_ = new Epetra_BlockMap( orig.RowMap() );

    Epetra_BlockMap * DomainMap = &(const_cast<Epetra_BlockMap&>(orig.DomainMap()));
    Epetra_BlockMap * RangeMap = &(const_cast<Epetra_BlockMap&>(orig.RangeMap()));

    for( int level = 0; level < levelOverlap_; ++level )
    {
      Epetra_BlockMap * OldRowMap = OverlapMap_;
      Epetra_CrsGraph * OldGraph = OverlapGraph;

      Epetra_Import & OverlapImporter = *(const_cast<Epetra_Import *>( OldGraph->Importer() ));
      OverlapMap_ = new Epetra_BlockMap( OverlapImporter.TargetMap() );

      //filter to local square block on last level if required
      if( squareLocalBlock_ && level==(levelOverlap_-1) )
        OverlapGraph = new Epetra_CrsGraph( Copy, *OverlapMap_, *OverlapMap_, 0 );
      else
        OverlapGraph = new Epetra_CrsGraph( Copy, *OverlapMap_, 0 );

      OverlapGraph->Import( *OldGraph, OverlapImporter, Insert );
      OverlapGraph->FillComplete( *DomainMap, *RangeMap );

      delete OldGraph;
      delete OldRowMap;
    }

    newObj_ = OverlapGraph;
  }
  else //just create a copy since this is not a InPlaceTransform
    newObj_ = new Epetra_CrsGraph( orig );

  return *newObj_;
}

} // namespace EpetraExt
