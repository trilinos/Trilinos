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

#include "EpetraExt_SubCopy_CrsMatrix.h"

#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_IntSerialDenseVector.h>

#include <vector>

namespace EpetraExt {

CrsMatrix_SubCopy::
~CrsMatrix_SubCopy()
{
  if( newObj_ ) delete newObj_;
}

CrsMatrix_SubCopy::NewTypeRef
CrsMatrix_SubCopy::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  //Error, must be local indices
  assert( orig.Filled() );

  //test maps, new map must be subset of old
  const Epetra_Map & oRowMap = orig.RowMap();
  const Epetra_Map & oColMap = orig.ColMap();

  int oNumRows = oRowMap.NumMyElements();
  (void) oNumRows; // Silence "unused variable" compiler warning.
  int oNumCols = oColMap.NumMyElements();
  int nNumRows = newRowMap_.NumMyElements();
  int nNumDomain = newDomainMap_.NumMyElements();

  bool matched = true;

  // Make sure all rows in newRowMap are already on this processor
  for( int i = 0; i < nNumRows; ++i )
    matched = matched && ( oRowMap.MyGID(newRowMap_.GID(i)) );
  if( !matched ) cerr << "EDT_CrsMatrix_SubCopy: Bad new_row_Map.  GIDs of new row map must be GIDs of the original row map on the same processor.\n";

  // Make sure all GIDs in the new domain map are GIDs in the old domain map
  if( !newRangeMap_.SameAs(newDomainMap_) ) {
    Epetra_IntSerialDenseVector pidList(nNumDomain);
    oColMap.RemoteIDList(newDomainMap_.NumMyElements(), newDomainMap_.MyGlobalElements(), pidList.Values(), 0);
    for( int i = 0; i < nNumDomain; ++i )
      matched = matched && ( pidList[i]>=0 );
  }

  if( !matched ) cout << "EDT_CrsMatrix_SubCopy: Bad newDomainMap.  One or more GIDs in new domain map are not part of original domain map.\n";
  assert( matched );


  // Next build new column map
  Epetra_IntSerialDenseVector pidList(oNumCols);
  Epetra_IntSerialDenseVector lidList(oNumCols);
  Epetra_IntSerialDenseVector sizeList(oNumCols);
  newDomainMap_.RemoteIDList(oColMap.NumMyElements(), oColMap.MyGlobalElements(), pidList.Values(), 0);
  int numNewCols = 0;
  Epetra_IntSerialDenseVector newColMapGidList(oNumCols);
  int * origColGidList = oColMap.MyGlobalElements();
  for( int i = 0; i < oNumCols; ++i )
    if (pidList[i] >=0)
      newColMapGidList[numNewCols++]= origColGidList[i];
  newColMap_ = Epetra_Map(-1, numNewCols, newColMapGidList.Values(), 0, oColMap.Comm());

  importer_ = new Epetra_Import(newRowMap_, oRowMap);

  Epetra_CrsMatrix * newMatrix = new Epetra_CrsMatrix(Copy, newRowMap_, newColMap_, 0);

  newObj_ = newMatrix;

  newObj_->Import(*origObj_, *importer_, Add);

  newObj_->FillComplete();

  return *newObj_;
}

//==============================================================================
bool CrsMatrix_SubCopy::fwd()
{

  if (newObj_->Filled()) newObj_->PutScalar(0.0); // zero contents

  newObj_->Import(*origObj_, *importer_, Add);

  newObj_->FillComplete();


  return (true);
}

//==============================================================================
bool CrsMatrix_SubCopy::rvs()
{
  if (!newObj_->Filled()) return(false); // Must have fillCompleted

  origObj_->Export(*newObj_, *importer_, Add);

  origObj_->FillComplete();

  return (true);
}

} // namespace EpetraExt

