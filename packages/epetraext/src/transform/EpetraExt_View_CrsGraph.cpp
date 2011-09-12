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

#include <EpetraExt_View_CrsGraph.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_BlockMap.h>

#include <vector>

namespace EpetraExt {

CrsGraph_View::
~CrsGraph_View()
{
  if( newObj_ ) delete newObj_;
}

CrsGraph_View::NewTypeRef
CrsGraph_View::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  //Error, must be local indices
  assert( !orig.IndicesAreGlobal() );

  //test maps, new std::map must be left subset of old
  const Epetra_BlockMap & oRowMap = orig.RowMap();
  const Epetra_BlockMap & oColMap = orig.ColMap();

  int nNumRows = NewRowMap_->NumMyElements();
  int nNumCols = 0;
  if( NewColMap_ ) nNumCols = NewColMap_->NumMyElements();

  bool matched = true;
  for( int i = 0; i < nNumRows; ++i )
    matched = matched && ( oRowMap.GID(i) == NewRowMap_->GID(i) );
  if( nNumCols )
    for( int i = 0; i < nNumCols; ++i )
      matched = matched && ( oColMap.GID(i) == NewColMap_->GID(i) );

  if( !matched ) std::cout << "EDT_CrsGraph_View: Bad Row or Col Mapping\n";
  assert( matched );

  //intial construction of graph
  std::vector<int> numIndices( nNumRows );
  std::vector<int*> indices( nNumRows );
  for( int i = 0; i < nNumRows; ++i )
  {
    orig.ExtractMyRowView( i, numIndices[i], indices[i] );
    int j = 0;
    if( nNumCols )
    {
      while( j < numIndices[i] && NewColMap_->GID(indices[i][j]) != -1 ) ++j;
      numIndices[i] = j;
    }
  }

  Epetra_CrsGraph * newGraph( new Epetra_CrsGraph( View,
                                                   *NewRowMap_,
                                                   *NewColMap_,
                                                   &numIndices[0] ) );

  //insert views of row indices
  for( int i = 0; i < nNumRows; ++i )
    newGraph->InsertMyIndices( i, numIndices[i], indices[i] );

  newGraph->FillComplete();

  newObj_ = newGraph;

  return *newGraph;
}

} // namespace EpetraExt

