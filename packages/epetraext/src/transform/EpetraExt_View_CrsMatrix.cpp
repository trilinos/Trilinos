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

#include <EpetraExt_View_CrsMatrix.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>

namespace EpetraExt {

CrsMatrix_View::
~CrsMatrix_View()
{
  if( newObj_ ) delete newObj_;
}

CrsMatrix_View::NewTypeRef
CrsMatrix_View::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  if( orig.IndicesAreGlobal() ) cout << "EDT_CrsMatrix_View: Indices must be LOCAL!\n";
  assert( !orig.IndicesAreGlobal() );

  //test graph, new graph must be contiguous subset of old

  //intial construction of matrix 
  Epetra_CrsMatrix * newMatrix( new Epetra_CrsMatrix( View, NewGraph_ ) );

  //insert views of row values
  int * myIndices;
  double * myValues;
  int indicesCnt;
  int numMyRows = newMatrix->NumMyRows();
  for( int i = 0; i < numMyRows; ++i )
  {
    orig.ExtractMyRowView( i, indicesCnt, myValues, myIndices );

    int newIndicesCnt = indicesCnt;
    bool done = false;
    for( int j = 0; j < indicesCnt; ++j )
      if( !done && NewGraph_.GCID( myIndices[j] ) == -1 )
      {
        newIndicesCnt = j;
        done = true;
      }

    newMatrix->InsertMyValues( i, newIndicesCnt, myValues, myIndices );
  }

  newMatrix->FillComplete(false);

  newObj_ = newMatrix;

  return *newMatrix;
}

} // namespace EpetraExt

