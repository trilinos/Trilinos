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

#include <EpetraExt_MapColoringIndex.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_MapColoring.h>
#include <Epetra_IntVector.h>
#include <Epetra_Map.h>

#include <vector>
#include <map>

using std::vector;
using std::map;

namespace EpetraExt {

CrsGraph_MapColoringIndex::NewTypeRef
CrsGraph_MapColoringIndex::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  const Epetra_BlockMap & RowMap = orig.RowMap();
  int nRows = RowMap.NumMyElements();

  int NumColors = ColorMap_.NumColors();
  int * ListOfColors = ColorMap_.ListOfColors();

  map<int,int> MapOfColors;
  for( int i = 0; i < NumColors; ++i ) MapOfColors[ ListOfColors[i] ] = i;

  //initial setup of stl vector of IntVectors for indexing
  vector<int> dummy( nRows, -1 );
  NewTypePtr IndexVec = new NewType( NumColors, Epetra_IntVector( Copy, RowMap, &dummy[0] ) );

  int MaxNumIndices = orig.MaxNumIndices();
  int NumIndices;
  vector<int> Indices( MaxNumIndices );

  for( int i = 0; i < nRows; ++i )
  {
    orig.ExtractGlobalRowCopy( orig.GRID(i), MaxNumIndices, NumIndices, &Indices[0] );

    for( int j = 0; j < NumIndices; ++j )
     (*IndexVec)[ MapOfColors[ColorMap_(Indices[j])] ][i] = Indices[j];
  }

  newObj_ = IndexVec;

  return *IndexVec;
}

} // namespace EpetraExt

