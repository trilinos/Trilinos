#if 0 // dead source
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

#include <EpetraExt_Dirichlet_CrsMatrix.h>

#include <Epetra_CrsMatrix.h>
#include <Epetra_IntVector.h>
#include <Epetra_Import.h>
#include <Epetra_Map.h>

namespace EpetraExt {

bool
CrsMatrix_Dirichlet::
fwd()
{
  Epetra_CrsMatrix & Matrix = *origObj_;

  const Epetra_Map & RowMap = Matrix.RowMap();
  const Epetra_Map & ColMap = Matrix.ColMap();

  int NumMyElements = RowMap.NumMyElements();
  int NumMyColElements = ColMap.NumMyElements();

  if( symmetric_ && colSet_.empty() ) // get non-local column info
  {
    if( Matrix.IndicesAreGlobal() )
    {
      Epetra_Import Importer( ColMap, RowMap );
      Epetra_IntVector colLocations( ColMap );
      colLocations.Import( locations_, Importer, Insert );
      for( int i = 0; i < NumMyColElements; ++ i )
        if( colLocations[i] ) colSet_.insert(i);
    }
    else
    {
      for( int i = 0; i < NumMyElements; ++i )
        if( locations_[i] ) colSet_.insert(i);
    }
  }

  for( int i = 0; i < NumMyElements; ++i ) 
  {
    int * Indices;
    double * Vals;
    int NumIndices;
    if( locations_[i] ) //this is a Dirichlet BC location
    {
      Matrix.ExtractMyRowView( i, NumIndices, Vals, Indices );
      for( int j = 0; j < NumIndices; ++j )
      {
        if( Indices[j] == i ) Vals[i] = 1.0;
        else                  Vals[i] = 0.0;
      }
    }
    else if( symmetric_ )
    {
      Matrix.ExtractMyRowView( i, NumIndices, Vals, Indices );
      for( int j = 0; j < NumIndices; ++j )
        if( colSet_.count( Indices[j] ) ) Vals[j] = 0.0;
    }
  }

  return true;
}

bool
CrsMatrix_Dirichlet::
rvs()
{
  return true;
}

} //namespace EpetraExt

#endif // if 0
