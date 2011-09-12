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

#include <EpetraExt_SolverMap_CrsMatrix.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>

#include <vector>

namespace EpetraExt {

CrsMatrix_SolverMap::
~CrsMatrix_SolverMap()
{
  if( newObj_ && newObj_ != origObj_ ) delete newObj_;
  if( NewGraph_ ) delete NewGraph_;
  if( NewColMap_ ) delete NewColMap_;
}

CrsMatrix_SolverMap::NewTypeRef
CrsMatrix_SolverMap::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  assert( !orig.IndicesAreGlobal() );

  //test if matrix has missing local columns in its col std::map
  const Epetra_Map & RowMap = orig.RowMap();
  const Epetra_Map & DomainMap = orig.DomainMap();
  const Epetra_Map & ColMap = orig.ColMap();
  const Epetra_Comm & Comm = RowMap.Comm();
  int NumMyRows = RowMap.NumMyElements();
  int NumCols = DomainMap.NumMyElements();
  int Match = 0;
  for( int i = 0; i < NumCols; ++i )
    if( DomainMap.GID(i) != ColMap.GID(i) )
    {
      Match = 1;
      break;
    }

  int MatchAll = 0;
  Comm.SumAll( &Match, &MatchAll, 1 );

  if( !MatchAll )
  {
    newObj_ = origObj_;
  }
  else
  {
    //create ColMap with all local rows included
    std::vector<int> Cols(NumCols);
    //fill Cols list with GIDs of all local columns 
    for( int i = 0; i < NumCols; ++i )
      Cols[i] = DomainMap.GID(i);

    //now append to Cols any ghost column entries
    int NumMyCols = ColMap.NumMyElements();
    for( int i = 0; i < NumMyCols; ++i )
      if( !DomainMap.MyGID( ColMap.GID(i) ) ) Cols.push_back( ColMap.GID(i) );
    
    int NewNumMyCols = Cols.size();
    int NewNumGlobalCols;
    Comm.SumAll( &NewNumMyCols, &NewNumGlobalCols, 1 );
    //create new column std::map
    NewColMap_ = new Epetra_Map( NewNumGlobalCols, NewNumMyCols,&Cols[0], DomainMap.IndexBase(), Comm );

    //New Graph
    std::vector<int> NumIndicesPerRow( NumMyRows );
    for( int i = 0; i < NumMyRows; ++i )
      NumIndicesPerRow[i] = orig.NumMyEntries(i);
    NewGraph_ = new Epetra_CrsGraph( Copy, RowMap, *NewColMap_, &NumIndicesPerRow[0] );

    int MaxNumEntries = orig.MaxNumEntries();
    int NumEntries;
    std::vector<int> Indices( MaxNumEntries );
    for( int i = 0; i < NumMyRows; ++i )
    {
      int RowGID = RowMap.GID(i);
      orig.Graph().ExtractGlobalRowCopy( RowGID, MaxNumEntries, NumEntries, &Indices[0] );
      NewGraph_->InsertGlobalIndices( RowGID, NumEntries, &Indices[0] );
    }
    const Epetra_Map & RangeMap = orig.RangeMap();
    NewGraph_->FillComplete(DomainMap,RangeMap);

    //intial construction of matrix 
    Epetra_CrsMatrix * NewMatrix = new Epetra_CrsMatrix( View, *NewGraph_ );

    //insert views of row values
    int * myIndices;
    double * myValues;
    int indicesCnt;
    int numMyRows = NewMatrix->NumMyRows();
    for( int i = 0; i < numMyRows; ++i )
    {
      orig.ExtractMyRowView( i, indicesCnt, myValues, myIndices );
      NewGraph_->ExtractMyRowView( i, indicesCnt, myIndices );

      NewMatrix->InsertMyValues( i, indicesCnt, myValues, myIndices );
    }

    NewMatrix->FillComplete(DomainMap,RangeMap);

    newObj_ = NewMatrix;
  }

  return *newObj_;
}

} // namespace EpetraExt

