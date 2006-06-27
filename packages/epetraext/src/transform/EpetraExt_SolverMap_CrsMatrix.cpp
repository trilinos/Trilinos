// @HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
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
// ***********************************************************************
// @HEADER

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

  //test if matrix has missing local columns in its col map
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
    vector<int> Cols(NumCols);
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
    //create new column map
    NewColMap_ = new Epetra_Map( NewNumGlobalCols, NewNumMyCols,&Cols[0], DomainMap.IndexBase(), Comm );

    //New Graph
    vector<int> NumIndicesPerRow( NumMyRows );
    for( int i = 0; i < NumMyRows; ++i )
      NumIndicesPerRow[i] = orig.NumMyEntries(i);
    NewGraph_ = new Epetra_CrsGraph( Copy, RowMap, *NewColMap_, &NumIndicesPerRow[0] );

    int MaxNumEntries = orig.MaxNumEntries();
    int NumEntries;
    vector<int> Indices( MaxNumEntries );
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

