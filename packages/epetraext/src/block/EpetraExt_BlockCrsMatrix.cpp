//@HEADER
// ************************************************************************
// 
//               EpetraExt: Extended Linear Algebra Services Package 
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

#include "EpetraExt_BlockCrsMatrix.h"
#include "EpetraExt_BlockUtility.h"
#include "Epetra_Map.h"

namespace EpetraExt {

using std::vector;

//==============================================================================
BlockCrsMatrix::BlockCrsMatrix(
        const Epetra_CrsGraph & BaseGraph,
        const vector<int> & RowStencil,
        int RowIndex,
        const Epetra_Comm & GlobalComm  ) 
  : Epetra_CrsMatrix( Copy, BlockUtility::GenerateBlockGraph( BaseGraph, RowStencil, RowIndex, GlobalComm ) ),
    RowIndex_( RowIndex ),
    RowStencil_( RowStencil ),
    BaseGraph_( BaseGraph )
{
  AllocateBlocks_();
}

//==============================================================================
BlockCrsMatrix::BlockCrsMatrix( const BlockCrsMatrix & Matrix ) 
  : Epetra_CrsMatrix( dynamic_cast<const Epetra_CrsMatrix &>( Matrix ) ),
    RowIndex_( Matrix.RowIndex_ ),
    RowStencil_( Matrix.RowStencil_ ),
    BaseGraph_( Matrix.BaseGraph_ )
{
  AllocateBlocks_();
}

//==============================================================================
BlockCrsMatrix::~BlockCrsMatrix()
{
  DeleteBlocks_();
}

//==============================================================================
void BlockCrsMatrix::AllocateBlocks_()
{
  const Epetra_BlockMap & BaseRowMap = BaseGraph_.RowMap();
  const Epetra_BlockMap & RowMap = Graph().RowMap();

  int NumMyRows = Graph().NumMyRows();

  vector<int> BaseNumIndices( NumMyRows );
  vector<int*> BaseIndices( NumMyRows );

  for( int i = 0; i < NumMyRows; ++i )
    BaseGraph_.ExtractMyRowView( i, BaseNumIndices[i], BaseIndices[i] );

  vector<double*> Values( NumMyRows );
  vector<int> NumValues;

  for( int i = 0; i < NumMyRows; ++i )
    ExtractMyRowView( i, NumValues[i], Values[i] );

  int NumBlockCols = RowStencil_.size();
  Blocks_.resize( NumBlockCols );

  for( int i = 0; i < NumBlockCols; ++ i )
  {
    Epetra_CrsMatrix * bMat = new Epetra_CrsMatrix( View, BaseGraph_ );

    for( int j = 0; j < NumMyRows; ++j )
      bMat->InsertMyValues( j, BaseNumIndices[j], Values[j]+i*BaseNumIndices[j], BaseIndices[j] );

    Blocks_[i] = bMat;
  }
}

//==============================================================================
void BlockCrsMatrix::DeleteBlocks_()
{
  int NumBlockCols = RowStencil_.size();
  for( int i = 0; i < NumBlockCols; ++i )
    delete Blocks_[i];

  Blocks_.clear();
}

} //namespace EpetraExt
