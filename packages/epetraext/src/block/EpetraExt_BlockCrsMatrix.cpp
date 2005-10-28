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
  : Epetra_CrsMatrix( Copy, *(BlockUtility::GenerateBlockGraph( BaseGraph, vector< vector<int> >(1,RowStencil), vector<int>(1,RowIndex), GlobalComm )) ),
    RowIndices_( vector<int>(1,RowIndex) ),
    RowStencil_( vector< vector<int> >(1,RowStencil) ),
    BaseGraph_( BaseGraph ),
    Offset_(BlockUtility::CalculateOffset(BaseGraph.RowMap()))
{
  AllocateBlocks_();
}

//==============================================================================
BlockCrsMatrix::BlockCrsMatrix(
        const Epetra_CrsGraph & BaseGraph,
        const vector< vector<int> > & RowStencil,
        const vector<int> & RowIndices,
        const Epetra_Comm & GlobalComm  ) 
  : Epetra_CrsMatrix( Copy, *(BlockUtility::GenerateBlockGraph( BaseGraph, RowStencil, RowIndices, GlobalComm )) ),
    RowIndices_( RowIndices ),
    RowStencil_( RowStencil ),
    BaseGraph_( BaseGraph ),
    Offset_(BlockUtility::CalculateOffset(BaseGraph.RowMap()))
{
  AllocateBlocks_();
}

//==============================================================================
BlockCrsMatrix::BlockCrsMatrix(
        const Epetra_RowMatrix & BaseMatrix,
        const vector< vector<int> > & RowStencil,
        const vector<int> & RowIndices,
        const Epetra_Comm & GlobalComm  ) 
  : Epetra_CrsMatrix( Copy, *(BlockUtility::GenerateBlockGraph( BaseMatrix, RowStencil, RowIndices, GlobalComm )) ),
    RowIndices_( RowIndices ),
    RowStencil_( RowStencil ),
    BaseGraph_( Copy, BaseMatrix.RowMatrixRowMap(), 1 ), //Junk to satisfy constructor
    Offset_(BlockUtility::CalculateOffset(BaseMatrix.RowMatrixRowMap()))
{
   AllocateBlocks_(); // This is now garbage because bogus baseGraph
}

//==============================================================================
BlockCrsMatrix::BlockCrsMatrix( const BlockCrsMatrix & Matrix ) 
  : Epetra_CrsMatrix( dynamic_cast<const Epetra_CrsMatrix &>( Matrix ) ),
    RowIndices_( Matrix.RowIndices_ ),
    RowStencil_( Matrix.RowStencil_ ),
    BaseGraph_( Matrix.BaseGraph_ ),
    Offset_( Matrix.Offset_ )
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

  int NumBaseRows = BaseGraph_.NumMyRows();
  int NumBlockRows = RowIndices_.size();
  int NumMyRows = Graph().NumMyRows();

  vector<int> BaseNumIndices( NumBaseRows );
  vector<int*> BaseIndices( NumBaseRows );

  for( int i = 0; i < NumBaseRows; ++i )
    BaseGraph_.ExtractMyRowView( i, BaseNumIndices[i], BaseIndices[i] );

  vector<double*> Values( NumMyRows );
  vector<int> NumValues( NumMyRows );

  Blocks_.resize( NumBlockRows );

  for( int i = 0; i < NumBlockRows; ++i )
  {
    for( int j = 0; j < NumBaseRows; ++j )
      ExtractMyRowView( i*NumBaseRows+j, NumValues[j], Values[j] );

    int NumBlockCols = RowStencil_[i].size();
    Blocks_[i].resize( NumBlockCols );

    for( int j = 0; j < NumBlockCols; ++j )
    {
      Epetra_CrsMatrix * bMat = new Epetra_CrsMatrix( View, BaseGraph_ );

      for( int k = 0; k < NumBaseRows; ++k )
        bMat->InsertMyValues( k, BaseNumIndices[k], Values[k]+j*BaseNumIndices[k], BaseIndices[k] );

      Blocks_[i][j] = bMat;
    }
  }
}

//==============================================================================
void BlockCrsMatrix::DeleteBlocks_()
{
  for( int i = 0; i < RowIndices_.size(); ++i )
  {
    int NumBlockCols = RowStencil_[i].size();
    for( int j = 0; j < NumBlockCols; ++j )
      delete Blocks_[i][j];
  }

  Blocks_.clear();
}

//==============================================================================
void BlockCrsMatrix::LoadBlock(const Epetra_RowMatrix & BaseMatrix, const int Row, const int Col)
{
  int RowOffset = RowIndices_[Row] * Offset_;
  int ColOffset = (RowIndices_[Row] + RowStencil_[Row][Col]) * Offset_;

//  const Epetra_CrsGraph & BaseGraph = BaseMatrix.Graph();
  const Epetra_BlockMap & BaseMap = BaseMatrix.RowMatrixRowMap();
  const Epetra_BlockMap & BaseColMap = BaseMatrix.RowMatrixColMap();

  // This routine copies entries of a BaseMatrix into big  BlockCrsMatrix
  // It performs the following operation on the global IDs row-by-row
  // this->val[i+rowOffset][j+ColOffset] = BaseMatrix.val[i][j]

  int MaxIndices = BaseMatrix.MaxNumEntries();
  vector<int> Indices(MaxIndices);
  vector<double> Values(MaxIndices);
  int NumIndices;
  int ierr=0;

  for (int i=0; i<BaseMap.NumMyElements(); i++) {
    BaseMatrix.ExtractMyRowCopy( i, MaxIndices, NumIndices, &Values[0], &Indices[0] );

    // Convert to BlockMatrix Global numbering scheme
    for( int l = 0; l < NumIndices; ++l )
       Indices[l] = ColOffset +  BaseColMap.GID(Indices[l]);

    int BaseRow = BaseMap.GID(i);
    ierr = this->ReplaceGlobalValues(BaseRow + RowOffset, NumIndices, &Values[0], &Indices[0]); 
    if (ierr != 0) cout << "WARNING BlockCrsMatrix::LoadBlock ReplaceGlobalValues err = " << ierr <<
	    "\n\t  Row " << BaseRow + RowOffset << "Col start" << Indices[0] << endl;

  }
}

//==============================================================================
void BlockCrsMatrix::ExtractBlock(Epetra_CrsMatrix & BaseMatrix, const int Row, const int Col)
{
  int RowOffset = RowIndices_[Row] * Offset_;
  int ColOffset = (RowIndices_[Row] + RowStencil_[Row][Col]) * Offset_;

//  const Epetra_CrsGraph & BaseGraph = BaseMatrix.Graph();
  const Epetra_BlockMap & BaseMap = BaseMatrix.RowMatrixRowMap();
  const Epetra_BlockMap & BaseColMap = BaseMatrix.RowMatrixColMap();

  // This routine extracts entries of a BaseMatrix from a big  BlockCrsMatrix
  // It performs the following operation on the global IDs row-by-row
  // BaseMatrix.val[i][j] = this->val[i+rowOffset][j+ColOffset] 

  int MaxIndices = BaseMatrix.MaxNumEntries();
  vector<int> Indices(MaxIndices);
  vector<double> Values(MaxIndices);
  int NumIndices;
  int indx,icol;
  double* BlkValues;
  int *BlkIndices;
  int BlkNumIndices;
  int ierr=0;

  for (int i=0; i<BaseMap.NumMyElements(); i++) {

    // Get pointers to values and indices of whole block matrix row
    int BaseRow = BaseMap.GID(i);
    int myBlkBaseRow = this->RowMatrixRowMap().LID(BaseRow + RowOffset);
    ierr = this->ExtractMyRowView(myBlkBaseRow, BlkNumIndices, BlkValues, BlkIndices); 

    NumIndices = 0;
    // Grab columns with global indices in correct range for this block
    for( int l = 0; l < BlkNumIndices; ++l ) {
       icol = this->RowMatrixColMap().GID(BlkIndices[l]);
       indx = icol - ColOffset;
       if (indx >= 0 && indx < Offset_) {
         Indices[NumIndices] = indx;
         Values[NumIndices] = BlkValues[l];
	 NumIndices++;
       }
    }

    //Load this row into base matrix
    BaseMatrix.ReplaceGlobalValues(BaseRow, NumIndices, &Values[0], &Indices[0] );

  }
}

} //namespace EpetraExt
