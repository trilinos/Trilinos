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
    BaseGraph_( BaseGraph ),
    RowStencil_( vector< vector<int> >(1,RowStencil) ),
    RowIndices_( vector<int>(1,RowIndex) ),
    ROffset_(BlockUtility::CalculateOffset(BaseGraph.RowMap())),
    COffset_(BlockUtility::CalculateOffset(BaseGraph.ColMap()))
{
}

//==============================================================================
BlockCrsMatrix::BlockCrsMatrix(
        const Epetra_CrsGraph & BaseGraph,
        const vector< vector<int> > & RowStencil,
        const vector<int> & RowIndices,
        const Epetra_Comm & GlobalComm  ) 
  : Epetra_CrsMatrix( Copy, *(BlockUtility::GenerateBlockGraph( BaseGraph, RowStencil, RowIndices, GlobalComm )) ),
    BaseGraph_( BaseGraph ),
    RowStencil_( RowStencil ),
    RowIndices_( RowIndices ),
    ROffset_(BlockUtility::CalculateOffset(BaseGraph.RowMap())),
    COffset_(BlockUtility::CalculateOffset(BaseGraph.ColMap()))
{
}

//==============================================================================
BlockCrsMatrix::BlockCrsMatrix(
        const Epetra_CrsGraph & BaseGraph,
        const Epetra_CrsGraph & LocalBlockGraph,
        const Epetra_Comm & GlobalComm  ) 
  : Epetra_CrsMatrix( Copy, *(BlockUtility::GenerateBlockGraph( BaseGraph, LocalBlockGraph, GlobalComm )) ),
    BaseGraph_( BaseGraph ),
    RowStencil_( ),
    RowIndices_( ),
    ROffset_(BlockUtility::CalculateOffset(BaseGraph.RowMap())),
    COffset_(BlockUtility::CalculateOffset(BaseGraph.ColMap()))
{
  BlockUtility::GenerateRowStencil(LocalBlockGraph, RowIndices_, RowStencil_);
}

//==============================================================================
BlockCrsMatrix::BlockCrsMatrix(
        const Epetra_RowMatrix & BaseMatrix,
        const vector< vector<int> > & RowStencil,
        const vector<int> & RowIndices,
        const Epetra_Comm & GlobalComm  ) 
  : Epetra_CrsMatrix( Copy, *(BlockUtility::GenerateBlockGraph( BaseMatrix, RowStencil, RowIndices, GlobalComm )) ),
    BaseGraph_( Copy, BaseMatrix.RowMatrixRowMap(), 1 ), //Junk to satisfy constructor
    RowStencil_( RowStencil ),
    RowIndices_( RowIndices ),
    ROffset_(BlockUtility::CalculateOffset(BaseMatrix.RowMatrixRowMap())),
    COffset_(BlockUtility::CalculateOffset(BaseMatrix.RowMatrixColMap()))
{
}

//==============================================================================
BlockCrsMatrix::BlockCrsMatrix( const BlockCrsMatrix & Matrix ) 
  : Epetra_CrsMatrix( dynamic_cast<const Epetra_CrsMatrix &>( Matrix ) ),
    BaseGraph_( Matrix.BaseGraph_ ),
    RowStencil_( Matrix.RowStencil_ ),
    RowIndices_( Matrix.RowIndices_ ),
    ROffset_( Matrix.ROffset_ ),
    COffset_( Matrix.COffset_ )
{
}

//==============================================================================
BlockCrsMatrix::~BlockCrsMatrix()
{
}

//==============================================================================
void BlockCrsMatrix::LoadBlock(const Epetra_RowMatrix & BaseMatrix, const int Row, const int Col)
{
  int RowOffset = RowIndices_[Row] * ROffset_;
  int ColOffset = (RowIndices_[Row] + RowStencil_[Row][Col]) * COffset_;

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
  void BlockCrsMatrix::SumIntoBlock(double alpha, const Epetra_RowMatrix & BaseMatrix, const int Row, const int Col)
{
  int RowOffset = RowIndices_[Row] * ROffset_;
  int ColOffset = (RowIndices_[Row] + RowStencil_[Row][Col]) * COffset_;

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
    for( int l = 0; l < NumIndices; ++l ) {
       Indices[l] = ColOffset +  BaseColMap.GID(Indices[l]);
       Values[l] *= alpha;
    }

    int BaseRow = BaseMap.GID(i);
    ierr = this->SumIntoGlobalValues(BaseRow + RowOffset, NumIndices, &Values[0], &Indices[0]); 
    if (ierr != 0) cout << "WARNING BlockCrsMatrix::SumIntoBlock SumIntoGlobalValues err = " << ierr <<
	    "\n\t  Row " << BaseRow + RowOffset << "Col start" << Indices[0] << endl;

  }
}

//==============================================================================
  void BlockCrsMatrix::SumIntoGlobalBlock(double alpha, const Epetra_RowMatrix & BaseMatrix, const int Row, const int Col)
{
  int RowOffset = Row * ROffset_;
  int ColOffset = Col * COffset_;

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
    for( int l = 0; l < NumIndices; ++l ) {
       Indices[l] = ColOffset +  BaseColMap.GID(Indices[l]);
       Values[l] *= alpha;
    }

    int BaseRow = BaseMap.GID(i);
    ierr = this->SumIntoGlobalValues(BaseRow + RowOffset, NumIndices, &Values[0], &Indices[0]); 
    if (ierr != 0) cout << "WARNING BlockCrsMatrix::SumIntoBlock SumIntoGlobalValues err = " << ierr <<
	    "\n\t  Row " << BaseRow + RowOffset << "Col start" << Indices[0] << endl;

  }
}

//==============================================================================
void BlockCrsMatrix::BlockSumIntoGlobalValues(const int BaseRow, int NumIndices,
     double* Values, const int* Indices, const int Row, const int Col)
//All arguments could be const, except some were not set as const in CrsMatrix
{
  int RowOffset = RowIndices_[Row] * ROffset_;
  int ColOffset = (RowIndices_[Row] + RowStencil_[Row][Col]) * COffset_;

  // Convert to BlockMatrix Global numbering scheme
  vector<int> OffsetIndices(NumIndices);
  for( int l = 0; l < NumIndices; ++l ) OffsetIndices[l] = Indices[l] + ColOffset;

  int ierr = this->SumIntoGlobalValues(BaseRow + RowOffset, NumIndices,
                                   Values, &OffsetIndices[0]); 

  if (ierr != 0) cout << "WARNING BlockCrsMatrix::BlockSumIntoGlobalValues err = "
     << ierr << "\n\t  Row " << BaseRow + RowOffset << "Col start" << Indices[0] << endl;
}

//==============================================================================
void BlockCrsMatrix::BlockReplaceGlobalValues(const int BaseRow, int NumIndices,
     double* Values, const int* Indices, const int Row, const int Col)
//All arguments could be const, except some were not set as const in CrsMatrix
{
  int RowOffset = RowIndices_[Row] * ROffset_;
  int ColOffset = (RowIndices_[Row] + RowStencil_[Row][Col]) * COffset_;

  // Convert to BlockMatrix Global numbering scheme
  vector<int> OffsetIndices(NumIndices);
  for( int l = 0; l < NumIndices; ++l ) OffsetIndices[l] = Indices[l] + ColOffset;

  int ierr = this->ReplaceGlobalValues(BaseRow + RowOffset, NumIndices,
                                   Values, &OffsetIndices[0]); 

  if (ierr != 0) cout << "WARNING BlockCrsMatrix::BlockReplaceGlobalValues err = "
     << ierr << "\n\t  Row " << BaseRow + RowOffset << "Col start" << Indices[0] << endl;
}

//==============================================================================
void BlockCrsMatrix::BlockExtractGlobalRowView(const int BaseRow, 
					       int& NumEntries, 
					       double*& Values, 
					       const int Row, 
					       const int Col)
//All arguments could be const, except some were not set as const in CrsMatrix
{
  int RowOffset = RowIndices_[Row] * ROffset_;
  int ColOffset = (RowIndices_[Row] + RowStencil_[Row][Col]) * COffset_;

  // Get the whole row
  int ierr = this->ExtractGlobalRowView(BaseRow + RowOffset, NumEntries,
					Values); 

  // Adjust for just this block column
  Values += ColOffset;
  NumEntries -= ColOffset;

  if (ierr != 0) cout << "WARNING BlockCrsMatrix::BlockExtractGlobalRowView err = "
     << ierr << "\n\t  Row " << BaseRow + RowOffset << "Col " << Col+ColOffset << endl;
}

//==============================================================================
void BlockCrsMatrix::ExtractBlock(Epetra_CrsMatrix & BaseMatrix, const int Row, const int Col)
{
  int RowOffset = RowIndices_[Row] * ROffset_;
  int ColOffset = (RowIndices_[Row] + RowStencil_[Row][Col]) * COffset_;

//  const Epetra_CrsGraph & BaseGraph = BaseMatrix.Graph();
  const Epetra_BlockMap & BaseMap = BaseMatrix.RowMatrixRowMap();
  //const Epetra_BlockMap & BaseColMap = BaseMatrix.RowMatrixColMap();

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
       if (indx >= 0 && indx < COffset_) {
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
