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
