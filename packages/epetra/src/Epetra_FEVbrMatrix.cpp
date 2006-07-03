
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
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

#include <Epetra_FEVbrMatrix.h>
#include <Epetra_BlockMap.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Export.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Comm.h>
#include <Epetra_Distributor.h>
#include <Epetra_Util.h>

//----------------------------------------------------------------------------
Epetra_FEVbrMatrix::Epetra_FEVbrMatrix(Epetra_DataAccess CV,
				       const Epetra_BlockMap& RowMap,
				       int *NumBlockEntriesPerRow,
				       bool ignoreNonLocalEntries) 
  : Epetra_VbrMatrix(CV, RowMap, NumBlockEntriesPerRow),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    numNonlocalBlockRows_(0),
    nonlocalBlockRows_(NULL),
    nonlocalBlockRowLengths_(NULL),
    nonlocalBlockRowAllocLengths_(NULL),
    nonlocalBlockCols_(NULL),
    nonlocalCoefs_(NULL),
    curRowOffset_(-1),
    curColOffset_(-1),
    curNumCols_(0),
    curCols_(NULL),
    curMode_(Add)
{
}

//----------------------------------------------------------------------------
Epetra_FEVbrMatrix::Epetra_FEVbrMatrix(Epetra_DataAccess CV,
				       const Epetra_BlockMap& RowMap,
				       int NumBlockEntriesPerRow,
				       bool ignoreNonLocalEntries) 
  : Epetra_VbrMatrix(CV, RowMap, NumBlockEntriesPerRow),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    numNonlocalBlockRows_(0),
    nonlocalBlockRows_(NULL),
    nonlocalBlockRowLengths_(NULL),
    nonlocalBlockRowAllocLengths_(NULL),
    nonlocalBlockCols_(NULL),
    nonlocalCoefs_(NULL),
    curRowOffset_(-1),
    curColOffset_(0),
    curNumCols_(0),
    curCols_(NULL),
    curMode_(Add)
{
}

//----------------------------------------------------------------------------
Epetra_FEVbrMatrix::Epetra_FEVbrMatrix(Epetra_DataAccess CV,
				       const Epetra_BlockMap& RowMap,
				       const Epetra_BlockMap& ColMap,
				       int *NumBlockEntriesPerRow,
				       bool ignoreNonLocalEntries) 
  : Epetra_VbrMatrix(CV, RowMap, ColMap, NumBlockEntriesPerRow),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    numNonlocalBlockRows_(0),
    nonlocalBlockRows_(NULL),
    nonlocalBlockRowLengths_(NULL),
    nonlocalBlockRowAllocLengths_(NULL),
    nonlocalBlockCols_(NULL),
    nonlocalCoefs_(NULL),
    curRowOffset_(-1),
    curColOffset_(-1),
    curNumCols_(0),
    curCols_(NULL),
    curMode_(Add)
{
}

//----------------------------------------------------------------------------
Epetra_FEVbrMatrix::Epetra_FEVbrMatrix(Epetra_DataAccess CV,
				       const Epetra_BlockMap& RowMap,
				       const Epetra_BlockMap& ColMap,
				       int NumBlockEntriesPerRow,
				       bool ignoreNonLocalEntries) 
  : Epetra_VbrMatrix(CV, RowMap, ColMap, NumBlockEntriesPerRow),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    numNonlocalBlockRows_(0),
    nonlocalBlockRows_(NULL),
    nonlocalBlockRowLengths_(NULL),
    nonlocalBlockRowAllocLengths_(NULL),
    nonlocalBlockCols_(NULL),
    nonlocalCoefs_(NULL),
    curRowOffset_(-1),
    curColOffset_(0),
    curNumCols_(0),
    curCols_(NULL),
    curMode_(Add)
{
}

//----------------------------------------------------------------------------
Epetra_FEVbrMatrix::Epetra_FEVbrMatrix(Epetra_DataAccess CV,
				       const Epetra_CrsGraph& Graph,
				       bool ignoreNonLocalEntries) 
  : Epetra_VbrMatrix(CV, Graph),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    numNonlocalBlockRows_(0),
    nonlocalBlockRows_(NULL),
    nonlocalBlockRowLengths_(NULL),
    nonlocalBlockRowAllocLengths_(NULL),
    nonlocalBlockCols_(NULL),
    nonlocalCoefs_(NULL),
    curRowOffset_(-1),
    curColOffset_(0),
    curNumCols_(0),
    curCols_(NULL),
    curMode_(Add)
{
}

//----------------------------------------------------------------------------
Epetra_FEVbrMatrix::Epetra_FEVbrMatrix(const Epetra_FEVbrMatrix& src)
  : Epetra_VbrMatrix(src),
    ignoreNonLocalEntries_(src.ignoreNonLocalEntries_),
    numNonlocalBlockRows_(0),
    nonlocalBlockRows_(NULL),
    nonlocalBlockRowLengths_(NULL),
    nonlocalBlockRowAllocLengths_(NULL),
    nonlocalBlockCols_(NULL),
    nonlocalCoefs_(NULL),
    curRowOffset_(-1),
    curColOffset_(0),
    curNumCols_(0),
    curCols_(NULL),
    curMode_(Add)
{
  operator=(src);
}

//----------------------------------------------------------------------------
Epetra_FEVbrMatrix& Epetra_FEVbrMatrix::operator=(const Epetra_FEVbrMatrix& src)
{
  if (this == &src) {
    return( *this );
  }

  Epetra_VbrMatrix::operator=(src);

  numNonlocalBlockRows_ = src.numNonlocalBlockRows_;

  nonlocalBlockRows_ = new int[numNonlocalBlockRows_];
  nonlocalBlockRowLengths_ = new int[numNonlocalBlockRows_];
  nonlocalBlockRowAllocLengths_ = new int[numNonlocalBlockRows_];
  nonlocalBlockCols_ = new int*[numNonlocalBlockRows_];
  nonlocalCoefs_ = new Epetra_SerialDenseMatrix**[numNonlocalBlockRows_];

  for(int i=0; i<numNonlocalBlockRows_; ++i) {
    nonlocalBlockRows_[i] = src.nonlocalBlockRows_[i];
    nonlocalBlockRowLengths_[i] = src.nonlocalBlockRowLengths_[i];
    nonlocalBlockRowAllocLengths_[i] = src.nonlocalBlockRowAllocLengths_[i];

    for(int j=0; j<nonlocalBlockRowLengths_[i]; ++j) {
      nonlocalBlockCols_[i][j] = src.nonlocalBlockCols_[i][j];

      nonlocalCoefs_[i][j] = new Epetra_SerialDenseMatrix(*(src.nonlocalCoefs_[i][j]));
    }
  }

  return( *this );
}

//----------------------------------------------------------------------------
Epetra_FEVbrMatrix::~Epetra_FEVbrMatrix()
{
  destroyNonlocalData();
}

//----------------------------------------------------------------------------
void Epetra_FEVbrMatrix::destroyNonlocalData()
{
  for(int i=0; i<numNonlocalBlockRows_; ++i) {
    delete [] nonlocalBlockCols_[i];
    for(int j=0; j<nonlocalBlockRowLengths_[i]; ++j) {
      delete nonlocalCoefs_[i][j];
    }
    delete [] nonlocalCoefs_[i];
  }

  delete [] nonlocalCoefs_;
  delete [] nonlocalBlockCols_;
  delete [] nonlocalBlockRowAllocLengths_;
  delete [] nonlocalBlockRowLengths_;
  delete [] nonlocalBlockRows_;

  numNonlocalBlockRows_ = 0;
  nonlocalBlockRows_ = NULL;
  nonlocalBlockRowLengths_ = NULL;
  nonlocalBlockRowAllocLengths_ = NULL;
  nonlocalBlockCols_ = NULL;
  nonlocalCoefs_ = NULL;
}

//----------------------------------------------------------------------------
int Epetra_FEVbrMatrix::PutScalar(double ScalarConstant) 
{
  for(int i=0; i<numNonlocalBlockRows_; ++i) {
    for(int j=0; j<nonlocalBlockRowLengths_[i]; ++j) {
      Epetra_SerialDenseMatrix& A = *(nonlocalCoefs_[i][j]);
      double* values = A.A();
      int lda = A.LDA();
      int n = A.N();
      for(int k=0; k<lda*n; ++k) {
	values[k] = ScalarConstant;
      }
    }
  }

  return( Epetra_VbrMatrix::PutScalar(ScalarConstant) );
}

//----------------------------------------------------------------------------
int Epetra_FEVbrMatrix::GlobalAssemble(bool callFillComplete) 
{
  if(Map().Comm().NumProc() < 2 || ignoreNonLocalEntries_) {
    if(callFillComplete) {
      EPETRA_CHK_ERR(FillComplete());
    }

    return(0);
  }

  int i;

  //In this method we need to gather all the non-local (overlapping) data
  //that's been input on each processor, into the
  //non-overlapping distribution defined by the map that 'this' matrix was
  //constructed with.

  //Need to build a map that describes our nonlocal data.

  //First, create a list of the sizes (point-rows per block-row) of the
  //nonlocal rows we're holding.
  int* pointRowsPerNonlocalBlockRow = numNonlocalBlockRows_>0 ?
    new int[numNonlocalBlockRows_] : NULL;

  for(i=0; i<numNonlocalBlockRows_; ++i) {
    pointRowsPerNonlocalBlockRow[i] = nonlocalCoefs_[i][0]->M();
  }

  //We'll use the arbitrary distribution constructor of BlockMap.

  Epetra_BlockMap sourceMap(-1, numNonlocalBlockRows_, nonlocalBlockRows_,
			    pointRowsPerNonlocalBlockRow,
			    RowMap().IndexBase(), RowMap().Comm());

  delete [] pointRowsPerNonlocalBlockRow;

  //If sourceMap has global size 0, then no nonlocal data exists and we can
  //skip most of this function.
  if(sourceMap.NumGlobalElements() < 1) {
    if(callFillComplete) {
      EPETRA_CHK_ERR(FillComplete());
    }
    return(0);
  }

  //We also need to build a column-map, containing the columns in our
  //nonlocal data. To do that, create a list of all column-indices that
  //occur in our nonlocal rows.

  int numCols = 0, allocLen = 0;
  int* cols = NULL;
  int* pointColsPerBlockCol = NULL;
  int ptColAllocLen = 0;
  int insertPoint = -1;

  for(i=0; i<numNonlocalBlockRows_; ++i) {
    for(int j=0; j<nonlocalBlockRowLengths_[i]; ++j) {
      int col = nonlocalBlockCols_[i][j];
      int offset = Epetra_Util_binary_search(col, cols, numCols, insertPoint);
      if (offset < 0) {
	EPETRA_CHK_ERR( Epetra_Util_insert(col, insertPoint, cols,
					   numCols, allocLen) );
	int tmpNumCols = numCols-1;
	EPETRA_CHK_ERR( Epetra_Util_insert(nonlocalCoefs_[i][j]->N(),
					   insertPoint,
					   pointColsPerBlockCol,
					   tmpNumCols, ptColAllocLen) );
      }
    }
  }

  Epetra_BlockMap colMap(-1, numCols, cols,
			 pointColsPerBlockCol,
			 RowMap().IndexBase(), RowMap().Comm());

  delete [] cols;
  delete [] pointColsPerBlockCol;
  numCols = 0;
  allocLen = 0;

  //now we need to create a matrix with sourceMap and colMap, and fill it with
  //our nonlocal data so we can then export it to the correct owning
  //processors.

  Epetra_VbrMatrix tempMat(Copy, sourceMap, colMap, nonlocalBlockRowLengths_);


  //Next we need to make sure the 'indices-are-global' attribute of tempMat's
  //graph is set to true, in case this processor doesn't end up calling the
  //InsertGlobalValues method...

  const Epetra_CrsGraph& graph = tempMat.Graph();
  Epetra_CrsGraph& nonconst_graph = const_cast<Epetra_CrsGraph&>(graph);
  nonconst_graph.SetIndicesAreGlobal(true);

  for(i=0; i<numNonlocalBlockRows_; ++i) {
    EPETRA_CHK_ERR( tempMat.BeginInsertGlobalValues(nonlocalBlockRows_[i],
						    nonlocalBlockRowLengths_[i],
						    nonlocalBlockCols_[i]) );

    for(int j=0; j<nonlocalBlockRowLengths_[i]; ++j) {
      Epetra_SerialDenseMatrix& subblock = *(nonlocalCoefs_[i][j]);

      EPETRA_CHK_ERR( tempMat.SubmitBlockEntry(subblock.A(),
					       subblock.LDA(),
					       subblock.M(),
					       subblock.N()) );
    }

    EPETRA_CHK_ERR( tempMat.EndSubmitEntries() );
  }

  //Now we need to call FillComplete on our temp matrix. We need to
  //pass a DomainMap and RangeMap, which are not the same as the RowMap
  //and ColMap that we constructed the matrix with.
  EPETRA_CHK_ERR(tempMat.FillComplete(RowMap(), sourceMap));

  //Finally, we're ready to create the exporter and export non-local data to
  //the appropriate owning processors.

  Epetra_Export exporter(sourceMap, RowMap());

  EPETRA_CHK_ERR( Export(tempMat, exporter, Add) );

  if(callFillComplete) {
    EPETRA_CHK_ERR(FillComplete());
  }

  destroyNonlocalData();

  return(0);
}

//----------------------------------------------------------------------------
int Epetra_FEVbrMatrix::InputNonlocalBlockEntry(double *Values, int LDA,
						int NumRows, int NumCols)
{
  if (curRowOffset_ < 0) {
    return(-1);
  }

  int insertPoint;
  int col = curCols_[curColOffset_++];
  int coloffset =
    Epetra_Util_binary_search(col, nonlocalBlockCols_[curRowOffset_],
			      nonlocalBlockRowLengths_[curRowOffset_],
			      insertPoint);
  if (coloffset < 0) return(-1);

  Epetra_SerialDenseMatrix*& subblock = nonlocalCoefs_[curRowOffset_][coloffset];

  if (subblock == NULL) {
    //For the construction of the serial dense matrix here, we choose Copy mode
    //in case the user deletes or reuses the Values array after this method is
    //called.
    subblock = new Epetra_SerialDenseMatrix(Copy, Values, LDA, NumRows, NumCols);

    if (subblock == NULL) {
      return(-1);
    }
  }
  else {
    int nrows = subblock->M();
    int ncols = subblock->N();
    int lda   = subblock->LDA();
    if (nrows != NumRows || ncols != NumCols || lda != LDA) {
      return(-1);
    }

    double* subblockvals = subblock->A();
    if (curMode_ == Add) {
      for(int i=0; i<lda*ncols; ++i) subblockvals[i] += Values[i];
    }
    else {
      for(int i=0; i<lda*ncols; ++i) subblockvals[i] = Values[i];
    }
  }

  return(0);
}

//----------------------------------------------------------------------------
int Epetra_FEVbrMatrix::InsertNonlocalRow(int row, int offset, int numCols)
{
  //insert a new row in our list of nonlocal rows.
  //also insert new arrays to hold block-column-index information

  int alloc_len = numNonlocalBlockRows_;
  EPETRA_CHK_ERR( Epetra_Util_insert(row, offset, nonlocalBlockRows_,
				     numNonlocalBlockRows_, alloc_len, 1) );

  int tmp1 = numNonlocalBlockRows_ - 1;
  int tmp2 = alloc_len - 1;

  EPETRA_CHK_ERR( Epetra_Util_insert(0, offset, nonlocalBlockRowLengths_,
				     tmp1, tmp2, 1) );

  --tmp1;
  --tmp2;
  int initialAllocLen = numCols*2;
  EPETRA_CHK_ERR( Epetra_Util_insert(initialAllocLen, offset,
				     nonlocalBlockRowAllocLengths_,
				     tmp1, tmp2, 1) );

  int** newCols = new int*[numNonlocalBlockRows_];
  Epetra_SerialDenseMatrix*** newCoefs =
    new Epetra_SerialDenseMatrix**[numNonlocalBlockRows_];

  if (newCols == NULL || newCoefs == NULL) {
    return(-1);
  }

  newCols[offset] = new int[initialAllocLen];
  newCoefs[offset] = new Epetra_SerialDenseMatrix*[initialAllocLen];

  for(int j=0; j<initialAllocLen; ++j) {
    newCols[offset][j] = 0;
    newCoefs[offset][j] = (Epetra_SerialDenseMatrix*)NULL;
  }

  int index = 0;
  for(int i=0; i<numNonlocalBlockRows_-1; ++i) {
    if (i == offset) {
      ++index;
    }

    newCols[index] = nonlocalBlockCols_[i];
    newCoefs[index++] = nonlocalCoefs_[i];
  }

  delete [] nonlocalBlockCols_;
  delete [] nonlocalCoefs_;

  nonlocalBlockCols_ = newCols;
  nonlocalCoefs_ = newCoefs;

  return(0);
}

//--------------------------------------------------------------------------
int Epetra_FEVbrMatrix::BeginInsertGlobalValues(int BlockRow,
						int NumBlockEntries,
						int * BlockIndices)
{
  int myRow = LRID(BlockRow);

  if (myRow > -1) {
    return( Epetra_VbrMatrix::BeginInsertGlobalValues(BlockRow,
						       NumBlockEntries,
						       BlockIndices) );
  }

  return( SetupForNonlocalSubmits(BlockRow, NumBlockEntries,
				  BlockIndices, false, Add) );
}

//--------------------------------------------------------------------------
int Epetra_FEVbrMatrix::BeginReplaceGlobalValues(int BlockRow,
						 int NumBlockEntries,
						 int *BlockIndices)
{
  int myRow = LRID(BlockRow);

  if (myRow > -1) {
    return( Epetra_VbrMatrix::BeginReplaceGlobalValues(BlockRow,
						       NumBlockEntries,
						       BlockIndices) );
  }

  return( SetupForNonlocalSubmits(BlockRow, NumBlockEntries,
				  BlockIndices, false, Insert) );
}

//--------------------------------------------------------------------------
int Epetra_FEVbrMatrix::BeginSumIntoGlobalValues(int BlockRow,
						 int NumBlockEntries,
						 int *BlockIndices)
{
  int myRow = LRID(BlockRow);

  if (myRow > -1) {
    return( Epetra_VbrMatrix::BeginSumIntoGlobalValues(BlockRow,
						       NumBlockEntries,
						       BlockIndices) );
  }

  return( SetupForNonlocalSubmits(BlockRow, NumBlockEntries,
				  BlockIndices, false, Add) );
}

//--------------------------------------------------------------------------
int Epetra_FEVbrMatrix::SetupForNonlocalSubmits(int BlockRow,
						int NumBlockEntries,
						int * BlockIndices, 
						bool IndicesAreLocal,
						Epetra_CombineMode SubmitMode)
{
  (void)IndicesAreLocal;
  if (ignoreNonLocalEntries_) {
    curRowOffset_ = 0;
    return(0);
  }

  int insertPoint = -1;

  //find offset of this row in our list of nonlocal rows
  int rowoffset = Epetra_Util_binary_search(BlockRow, nonlocalBlockRows_,
					    numNonlocalBlockRows_, insertPoint);

  //if this row is not already present, insert it
  if (rowoffset < 0) {
    EPETRA_CHK_ERR( InsertNonlocalRow(BlockRow, insertPoint, NumBlockEntries) );
    rowoffset = insertPoint;
  }

  //now insert each incoming block-column-index in this list of column-indices,
  //maintaining sortedness.
  for(int i=0; i<NumBlockEntries; ++i) {
    int col = BlockIndices[i];
    int coloffset = Epetra_Util_binary_search(col, nonlocalBlockCols_[rowoffset],
					     nonlocalBlockRowLengths_[rowoffset],
					      insertPoint);
    if (coloffset < 0) {
      int tmp1 = nonlocalBlockRowLengths_[rowoffset];
      int tmp2 = nonlocalBlockRowAllocLengths_[rowoffset];

      EPETRA_CHK_ERR( Epetra_Util_insert(col, insertPoint,
					 nonlocalBlockCols_[rowoffset],
				  nonlocalBlockRowLengths_[rowoffset],
	       			 nonlocalBlockRowAllocLengths_[rowoffset]));

      EPETRA_CHK_ERR( Epetra_Util_insert((Epetra_SerialDenseMatrix*)NULL,
					 insertPoint,
					 nonlocalCoefs_[rowoffset],
					 tmp1, tmp2) );
    }
  }

  curRowOffset_ = rowoffset;
  curColOffset_ = 0;
  curNumCols_ = NumBlockEntries;
  curCols_ = new int[NumBlockEntries];
  for(int j=0; j<NumBlockEntries; ++j) {
    curCols_[j] = BlockIndices[j];
  }

  curMode_ = SubmitMode;

  return(0);
}

//--------------------------------------------------------------------------
int Epetra_FEVbrMatrix::SubmitBlockEntry(double *Values,
					 int LDA,
					 int NumRows,
					 int NumCols)
{
  if (curRowOffset_ < 0) {
    EPETRA_CHK_ERR( Epetra_VbrMatrix::SubmitBlockEntry(Values, LDA,
							NumRows, NumCols) );
  }
  else {
    if (ignoreNonLocalEntries_) {
      return(0);
    }

    EPETRA_CHK_ERR( InputNonlocalBlockEntry(Values, LDA,
					    NumRows, NumCols) );
  }

  return(0);
}
//--------------------------------------------------------------------------
int Epetra_FEVbrMatrix::EndSubmitEntries()
{
  if (curRowOffset_ < 0) {
    EPETRA_CHK_ERR( Epetra_VbrMatrix::EndSubmitEntries() );
  }
  else {
    curRowOffset_ = -1;
    curColOffset_ = -1;
    curNumCols_ = 0;
    delete [] curCols_;
  }

  return(0);
}

