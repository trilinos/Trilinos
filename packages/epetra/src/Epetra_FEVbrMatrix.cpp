
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

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
    curNumCols_(0)
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
    curNumCols_(0)
{
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
int Epetra_FEVbrMatrix::GlobalAssemble(bool callTransformToLocal) 
{
  if (Map().Comm().NumProc() < 2) {
    return(0);
  }

  if (ignoreNonLocalEntries_) {
    return(0);
  }

  int i, j;

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
  if (sourceMap.NumGlobalElements() < 1) {
    if (callTransformToLocal) {
      EPETRA_CHK_ERR( TransformToLocal() );
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
    for(j=0; j<nonlocalBlockRowLengths_[i]; ++j) {
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

  //Now we need to call TransformToLocal on our temp matrix. We need to
  //pass a DomainMap and RangeMap, which are not the same as the RowMap
  //and ColMap that we constructed the matrix with.
  EPETRA_CHK_ERR( tempMat.TransformToLocal(&(RowMap()), &sourceMap ) );

  Epetra_Export exporter(sourceMap, RowMap());

  EPETRA_CHK_ERR( Export(tempMat, exporter, Add) );

  if (callTransformToLocal) {
    EPETRA_CHK_ERR( TransformToLocal() );
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

  if (curColOffset_ >= nonlocalBlockRowAllocLengths_[curRowOffset_]) {
    return(-1);
  }

  int col = curColOffset_++;

  //For the construction of the serial dense matrix here, we choose Copy mode
  //in case the user deletes or reuses the Values array after this method is
  //called.
  nonlocalCoefs_[curRowOffset_][col] = 
    new Epetra_SerialDenseMatrix(Copy, Values, LDA, NumRows, NumCols);

  if (nonlocalCoefs_[curRowOffset_][col] == NULL) {
    return(-1);
  }

  nonlocalBlockRowLengths_[curRowOffset_] = curColOffset_;

  return(0);
}

//----------------------------------------------------------------------------
int Epetra_FEVbrMatrix::InsertNonlocalRow(int row, int offset)
{
  int alloc_len = numNonlocalBlockRows_;
  EPETRA_CHK_ERR( Epetra_Util_insert(row, offset, nonlocalBlockRows_,
				     numNonlocalBlockRows_, alloc_len, 1) );

  int tmp1 = numNonlocalBlockRows_ - 1;
  int tmp2 = alloc_len - 1;

  EPETRA_CHK_ERR( Epetra_Util_insert(0, offset, nonlocalBlockRowLengths_,
				     tmp1, tmp2, 1) );

  --tmp1;
  --tmp2;
  int initialAllocLen = 16;
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
  return(-1);
}

//--------------------------------------------------------------------------
int Epetra_FEVbrMatrix::BeginSumIntoGlobalValues(int BlockRow,
						 int NumBlockEntries,
						 int *BlockIndices)
{
  return(-1);
}

//--------------------------------------------------------------------------
int Epetra_FEVbrMatrix::SetupForNonlocalSubmits(int BlockRow,
						int NumBlockEntries,
						int * BlockIndices, 
						bool IndicesAreLocal,
						Epetra_CombineMode SubmitMode)
{
  if (ignoreNonLocalEntries_) {
    curRowOffset_ = 0;
    return(0);
  }

  int insertPoint = -1;

  //find offset of this row in our list of nonlocal rows
  int rowoffset = Epetra_Util_binary_search(BlockRow, nonlocalBlockRows_,
					    numNonlocalBlockRows_, insertPoint);

  if (rowoffset < 0) {
    EPETRA_CHK_ERR( InsertNonlocalRow(BlockRow, insertPoint) );
    rowoffset = insertPoint;
  }

  for(int i=0; i<NumBlockEntries; ++i) {
    nonlocalBlockCols_[rowoffset][i] = BlockIndices[i];
  }

  curRowOffset_ = rowoffset;
  curColOffset_ = nonlocalBlockRowLengths_[rowoffset];
  curNumCols_ = NumBlockEntries;

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
  }

  return(0);
}

