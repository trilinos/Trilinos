
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

#include "Epetra_FECrsGraph.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Distributor.h"
#include "Epetra_Util.h"
#include "Epetra_Comm.h"
#include "Epetra_HashTable.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_OffsetIndex.h"

//-------------------------------------------------------------------------------
Epetra_FECrsGraph::Epetra_FECrsGraph(Epetra_DataAccess CV,
				     const Epetra_BlockMap& RowMap,
				     int* NumIndicesPerRow,
				     bool ignoreNonLocalEntries)
  : Epetra_CrsGraph(CV, RowMap, NumIndicesPerRow),
    myFirstRow_(0),
    myNumRows_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    numNonlocalRows_(0),
    nonlocalRows_(NULL),
    nonlocalRowLengths_(NULL),
    nonlocalRowAllocLengths_(NULL),
    nonlocalCols_(NULL)
{
}

//-------------------------------------------------------------------------------
Epetra_FECrsGraph::Epetra_FECrsGraph(Epetra_DataAccess CV,
				     const Epetra_BlockMap& RowMap,
				     int NumIndicesPerRow,
				     bool ignoreNonLocalEntries)
  : Epetra_CrsGraph(CV, RowMap, NumIndicesPerRow),
    myFirstRow_(0),
    myNumRows_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    numNonlocalRows_(0),
    nonlocalRows_(NULL),
    nonlocalRowLengths_(NULL),
    nonlocalRowAllocLengths_(NULL),
    nonlocalCols_(NULL)
{
  myFirstRow_ = RowMap.MinMyGID();
  myNumRows_ = RowMap.NumMyElements();
}

//-------------------------------------------------------------------------------
Epetra_FECrsGraph::Epetra_FECrsGraph(Epetra_DataAccess CV,
				     const Epetra_BlockMap& RowMap,
				     const Epetra_BlockMap& ColMap,
				     int* NumIndicesPerRow,
				     bool ignoreNonLocalEntries)
  : Epetra_CrsGraph(CV, RowMap, ColMap, NumIndicesPerRow),
    myFirstRow_(0),
    myNumRows_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    numNonlocalRows_(0),
    nonlocalRows_(NULL),
    nonlocalRowLengths_(NULL),
    nonlocalRowAllocLengths_(NULL),
    nonlocalCols_(NULL)
{
  myFirstRow_ = RowMap.MinMyGID();
  myNumRows_ = RowMap.NumMyElements();
}

//-------------------------------------------------------------------------------
Epetra_FECrsGraph::Epetra_FECrsGraph(Epetra_DataAccess CV,
				     const Epetra_BlockMap& RowMap,
				     const Epetra_BlockMap& ColMap,
				     int NumIndicesPerRow,
				     bool ignoreNonLocalEntries)
  : Epetra_CrsGraph(CV, RowMap, ColMap, NumIndicesPerRow),
    myFirstRow_(0),
    myNumRows_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    numNonlocalRows_(0),
    nonlocalRows_(NULL),
    nonlocalRowLengths_(NULL),
    nonlocalRowAllocLengths_(NULL),
    nonlocalCols_(NULL)
{
  myFirstRow_ = RowMap.MinMyGID();
  myNumRows_ = RowMap.NumMyElements();
}

//-------------------------------------------------------------------------------
Epetra_FECrsGraph::~Epetra_FECrsGraph()
{
  DeleteMemory();
}

//----------------------------------------------------------------------------
void Epetra_FECrsGraph::DeleteMemory()
{
  if (numNonlocalRows_ > 0) {
    for(int i=0; i<numNonlocalRows_; ++i) {
      delete [] nonlocalCols_[i];
    }
    delete [] nonlocalCols_;
    delete [] nonlocalRows_;
    delete [] nonlocalRowLengths_;
    delete [] nonlocalRowAllocLengths_;
    numNonlocalRows_ = 0;
  }
}

//----------------------------------------------------------------------------
int Epetra_FECrsGraph::InsertGlobalIndices(int numRows, const int* rows,
					   int numCols, const int* cols)
{
  int returncode = 0;
  int err = 0;

  Epetra_CrsGraph* thisgraph = static_cast<Epetra_CrsGraph*>(this);

  for(int i=0; i<numRows; ++i) {
    if (Map().MyGID(rows[i])) {
      err = thisgraph->InsertGlobalIndices(rows[i], numCols, (int*)cols);
    }
    else {
      err = InputNonlocalIndices(rows[i], numCols, cols);
    }

    if (err<0) return(err);
    if (err>0) returncode = err;
  }

  return(returncode);
}

//----------------------------------------------------------------------------
int Epetra_FECrsGraph::InsertNonlocalRow(int row, int offset)
{
  int alloc_len = numNonlocalRows_;
  EPETRA_CHK_ERR( Epetra_Util_insert(row, offset, nonlocalRows_, numNonlocalRows_,
				     alloc_len, 1) );

  int tmp1 = numNonlocalRows_-1;
  int tmp2 = alloc_len-1;

  EPETRA_CHK_ERR( Epetra_Util_insert(0, offset, nonlocalRowLengths_,
				     tmp1, tmp2, 1) );

  --tmp1;
  --tmp2;
  int initialAllocLen = 16;
  EPETRA_CHK_ERR( Epetra_Util_insert(initialAllocLen, offset,
				     nonlocalRowAllocLengths_,
				     tmp1, tmp2, 1) );

  int** newCols = new int*[numNonlocalRows_];

  newCols[offset] = new int[initialAllocLen];

  int index = 0;
  for(int i=0; i<numNonlocalRows_-1; ++i) {
    if (i == offset) {
      ++index;
    }

    newCols[index++] = nonlocalCols_[i];
  }

  delete [] nonlocalCols_;

  nonlocalCols_ = newCols;

  return(0);
}

//----------------------------------------------------------------------------
int Epetra_FECrsGraph::InputNonlocalIndices(int row,
					    int numCols,
					    const int* cols)
{
  int insertPoint = -1;

  //find offset of this row in our list of nonlocal rows.
  int rowoffset = Epetra_Util_binary_search(row, nonlocalRows_, numNonlocalRows_,
					    insertPoint);

  if (rowoffset < 0) {
    EPETRA_CHK_ERR( InsertNonlocalRow(row, insertPoint) );
    rowoffset = insertPoint;
  }

  for(int i=0; i<numCols; ++i) {
    EPETRA_CHK_ERR( InputNonlocalIndex(rowoffset, cols[i]) );
  }

  return(0); 
}

//----------------------------------------------------------------------------
int Epetra_FECrsGraph::InputNonlocalIndex(int rowoffset,
					  int col)
{
  int*& colIndices = nonlocalCols_[rowoffset];

  int insertPoint = -1;
  int coloffset = Epetra_Util_binary_search(col, colIndices,
					    nonlocalRowLengths_[rowoffset],
					    insertPoint);

  if (coloffset < 0) {
    //  insert col in colIndices
    EPETRA_CHK_ERR( Epetra_Util_insert(col, insertPoint, colIndices,
				       nonlocalRowLengths_[rowoffset],
				       nonlocalRowAllocLengths_[rowoffset]));
  }

  return(0);
}

//----------------------------------------------------------------------------
int Epetra_FECrsGraph::GlobalAssemble(bool callFillComplete)
{
  if (Map().Comm().NumProc() < 2 || ignoreNonLocalEntries_) {
    if (callFillComplete) {
      EPETRA_CHK_ERR( FillComplete() );
    }
    return(0);
  }

  int i, j;

  //In this method we need to gather all the non-local (overlapping) data
  //that's been input on each processor, into the
  //non-overlapping distribution defined by the map that 'this' graph was
  //constructed with.

  //First build a map that describes our nonlocal data.
  //We'll use the arbitrary distribution constructor of Map.

  Epetra_Map* sourceMap = new Epetra_Map(-1, numNonlocalRows_, nonlocalRows_,
					 Map().IndexBase(), Map().Comm());

  //If sourceMap has global size 0, then no nonlocal data exists and we can
  //skip most of this function.
  if (sourceMap->NumGlobalElements() < 1) {
    if (callFillComplete) {
      EPETRA_CHK_ERR( FillComplete() );
    }
    delete sourceMap;
    return(0);
  }

  //We also need to build a column-map, containing the columns in our
  //nonlocal data. To do that, create a list of all column-indices that
  //occur in our nonlocal rows.

  int numCols = 0, allocLen = 0;
  int* cols = NULL;
  int insertPoint = -1;

  for(i=0; i<numNonlocalRows_; ++i) {
    for(j=0; j<nonlocalRowLengths_[i]; ++j) {
      int col = nonlocalCols_[i][j];
      int offset = Epetra_Util_binary_search(col, cols, numCols, insertPoint);
      if (offset < 0) {
	EPETRA_CHK_ERR( Epetra_Util_insert(col, insertPoint, cols,
					   numCols, allocLen) );
      }
    }
  }

  Epetra_Map* colMap = new Epetra_Map(-1, numCols, cols,
				      Map().IndexBase(), Map().Comm());

  delete [] cols;
  numCols = 0;
  allocLen = 0;

  //now we need to create a graph with sourceMap and colMap, and fill it with
  //our nonlocal data so we can then export it to the correct owning processors.

  Epetra_CrsGraph* tempGrph = new Epetra_CrsGraph(Copy, *sourceMap, *colMap,
						  nonlocalRowLengths_);


  //Next we need to make sure the 'indices-are-global' attribute of tempGrph
  //is set to true, in case this processor doesn't end up calling the
  //InsertGlobalIndices method...

  tempGrph->SetIndicesAreGlobal(true);

  for(i=0; i<numNonlocalRows_; ++i) {
    EPETRA_CHK_ERR( tempGrph->InsertGlobalIndices(nonlocalRows_[i],
						  nonlocalRowLengths_[i],
						  nonlocalCols_[i]) );
  }

  //Now we need to call FillComplete on our temp graph. We need to
  //pass a DomainMap and RangeMap, which are not the same as the RowMap
  //and ColMap that we constructed the graph with.
  EPETRA_CHK_ERR(tempGrph->FillComplete(RowMap(), *sourceMap));

  Epetra_Export* exporter = new Epetra_Export(*sourceMap, RowMap());

  EPETRA_CHK_ERR(Export(*tempGrph, *exporter, Add));

  if(callFillComplete) {
    EPETRA_CHK_ERR(FillComplete());
  }

  //now reset the values in our nonlocal data
  for(i=0; i<numNonlocalRows_; ++i) {
    for(j=0; j<nonlocalRowLengths_[i]; ++j) {
      nonlocalCols_[i][j] = 0;
    }
    nonlocalRowLengths_[i] = 0;
  }

  delete exporter;
  delete tempGrph;
  delete sourceMap;
  delete colMap;

  return(0);
}
