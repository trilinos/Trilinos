
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

#include <Epetra_FECrsMatrix.h>
#include <Epetra_IntSerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_Export.h>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Util.h>

//----------------------------------------------------------------------------
Epetra_FECrsMatrix::Epetra_FECrsMatrix(Epetra_DataAccess CV,
				       const Epetra_Map& RowMap,
				       int* NumEntriesPerRow,
				       bool ignoreNonLocalEntries)
  : Epetra_CrsMatrix(CV, RowMap, NumEntriesPerRow),
    myFirstRow_(0),
    myNumRows_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    numNonlocalRows_(0),
    nonlocalRows_(NULL),
    nonlocalRowLengths_(NULL),
    nonlocalRowAllocLengths_(NULL),
    nonlocalCols_(NULL),
    nonlocalCoefs_(NULL),
    workData_(NULL),
    workDataLength_(0)
{
  myFirstRow_ = RowMap.MinMyGID();
  myNumRows_ = RowMap.NumMyElements();

  workData_ = new double[128];
  workDataLength_ = 128;
}

//----------------------------------------------------------------------------
Epetra_FECrsMatrix::Epetra_FECrsMatrix(Epetra_DataAccess CV,
				       const Epetra_Map& RowMap,
				       int NumEntriesPerRow,
				       bool ignoreNonLocalEntries)
  : Epetra_CrsMatrix(CV, RowMap, NumEntriesPerRow),
    myFirstRow_(0),
    myNumRows_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    numNonlocalRows_(0),
    nonlocalRows_(NULL),
    nonlocalRowLengths_(NULL),
    nonlocalRowAllocLengths_(NULL),
    nonlocalCols_(NULL),
    nonlocalCoefs_(NULL),
    workData_(NULL),
    workDataLength_(0)
{
  myFirstRow_ = RowMap.MinMyGID();
  myNumRows_ = RowMap.NumMyElements();

  workData_ = new double[128];
  workDataLength_ = 128;
}

//----------------------------------------------------------------------------
Epetra_FECrsMatrix::Epetra_FECrsMatrix(Epetra_DataAccess CV,
				       const Epetra_Map& RowMap,
				       const Epetra_Map& ColMap,
				       int* NumEntriesPerRow,
				       bool ignoreNonLocalEntries)
  : Epetra_CrsMatrix(CV, RowMap, ColMap, NumEntriesPerRow),
    myFirstRow_(0),
    myNumRows_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    numNonlocalRows_(0),
    nonlocalRows_(NULL),
    nonlocalRowLengths_(NULL),
    nonlocalRowAllocLengths_(NULL),
    nonlocalCols_(NULL),
    nonlocalCoefs_(NULL),
    workData_(NULL),
    workDataLength_(0)
{
  myFirstRow_ = RowMap.MinMyGID();
  myNumRows_ = RowMap.NumMyElements();

  workData_ = new double[128];
  workDataLength_ = 128;
}

//----------------------------------------------------------------------------
Epetra_FECrsMatrix::Epetra_FECrsMatrix(Epetra_DataAccess CV,
				       const Epetra_Map& RowMap,
				       const Epetra_Map& ColMap,
				       int NumEntriesPerRow,
				       bool ignoreNonLocalEntries)
  : Epetra_CrsMatrix(CV, RowMap, ColMap, NumEntriesPerRow),
    myFirstRow_(0),
    myNumRows_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    numNonlocalRows_(0),
    nonlocalRows_(NULL),
    nonlocalRowLengths_(NULL),
    nonlocalRowAllocLengths_(NULL),
    nonlocalCols_(NULL),
    nonlocalCoefs_(NULL),
    workData_(NULL),
    workDataLength_(0)
{
  myFirstRow_ = RowMap.MinMyGID();
  myNumRows_ = RowMap.NumMyElements();

  workData_ = new double[128];
  workDataLength_ = 128;
}

//----------------------------------------------------------------------------
Epetra_FECrsMatrix::Epetra_FECrsMatrix(Epetra_DataAccess CV,
				       const Epetra_CrsGraph& Graph,
				       bool ignoreNonLocalEntries)
  : Epetra_CrsMatrix(CV, Graph),
    myFirstRow_(0),
    myNumRows_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    numNonlocalRows_(0),
    nonlocalRows_(NULL),
    nonlocalRowLengths_(NULL),
    nonlocalRowAllocLengths_(NULL),
    nonlocalCols_(NULL),
    nonlocalCoefs_(NULL),
    workData_(NULL),
    workDataLength_(0)
{
  myFirstRow_ = RowMap().MinMyGID();
  myNumRows_ = RowMap().NumMyElements();

  workData_ = new double[128];
  workDataLength_ = 128;
}
   
//----------------------------------------------------------------------------
Epetra_FECrsMatrix::Epetra_FECrsMatrix(const Epetra_FECrsMatrix& src)
 : Epetra_CrsMatrix(src),
   myFirstRow_(0),
   myNumRows_(0),
   ignoreNonLocalEntries_(false),
   numNonlocalRows_(0),
   nonlocalRows_(NULL),
   nonlocalRowLengths_(NULL),
   nonlocalRowAllocLengths_(NULL),
   nonlocalCols_(NULL),
   nonlocalCoefs_(NULL),
   workData_(NULL),
   workDataLength_(0)
{
  operator=(src);
}

//----------------------------------------------------------------------------
Epetra_FECrsMatrix& Epetra_FECrsMatrix::operator=(const Epetra_FECrsMatrix& src)
{
  if (this == &src) {
    return( *this );
  }

  DeleteMemory();

  Epetra_CrsMatrix::operator=(src);

  myFirstRow_ = src.myFirstRow_;
  myNumRows_ = src.myNumRows_;
  ignoreNonLocalEntries_ = src.ignoreNonLocalEntries_;
  numNonlocalRows_ = src.numNonlocalRows_;

  workDataLength_ = 128;
  workData_ = new double[workDataLength_];

  if (numNonlocalRows_ < 1) {
    return( *this );
  }

  nonlocalRows_ = new int[numNonlocalRows_];
  nonlocalRowLengths_ = new int[numNonlocalRows_];
  nonlocalRowAllocLengths_ = new int[numNonlocalRows_];
  nonlocalCols_ = new int*[numNonlocalRows_];
  nonlocalCoefs_= new double*[numNonlocalRows_];

  for(int i=0; i<numNonlocalRows_; ++i) {
    nonlocalRows_[i] = src.nonlocalRows_[i];
    nonlocalRowLengths_[i] = src.nonlocalRowLengths_[i];
    nonlocalRowAllocLengths_[i] = src.nonlocalRowAllocLengths_[i];

    nonlocalCols_[i] = new int[nonlocalRowAllocLengths_[i]];
    nonlocalCoefs_[i] = new double[nonlocalRowAllocLengths_[i]];

    for(int j=0; j<nonlocalRowLengths_[i]; ++j) {
      nonlocalCols_[i][j] = src.nonlocalCols_[i][j];
      nonlocalCoefs_[i][j] = src.nonlocalCoefs_[i][j];
    }
  }

  return( *this );
}

//----------------------------------------------------------------------------
Epetra_FECrsMatrix::~Epetra_FECrsMatrix()
{
  DeleteMemory();
}

//----------------------------------------------------------------------------
void Epetra_FECrsMatrix::DeleteMemory()
{
  if (numNonlocalRows_ > 0) {
    for(int i=0; i<numNonlocalRows_; ++i) {
      delete [] nonlocalCols_[i];
      delete [] nonlocalCoefs_[i];
    }
    delete [] nonlocalCols_;
    delete [] nonlocalCoefs_;
    delete [] nonlocalRows_;
    delete [] nonlocalRowLengths_;
    delete [] nonlocalRowAllocLengths_;
    numNonlocalRows_ = 0;
  }

  delete [] workData_;
  workDataLength_ = 0;
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::SumIntoGlobalValues(int numIndices, const int* indices,
					    const double* const* values,
					    int format)
{
  return(InputGlobalValues(numIndices, indices,
                           numIndices, indices,
                           values, format, SUMINTO));
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::SumIntoGlobalValues(int numRows, const int* rows,
					    int numCols, const int* cols,
					    const double* const* values,
					    int format)
{
  return(InputGlobalValues(numRows, rows,
                           numCols, cols,
                           values, format, SUMINTO));
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::SumIntoGlobalValues(int numIndices, const int* indices,
					    const double* values,
					    int format)
{
  return(InputGlobalValues(numIndices, indices,
                           numIndices, indices,
                           values, format, SUMINTO));
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::SumIntoGlobalValues(int numRows, const int* rows,
					    int numCols, const int* cols,
					    const double* values,
					    int format)
{
  return(InputGlobalValues(numRows, rows,
                           numCols, cols,
                           values, format, SUMINTO));
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::SumIntoGlobalValues(const Epetra_IntSerialDenseVector& indices,
					    const Epetra_SerialDenseMatrix& values,
					    int format)
{
  if (indices.Length() != values.M() || indices.Length() != values.N()) {
    return(-1);
  }

  return( SumIntoGlobalValues(indices.Length(), indices.Values(),
			      values.A(), format) );
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::InsertGlobalValues(const Epetra_IntSerialDenseVector& indices,
					    const Epetra_SerialDenseMatrix& values,
					    int format)
{
  if (indices.Length() != values.M() || indices.Length() != values.N()) {
    return(-1);
  }

  return( InsertGlobalValues(indices.Length(), indices.Values(),
			      values.A(), format) );
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::ReplaceGlobalValues(const Epetra_IntSerialDenseVector& indices,
					    const Epetra_SerialDenseMatrix& values,
					    int format)
{
  if (indices.Length() != values.M() || indices.Length() != values.N()) {
    return(-1);
  }

  return( ReplaceGlobalValues(indices.Length(), indices.Values(),
			      values.A(), format) );
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::SumIntoGlobalValues(const Epetra_IntSerialDenseVector& rows,
					    const Epetra_IntSerialDenseVector& cols,
					    const Epetra_SerialDenseMatrix& values,
					    int format)
{
  if (rows.Length() != values.M() || cols.Length() != values.N()) {
    return(-1);
  }

  return( SumIntoGlobalValues(rows.Length(), rows.Values(),
			      cols.Length(), cols.Values(),
			      values.A(), format) );
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::InsertGlobalValues(const Epetra_IntSerialDenseVector& rows,
					   const Epetra_IntSerialDenseVector& cols,
					   const Epetra_SerialDenseMatrix& values,
					   int format)
{
  if (rows.Length() != values.M() || cols.Length() != values.N()) {
    return(-1);
  }

  return( InsertGlobalValues(rows.Length(), rows.Values(),
			     cols.Length(), cols.Values(),
			      values.A(), format) );
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::ReplaceGlobalValues(const Epetra_IntSerialDenseVector& rows,
					    const Epetra_IntSerialDenseVector& cols,
					    const Epetra_SerialDenseMatrix& values,
					    int format)
{
  if (rows.Length() != values.M() || cols.Length() != values.N()) {
    return(-1);
  }

  return( ReplaceGlobalValues(rows.Length(), rows.Values(),
			      cols.Length(), cols.Values(),
			      values.A(), format) );
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::InsertGlobalValues(int numIndices, const int* indices,
					    const double* const* values,
					    int format)
{
  return(InputGlobalValues(numIndices, indices,
                           numIndices, indices,
                           values, format, INSERT));
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::InsertGlobalValues(int numRows, const int* rows,
					    int numCols, const int* cols,
					    const double* const* values,
					    int format)
{
  return(InputGlobalValues(numRows, rows,
                           numCols, cols,
                           values, format, INSERT));
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::InsertGlobalValues(int numIndices, const int* indices,
					    const double* values,
					    int format)
{
  return(InputGlobalValues(numIndices, indices,
                           numIndices, indices,
                           values, format, INSERT));
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::InsertGlobalValues(int numRows, const int* rows,
					    int numCols, const int* cols,
					    const double* values,
					    int format)
{
  return(InputGlobalValues(numRows, rows,
                           numCols, cols,
                           values, format, INSERT));
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::ReplaceGlobalValues(int numIndices, const int* indices,
					    const double* const* values,
					    int format)
{
  return(InputGlobalValues(numIndices, indices,
                           numIndices, indices,
                           values, format, REPLACE));
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::ReplaceGlobalValues(int numRows, const int* rows,
					    int numCols, const int* cols,
					    const double* const* values,
					    int format)
{
  return(InputGlobalValues(numRows, rows,
                           numCols, cols,
                           values, format, REPLACE));
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::ReplaceGlobalValues(int numIndices, const int* indices,
					    const double* values,
					    int format)
{
  return(InputGlobalValues(numIndices, indices,
                           numIndices, indices,
                           values, format, REPLACE));
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::ReplaceGlobalValues(int numRows, const int* rows,
					    int numCols, const int* cols,
					    const double* values,
					    int format)
{
  return(InputGlobalValues(numRows, rows,
                           numCols, cols,
                           values, format, REPLACE));
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::GlobalAssemble(bool callFillComplete)
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
  //non-overlapping distribution defined by the map that 'this' matrix was
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

  //now we need to create a matrix with sourceMap and colMap, and fill it with
  //our nonlocal data so we can then export it to the correct owning processors.

  Epetra_CrsMatrix* tempMat = new Epetra_CrsMatrix(Copy, *sourceMap, *colMap,
						   nonlocalRowLengths_);


  //Next we need to make sure the 'indices-are-global' attribute of tempMat's
  //graph is set to true, in case this processor doesn't end up calling the
  //InsertGlobalValues method...

  const Epetra_CrsGraph& graph = tempMat->Graph();
  Epetra_CrsGraph& nonconst_graph = const_cast<Epetra_CrsGraph&>(graph);
  nonconst_graph.SetIndicesAreGlobal(true);

  for(i=0; i<numNonlocalRows_; ++i) {
    EPETRA_CHK_ERR( tempMat->InsertGlobalValues(nonlocalRows_[i],
					       nonlocalRowLengths_[i],
					       nonlocalCoefs_[i],
					       nonlocalCols_[i]) );
  }

  //Now we need to call FillComplete on our temp matrix. We need to
  //pass a DomainMap and RangeMap, which are not the same as the RowMap
  //and ColMap that we constructed the matrix with.
  EPETRA_CHK_ERR(tempMat->FillComplete(RowMap(), *sourceMap));

  Epetra_Export* exporter = new Epetra_Export(*sourceMap, RowMap());

  EPETRA_CHK_ERR(Export(*tempMat, *exporter, Add));

  if(callFillComplete) {
    EPETRA_CHK_ERR(FillComplete());
  }

  //now reset the values in our nonlocal data
  for(i=0; i<numNonlocalRows_; ++i) {
    for(j=0; j<nonlocalRowLengths_[i]; ++j) {
      nonlocalCols_[i][j] = 0;
      nonlocalCoefs_[i][j] = 0.0;
    }
    nonlocalRowLengths_[i] = 0;
  }

  delete exporter;
  delete tempMat;
  delete sourceMap;
  delete colMap;

  return(0);
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::InputGlobalValues(int numRows, const int* rows,
					  int numCols, const int* cols,
					  const double*const* values,
					  int format, int mode)
{
  if (format != Epetra_FECrsMatrix::ROW_MAJOR &&
      format != Epetra_FECrsMatrix::COLUMN_MAJOR) {
    cerr << "Epetra_FECrsMatrix: unrecognized format specifier."<<endl;
    return(-1);
  }

  Epetra_CrsMatrix* crsthis = static_cast<Epetra_CrsMatrix*>(this);
  if (crsthis == NULL) {
    cerr << "Epetra_FECrsMatrix::InputGlobalValues: static_cast failed."<<endl;
    return(-1);
  }

  if (format == Epetra_FECrsMatrix::COLUMN_MAJOR) {
    if (numCols > workDataLength_) {
      delete [] workData_;
      workDataLength_ = numCols*2;
      workData_ = new double[workDataLength_];
    }
  }

  int returncode = 0;
  int err = 0;

  for(int i=0; i<numRows; ++i) {
    double* valuesptr = (double*)values[i];

    if (format == Epetra_FECrsMatrix::COLUMN_MAJOR) {
      //if the data is in column-major order, then we need to copy the i-th
      //column of the values table into workData_, and that will be the i-th
      //row. ... Is that clear?

      for(int j=0; j<numCols; ++j) {
	workData_[j] = values[j][i];
      }
      valuesptr = workData_;
    }

    if (Map().MyGID(rows[i])) {
      switch(mode) {
      case Epetra_FECrsMatrix::SUMINTO:
        err = crsthis->SumIntoGlobalValues(rows[i], numCols,
					   valuesptr, (int*)cols);
	if (err<0) return(err);
	if (err>0) returncode = err;
	break;
      case Epetra_FECrsMatrix::REPLACE:
	err = crsthis->ReplaceGlobalValues(rows[i], numCols,
					   valuesptr, (int*)cols);
	if (err<0) return(err);
	if (err>0) returncode = err;
	break;
      case Epetra_FECrsMatrix::INSERT:
	err = crsthis->InsertGlobalValues(rows[i], numCols,
					  valuesptr, (int*)cols);
	if (err<0) return(err);
	if (err>0) returncode = err;
	break;
      default:
	cerr << "Epetra_FECrsMatrix: internal error, bad input mode."<<endl;
	return(-1);
      }
    }
    else {
      err = InputNonlocalGlobalValues(rows[i],
				      numCols, cols,
				      valuesptr, mode);
      if (err<0) return(err);
      if (err>0) returncode = err;
    }
  }

  return(returncode);
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::InputGlobalValues(int numRows, const int* rows,
					  int numCols, const int* cols,
					  const double* values,
					  int format, int mode)
{
  const double** values_2d = new const double*[numRows];
  int offset = 0;
  for(int i=0; i<numRows; ++i) {
    values_2d[i] = &(values[offset]);
    offset += numCols;
  }

  int err = InputGlobalValues(numRows, rows, numCols, cols,
			      values_2d, format, mode);
  delete [] values_2d;

  return(err);
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::InputNonlocalGlobalValues(int row,
						  int numCols, const int* cols,
						  const double* values,
						  int mode)
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
    EPETRA_CHK_ERR( InputNonlocalValue(rowoffset, cols[i], values[i],
				       mode) );
  }

  return(0);
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::InsertNonlocalRow(int row, int offset)
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
  double** newCoefs = new double*[numNonlocalRows_];

  if (newCols == NULL || newCoefs == NULL) {
    return(-1);
  }

  newCols[offset] = new int[initialAllocLen];
  newCoefs[offset] = new double[initialAllocLen];

  int index = 0;
  for(int i=0; i<numNonlocalRows_-1; ++i) {
    if (i == offset) {
      ++index;
    }

    newCols[index] = nonlocalCols_[i];
    newCoefs[index++] = nonlocalCoefs_[i];
  }

  delete [] nonlocalCols_;
  delete [] nonlocalCoefs_;

  nonlocalCols_ = newCols;
  nonlocalCoefs_ = newCoefs;

  return(0);
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::InputNonlocalValue(int rowoffset,
					   int col, double value,
					   int mode)
{
  int*& colIndices = nonlocalCols_[rowoffset];
  double*& coefs = nonlocalCoefs_[rowoffset];

  int insertPoint = -1;
  int coloffset = Epetra_Util_binary_search(col, colIndices,
					    nonlocalRowLengths_[rowoffset],
					    insertPoint);

  if (coloffset >= 0) {
    if (mode == SUMINTO || mode == INSERT) {
      coefs[coloffset] += value;
    }
    else {
      coefs[coloffset] = value;
    }
  }
  else {
    //else
    //  insert col in colIndices
    //  insert value in coefs

    int tmp1 = nonlocalRowLengths_[rowoffset];
    int tmp2 = nonlocalRowAllocLengths_[rowoffset];
    EPETRA_CHK_ERR( Epetra_Util_insert(col, insertPoint, colIndices, tmp1, tmp2));
    EPETRA_CHK_ERR( Epetra_Util_insert(value, insertPoint, coefs,
				       nonlocalRowLengths_[rowoffset],
				       nonlocalRowAllocLengths_[rowoffset]));
  }

  return(0);
}
