
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

#include <Epetra_FECrsMatrix.h>

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
    numNonlocalRows_(0),
    nonlocalRows_(NULL),
    nonlocalRowLengths_(NULL),
    nonlocalRowAllocLengths_(NULL),
    nonlocalCols_(NULL),
    nonlocalCoefs_(NULL),
    workData_(NULL),
    workDataLength_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries)
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
    numNonlocalRows_(0),
    nonlocalRows_(NULL),
    nonlocalRowLengths_(NULL),
    nonlocalRowAllocLengths_(NULL),
    nonlocalCols_(NULL),
    nonlocalCoefs_(NULL),
    workData_(NULL),
    workDataLength_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries)
{
  myFirstRow_ = RowMap.MinMyGID();
  myNumRows_ = RowMap.NumMyElements();

  workData_ = new double[128];
  workDataLength_ = 128;
}

//----------------------------------------------------------------------------
Epetra_FECrsMatrix::~Epetra_FECrsMatrix()
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
                           values, format, true));
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::SumIntoGlobalValues(int numRows, const int* rows,
					    int numCols, const int* cols,
					    const double* const* values,
					    int format)
{
  return(InputGlobalValues(numRows, rows,
                           numCols, cols,
                           values, format, true));
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::SumIntoGlobalValues(int numIndices, const int* indices,
					    const double* values,
					    int format)
{
  return(InputGlobalValues(numIndices, indices,
                           numIndices, indices,
                           values, format, true));
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::SumIntoGlobalValues(int numRows, const int* rows,
					    int numCols, const int* cols,
					    const double* values,
					    int format)
{
  return(InputGlobalValues(numRows, rows,
                           numCols, cols,
                           values, format, true));
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::ReplaceGlobalValues(int numIndices, const int* indices,
					    const double* const* values,
					    int format)
{
  return(InputGlobalValues(numIndices, indices,
                           numIndices, indices,
                           values, format, false));
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::ReplaceGlobalValues(int numRows, const int* rows,
					    int numCols, const int* cols,
					    const double* const* values,
					    int format)
{
  return(InputGlobalValues(numRows, rows,
                           numCols, cols,
                           values, format, false));
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::ReplaceGlobalValues(int numIndices, const int* indices,
					    const double* values,
					    int format)
{
  return(InputGlobalValues(numIndices, indices,
                           numIndices, indices,
                           values, format, false));
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::ReplaceGlobalValues(int numRows, const int* rows,
					    int numCols, const int* cols,
					    const double* values,
					    int format)
{
  return(InputGlobalValues(numRows, rows,
                           numCols, cols,
                           values, format, false));
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::GlobalAssemble(bool callTransformToLocal)
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

  //First build a map that describes our nonlocal data.
  //We'll use the arbitrary distribution constructor of Map.

  Epetra_Map sourceMap(-1, numNonlocalRows_, nonlocalRows_,
		       RowMap().IndexBase(), RowMap().Comm());

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

  Epetra_Map colMap(-1, numCols, cols, RowMap().IndexBase(), RowMap().Comm());

  delete [] cols;
  numCols = 0;
  allocLen = 0;

  //now we need to create a matrix with sourceMap and colMap, and fill it with
  //our nonlocal data so we can then export it to the correct owning processors.

  Epetra_CrsMatrix tempMat(Copy, sourceMap, colMap, nonlocalRowLengths_);


  //Next we need to make sure the 'indices-are-global' attribute of tempMat's
  //graph is set to true, in case this processor doesn't end up calling the
  //InsertGlobalValues method...

  const Epetra_CrsGraph& graph = tempMat.Graph();
  Epetra_CrsGraph& nonconst_graph = const_cast<Epetra_CrsGraph&>(graph);
  nonconst_graph.SetIndicesAreGlobal(true);

  for(i=0; i<numNonlocalRows_; ++i) {
    EPETRA_CHK_ERR( tempMat.InsertGlobalValues(nonlocalRows_[i],
					       nonlocalRowLengths_[i],
					       nonlocalCoefs_[i],
					       nonlocalCols_[i]) );
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

  //now reset the values in our nonlocal data
  for(i=0; i<numNonlocalRows_; ++i) {
    for(j=0; j<nonlocalRowLengths_[i]; ++j) {
      nonlocalCols_[i][j] = 0;
      nonlocalCoefs_[i][j] = 0.0;
    }
    nonlocalRowLengths_[i] = 0;
  }

  return(0);
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::InputGlobalValues(int numRows, const int* rows,
					  int numCols, const int* cols,
					  const double*const* values,
					  int format, bool accumulate)
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
      if (IndicesAreLocal()) {
        if (accumulate) {
          EPETRA_CHK_ERR( crsthis->SumIntoGlobalValues(rows[i], numCols,
						       valuesptr, (int*)cols) );
        }
        else {
          EPETRA_CHK_ERR( crsthis->ReplaceGlobalValues(rows[i], numCols,
						       valuesptr, (int*)cols) );
        }
      }
      else {
        int err = crsthis->InsertGlobalValues(rows[i], numCols,
					      valuesptr, (int*)cols);
	if (err < 0) return(err);
      }
    }
    else {
      EPETRA_CHK_ERR( InputNonlocalGlobalValues(rows[i],
						numCols, cols,
						valuesptr, accumulate) );
    }
  }

  return(0);
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::InputGlobalValues(int numRows, const int* rows,
					  int numCols, const int* cols,
					  const double* values,
					  int format, bool accumulate)
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

  int offset = 0;
  for(int i=0; i<numRows; ++i) {
    double* valuesptr = (double*)(&(values[offset]));

    if (format == Epetra_FECrsMatrix::COLUMN_MAJOR) {
      //if the data is in column-major order, then we need to copy the i-th
      //column of the values table into workData_, and that will be the i-th
      //row. ... Is that clear?

      for(int j=0; j<numCols; ++j) {
	workData_[j] = values[i+j*numRows];
      }
      valuesptr = workData_;
    }

    if (Map().MyGID(rows[i])) {
      if (IndicesAreLocal()) {
        if (accumulate) {
          EPETRA_CHK_ERR( crsthis->SumIntoGlobalValues(rows[i], numCols,
						       valuesptr, (int*)cols) );
        }
        else {
          EPETRA_CHK_ERR( crsthis->ReplaceGlobalValues(rows[i], numCols,
						       valuesptr, (int*)cols) );
        }
      }
      else {
        int err = crsthis->InsertGlobalValues(rows[i], numCols,
					      valuesptr, (int*)cols);
	if (err < 0) return(err);
      }
    }
    else {
      if (!ignoreNonLocalEntries_) {
	EPETRA_CHK_ERR( InputNonlocalGlobalValues(rows[i],
						  numCols, cols,
						  valuesptr, accumulate) );
      }
    }
    offset += numCols;
  }

  return(0);
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::InputNonlocalGlobalValues(int row,
						  int numCols, const int* cols,
						  const double* values,
						  bool accumulate)
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
				       accumulate) );
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
					   bool accumulate)
{
  int* colIndices = nonlocalCols_[rowoffset];
  double* coefs = nonlocalCoefs_[rowoffset];

  int insertPoint = -1;
  int coloffset = Epetra_Util_binary_search(col, colIndices,
					    nonlocalRowLengths_[rowoffset],
					    insertPoint);

  if (coloffset >= 0) {
    if (accumulate) {
      nonlocalCoefs_[rowoffset][coloffset] += value;
    }
    else {
      nonlocalCoefs_[rowoffset][coloffset] = value;
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
