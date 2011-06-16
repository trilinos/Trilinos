
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2001 Sandia Corporation
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
// ************************************************************************
//@HEADER

#include <Epetra_FECrsMatrix.h>
#include <Epetra_IntSerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_FECrsGraph.h>
#include <Epetra_Export.h>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Util.h>

//----------------------------------------------------------------------------
Epetra_FECrsMatrix::Epetra_FECrsMatrix(Epetra_DataAccess CV,
				       const Epetra_Map& rowMap,
				       int* NumEntriesPerRow,
				       bool ignoreNonLocalEntries)
  : Epetra_CrsMatrix(CV, rowMap, NumEntriesPerRow),
    myFirstRow_(0),
    myNumRows_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    numNonlocalRows_(0),
    nonlocalRows_(NULL),
    nonlocalRowLengths_(NULL),
    nonlocalRowAllocLengths_(NULL),
    nonlocalCols_(NULL),
    nonlocalCoefs_(NULL),
    workData_(128),
    useNonlocalMatrix_ (false),
    nonlocalMatrix_ (NULL)
{
  myFirstRow_ = rowMap.MinMyGID();
  myNumRows_ = rowMap.NumMyElements();
}

//----------------------------------------------------------------------------
Epetra_FECrsMatrix::Epetra_FECrsMatrix(Epetra_DataAccess CV,
				       const Epetra_Map& rowMap,
				       int NumEntriesPerRow,
				       bool ignoreNonLocalEntries)
  : Epetra_CrsMatrix(CV, rowMap, NumEntriesPerRow),
    myFirstRow_(0),
    myNumRows_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    numNonlocalRows_(0),
    nonlocalRows_(NULL),
    nonlocalRowLengths_(NULL),
    nonlocalRowAllocLengths_(NULL),
    nonlocalCols_(NULL),
    nonlocalCoefs_(NULL),
    workData_(128),
    useNonlocalMatrix_ (false),
    nonlocalMatrix_ (NULL)
{
  myFirstRow_ = rowMap.MinMyGID();
  myNumRows_ = rowMap.NumMyElements();
}

//----------------------------------------------------------------------------
Epetra_FECrsMatrix::Epetra_FECrsMatrix(Epetra_DataAccess CV,
				       const Epetra_Map& rowMap,
				       const Epetra_Map& colMap,
				       int* NumEntriesPerRow,
				       bool ignoreNonLocalEntries)
  : Epetra_CrsMatrix(CV, rowMap, colMap, NumEntriesPerRow),
    myFirstRow_(0),
    myNumRows_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    numNonlocalRows_(0),
    nonlocalRows_(NULL),
    nonlocalRowLengths_(NULL),
    nonlocalRowAllocLengths_(NULL),
    nonlocalCols_(NULL),
    nonlocalCoefs_(NULL),
    workData_(128),
    useNonlocalMatrix_ (false),
    nonlocalMatrix_ (NULL)
{
  myFirstRow_ = rowMap.MinMyGID();
  myNumRows_ = rowMap.NumMyElements();
}

//----------------------------------------------------------------------------
Epetra_FECrsMatrix::Epetra_FECrsMatrix(Epetra_DataAccess CV,
				       const Epetra_Map& rowMap,
				       const Epetra_Map& colMap,
				       int NumEntriesPerRow,
				       bool ignoreNonLocalEntries)
  : Epetra_CrsMatrix(CV, rowMap, colMap, NumEntriesPerRow),
    myFirstRow_(0),
    myNumRows_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    numNonlocalRows_(0),
    nonlocalRows_(NULL),
    nonlocalRowLengths_(NULL),
    nonlocalRowAllocLengths_(NULL),
    nonlocalCols_(NULL),
    nonlocalCoefs_(NULL),
    workData_(128),
    useNonlocalMatrix_ (false),
    nonlocalMatrix_ (NULL)
{
  myFirstRow_ = rowMap.MinMyGID();
  myNumRows_ = rowMap.NumMyElements();
}

//----------------------------------------------------------------------------
Epetra_FECrsMatrix::Epetra_FECrsMatrix(Epetra_DataAccess CV,
				       const Epetra_CrsGraph& graph,
				       bool ignoreNonLocalEntries)
  : Epetra_CrsMatrix(CV, graph),
    myFirstRow_(0),
    myNumRows_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    numNonlocalRows_(0),
    nonlocalRows_(NULL),
    nonlocalRowLengths_(NULL),
    nonlocalRowAllocLengths_(NULL),
    nonlocalCols_(NULL),
    nonlocalCoefs_(NULL),
    workData_(128),
    useNonlocalMatrix_ (false),
    nonlocalMatrix_ (NULL)
{
  myFirstRow_ = RowMap().MinMyGID();
  myNumRows_ = RowMap().NumMyElements();
}
   
//----------------------------------------------------------------------------
Epetra_FECrsMatrix::Epetra_FECrsMatrix(Epetra_DataAccess CV,
              const Epetra_FECrsGraph& graph,
              bool ignoreNonLocalEntries)
  : Epetra_CrsMatrix(CV, graph),
    myFirstRow_(0),
    myNumRows_(0),
    ignoreNonLocalEntries_(ignoreNonLocalEntries),
    numNonlocalRows_(0),
    nonlocalRows_(NULL),
    nonlocalRowLengths_(NULL),
    nonlocalRowAllocLengths_(NULL),
    nonlocalCols_(NULL),
    nonlocalCoefs_(NULL),
    workData_(128),
    useNonlocalMatrix_ (graph.UseNonlocalGraph() && graph.nonlocalGraph_ != 0),
    nonlocalMatrix_ (useNonlocalMatrix_ ? 
        new Epetra_CrsMatrix(Copy,*graph.nonlocalGraph_) : NULL)
{
  myFirstRow_ = RowMap().MinMyGID();
  myNumRows_ = RowMap().NumMyElements();
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
   workData_(128),
   nonlocalMatrix_ (NULL)
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

  useNonlocalMatrix_ = src.useNonlocalMatrix_;

  myFirstRow_ = src.myFirstRow_;
  myNumRows_ = src.myNumRows_;
  ignoreNonLocalEntries_ = src.ignoreNonLocalEntries_;
  numNonlocalRows_ = src.numNonlocalRows_;

  if (numNonlocalRows_ < 1) {
    return( *this );
  }

  if (useNonlocalMatrix_ && src.nonlocalMatrix_ != 0) {
    *nonlocalMatrix_ = *src.nonlocalMatrix_;
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

  if (nonlocalMatrix_ != 0)
    delete nonlocalMatrix_;
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
int Epetra_FECrsMatrix::SumIntoGlobalValues(int GlobalRow, int NumEntries,
                                            const double* values, const int* Indices)
{
  if (Map().MyGID(GlobalRow))
    return Epetra_CrsMatrix::SumIntoGlobalValues(GlobalRow, NumEntries,
            values, Indices);
  else if (useNonlocalMatrix_)
    return nonlocalMatrix_->SumIntoGlobalValues(GlobalRow,
           NumEntries, values, Indices);
  else
    return InputNonlocalGlobalValues(GlobalRow, NumEntries, Indices, values, SUMINTO);
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::InsertGlobalValues(int GlobalRow, int NumEntries,
                                            const double* values, const int* Indices)
{
  return(InputGlobalValues(1, &GlobalRow,
                           NumEntries, Indices, values,
                           Epetra_FECrsMatrix::ROW_MAJOR, INSERT));
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::ReplaceGlobalValues(int GlobalRow, int NumEntries,
                                            const double* values, const int* Indices)
{
  return(InputGlobalValues(1, &GlobalRow,
                           NumEntries, Indices, values,
                           Epetra_FECrsMatrix::ROW_MAJOR, REPLACE));
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
int Epetra_FECrsMatrix::GlobalAssemble(bool callFillComplete, Epetra_CombineMode combineMode)
{
  return( GlobalAssemble(DomainMap(), RangeMap(), callFillComplete), combineMode );
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::GlobalAssemble(const Epetra_Map& domain_map,
                                       const Epetra_Map& range_map,
                                       bool callFillComplete,
                                       Epetra_CombineMode combineMode)
{
  if (Map().Comm().NumProc() < 2 || ignoreNonLocalEntries_) {
    if (callFillComplete) {
      EPETRA_CHK_ERR( FillComplete(domain_map, range_map) );
    }
    return(0);
  }

  Epetra_CrsMatrix* tempMat;
  if (useNonlocalMatrix_) {
    tempMat = nonlocalMatrix_;
  }
  else {
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
        EPETRA_CHK_ERR( FillComplete(domain_map, range_map) );
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

    for(int i=0; i<numNonlocalRows_; ++i) {
      for(int j=0; j<nonlocalRowLengths_[i]; ++j) {
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

    tempMat = new Epetra_CrsMatrix(Copy, *sourceMap, *colMap,
          nonlocalRowLengths_);


    for(int i=0; i<numNonlocalRows_; ++i) {
      EPETRA_CHK_ERR( tempMat->InsertGlobalValues(nonlocalRows_[i],
             nonlocalRowLengths_[i],
             nonlocalCoefs_[i],
             nonlocalCols_[i]) );
    }

    delete sourceMap;
    delete colMap;
  }

  //Next we need to make sure the 'indices-are-global' attribute of tempMat's
  //graph is set to true, in case this processor doesn't end up calling the
  //InsertGlobalValues method...

  const Epetra_CrsGraph& graph = tempMat->Graph();
  Epetra_CrsGraph& nonconst_graph = const_cast<Epetra_CrsGraph&>(graph);
  nonconst_graph.SetIndicesAreGlobal(true);

  //Now we need to call FillComplete on our temp matrix. We need to
  //pass a DomainMap and RangeMap, which are not the same as the RowMap
  //and ColMap that we constructed the matrix with.
  EPETRA_CHK_ERR(tempMat->FillComplete(domain_map, range_map));

  Epetra_Export* exporter = new Epetra_Export(tempMat->RowMap(), RowMap());

  EPETRA_CHK_ERR(Export(*tempMat, *exporter, combineMode));

  if(callFillComplete) {
    EPETRA_CHK_ERR(FillComplete(domain_map, range_map));
  }

  //now reset the values in our nonlocal data
  if (!useNonlocalMatrix_) {
    for(int i=0; i<numNonlocalRows_; ++i) {
      for(int j=0; j<nonlocalRowLengths_[i]; ++j) {
        nonlocalCols_[i][j] = 0;
        nonlocalCoefs_[i][j] = 0.0;
      }
      nonlocalRowLengths_[i] = 0;
    }

    delete tempMat;
  }

  delete exporter;

  return(0);
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::InputGlobalValues_RowMajor(
            int numRows, const int* rows,
					  int numCols, const int* cols,
					  const double* values,
					  int mode)
{
  int returncode = 0;
  int err = 0;

  for(int i=0; i<numRows; ++i) {
    double* valuesptr = (double*)values + i*numRows;

    int local_row_id = Map().LID(rows[i]);
    if (local_row_id >= 0) {
      switch(mode) {
        case Epetra_FECrsMatrix::SUMINTO:
          err = this->Epetra_CrsMatrix::SumIntoGlobalValues(rows[i], numCols,
              valuesptr, (int*)cols);
          if (err<0) return(err);
          if (err>0) returncode = err;
          break;
        case Epetra_FECrsMatrix::REPLACE:
          err = this->Epetra_CrsMatrix::ReplaceGlobalValues(rows[i], numCols,
              valuesptr, (int*)cols);
          if (err<0) return(err);
          if (err>0) returncode = err;
          break;
        case Epetra_FECrsMatrix::INSERT:
          err = this->Epetra_CrsMatrix::InsertGlobalValues(rows[i], numCols,
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
      err = InputNonlocalGlobalValues(rows[i], numCols, cols,
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
					  const double*const* values,
					  int format, int mode)
{
  if (format != Epetra_FECrsMatrix::ROW_MAJOR &&
      format != Epetra_FECrsMatrix::COLUMN_MAJOR) {
    cerr << "Epetra_FECrsMatrix: unrecognized format specifier."<<endl;
    return(-1);
  }

  if (format == Epetra_FECrsMatrix::COLUMN_MAJOR) {
    workData_.resize(numCols);
  }

  int returncode = 0;

  for(int i=0; i<numRows; ++i) {
    if (format == Epetra_FECrsMatrix::ROW_MAJOR) {
      returncode += InputGlobalValues_RowMajor(1, &rows[i], numCols, cols,
                                          values[i], mode);
      if (returncode < 0) return returncode;
      continue;
    }

    //If we get to here, the data is in column-major order.

    double* valuesptr = &workData_[0];

    //Since the data is in column-major order, then we copy the i-th row
    //of the values table into workData_, in order to have the row in
    //contiguous memory.
    //This is slow and not thread-safe.

    for(int j=0; j<numCols; ++j) {
      valuesptr[j] = values[i][j];
    }

    returncode += InputGlobalValues_RowMajor(1, &rows[i], numCols, cols, valuesptr, mode);
    if (returncode < 0) return returncode;
  }

  return(returncode);
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::InputGlobalValues(int numRows, const int* rows,
					  int numCols, const int* cols,
					  const double* values,
					  int format, int mode)
{
  if (format == Epetra_FECrsMatrix::ROW_MAJOR) {
    return InputGlobalValues_RowMajor(numRows, rows, numCols, cols, values, mode);
  }

  workData_.resize(numCols);

  int returncode = 0;
  for(int i=0; i<numRows; ++i) {
    //copy each row out of the column-major values array, so we can pass it
    //to a row-major input function.
    for(int j=0; j<numCols; ++j) {
      workData_[j] = values[i+j*numRows];
    }
    int err = InputGlobalValues_RowMajor(1, &rows[i], numCols, cols, &workData_[0], mode);
    if (err < 0) return err;
    returncode += err;
  }

  return(returncode);
}

//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::InputNonlocalGlobalValues(int row,
						  int numCols, const int* cols,
						  const double* values,
						  int mode)
{
  // if we already have a nonlocal matrix object, this is easier...
  if (useNonlocalMatrix_) {
    int err, returncode = 0;
    double* valuesptr = (double*)values;
    switch(mode) {
    case Epetra_FECrsMatrix::SUMINTO:
      err = nonlocalMatrix_->SumIntoGlobalValues(row, numCols,
            valuesptr, (int*)cols);
      if (err<0) return(err);
      if (err>0) returncode = err;
      break;
    case Epetra_FECrsMatrix::REPLACE:
      err = nonlocalMatrix_->ReplaceGlobalValues(row, numCols,
            valuesptr, (int*)cols);
      if (err<0) return(err);
      if (err>0) returncode = err;
      break;
    case Epetra_FECrsMatrix::INSERT:
      err = nonlocalMatrix_->InsertGlobalValues(row, numCols,
           valuesptr, (int*)cols);
      if (err<0) return(err);
      if (err>0) returncode = err;
      break;
    default:
      cerr << "Epetra_FECrsMatrix: internal error, bad input mode."<<endl;
      return(-1);
    }
    return (returncode);
  }

  int insertPoint = -1;

  //find offset of this row in our list of nonlocal rows.
  int rowoffset = Epetra_Util_binary_search(row, nonlocalRows_,
             numNonlocalRows_, insertPoint);

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
  EPETRA_CHK_ERR( Epetra_Util_insert(row, offset, nonlocalRows_,
                                     numNonlocalRows_, alloc_len, 1) );

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
  int len = nonlocalRowLengths_[rowoffset];

  int insertPoint = -1;
  int coloffset = Epetra_Util_binary_search(col, colIndices,
					    len, insertPoint);

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

