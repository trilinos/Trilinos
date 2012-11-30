
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
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

#include <Epetra_ConfigDefs.h>
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
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    nonlocalRows_int_(),
    nonlocalCols_int_(),
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    nonlocalRows_LL_(),
    nonlocalCols_LL_(),
#endif
    nonlocalCoefs_(),
    workData_(128),
    useNonlocalMatrix_ (false),
    nonlocalMatrix_ (NULL),
    sourceMap_(NULL),
    colMap_(NULL),
    exporter_(NULL)
{
  myFirstRow_ = rowMap.MinMyGID64();
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
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    nonlocalRows_int_(),
    nonlocalCols_int_(),
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    nonlocalRows_LL_(),
    nonlocalCols_LL_(),
#endif
    nonlocalCoefs_(),
    workData_(128),
    useNonlocalMatrix_ (false),
    nonlocalMatrix_ (NULL),
    sourceMap_(NULL),
    colMap_(NULL),
    exporter_(NULL)
{
  myFirstRow_ = rowMap.MinMyGID64();
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
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    nonlocalRows_int_(),
    nonlocalCols_int_(),
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    nonlocalRows_LL_(),
    nonlocalCols_LL_(),
#endif
    nonlocalCoefs_(),
    workData_(128),
    useNonlocalMatrix_ (false),
    nonlocalMatrix_ (NULL),
    sourceMap_(NULL),
    colMap_(NULL),
    exporter_(NULL)
{
  myFirstRow_ = rowMap.MinMyGID64();
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
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    nonlocalRows_int_(),
    nonlocalCols_int_(),
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    nonlocalRows_LL_(),
    nonlocalCols_LL_(),
#endif
    nonlocalCoefs_(),
    workData_(128),
    useNonlocalMatrix_ (false),
    nonlocalMatrix_ (NULL),
    sourceMap_(NULL),
    colMap_(NULL),
    exporter_(NULL)
{
  myFirstRow_ = rowMap.MinMyGID64();
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
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    nonlocalRows_int_(),
    nonlocalCols_int_(),
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    nonlocalRows_LL_(),
    nonlocalCols_LL_(),
#endif
    nonlocalCoefs_(),
    workData_(128),
    useNonlocalMatrix_ (false),
    nonlocalMatrix_ (NULL),
    sourceMap_(NULL),
    colMap_(NULL),
    exporter_(NULL)
{
  myFirstRow_ = RowMap().MinMyGID64();
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
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    nonlocalRows_int_(),
    nonlocalCols_int_(),
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    nonlocalRows_LL_(),
    nonlocalCols_LL_(),
#endif
    nonlocalCoefs_(),
    workData_(128),
    useNonlocalMatrix_ (graph.UseNonlocalGraph() && graph.nonlocalGraph_ != 0),
    nonlocalMatrix_ (useNonlocalMatrix_ ? 
    		new Epetra_CrsMatrix(Copy,*graph.nonlocalGraph_) : NULL),
    sourceMap_(NULL),
    colMap_(NULL),
    exporter_(NULL)
{
  myFirstRow_ = RowMap().MinMyGID64();
  myNumRows_ = RowMap().NumMyElements();
}
   
//----------------------------------------------------------------------------
Epetra_FECrsMatrix::Epetra_FECrsMatrix(const Epetra_FECrsMatrix& src)
 : Epetra_CrsMatrix(src),
   myFirstRow_(0),
   myNumRows_(0),
   ignoreNonLocalEntries_(false),
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   nonlocalRows_int_(),
   nonlocalCols_int_(),
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   nonlocalRows_LL_(),
   nonlocalCols_LL_(),
#endif
   nonlocalCoefs_(),
   workData_(128),
   nonlocalMatrix_ (NULL),
   sourceMap_(NULL),
   colMap_(NULL),
   exporter_(NULL)
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

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if (src.RowMap().GlobalIndicesInt() && nonlocalRows_int_.size() < 1) {
    return( *this );
  }
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if (src.RowMap().GlobalIndicesLongLong() && nonlocalRows_LL_.size() < 1) {
    return( *this );
  }
#endif

  if (useNonlocalMatrix_ && src.nonlocalMatrix_ != 0) {
    *nonlocalMatrix_ = *src.nonlocalMatrix_;
    return( *this );
  }

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if (src.RowMap().GlobalIndicesInt()) {
    nonlocalRows_int_ = src.nonlocalRows_int_;
    nonlocalCols_int_ = src.nonlocalCols_int_;
  }
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if (src.RowMap().GlobalIndicesLongLong()) {
    nonlocalRows_LL_ = src.nonlocalRows_LL_;
    nonlocalCols_LL_ = src.nonlocalCols_LL_;
  }
#endif

  nonlocalCoefs_= src.nonlocalCoefs_;

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
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if (RowMap().GlobalIndicesInt()) {
    nonlocalRows_int_.clear();
    nonlocalCols_int_.clear();
  }
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if (RowMap().GlobalIndicesLongLong()) {
    nonlocalRows_LL_.clear();
    nonlocalCols_LL_.clear();
  }
#endif

  nonlocalCoefs_.clear();

  if (nonlocalMatrix_ != 0)
    delete nonlocalMatrix_;

  if ( sourceMap_ )
	delete sourceMap_;
  if ( colMap_ )
	delete colMap_;
  if ( exporter_ )
	delete exporter_;

}

//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::SumIntoGlobalValues(int numIndices, const int* indices,
              const double* const* values,
              int format)
{
  return(InputGlobalValues(numIndices, indices,
                           numIndices, indices,
                           values, format, SUMINTO));
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::SumIntoGlobalValues(int numIndices, const long long* indices,
              const double* const* values,
              int format)
{
  return(InputGlobalValues(numIndices, indices,
                           numIndices, indices,
                           values, format, SUMINTO));
}
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::SumIntoGlobalValues(int numRows, const int* rows,
              int numCols, const int* cols,
              const double* const* values,
              int format)
{
  return(InputGlobalValues(numRows, rows,
                           numCols, cols,
                           values, format, SUMINTO));
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::SumIntoGlobalValues(int numRows, const long long* rows,
              int numCols, const long long* cols,
              const double* const* values,
              int format)
{
  return(InputGlobalValues(numRows, rows,
                           numCols, cols,
                           values, format, SUMINTO));
}
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::SumIntoGlobalValues(int numIndices, const int* indices,
              const double* values,
              int format)
{
  return(InputGlobalValues(numIndices, indices,
                           numIndices, indices,
                           values, format, SUMINTO));
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::SumIntoGlobalValues(int numIndices, const long long* indices,
              const double* values,
              int format)
{
  return(InputGlobalValues(numIndices, indices,
                           numIndices, indices,
                           values, format, SUMINTO));
}
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::SumIntoGlobalValues(int numRows, const int* rows,
              int numCols, const int* cols,
              const double* values,
              int format)
{
  return(InputGlobalValues(numRows, rows,
                           numCols, cols,
                           values, format, SUMINTO));
}

#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::SumIntoGlobalValues(int numRows, const long long* rows,
              int numCols, const long long* cols,
              const double* values,
              int format)
{
  return(InputGlobalValues(numRows, rows,
                           numCols, cols,
                           values, format, SUMINTO));
}
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
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

#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::SumIntoGlobalValues(const Epetra_LongLongSerialDenseVector& indices,
              const Epetra_SerialDenseMatrix& values,
              int format)
{
  if (indices.Length() != values.M() || indices.Length() != values.N()) {
    return(-1);
  }

  return( SumIntoGlobalValues(indices.Length(), indices.Values(),
            values.A(), format) );
}
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
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
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::InsertGlobalValues(const Epetra_LongLongSerialDenseVector& indices,
              const Epetra_SerialDenseMatrix& values,
              int format)
{
  if (indices.Length() != values.M() || indices.Length() != values.N()) {
    return(-1);
  }

  return( InsertGlobalValues(indices.Length(), indices.Values(),
            values.A(), format) );
}
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
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
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::ReplaceGlobalValues(const Epetra_LongLongSerialDenseVector& indices,
              const Epetra_SerialDenseMatrix& values,
              int format)
{
  if (indices.Length() != values.M() || indices.Length() != values.N()) {
    return(-1);
  }

  return( ReplaceGlobalValues(indices.Length(), indices.Values(),
            values.A(), format) );
}
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
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
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::SumIntoGlobalValues(const Epetra_LongLongSerialDenseVector& rows,
              const Epetra_LongLongSerialDenseVector& cols,
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
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
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
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::InsertGlobalValues(const Epetra_LongLongSerialDenseVector& rows,
             const Epetra_LongLongSerialDenseVector& cols,
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
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
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
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES

int Epetra_FECrsMatrix::ReplaceGlobalValues(const Epetra_LongLongSerialDenseVector& rows,
              const Epetra_LongLongSerialDenseVector& cols,
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
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::InsertGlobalValues(int numIndices, const int* indices,
            const double* const* values,
            int format)
{
  return(InputGlobalValues(numIndices, indices,
                           numIndices, indices,
                           values, format, INSERT));
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES

int Epetra_FECrsMatrix::InsertGlobalValues(int numIndices, const long long* indices,
            const double* const* values,
            int format)
{
  return(InputGlobalValues(numIndices, indices,
                           numIndices, indices,
                           values, format, INSERT));
}
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::InsertGlobalValues(int numRows, const int* rows,
            int numCols, const int* cols,
            const double* const* values,
            int format)
{
  return(InputGlobalValues(numRows, rows,
                           numCols, cols,
                           values, format, INSERT));
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES

int Epetra_FECrsMatrix::InsertGlobalValues(int numRows, const long long* rows,
            int numCols, const long long* cols,
            const double* const* values,
            int format)
{
  return(InputGlobalValues(numRows, rows,
                           numCols, cols,
                           values, format, INSERT));
}
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::InsertGlobalValues(int numIndices, const int* indices,
              const double* values,
              int format)
{
  return(InputGlobalValues(numIndices, indices,
                           numIndices, indices,
                           values, format, INSERT));
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES

int Epetra_FECrsMatrix::InsertGlobalValues(int numIndices, const long long* indices,
              const double* values,
              int format)
{
  return(InputGlobalValues(numIndices, indices,
                           numIndices, indices,
                           values, format, INSERT));
}
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::InsertGlobalValues(int numRows, const int* rows,
              int numCols, const int* cols,
              const double* values,
              int format)
{
  return(InputGlobalValues(numRows, rows,
                           numCols, cols,
                           values, format, INSERT));
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES

int Epetra_FECrsMatrix::InsertGlobalValues(int numRows, const long long* rows,
              int numCols, const long long* cols,
              const double* values,
              int format)
{
  return(InputGlobalValues(numRows, rows,
                           numCols, cols,
                           values, format, INSERT));
}
#endif
//----------------------------------------------------------------------------
template<typename int_type>
int Epetra_FECrsMatrix::SumIntoGlobalValues(int_type GlobalRow, int NumEntries,
                                            const double* values, const int_type* Indices)
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

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::SumIntoGlobalValues(int GlobalRow, int NumEntries,
                                            const double* values, const int* Indices)
{
  if(RowMap().GlobalIndicesInt())
  return SumIntoGlobalValues<int>(GlobalRow, NumEntries, values, Indices);
  else
  throw ReportError("Epetra_FECrsMatrix::SumIntoGlobalValues int version called for a matrix that is not int.", -1);
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES


int Epetra_FECrsMatrix::SumIntoGlobalValues(long long GlobalRow, int NumEntries,
                                            const double* values, const long long* Indices)
{
  if(RowMap().GlobalIndicesLongLong())
  return SumIntoGlobalValues<long long>(GlobalRow, NumEntries, values, Indices);
  else
  throw ReportError("Epetra_FECrsMatrix::SumIntoGlobalValues long long version called for a matrix that is not long long.", -1);
}
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::InsertGlobalValues(int GlobalRow, int NumEntries,
                                            const double* values, const int* Indices)
{
  return(InputGlobalValues(1, &GlobalRow,
                           NumEntries, Indices, values,
                           Epetra_FECrsMatrix::ROW_MAJOR, INSERT));
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES

int Epetra_FECrsMatrix::InsertGlobalValues(long long GlobalRow, int NumEntries,
                                            const double* values, const long long* Indices)
{
  return(InputGlobalValues(1, &GlobalRow,
                           NumEntries, Indices, values,
                           Epetra_FECrsMatrix::ROW_MAJOR, INSERT));
}
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::InsertGlobalValues(int GlobalRow, int NumEntries,
                                            double* values, int* Indices)
{
  return(InputGlobalValues(1, &GlobalRow,
                           NumEntries, Indices, values,
                           Epetra_FECrsMatrix::ROW_MAJOR, INSERT));
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::InsertGlobalValues(long long GlobalRow, int NumEntries,
                                            double* values, long long* Indices)
{
  return(InputGlobalValues(1, &GlobalRow,
                           NumEntries, Indices, values,
                           Epetra_FECrsMatrix::ROW_MAJOR, INSERT));
}
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::ReplaceGlobalValues(int GlobalRow, int NumEntries,
                                            const double* values, const int* Indices)
{
  return(InputGlobalValues(1, &GlobalRow,
                           NumEntries, Indices, values,
                           Epetra_FECrsMatrix::ROW_MAJOR, REPLACE));
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::ReplaceGlobalValues(long long GlobalRow, int NumEntries,
                                            const double* values, const long long* Indices)
{
  return(InputGlobalValues(1, &GlobalRow,
                           NumEntries, Indices, values,
                           Epetra_FECrsMatrix::ROW_MAJOR, REPLACE));
}
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::ReplaceGlobalValues(int numIndices, const int* indices,
              const double* const* values,
              int format)
{
  return(InputGlobalValues(numIndices, indices,
                           numIndices, indices,
                           values, format, REPLACE));
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::ReplaceGlobalValues(int numIndices, const long long* indices,
              const double* const* values,
              int format)
{
  return(InputGlobalValues(numIndices, indices,
                           numIndices, indices,
                           values, format, REPLACE));
}
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::ReplaceGlobalValues(int numRows, const int* rows,
              int numCols, const int* cols,
              const double* const* values,
              int format)
{
  return(InputGlobalValues(numRows, rows,
                           numCols, cols,
                           values, format, REPLACE));
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::ReplaceGlobalValues(int numRows, const long long* rows,
              int numCols, const long long* cols,
              const double* const* values,
              int format)
{
  return(InputGlobalValues(numRows, rows,
                           numCols, cols,
                           values, format, REPLACE));
}
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::ReplaceGlobalValues(int numIndices, const int* indices,
              const double* values,
              int format)
{
  return(InputGlobalValues(numIndices, indices,
                           numIndices, indices,
                           values, format, REPLACE));
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::ReplaceGlobalValues(int numIndices, const long long* indices,
              const double* values,
              int format)
{
  return(InputGlobalValues(numIndices, indices,
                           numIndices, indices,
                           values, format, REPLACE));
}
#endif
//----------------------------------------------------------------------------
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::ReplaceGlobalValues(int numRows, const int* rows,
              int numCols, const int* cols,
              const double* values,
              int format)
{
  return(InputGlobalValues(numRows, rows,
                           numCols, cols,
                           values, format, REPLACE));
}
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
int Epetra_FECrsMatrix::ReplaceGlobalValues(int numRows, const long long* rows,
              int numCols, const long long* cols,
              const double* values,
              int format)
{
  return(InputGlobalValues(numRows, rows,
                           numCols, cols,
                           values, format, REPLACE));
}
#endif
//----------------------------------------------------------------------------
int Epetra_FECrsMatrix::GlobalAssemble(bool callFillComplete, Epetra_CombineMode combineMode,
                                       bool save_off_and_reuse_map_exporter)
{
  return( GlobalAssemble(DomainMap(), RangeMap(), callFillComplete, combineMode, save_off_and_reuse_map_exporter));
}

//----------------------------------------------------------------------------
template<typename int_type>
int Epetra_FECrsMatrix::GlobalAssemble(const Epetra_Map& domain_map,
                                       const Epetra_Map& range_map,
                                       bool callFillComplete,
                                       Epetra_CombineMode combineMode,
                                       bool save_off_and_reuse_map_exporter)
{
  if (Map().Comm().NumProc() < 2 || ignoreNonLocalEntries_) {
    if (callFillComplete) {
      EPETRA_CHK_ERR( FillComplete(domain_map, range_map) );
    }
    return(0);
  }

  std::vector<int_type>& nonlocalRows_var = nonlocalRows<int_type>();
  std::vector<std::vector<int_type> >& nonlocalCols_var = nonlocalCols<int_type>();

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

    int_type* nlr_ptr = nonlocalRows_var.size() > 0 ? &nonlocalRows_var[0] : 0;
    if (sourceMap_ == NULL)
      sourceMap_ = new Epetra_Map((int_type) -1, (int) nonlocalRows_var.size(), nlr_ptr,
            Map().IndexBase(), Map().Comm());

    //If sourceMap has global size 0, then no nonlocal data exists and we can
    //skip most of this function.
    if (sourceMap_->NumGlobalElements64() < 1) {
      if (callFillComplete) {
        EPETRA_CHK_ERR( FillComplete(domain_map, range_map) );
      }
      if (!save_off_and_reuse_map_exporter) {
        delete sourceMap_;
        sourceMap_ = NULL;
      }
      return(0);
    }

    //We also need to build a column-map, containing the columns in our
    //nonlocal data. To do that, create a list of all column-indices that
    //occur in our nonlocal rows.

    if ( colMap_ == NULL ) {
      std::vector<int_type> cols;

      for(size_t i=0; i<nonlocalRows_var.size(); ++i) {
        for(size_t j=0; j<nonlocalCols_var[i].size(); ++j) {
          int_type col = nonlocalCols_var[i][j];
          typename std::vector<int_type>::iterator it =
            std::lower_bound(cols.begin(), cols.end(), col);
          if (it == cols.end() || *it != col) {
            cols.insert(it, col);
          }
        }
      }

      int_type* cols_ptr = cols.size() > 0 ? &cols[0] : 0;

      colMap_ = new Epetra_Map((int_type) -1, (int) cols.size(), cols_ptr,
                               Map().IndexBase(), Map().Comm());
    }
    //now we need to create a matrix with sourceMap and colMap, and fill it with
    //our nonlocal data so we can then export it to the correct owning processors.

    std::vector<int> nonlocalRowLengths(nonlocalRows_var.size());
    for(size_t i=0; i<nonlocalRows_var.size(); ++i) {
      nonlocalRowLengths[i] = (int) nonlocalCols_var[i].size();
    }

    int* nlRLptr = nonlocalRowLengths.size()>0 ? &nonlocalRowLengths[0] : NULL;
    tempMat = new Epetra_CrsMatrix(Copy, *sourceMap_, *colMap_, nlRLptr);


    for(size_t i=0; i<nonlocalRows_var.size(); ++i) {
      EPETRA_CHK_ERR( tempMat->InsertGlobalValues(nonlocalRows_var[i],
                                                  (int) nonlocalCols_var[i].size(),
                                                  &nonlocalCoefs_[i][0],
                                                  &nonlocalCols_var[i][0]) );
    }

    if (!save_off_and_reuse_map_exporter) {
      delete sourceMap_;
      delete colMap_;
      sourceMap_ = colMap_ = NULL;
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

  if (exporter_ == NULL)
    exporter_ = new Epetra_Export(tempMat->RowMap(), RowMap());

  EPETRA_CHK_ERR(Export(*tempMat, *exporter_, combineMode));

  if(callFillComplete) {
    EPETRA_CHK_ERR(FillComplete(domain_map, range_map));
  }

  //now reset the values in our nonlocal data
  if (!useNonlocalMatrix_) {
    for(size_t i=0; i<nonlocalRows_var.size(); ++i) {
        nonlocalCols_var[i].resize(0); 
        nonlocalCoefs_[i].resize(0);
      }
    }

    delete tempMat;
  }

  if (!save_off_and_reuse_map_exporter) {
    delete exporter_;
    exporter_ = NULL;
  }
  return(0);
}

int Epetra_FECrsMatrix::GlobalAssemble(const Epetra_Map& domain_map,
                                       const Epetra_Map& range_map,
                                       bool callFillComplete,
                                       Epetra_CombineMode combineMode,
                                       bool save_off_and_reuse_map_exporter)
{
  if(!domain_map.GlobalIndicesTypeMatch(range_map))
    throw ReportError("Epetra_FECrsMatrix::GlobalAssemble: cannot be called with different indices types for domainMap and rangeMap", -1);

  if(!RowMap().GlobalIndicesTypeMatch(domain_map))
    throw ReportError("Epetra_FECrsMatrix::GlobalAssemble: cannot be called with different indices types for row map and incoming rangeMap", -1);

  if(RowMap().GlobalIndicesInt())
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    return GlobalAssemble<int>(domain_map, range_map, callFillComplete, combineMode, save_off_and_reuse_map_exporter);
#else
    throw ReportError("Epetra_FECrsMatrix::GlobalAssemble: ERROR, GlobalIndicesInt but no API for it.",-1);
#endif

  if(RowMap().GlobalIndicesLongLong())
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    return GlobalAssemble<long long>(domain_map, range_map, callFillComplete, combineMode, save_off_and_reuse_map_exporter);
#else
    throw ReportError("Epetra_FECrsMatrix::GlobalAssemble: ERROR, GlobalIndicesLongLong but no API for it.",-1);
#endif

  throw ReportError("Epetra_FECrsMatrix::GlobalAssemble: Internal error, unable to determine global index type of maps", -1);
}

//----------------------------------------------------------------------------
template<typename int_type>
int Epetra_FECrsMatrix::InputGlobalValues_RowMajor(
            int numRows, const int_type* rows,
            int numCols, const int_type* cols,
            const double* values,
            int mode)
{
  if(!RowMap().template GlobalIndicesIsType<int_type>())
  throw ReportError("Epetra_FECrsMatrix::InputGlobalValues_RowMajor mismatch between argument types (int/long long) and map type.", -1);

  int returncode = 0;
  int err = 0;

  for(int i=0; i<numRows; ++i) {
    double* valuesptr = (double*)values + i*numCols;

    int local_row_id = Map().LID(rows[i]);
    if (local_row_id >= 0) {
      switch(mode) {
        case Epetra_FECrsMatrix::SUMINTO:
          err = this->Epetra_CrsMatrix::SumIntoGlobalValues(rows[i], numCols,
              valuesptr, (int_type*)cols);
          if (err<0) return(err);
          if (err>0) returncode = err;
          break;
        case Epetra_FECrsMatrix::REPLACE:
          err = this->Epetra_CrsMatrix::ReplaceGlobalValues(rows[i], numCols,
              valuesptr, (int_type*)cols);
          if (err<0) return(err);
          if (err>0) returncode = err;
          break;
        case Epetra_FECrsMatrix::INSERT:
          err = this->Epetra_CrsMatrix::InsertGlobalValues(rows[i], numCols,
              valuesptr, (int_type*)cols);
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
template<typename int_type>
int Epetra_FECrsMatrix::InputGlobalValues(int numRows, const int_type* rows,
            int numCols, const int_type* cols,
            const double*const* values,
            int format, int mode)
{
  if(!RowMap().template GlobalIndicesIsType<int_type>())
  throw ReportError("Epetra_FECrsMatrix::InputGlobalValues mismatch between argument types (int/long long) and map type.", -1);

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
      valuesptr[j] = values[j][i];
    }

    returncode += InputGlobalValues_RowMajor(1, &rows[i], numCols, cols, valuesptr, mode);
    if (returncode < 0) return returncode;
  }

  return(returncode);
}

//----------------------------------------------------------------------------
template<typename int_type>
int Epetra_FECrsMatrix::InputGlobalValues(int numRows, const int_type* rows,
            int numCols, const int_type* cols,
            const double* values,
            int format, int mode)
{
  if(!RowMap().template GlobalIndicesIsType<int_type>())
  throw ReportError("Epetra_FECrsMatrix::InputGlobalValues mismatch between argument types (int/long long) and map type.", -1);

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
template<typename int_type>
int Epetra_FECrsMatrix::InputNonlocalGlobalValues(int_type row,
              int numCols, const int_type* cols,
              const double* values,
              int mode)
{
  if(!RowMap().template GlobalIndicesIsType<int_type>())
  throw ReportError("Epetra_FECrsMatrix::InputNonlocalGlobalValues mismatch between argument types (int/long long) and map type.", -1);

  // if we already have a nonlocal matrix object, this is easier...
  if (useNonlocalMatrix_) {
    int err, returncode = 0;
    double* valuesptr = (double*)values;
    switch(mode) {
    case Epetra_FECrsMatrix::SUMINTO:
      err = nonlocalMatrix_->SumIntoGlobalValues(row, numCols,
            valuesptr, (int_type*)cols);
      if (err<0) return(err);
      if (err>0) returncode = err;
      break;
    case Epetra_FECrsMatrix::REPLACE:
      err = nonlocalMatrix_->ReplaceGlobalValues(row, numCols,
            valuesptr, (int_type*)cols);
      if (err<0) return(err);
      if (err>0) returncode = err;
      break;
    case Epetra_FECrsMatrix::INSERT:
      err = nonlocalMatrix_->InsertGlobalValues(row, numCols,
           valuesptr, (int_type*)cols);
      if (err<0) return(err);
      if (err>0) returncode = err;
      break;
    default:
      cerr << "Epetra_FECrsMatrix: internal error, bad input mode."<<endl;
      return(-1);
    }
    return (returncode);
  }

  std::vector<int_type>& nonlocalRows_var = nonlocalRows<int_type>();

  //find offset of this row in our list of nonlocal rows.
  typename std::vector<int_type>::iterator it =
      std::lower_bound(nonlocalRows_var.begin(), nonlocalRows_var.end(), row);

  int rowoffset = (int) (it - nonlocalRows_var.begin());
  if (it == nonlocalRows_var.end() || *it != row) {
    EPETRA_CHK_ERR( InsertNonlocalRow(row, it) );
  }

  for(int i=0; i<numCols; ++i) {
    EPETRA_CHK_ERR( InputNonlocalValue(rowoffset, cols[i], values[i], mode) );
  }

  return(0);
}

//----------------------------------------------------------------------------
template<typename int_type>
int Epetra_FECrsMatrix::InsertNonlocalRow(int_type row, typename std::vector<int_type>::iterator iter)
{
  if(!RowMap().template GlobalIndicesIsType<int_type>())
  throw ReportError("Epetra_FECrsMatrix::InsertNonlocalRow mismatch between argument types (int/long long) and map type.", -1);

  std::vector<int_type>& nonlocalRows_var = nonlocalRows<int_type>();
  std::vector<std::vector<int_type> >& nonlocalCols_var = nonlocalCols<int_type>();

  int offset = (int) (iter - nonlocalRows_var.begin());
  nonlocalRows_var.insert(iter, row);
  typename std::vector<std::vector<int_type> >::iterator cols_iter = nonlocalCols_var.begin() + offset;
  nonlocalCols_var.insert(cols_iter, std::vector<int_type>());
  std::vector<std::vector<double> >::iterator coefs_iter = nonlocalCoefs_.begin() + offset;
  nonlocalCoefs_.insert(coefs_iter, std::vector<double>());

  return(0);
}

//----------------------------------------------------------------------------
template<typename int_type>
int Epetra_FECrsMatrix::InputNonlocalValue(int rowoffset,
             int_type col, double value,
             int mode)
{
  if(!RowMap().template GlobalIndicesIsType<int_type>())
  throw ReportError("Epetra_FECrsMatrix::InputNonlocalValue mismatch between argument types (int/long long) and map type.", -1);

  std::vector<int_type>& colIndices = nonlocalCols<int_type>()[rowoffset];
  std::vector<double>& coefs = nonlocalCoefs_[rowoffset];

  typename std::vector<int_type>::iterator it =
     std::lower_bound(colIndices.begin(), colIndices.end(), col);

  if (it == colIndices.end() || *it != col) {
    int offset = (int) (it - colIndices.begin());
    colIndices.insert(it, col);
    std::vector<double>::iterator dit = coefs.begin()+offset;
    coefs.insert(dit, value);
    return 0;
  }

  int coloffset = (int) (it - colIndices.begin());
  if (mode == SUMINTO || mode == INSERT) {
    coefs[coloffset] += value;
  }
  else {
    coefs[coloffset] = value;
  }

  return(0);
}

