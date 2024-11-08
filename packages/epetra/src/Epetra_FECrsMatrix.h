/*
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
*/

#ifndef EPETRA_FECRSMATRIX_H
#define EPETRA_FECRSMATRIX_H

#if defined(Epetra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Epetra package is deprecated"
#endif
#endif



#include <Epetra_ConfigDefs.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_CombineMode.h>

#include <vector>

class Epetra_Map;
class Epetra_IntSerialDenseVector;

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
class Epetra_LongLongSerialDenseVector;
#endif

class Epetra_SerialDenseMatrix;
class Epetra_FECrsGraph;

/** Epetra Finite-Element CrsMatrix. This class provides the ability to
    input finite-element style sub-matrix data, including sub-matrices with
    non-local rows (which could correspond to shared finite-element nodes for
    example). This class inherits Epetra_CrsMatrix, and so all Epetra_CrsMatrix
    functionality is also available.

    It is intended that this class will be used as follows:
    <ul>
    <li> Construct with either a map or graph that describes a (non-overlapping)
    data distribution.
    <li> Input data, including non-local data, using the methods
    InsertGlobalValues(), SumIntoGlobalValues() and/or ReplaceGlobalValues().
    <li> Call the method GlobalAssemble(), which gathers all non-local data
    onto the owning processors as determined by the map provided at
    construction. Users should note that the GlobalAssemble() method has an
    optional argument which determines whether GlobalAssemble() in turn calls
    FillComplete() after the data-exchange has occurred. If not explicitly
    supplied, this argument defaults to true.
    ***NOTE***: When GlobalAssemble() calls FillComplete(), it passes the
    arguments 'DomainMap()' and 'RangeMap()', which are the map attributes
    held by the base-class CrsMatrix and its graph. If a rectangular matrix
    is being assembled, the correct domain-map and range-map must be passed to
    GlobalAssemble (there are two overloadings of this method) -- otherwise, it
    has no way of knowing what these maps should really be.
    </ul>

    Sub-matrix data, which is assumed to be a rectangular 'table' of
    coefficients accompanied by 'scatter-indices', can be provided in three
    forms:
    <ul>
    <li>Fortran-style packed 1-D array.
    <li>C-style double-pointer, or list-of-rows.
    <li>Epetra_SerialDenseMatrix object.
    </ul>
    In all cases, a "format" parameter specifies whether the data is laid out
    in row-major or column-major order (i.e., whether coefficients for a row
    lie contiguously or whether coefficients for a column lie contiguously).
    See the documentation for the methods SumIntoGlobalValues() and
    ReplaceGlobalValues().

    Important notes:
    <ol>
    <li> Since Epetra_FECrsMatrix inherits Epetra_CrsMatrix, the semantics of
    the Insert/SumInto/Replace methods are the same as they are on
    Epetra_CrsMatrix, which is:
    <ul>
    <li>InsertGlobalValues() inserts values into the matrix only if the graph
    has not yet been finalized (FillComplete() has not yet been called). For
    non-local values, the call to InsertGlobalValues() may succeed but the
    GlobalAssemble() method may then fail because the non-local data is not
    actually inserted in the underlying matrix until GlobalAssemble() is called.
    <li>SumIntoGlobalValues() and ReplaceGlobalValues() only work for values
    that already exist in the matrix. In other words, these methods can not be
    used to put new values into the matrix.
    </ul>
    </ol>
*/
class EPETRA_LIB_DLL_EXPORT Epetra_FECrsMatrix : public Epetra_CrsMatrix {
  public:
  /** Constructor. */
   Epetra_FECrsMatrix(Epetra_DataAccess CV,
          const Epetra_Map& RowMap,
          int* NumEntriesPerRow,
          bool ignoreNonLocalEntries=false);

   /** Constructor. */
   Epetra_FECrsMatrix(Epetra_DataAccess CV,
          const Epetra_Map& RowMap,
          int NumEntriesPerRow,
          bool ignoreNonLocalEntries=false);

  /** Constructor. */
   Epetra_FECrsMatrix(Epetra_DataAccess CV,
          const Epetra_Map& RowMap,
          const Epetra_Map& ColMap,
          int* NumEntriesPerRow,
          bool ignoreNonLocalEntries=false);

   /** Constructor. */
   Epetra_FECrsMatrix(Epetra_DataAccess CV,
          const Epetra_Map& RowMap,
          const Epetra_Map& ColMap,
          int NumEntriesPerRow,
          bool ignoreNonLocalEntries=false);

   /** Constructor. */
   Epetra_FECrsMatrix(Epetra_DataAccess CV,
          const Epetra_CrsGraph& Graph,
          bool ignoreNonLocalEntries=false);

   /** Constructor. */
   Epetra_FECrsMatrix(Epetra_DataAccess CV,
         const Epetra_FECrsGraph& Graph,
         bool ignoreNonLocalEntries=false);

   /** Copy Constructor. */
   Epetra_FECrsMatrix(const Epetra_FECrsMatrix& src);

   /** Destructor. */
   virtual ~Epetra_FECrsMatrix();

   /** Assignment operator */
   Epetra_FECrsMatrix& operator=(const Epetra_FECrsMatrix& src);

   enum { ROW_MAJOR = 0, COLUMN_MAJOR = 3 };

   void Print(std::ostream& os) const;
#if !defined(EPETRA_NO_32BIT_GLOBAL_INDICES) || !defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
   using Epetra_CrsMatrix::SumIntoGlobalValues;
   using Epetra_CrsMatrix::InsertGlobalValues;
   using Epetra_CrsMatrix::ReplaceGlobalValues;
#endif

   /** override base-class Epetra_CrsMatrix::SumIntoGlobalValues method */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   int SumIntoGlobalValues(int GlobalRow, int NumEntries,
                           const double* Values, const int* Indices);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   int SumIntoGlobalValues(long long GlobalRow, int NumEntries,
                           const double* Values, const long long* Indices);
#endif

   /** override base-class Epetra_CrsMatrix::InsertGlobalValues method */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   int InsertGlobalValues(int GlobalRow, int NumEntries,
                           const double* Values, const int* Indices);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   int InsertGlobalValues(long long GlobalRow, int NumEntries,
                           const double* Values, const long long* Indices);
#endif

   /** override base-class Epetra_CrsMatrix::InsertGlobalValues method */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   int InsertGlobalValues(int GlobalRow, int NumEntries,
                           double* Values, int* Indices);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   int InsertGlobalValues(long long GlobalRow, int NumEntries,
                           double* Values, long long* Indices);
#endif

   /** override base-class Epetra_CrsMatrix::ReplaceGlobalValues method */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   int ReplaceGlobalValues(int GlobalRow, int NumEntries,
                           const double* Values, const int* Indices);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   int ReplaceGlobalValues(long long GlobalRow, int NumEntries,
                           const double* Values, const long long* Indices);
#endif

   /** Sum a Fortran-style table (single-dimensional packed-list) of
       coefficients into the matrix, adding them to any coefficients that
       may already exist at the specified row/column locations.

       @param numIndices Number of rows (and columns) in the sub-matrix.
       @param indices List of scatter-indices (rows and columns) for the
       sub-matrix.
       @param values List, length numIndices*numIndices. Square sub-matrix of
       coefficients, packed in a 1-D array. Data is packed either contiguously
       by row or by column, specified by the final parameter 'format'.
       @param format Specifies whether the data in 'values' is packed in
       column-major or row-major order. Valid values are
       Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR. This
       is an optional parameter, default value is COLUMN_MAJOR.
   */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   int SumIntoGlobalValues(int numIndices, const int* indices,
                           const double* values,
                           int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   int SumIntoGlobalValues(int numIndices, const long long* indices,
                           const double* values,
                           int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif

   /** Sum a Fortran-style table (single-dimensional packed-list) of
       coefficients into the matrix, adding them to any coefficients that
       may already exist at the specified row/column locations.

       @param numRows Number of rows in the sub-matrix.
       @param rows List of row-numbers (scatter-indices) for the sub-matrix.
       @param numCols Number of columns in the sub-matrix.
       @param cols List of column-numbers (scatter-indices) for the sub-matrix.
       @param values List, length numRows*numCols. Rectangular sub-matrix of
       coefficients, packed in a 1-D array. Data is packed either contiguously
       by row or by column, specified by the final parameter 'format'.
       @param format Specifies whether the data in 'values' is packed in
       column-major or row-major order. Valid values are
       Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR. This
       is an optional parameter, default value is COLUMN_MAJOR.
   */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   int SumIntoGlobalValues(int numRows, const int* rows,
                           int numCols, const int* cols,
                           const double* values,
                           int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   int SumIntoGlobalValues(int numRows, const long long* rows,
                           int numCols, const long long* cols,
                           const double* values,
                           int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif

   /** Sum C-style table (double-pointer, or list of lists) of coefficients
       into the matrix, adding them to any coefficients that
       may already exist at the specified row/column locations.

       @param numIndices Number of rows (and columns) in the sub-matrix.
       @param indices List of scatter-indices (rows and columns) for the
       sub-matrix.
       @param values Square sub-matrix of coefficients, provided in a 2-D
       array, or double-pointer.
       @param format Specifies whether the data in 'values' is packed in
       column-major or row-major order. Valid values are
       Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR. This
       is an optional parameter, default value is ROW_MAJOR.
   */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   int SumIntoGlobalValues(int numIndices, const int* indices,
                           const double* const* values,
                           int format=Epetra_FECrsMatrix::ROW_MAJOR);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   int SumIntoGlobalValues(int numIndices, const long long* indices,
                           const double* const* values,
                           int format=Epetra_FECrsMatrix::ROW_MAJOR);
#endif

   /** Sum C-style table (double-pointer, or list of lists) of coefficients
       into the matrix, adding them to any coefficients that
       may already exist at the specified row/column locations.

       @param numRows Number of rows in the sub-matrix.
       @param rows List of row-numbers (scatter-indices) for the sub-matrix.
       @param numCols Number of columns in the sub-matrix.
       @param cols List of column-numbers (scatter-indices) for the sub-matrix.
       @param values Rectangular sub-matrix of coefficients, provided in a 2-D
       array, or double-pointer.
       @param format Specifies whether the data in 'values' is packed in
       column-major or row-major order. Valid values are
       Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR. This
       is an optional parameter, default value is ROW_MAJOR.
   */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   int SumIntoGlobalValues(int numRows, const int* rows,
                     int numCols, const int* cols,
                           const double* const* values,
                           int format=Epetra_FECrsMatrix::ROW_MAJOR);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   int SumIntoGlobalValues(int numRows, const long long* rows,
                     int numCols, const long long* cols,
                           const double* const* values,
                           int format=Epetra_FECrsMatrix::ROW_MAJOR);
#endif

   /** Insert a Fortran-style table (single-dimensional packed-list) of
       coefficients into the matrix.

       @param numIndices Number of rows (and columns) in the sub-matrix.
       @param indices List of scatter-indices (rows and columns) for the
       sub-matrix.
       @param values List, length numIndices*numIndices. Square sub-matrix of
       coefficients, packed in a 1-D array. Data is packed either contiguously
       by row or by column, specified by the final parameter 'format'.
       @param format Specifies whether the data in 'values' is packed in
       column-major or row-major order. Valid values are
       Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR. This
       is an optional parameter, default value is COLUMN_MAJOR.
   */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   int InsertGlobalValues(int numIndices, const int* indices,
                           const double* values,
                           int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   int InsertGlobalValues(int numIndices, const long long* indices,
                           const double* values,
                           int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif

   /** Insert a Fortran-style table (single-dimensional packed-list) of
       coefficients into the matrix.

       @param numRows Number of rows in the sub-matrix.
       @param rows List of row-numbers (scatter-indices) for the sub-matrix.
       @param numCols Number of columns in the sub-matrix.
       @param cols List of column-numbers (scatter-indices) for the sub-matrix.
       @param values List, length numRows*numCols. Rectangular sub-matrix of
       coefficients, packed in a 1-D array. Data is packed either contiguously
       by row or by column, specified by the final parameter 'format'.
       @param format Specifies whether the data in 'values' is packed in
       column-major or row-major order. Valid values are
       Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR. This
       is an optional parameter, default value is COLUMN_MAJOR.
   */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   int InsertGlobalValues(int numRows, const int* rows,
                           int numCols, const int* cols,
                           const double* values,
                           int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   int InsertGlobalValues(int numRows, const long long* rows,
                           int numCols, const long long* cols,
                           const double* values,
                           int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif

   /** Insert a C-style table (double-pointer, or list of lists) of coefficients
       into the matrix.

       @param numIndices Number of rows (and columns) in the sub-matrix.
       @param indices List of scatter-indices (rows and columns) for the
       sub-matrix.
       @param values Square sub-matrix of coefficients, provided in a 2-D
       array, or double-pointer.
       @param format Specifies whether the data in 'values' is packed in
       column-major or row-major order. Valid values are
       Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR. This
       is an optional parameter, default value is ROW_MAJOR.
   */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   int InsertGlobalValues(int numIndices, const int* indices,
                           const double* const* values,
                           int format=Epetra_FECrsMatrix::ROW_MAJOR);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   int InsertGlobalValues(int numIndices, const long long* indices,
                           const double* const* values,
                           int format=Epetra_FECrsMatrix::ROW_MAJOR);
#endif

   /** Insert a C-style table (double-pointer, or list of lists) of coefficients
       into the matrix.

       @param numRows Number of rows in the sub-matrix.
       @param rows List of row-numbers (scatter-indices) for the sub-matrix.
       @param numCols Number of columns in the sub-matrix.
       @param cols List of column-numbers (scatter-indices) for the sub-matrix.
       @param values Rectangular sub-matrix of coefficients, provided in a 2-D
       array, or double-pointer.
       @param format Specifies whether the data in 'values' is packed in
       column-major or row-major order. Valid values are
       Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR. This
       is an optional parameter, default value is ROW_MAJOR.
   */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   int InsertGlobalValues(int numRows, const int* rows,
                     int numCols, const int* cols,
                           const double* const* values,
                           int format=Epetra_FECrsMatrix::ROW_MAJOR);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   int InsertGlobalValues(int numRows, const long long* rows,
                     int numCols, const long long* cols,
                           const double* const* values,
                           int format=Epetra_FECrsMatrix::ROW_MAJOR);
#endif

   /** Copy a Fortran-style table (single-dimensional packed-list) of
       coefficients into the matrix, replacing any coefficients that
       may already exist at the specified row/column locations.

       @param numIndices Number of rows (and columns) in the sub-matrix.
       @param indices List of scatter-indices (rows and columns) for the
       sub-matrix.
       @param values List, length numIndices*numIndices. Square sub-matrix of
       coefficients, packed in a 1-D array. Data is packed either contiguously
       by row or by column, specified by the final parameter 'format'.
       @param format Specifies whether the data in 'values' is packed in
       column-major or row-major order. Valid values are
       Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR. This
       is an optional parameter, default value is COLUMN_MAJOR.
   */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   int ReplaceGlobalValues(int numIndices, const int* indices,
                           const double* values,
                           int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   int ReplaceGlobalValues(int numIndices, const long long* indices,
                           const double* values,
                           int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif

   /** Copy Fortran-style table (single-dimensional packed-list) of coefficients
       into the matrix, replacing any coefficients that
       may already exist at the specified row/column locations.

       @param numRows Number of rows in the sub-matrix.
       @param rows List of row-numbers (scatter-indices) for the sub-matrix.
       @param numCols Number of columns in the sub-matrix.
       @param cols List, of column-numbers
       (scatter-indices) for the sub-matrix.
       @param values List, length numRows*numCols. Rectangular sub-matrix of
       coefficients, packed in a 1-D array. Data is packed either contiguously
       by row or by column, specified by the final parameter 'format'.
       @param format Specifies whether the data in 'values' is packed in
       column-major or row-major order. Valid values are
       Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR. This
       is an optional parameter, default value is COLUMN_MAJOR.
   */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   int ReplaceGlobalValues(int numRows, const int* rows,
                           int numCols, const int* cols,
                           const double* values,
                           int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   int ReplaceGlobalValues(int numRows, const long long* rows,
                           int numCols, const long long* cols,
                           const double* values,
                           int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif

   /** Copy C-style table (double-pointer, or list of lists) of coefficients
       into the matrix, replacing any coefficients that
       may already exist at the specified row/column locations.

       @param numIndices Number of rows (and columns) in the sub-matrix.
       @param indices List of scatter-indices (rows and columns) for the
       sub-matrix.
       @param values Square sub-matrix of coefficients, provided in a 2-D
       array, or double-pointer.
       @param format Specifies whether the data in 'values' is packed in
       column-major or row-major order. Valid values are
       Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR. This
       is an optional parameter, default value is ROW_MAJOR.
   */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   int ReplaceGlobalValues(int numIndices, const int* indices,
                           const double* const* values,
                           int format=Epetra_FECrsMatrix::ROW_MAJOR);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   int ReplaceGlobalValues(int numIndices, const long long* indices,
                           const double* const* values,
                           int format=Epetra_FECrsMatrix::ROW_MAJOR);
#endif

   /** Copy C-style table (double-pointer, or list of lists) of coefficients
       into the matrix, replacing any coefficients that
       may already exist at the specified row/column locations.

       @param numRows Number of rows in the sub-matrix.
       @param rows List of row-numbers (scatter-indices) for the sub-matrix.
       @param numCols Number of columns in the sub-matrix.
       @param cols List of column-numbers (scatter-indices) for the sub-matrix.
       @param values Rectangular sub-matrix of coefficients, provided in a 2-D
       array, or double-pointer.
       @param format Specifies whether the data in 'values' is packed in
       column-major or row-major order. Valid values are
       Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR. This
       is an optional parameter, default value is ROW_MAJOR.
   */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   int ReplaceGlobalValues(int numRows, const int* rows,
                           int numCols, const int* cols,
                           const double* const* values,
                           int format=Epetra_FECrsMatrix::ROW_MAJOR);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   int ReplaceGlobalValues(int numRows, const long long* rows,
                           int numCols, const long long* cols,
                           const double* const* values,
                           int format=Epetra_FECrsMatrix::ROW_MAJOR);
#endif

   /** Sum a square structurally-symmetric sub-matrix into the global matrix.
       For non-square sub-matrices, see the other overloading of this method.

       @param indices List of scatter-indices. indices.Length() must be the same
       as values.M() and values.N().

       @param values Sub-matrix of coefficients. Must be square.

       @param format Optional format specifier, defaults to COLUMN_MAJOR.
   */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   int SumIntoGlobalValues(const Epetra_IntSerialDenseVector& indices,
         const Epetra_SerialDenseMatrix& values,
         int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   int SumIntoGlobalValues(const Epetra_LongLongSerialDenseVector& indices,
         const Epetra_SerialDenseMatrix& values,
         int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif

   /** Sum a general sub-matrix into the global matrix.
       For square structurally-symmetric sub-matrices, see the other
       overloading of this method.

       @param rows List of row-indices. rows.Length() must be the same
       as values.M().

       @param cols List of column-indices. cols.Length() must be the same
       as values.N().

       @param values Sub-matrix of coefficients.

       @param format Optional format specifier, defaults to COLUMN_MAJOR.
   */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   int SumIntoGlobalValues(const Epetra_IntSerialDenseVector& rows,
         const Epetra_IntSerialDenseVector& cols,
         const Epetra_SerialDenseMatrix& values,
         int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   int SumIntoGlobalValues(const Epetra_LongLongSerialDenseVector& rows,
         const Epetra_LongLongSerialDenseVector& cols,
         const Epetra_SerialDenseMatrix& values,
         int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif

   /** Insert a square structurally-symmetric sub-matrix into the global matrix.
       For non-square sub-matrices, see the other overloading of this method.

       @param indices List of scatter-indices. indices.Length() must be the same
       as values.M() and values.N().

       @param values Sub-matrix of coefficients. Must be square.

       @param format Optional format specifier, defaults to COLUMN_MAJOR.
   */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   int InsertGlobalValues(const Epetra_IntSerialDenseVector& indices,
         const Epetra_SerialDenseMatrix& values,
         int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   int InsertGlobalValues(const Epetra_LongLongSerialDenseVector& indices,
         const Epetra_SerialDenseMatrix& values,
         int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif

   /** Insert a general sub-matrix into the global matrix.
       For square structurally-symmetric sub-matrices, see the other
       overloading of this method.

       @param rows List of row-indices. rows.Length() must be the same
       as values.M().

       @param cols List of column-indices. cols.Length() must be the same
       as values.N().

       @param values Sub-matrix of coefficients.

       @param format Optional format specifier, defaults to COLUMN_MAJOR.
   */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   int InsertGlobalValues(const Epetra_IntSerialDenseVector& rows,
         const Epetra_IntSerialDenseVector& cols,
         const Epetra_SerialDenseMatrix& values,
         int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   int InsertGlobalValues(const Epetra_LongLongSerialDenseVector& rows,
         const Epetra_LongLongSerialDenseVector& cols,
         const Epetra_SerialDenseMatrix& values,
         int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif

   /** Use a square structurally-symmetric sub-matrix to replace existing
       values in the global matrix.
       For non-square sub-matrices, see the other overloading of this method.

       @param indices List of scatter-indices. indices.Length() must be the same
       as values.M() and values.N().

       @param values Sub-matrix of coefficients. Must be square.

       @param format Optional format specifier, defaults to COLUMN_MAJOR.
   */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   int ReplaceGlobalValues(const Epetra_IntSerialDenseVector& indices,
         const Epetra_SerialDenseMatrix& values,
         int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   int ReplaceGlobalValues(const Epetra_LongLongSerialDenseVector& indices,
         const Epetra_SerialDenseMatrix& values,
         int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif

   /** Use a general sub-matrix to replace existing values.
       For square structurally-symmetric sub-matrices, see the other
       overloading of this method.

       @param rows List of row-indices. rows.Length() must be the same
       as values.M().

       @param cols List of column-indices. cols.Length() must be the same
       as values.N().

       @param values Sub-matrix of coefficients.

       @param format Optional format specifier, defaults to COLUMN_MAJOR.
   */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   int ReplaceGlobalValues(const Epetra_IntSerialDenseVector& rows,
         const Epetra_IntSerialDenseVector& cols,
         const Epetra_SerialDenseMatrix& values,
         int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   int ReplaceGlobalValues(const Epetra_LongLongSerialDenseVector& rows,
         const Epetra_LongLongSerialDenseVector& cols,
         const Epetra_SerialDenseMatrix& values,
         int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
#endif

   /** Gather any overlapping/shared data into the non-overlapping partitioning
      defined by the Map that was passed to this matrix at construction time.
      Data imported from other processors is stored on the owning processor
      with a "sumInto" or accumulate operation.
      This is a collective method -- every processor must enter it before any
      will complete it.

      ***NOTE***: When GlobalAssemble() calls FillComplete(), it passes the
      arguments 'DomainMap()' and 'RangeMap()', which are the map attributes
      held by the base-class CrsMatrix and its graph. If a rectangular matrix
      is being assembled, the domain-map and range-map must be specified by
      calling the other overloading of this method. Otherwise, GlobalAssemble()
      has no way of knowing what these maps should really be.


      @param callFillComplete option argument, defaults to true.
        Determines whether GlobalAssemble() internally calls the
        FillComplete() method on this matrix.

      @return error-code 0 if successful, non-zero if some error occurs
   */
   int GlobalAssemble(bool callFillComplete=true,
                      Epetra_CombineMode combineMode=Add,
                      bool save_off_and_reuse_map_exporter=false);

   /** Gather any overlapping/shared data into the non-overlapping partitioning
      defined by the Map that was passed to this matrix at construction time.
      Data imported from other processors is stored on the owning processor
      with a "sumInto" or accumulate operation.
      This is a collective method -- every processor must enter it before any
      will complete it.

      ***NOTE***: When GlobalAssemble() (the other overloading of this method)
      calls FillComplete(), it passes the arguments 'DomainMap()' and
      'RangeMap()', which are the map attributes already held by the base-class
      CrsMatrix and its graph. If a rectangular matrix is being assembled, the
      domain-map and range-map must be specified. Otherwise, GlobalAssemble()
      has no way of knowing what these maps should really be.


      @param domain_map user-supplied domain map for this matrix

      @param range_map user-supplied range map for this matrix

      @param callFillComplete option argument, defaults to true.
        Determines whether GlobalAssemble() internally calls the
        FillComplete() method on this matrix.

      @return error-code 0 if successful, non-zero if some error occurs
   */
   int GlobalAssemble(const Epetra_Map& domain_map,
                      const Epetra_Map& range_map,
                      bool callFillComplete=true,
                      Epetra_CombineMode combineMode=Add,
                      bool save_off_and_reuse_map_exporter=false);

   /** Set whether or not non-local data values should be ignored. By default,
       non-local data values are NOT ignored.
    */
   void setIgnoreNonLocalEntries(bool flag) {
     ignoreNonLocalEntries_ = flag;
   }

  private:
   void DeleteMemory();

   enum {SUMINTO = 0, REPLACE = 1, INSERT = 2};

   template<typename int_type>
   int InputGlobalValues(int numRows, const int_type* rows,
                         int numCols, const int_type* cols,
                         const double* const* values,
                         int format,
                         int mode);

   template<typename int_type>
   int InputGlobalValues(int numRows, const int_type* rows,
                         int numCols, const int_type* cols,
                         const double* values,
                         int format,
                         int mode);

   template<typename int_type>
   int InputNonlocalGlobalValues(int_type row,
         int numCols, const int_type* cols,
         const double* values,
         int mode);

  template<typename int_type>
  int InputGlobalValues_RowMajor(
            int numRows, const int_type* rows,
            int numCols, const int_type* cols,
            const double* values,
            int mode);

   template<typename int_type>
   int InsertNonlocalRow(int_type row, typename std::vector<int_type>::iterator offset);

   template<typename int_type>
   int InputNonlocalValue(int rowoffset,
        int_type col, double value,
        int mode);

   long long myFirstRow_;
   int myNumRows_;

   bool ignoreNonLocalEntries_;

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   std::vector<int> nonlocalRows_int_;
   std::vector<std::vector<int> > nonlocalCols_int_;
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   std::vector<long long> nonlocalRows_LL_;
   std::vector<std::vector<long long> > nonlocalCols_LL_;
#endif

   template<typename int_type> std::vector<int_type>& nonlocalRows();
   template<typename int_type> std::vector<std::vector<int_type> >& nonlocalCols();

   std::vector<std::vector<double> > nonlocalCoefs_;

   //IMPORTANT NOTE: The use of class-member work-data arrays is
   //**NOT** thread-safe.
   std::vector<double> workData_;
   std::vector<const double*> workData2d_;
   int workDataLength_;

   bool useNonlocalMatrix_;
   Epetra_CrsMatrix* nonlocalMatrix_;

   Epetra_Map* sourceMap_;
   Epetra_Map* colMap_;
   Epetra_Export* exporter_;
   Epetra_CrsMatrix* tempMat_;

   template<typename int_type>
   int SumIntoGlobalValues(int_type GlobalRow, int NumEntries, const double* values, const int_type* Indices);

   template<typename int_type>
   int GlobalAssemble(const Epetra_Map& domain_map,
                      const Epetra_Map& range_map,
                      bool callFillComplete=true,
                      Epetra_CombineMode combineMode=Add,
                      bool save_off_and_reuse_map_exporter=false);

};//class Epetra_FECrsMatrix

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
template<> inline std::vector<int>& Epetra_FECrsMatrix::nonlocalRows<int>()
{
  return nonlocalRows_int_;
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
template<> inline std::vector<long long>& Epetra_FECrsMatrix::nonlocalRows<long long>()
{
  return nonlocalRows_LL_;
}
#endif

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
template<> inline std::vector<std::vector<int> >& Epetra_FECrsMatrix::nonlocalCols<int>()
{
  return nonlocalCols_int_;
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
template<> inline std::vector<std::vector<long long> >& Epetra_FECrsMatrix::nonlocalCols<long long>()
{
  return nonlocalCols_LL_;
}
#endif


#endif /* EPETRA_FECRSMATRIX_H */
