
//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

#ifndef EPETRA_FECRSMATRIX_H
#define EPETRA_FECRSMATRIX_H

#include <Epetra_CrsMatrix.h>
class Epetra_Map;

/** Epetra Finite-Element CrsMatrix. This class provides the ability to
    input finite-element style sub-matrix data, including sub-matrices with
    non-local rows (which could correspond to shared finite-element nodes for
    example). This class inherits Epetra_CrsMatrix, and so all Epetra_CrsMatrix
    functionality is also available.

    It is intended that this class will be used as follows:
    <ul>
    <li> Construct with a map that describes a (non-overlapping) data
    distribution.
    <li> Input data, including non-local data, using the methods
    SumIntoGlobalValues() and/or ReplaceGlobalValues().
    <li> Call the method GlobalAssemble(), which gathers all non-local data
    onto the owning processors as determined by the map provided at
    construction. Users should note that the GlobalAssemble() method has an
    optional argument which determines whether GlobalAssemble() in turn calls
		FillComplete() after the data-exchange has occurred. If not explicitly
    supplied, this argument defaults to true.
    </ul>

    Sub-matrix data, which is assumed to be a rectangular 'table' of
    coefficients accompanied by 'scatter-indices', can be provided in two
    forms:
    <ul>
    <li>Fortran-style packed 1-D array
    <li>C-style double-pointer, or list-of-rows.
    </ul>
    In both cases, a "format" parameter specifies whether the data is laid out
    in row-major or column-major order (i.e., whether coefficients for a row
    lie contiguously or whether coefficients for a column lie contiguously).
    See the documentation for the methods SumIntoGlobalValues() and
    ReplaceGlobalValues().
*/
class Epetra_FECrsMatrix : public Epetra_CrsMatrix {
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

   /** Destructor. */
   virtual ~Epetra_FECrsMatrix();

   enum { ROW_MAJOR = 0, COLUMN_MAJOR = 3 };

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
       is an optional parameter, default value is ROW_MAJOR.
   */
   int SumIntoGlobalValues(int numIndices, const int* indices,
                           const double* values,
                           int format=Epetra_FECrsMatrix::ROW_MAJOR);

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
       is an optional parameter, default value is ROW_MAJOR.
   */
   int SumIntoGlobalValues(int numRows, const int* rows,
                           int numCols, const int* cols,
                           const double* values,
                           int format=Epetra_FECrsMatrix::ROW_MAJOR);

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
   int SumIntoGlobalValues(int numIndices, const int* indices,
                           const double* const* values,
                           int format=Epetra_FECrsMatrix::ROW_MAJOR);

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
   int SumIntoGlobalValues(int numRows, const int* rows,
	                   int numCols, const int* cols,
                           const double* const* values,
                           int format=Epetra_FECrsMatrix::ROW_MAJOR);

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
       is an optional parameter, default value is ROW_MAJOR.
   */
   int ReplaceGlobalValues(int numIndices, const int* indices,
                           const double* values,
                           int format=Epetra_FECrsMatrix::ROW_MAJOR);

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
       is an optional parameter, default value is ROW_MAJOR.
   */
   int ReplaceGlobalValues(int numRows, const int* rows,
                           int numCols, const int* cols,
                           const double* values,
                           int format=Epetra_FECrsMatrix::ROW_MAJOR);

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
   int ReplaceGlobalValues(int numIndices, const int* indices,
                           const double* const* values,
                           int format=Epetra_FECrsMatrix::ROW_MAJOR);

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
   int ReplaceGlobalValues(int numRows, const int* rows,
                           int numCols, const int* cols,
                           const double* const* values,
                           int format=Epetra_FECrsMatrix::ROW_MAJOR);

   /** Gather any overlapping/shared data into the non-overlapping partitioning
      defined by the Map that was passed to this matrix at construction time.
      Data imported from other processors is stored on the owning processor
      with a "sumInto" or accumulate operation.
      This is a collective method -- every processor must enter it before any
      will complete it.

      @param callFillComplete option argument, defaults to true.
        Determines whether GlobalAssemble() internally calls the
        FillComplete() method on this matrix.

      @return error-code 0 if successful, non-zero if some error occurs
   */
   int GlobalAssemble(bool callFillComplete=true);

   /** Set whether or not non-local data values should be ignored.
    */
   void setIgnoreNonLocalEntries(bool flag) {
     ignoreNonLocalEntries_ = flag;
   }

  private:
   int InputGlobalValues(int numRows, const int* rows,
                         int numCols, const int* cols,
                         const double* const* values,
                         int format,
                         bool accumulate);

   int InputGlobalValues(int numRows, const int* rows,
                         int numCols, const int* cols,
                         const double* values,
                         int format,
                         bool accumulate);

   int InputNonlocalGlobalValues(int row,
				 int numCols, const int* cols,
				 const double* values,
				 bool accumulate);

   int InsertNonlocalRow(int row, int offset);

   int InputNonlocalValue(int rowoffset,
			  int col, double value,
			  bool accumulate);

   int myFirstRow_;
   int myNumRows_;

   int numNonlocalRows_;
   int* nonlocalRows_;
   int* nonlocalRowLengths_;
   int* nonlocalRowAllocLengths_;
   int** nonlocalCols_;
   double** nonlocalCoefs_;

   double* workData_;
   int workDataLength_;
   bool ignoreNonLocalEntries_;
};//class Epetra_FECrsMatrix

#endif /* EPETRA_FECRSMATRIX_H */
