
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

#ifndef _EPETRA_FECRSMATRIX_H_
#define _EPETRA_FECRSMATRIX_H_

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
    <li> Call the method GlobalAssemble(), which gathers all non-local data onto
    the owning processors as determined by the map provided at construction.
    </ul>

    Sub-matrix data, which is assumed to be a rectangular 'table' of
    coefficients accompanied by 'scatter-indices', can be provided in two forms:
    <ul>
    <li>Fortran-style packed 1-D array
    <li>C-style double-pointer, or list-of-rows.
    </ul>
    In both cases, a "format" parameter specifies whether the data is laid out
    in row-major or column-major order (i.e., whether coefficients for a row lie
    contiguously or whether coefficients for a column lie contiguously). See
    the documentation for the methods SumIntoGlobalValues() and
    ReplaceGlobalValues().
*/
class Epetra_FECrsMatrix : public Epetra_CrsMatrix {
  public:
  /** Constructor. */
   Epetra_FECrsMatrix(Epetra_DataAccess CV,
                   const Epetra_Map& RowMap,
                   int* NumEntriesPerRow);

   /** Constructor. */
   Epetra_FECrsMatrix(Epetra_DataAccess CV,
                   const Epetra_Map& RowMap,
                   int NumEntriesPerRow);

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
   */
   int GlobalAssemble();

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
};//class Epetra_FECrsMatrix

#endif //_EPETRA_FECRSMATRIX_H_
