/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#ifndef _fei_Matrix_Local_hpp_
#define _fei_Matrix_Local_hpp_

#include <fei_SharedPtr.hpp>
#include <fei_MatrixGraph.hpp>
#include <fei_Matrix.hpp>
#include <fei_SparseRowGraph.hpp>

#include <vector>

namespace fei {
class Matrix_Local : public fei::Matrix {
 public:
  Matrix_Local(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
               fei::SharedPtr<fei::SparseRowGraph> sparseRowGraph);

  virtual ~Matrix_Local();

  static fei::SharedPtr<fei::Matrix>
    create_Matrix_Local(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
                        bool blockEntry);

  const char* typeName();

  /** Method for supplying parameters
  */
  int parameters(const fei::ParameterSet& paramset);

  /** Method for supplying parameters
  */
  int parameters(int numParams, const char* const* paramStrings);

    /** Obtain the fei::MatrixGraph associated with this matrix */
    fei::SharedPtr<fei::MatrixGraph> getMatrixGraph() const;

    /** Set the fei::MatrixGraph associated with this matrix */
    void setMatrixGraph(fei::SharedPtr<fei::MatrixGraph> matrixGraph);

    /** Get the global number of rows in the matrix.
     */
    int getGlobalNumRows() const;

    /** Get the local number of rows in the matrix.
     */
    int getLocalNumRows() const;

   /** Get the length of a row of the matrix.
       @param row Global 0-based equation number
       @param length Output. Length of the row.
       @return error-code non-zero if any error occurs.
   */
    int getRowLength(int row, int& length) const;

    /** Set a specified scalar throughout the matrix. */
    int putScalar(double scalar);

   /** Obtain a copy of the coefficients and indices for a row of the matrix.
       @param row Global 0-based equation number
       @param coefs Caller-allocated array, length 'len', to be filled with
       coefficients
       @param indices Caller-allocated array, length 'len', to be filled with
       indices. (These indices will be global 0-based equation numbers.)
       @param len Length of the caller-allocated coefs and indices arrays
       @return error-code non-zero if any error occurs.
   */
    int copyOutRow(int row, int len, double* coefs, int* indices) const;

    /** Sum coefficients into the matrix, adding them to any coefficients that
        may already exist at the specified row/column locations.

        @param numRows
        @param rows
        @param numCols
        @param cols
        @param values
        @param format For compatibility with old FEI elemFormat...
        0 means row-wise or row-major, 3 means column-major.
        Others not recognized
     */
    int sumIn(int numRows, const int* rows,
                      int numCols, const int* cols,
                      const double* const* values,
                      int format=0);

    /** Copy coefficients into the matrix, overwriting any coefficients that
        may already exist at the specified row/column locations.

        @param numRows
        @param rows
        @param numCols
        @param cols
        @param values
        @param format For compatibility with old FEI elemFormat...
        0 means row-wise or row-major, 3 means column-major.
        Others not recognized
    */
    int copyIn(int numRows, const int* rows,
                       int numCols, const int* cols,
                       const double* const* values,
                      int format=0);

    /** Sum coefficients into the matrix, specifying row/column locations by
        identifier/fieldID pairs.
        @param fieldID Input. field-identifier for which data is being input.
        @param idType Input. The identifier-type of the identifiers.
        @param rowID Input. Identifier in row-space, for which data is being
        input.
        @param colID Input. Identifier in column-space, for which data is being
        input.
        @param data Input. C-style table of data. num-rows is the field-size
        (i.e., number of scalar components that make up the field) of 'fieldID',
        as is num-columns.
        @param format For compatibility with old FEI elemFormat...
        0 means row-wise or row-major, 3 means column-major. Others not recognized
        @return error-code 0 if successful
    */
    int sumInFieldData(int fieldID,
                               int idType,
                               int rowID,
                               int colID,
                               const double* const* data,
                               int format=0);

    /** Sum coefficients into the matrix, specifying row/column locations by
        identifier/fieldID pairs.
        @param fieldID Input. field-identifier for which data is being input.
        @param idType Input. The identifier-type of the identifiers.
        @param rowID Input. Identifier in row-space, for which data is being
        input.
        @param colID Input. Identifier in column-space, for which data is being
        input.
        @param data Input. 1-D list representing a packed table of data. Data may
        be backed in row-major or column-major order and this may be specified with
        the 'format' argument. The "table" of data is of size num-rows X num-columns
        and num-rows is the field-size (i.e., number of scalar components that
        make up the field) of 'fieldID', as is num-columns.
        @param format For compatibility with old FEI elemFormat...
        0 means row-wise or row-major, 3 means column-major. Others not recognized
        @return error-code 0 if successful
    */
    int sumInFieldData(int fieldID,
                               int idType,
                               int rowID,
                               int colID,
                               const double* data,
                               int format=0);

    /** Sum coefficients, associated with a connectivity-block that was
        initialized on the MatrixGraph object, into this matrix.

        @param blockID
        @param connectivityID
        @param values
        @param format For compatibility with old FEI elemFormat...
        0 means row-wise or row-major, 3 means column-major. Others not recognized
     */
    int sumIn(int blockID, int connectivityID,
                      const double* const* values,
                      int format=0);

    /** Perform any necessary internal communications/synchronizations or other
        operations appropriate at end of data input. For some implementations this
        will be a no-op.
    */
    int globalAssemble();

    /** Form a matrix-vector product y = 'this' * x
     */
    int multiply(fei::Vector* x,
                         fei::Vector* y);

    void setCommSizes();

    /** After local overlapping data has been input, (e.g., element-data for a
        finite-element application) call this method to have data that
        corresponds to shared identifiers be communicated from sharing-but-not-
        owning processors, to owning processors.
    */
    int gatherFromOverlap(bool accumulate = true);

    /** Write the matrix contents into the specified file.
        @param filename Text name of the file to be created or overwritten.
        If in a parallel environment, each processor will take turns writing
        into the file.
        @param matrixMarketFormat Optional argument, defaults to true. If true
        the contents of the file will be MatrixMarket real array format. If not
        true, the contents of the file will contain the matrix global
        dimensions on the first line, and all following lines will contain a
        space-separated triple with global row index first, global column index
        second and coefficient value third.
        Note also that if matrixMarketFormat is true, indices will be output in
        1-based form, but if not true, indices will be 0-based.
        @return error-code 0 if successful, -1 if some error occurs such as
        failure to open file.
     */
    int writeToFile(const char* filename,
                            bool matrixMarketFormat=true);

    /** Write the matrix contents into the specified ostream.
        @param ostrm ostream to be written to.
        @param matrixMarketFormat Optional argument, defaults to true. If true
        the data will be written in MatrixMarket real array format. If not
        true, the stream will receive the matrix global
        dimensions on the first line, and all following lines will contain a
        space-separated triple with global row index first, global column index
        second and coefficient value third.
        Note also that if matrixMarketFormat is true, indices will be output in
        1-based form, but if not true, indices will be 0-based.
        @return error-code 0 if successful, -1 if some error occurs.
     */
    int writeToStream(FEI_OSTREAM& ostrm,
                              bool matrixMarketFormat=true);

    /** Query whether the underlying matrix object is a block-entry matrix.
     */
    bool usingBlockEntryStorage();

    /** Set a "mark" point on the current state of the matrix, so that later
        a query can be made to see if the matrix has changed since this mark
        was set.
    */
    void markState();

    /** Query whether the matrix has changed since markState() was called. If
        markState() hasn't been called since the matrix was constructed, then
        this query will return true.
    */
    bool changedSinceMark();

    const std::vector<int>& getRowNumbers() const;

    const std::vector<int>& getRowOffsets() const;

    const std::vector<int>& getColumnIndices() const;

    const std::vector<double>& getCoefs() const;

 private:
  int getRowIndex(int rowNumber) const;

  int giveToMatrix(int numRows, const int* rows,
                      int numCols, const int* cols,
                      const double* const* values,
                      bool sumInto, int format);

  fei::SharedPtr<fei::MatrixGraph> matrixGraph_;
  fei::SharedPtr<fei::SparseRowGraph> sparseRowGraph_;
  std::vector<double> coefs_;
  bool stateChanged_;
  std::vector<double> work_data1D_;
  std::vector<const double*> work_data2D_;
};//class Matrix_Local
}//namespace fei

#endif

