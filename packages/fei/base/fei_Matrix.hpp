/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_Matrix_hpp_
#define _fei_Matrix_hpp_

#include "fei_iosfwd.hpp"
#include "fei_SharedPtr.hpp"
#include "fei_MatrixGraph.hpp"
#include "fei_defs.h"

namespace fei {
  /** Abstract representation of an algebraic matrix. This representation does
      not require that data be accessed only on the 'owning' processor. In
      other words, this representation may be used with an overlapping data
      decomposition. In most cases the underlying library-specific matrix will
      have a non-overlapping data decomposition (each equation uniquely owned
      by a single processor). Overlapping data may be assembled into this
      abstract matrix locally, and will be funneled into the underlying non-
      overlapping matrix on the correct processor when the gatherFromOverlap()
      method is called. Conversely, if the user wants to retrieve overlapping
      data from the matrix locally, that data is not guaranteed to be available
      until the scatterToOverlap() method is called.
  */
  class Matrix {
  public:
    /** Matrix Factory interface */
    class Factory {
    public:
      /** Usual virtual destructor. */
      virtual ~Factory(){}

      /** Produce an instance of a Matrix. */
      virtual fei::SharedPtr<fei::Matrix>
        createMatrix(fei::SharedPtr<fei::MatrixGraph> matrixGraph) = 0;
    };

    /** Virtual destructor. */
    virtual ~Matrix(){}

    /** Return an implementation-dependent name describing the run-time type
        of this object.
    */
    virtual const char* typeName() = 0;

    /** Method for supplying parameters
     */
    virtual int parameters(const fei::ParameterSet& paramset) = 0;

    /** Obtain the fei::MatrixGraph associated with this matrix */
    virtual fei::SharedPtr<fei::MatrixGraph> getMatrixGraph() const = 0;

    /** Set the fei::MatrixGraph associated with this matrix */
    virtual void setMatrixGraph(fei::SharedPtr<fei::MatrixGraph> matrixGraph) = 0;

    /** Get the global number of rows in the matrix.
     */
    virtual int getGlobalNumRows() const = 0;

    /** Get the local number of rows in the matrix.
     */
    virtual int getLocalNumRows() const = 0;

    /** Get the length of a row of the matrix.
        @param row Global 0-based equation number
        @param length Output. Length of the row.
        @return error-code non-zero if any error occurs.
    */
    virtual int getRowLength(int row, int& length) const = 0;

    /** Set a specified scalar throughout the matrix. */
    virtual int putScalar(double scalar) = 0;

   /** Obtain a copy of the coefficients and indices for a row of the matrix.
       @param row Global 0-based equation number
       @param coefs Caller-allocated array, length 'len', to be filled with
       coefficients
       @param indices Caller-allocated array, length 'len', to be filled with
       indices. (These indices will be global 0-based equation numbers.)
       @param len Length of the caller-allocated coefs and indices arrays
       @return error-code non-zero if any error occurs.
   */
    virtual int copyOutRow(int row, int len, double* coefs, int* indices) const = 0;

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
    virtual int sumIn(int numRows, const int* rows,
                      int numCols, const int* cols,
                      const double* const* values,
                      int format=0) = 0;

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
    virtual int copyIn(int numRows, const int* rows,
                       int numCols, const int* cols,
                       const double* const* values,
                      int format=0) = 0;

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
    virtual int sumInFieldData(int fieldID,
                               int idType,
                               int rowID,
                               int colID,
                               const double* const* data,
                               int format=0) = 0;

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
    virtual int sumInFieldData(int fieldID,
                               int idType,
                               int rowID,
                               int colID,
                               const double* data,
                               int format=0) = 0;

    /** Sum coefficients, associated with a connectivity-block that was
        initialized on the MatrixGraph object, into this matrix.

        @param blockID
        @param connectivityID
        @param values
        @param format For compatibility with old FEI elemFormat...
        0 means row-wise or row-major, 3 means column-major. Others not recognized
     */
    virtual int sumIn(int blockID, int connectivityID,
                      const double* const* values,
                      int format=0) = 0;

    /** Perform any necessary internal communications/synchronizations or other
        operations appropriate at end of data input. For some implementations this
        will be a no-op.
    */
    virtual int globalAssemble() = 0;

    /** Form a matrix-vector product y = 'this' * x
     */
    virtual int multiply(fei::Vector* x,
                         fei::Vector* y) = 0;

    /** perform initial communication to establish message sizes that will
      be needed for exchanging shared-node data.
      Called from within gatherFromOverlap usually, doesn't usually need to
      be explicitly called by client code. (Power users only...)
    */
    virtual void setCommSizes() = 0;

    /** After local overlapping data has been input, (e.g., element-data for a
        finite-element application) call this method to have data that 
        corresponds to shared identifiers be communicated from sharing-but-not-
        owning processors, to owning processors.
    */
    virtual int gatherFromOverlap(bool accumulate = true) = 0;

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
    virtual int writeToFile(const char* filename,
                            bool matrixMarketFormat=true) = 0;

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
    virtual int writeToStream(FEI_OSTREAM& ostrm,
                              bool matrixMarketFormat=true) = 0;

    /** Query whether the underlying matrix object is a block-entry matrix.
     */
    virtual bool usingBlockEntryStorage() = 0;

    /** Set a "mark" point on the current state of the matrix, so that later
        a query can be made to see if the matrix has changed since this mark
        was set.
    */
    virtual void markState() = 0;

    /** Query whether the matrix has changed since markState() was called. If
        markState() hasn't been called since the matrix was constructed, then
        this query will return true.
    */
    virtual bool changedSinceMark() = 0;

    virtual double* getBeginPointer() { return NULL; }
    virtual int getOffset(int row, int col) { return -1; }

  };//class Matrix
}//namespace fei

#ifndef _fei_ostream_ops_hpp_
#include <fei_ostream_ops.hpp>
#endif

#endif // _fei_Matrix_hpp_
