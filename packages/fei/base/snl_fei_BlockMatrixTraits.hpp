/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _snl_fei_BlockMatrixTraits_hpp_
#define _snl_fei_BlockMatrixTraits_hpp_

#include <fei_macros.hpp>

namespace snl_fei {
  /** Internal implementation block-entry matrix traits. Define a "template"
      for accessing matrix data.
      Provide function stubs for default type "T", which will catch the
      use of any matrix type for which specialized traits have not been
      defined.
  */
  template<typename T>
  struct BlockMatrixTraits {

    /** Return a string type-name for the underlying matrix */
    static const char* typeName()
      { static const char name[] = "unsupported"; return(name); }

    /** Set a specified scalar value throughout the matrix.
     */
    static int putScalar(T* /*mat*/, double /*scalar*/)
      { return(-1); }

    /** Given a global (zero-based) row number, query the length of that row.
     */
    static int getRowLength(T* /*mat*/, int /*row*/, int& /*length*/)
    { return(-1); }

    /** Given a global (zero-based) point-row number, query the length
        of that row.
    */
    static int getPointRowLength(T* /*mat*/, int /*row*/, int& /*length*/)
    { return(-1); }

    /** Given a global (zero-based) row number, pass out a copy of the contents
        of that row.
        @param mat
        @param row
        @param numBlkCols Length of the user-allocated arrays indices and colDims.
        @param rowDim Number of point-equations associated with this block-row.
        @param blkCols User-allocated array which will hold column-indices on
        output.
        @param colDims User-allocated array which will hold the number of point-
        indices per block-column.
        @param coefs User-allocated array of arrays. First dimension (number of
        arrays) must be 'len'. i-th array will hold the coefficients from the i-th
        block-entry in this matrix row, packed in column-major order as a 1D list.
        @param coefsLen Length of the user-allocated coefs array.
        @param blkRowLength Output value, will be the length of the matrix row,
        (number of blk-cols) which may be more or less than the length of the
        input 'numBlkCols'. If blkRowLength is less than numBlkCols, then only
        'blkRowLength' positions in the above array arguments will be referenced.
        If blkRowLength is greater than numBlkCols, then only numBlkCols positions
        will be referenced.
        @return error-code 0 if successful. Non-zero return-value may indicate
        that the specified row is not locally owned.
    */
    static int copyOutRow(T* /*mat*/,
                          int /*row*/,
                          int /*numBlkCols*/,
                          int /*rowDim*/,
                          int* /*blkCols*/,
                          int* /*colDims*/,
                          double* /*coefs*/,
                          int /*coefsLen*/,
                          int& /*blkRowLength*/)
      { return(-1); }

    /** Given a global (zero-based) point-row number, pass out a copy of the
        contents of that row.
        @param mat
        @param firstLocalOffset First point-equation that is owned by the local
        processor.
        @param row Global equation-number of the point-row being requested.
        @param len Length of the user-allocated arrays coefs and indices.
        @param coefs User-allocated array which will hold matrix coefficients
        on output.
        @param indices User-allocated array which will hold column-indices on
        output.
        @param rowLength Output value, will be the length of the matrix row,
        which may be more or less than the length of the above user-allocated
        arrays.
        @return error-code 0 if successful. Non-zero return-value may indicate
        that the specified row is not locally owned.
    */
    static int copyOutPointRow(T* /*mat*/,
                               int /*firstLocalOffset*/,
                               int /*row*/, int /*len*/,
                               double* /*coefs*/, int* /*indices*/,
                               int& /*rowLength*/)
    { return(-1); }

    /** Sum a flat Fortran-style array of coefficient data into the
        underlying matrix.
    */
    static int sumIn(T* /*mat*/,
                     int /*blockRow*/,
                     int /*rowDim*/,
                     int /*numBlockCols*/,
                     const int* /*blockCols*/,
                     const int* /*colDims*/,
                     int /*LDA*/,
                     const double* /*values*/)
    { return(-1); }

    /** Copy a flat Fortran-style array of coefficient data into the
        underlying matrix.
    */
    static int copyIn(T* /*mat*/,
                     int /*blockRow*/,
                     int /*rowDim*/,
                     int /*numBlockCols*/,
                     const int* /*blockCols*/,
                     const int* /*colDims*/,
                     int /*LDA*/,
                     const double* /*values*/)
    { return(-1); }

    /** Sum a C-style table of coefficient data into the underlying matrix.
     */
    static int sumIn(T* /*mat*/,
                     int /*row*/, int /*rowDim*/,
                     int /*numCols*/, const int* /*cols*/,
                     const int* /*LDAs*/,
                     const int* /*colDims*/,
                     const double* const* /*values*/)
    { return(-1); }

    /** Copy (replacing any already-existing values at the specified locations)
        a C-style table of coefficient data into the underlying matrix.
    */
    static int copyIn(T* /*mat*/,
                      int /*row*/, int /*rowDim*/,
                      int /*numCols*/, const int* /*cols*/,
                      const int* /*LDAs*/,
                      const int* /*colDims*/,
                      const double* const* /*values*/)
      { return(-1); }

    /** Have the underlying matrix perform any global synchronization or
        assembly that needs to be done after all data has been input.
    */
    static int globalAssemble(T* /*mat*/)
    { return(-1); }

  };//struct BlockMatrixTraits
}//namespace snl_fei

#endif // _snl_fei_BlockMatrixTraits_hpp_
