/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_MatrixTraits_LinProbMgr_hpp_
#define _fei_MatrixTraits_LinProbMgr_hpp_

//This file defines matrix traits for LinearProblemManager matrix
//representations.
//

#include <fei_LinearProblemManager.hpp>

namespace fei {

  /** Specialization for LinearProblemManager. */
  template<>
  struct MatrixTraits<fei::LinearProblemManager> {

    /** Return a string type-name for the underlying matrix */
    static const char* typeName()
      { return("fei::LinearProblemManager"); }

    static double* getBeginPointer(fei::LinearProblemManager* /*mat*/)
      {
        return NULL;
      }

    static int getOffset(fei::LinearProblemManager* /*mat*/, int /*row*/, int /*col*/)
      {
        return -1;
      }

    /** Set a specified scalar value throughout the matrix.
     */
    static int setValues(fei::LinearProblemManager* mat, double scalar)
      {
	mat->setMatrixValues(scalar);
        return(0);
      }

    /** Query the number of rows. This is expected to be the number of rows
        on the local processor.
    */
    static int getNumLocalRows(fei::LinearProblemManager* mat, int& numRows)
    {
      numRows = mat->getLocalNumRows();
      return(0);
    }

    /** Given a global (zero-based) row number, query the length of that row.
     */
    static int getRowLength(fei::LinearProblemManager* mat, int row, int& length)
      {
	length = mat->getRowLength(row);
        if (length < 0) return(length);
        return(0);
      }

    /** Given a global (zero-based) row number, pass out a copy of the contents
        of that row.
        @param mat
        @param row
        @param len Length of the user-allocated arrays coefs and indices.
        @param coefs User-allocated array which will hold matrix coefficients
        on output.
        @param indices User-allocated array which will hold column-indices on
        output.
        @return error-code 0 if successful. Non-zero return-value may indicate
        that the specified row is not locally owned.
    */
    static int copyOutRow(fei::LinearProblemManager* mat,
		      int row, int len, double* coefs, int* indices)
      {
	return( mat->copyOutMatrixRow(row, len, coefs, indices) );
      }

    /** Sum a C-style table of coefficient data into the underlying matrix.
     */
    static int putValuesIn(fei::LinearProblemManager* mat,
                           int numRows, const int* rows,
                           int numCols, const int* cols,
                           const double* const* values,
                           bool sum_into)
      {
        return( mat->insertMatrixValues(numRows, rows,
                                        numCols, cols,
                                        values, sum_into) );
      }

    /** Perform any necessary internal communications/synchronizations or other
        operations appropriate at end of data input.
    */
    static int globalAssemble(fei::LinearProblemManager* mat)
    {
      return( mat->globalAssemble() );
    }

    /** Compute the matrix-vector product y = A*x */
    static int matvec(fei::LinearProblemManager* mat,
		      fei::Vector* x,
		      fei::Vector* y)
    {
      return( -1 );
    }
  };//struct MatrixTraits
}//namespace fei

#endif // _fei_MatrixTraits_LinProbMgr_hpp_

