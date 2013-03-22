/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_MatrixTraits_FEData_hpp_
#define _fei_MatrixTraits_FEData_hpp_

//This file defines matrix traits for FiniteElementData matrices
//

#include <fei_FiniteElementData.hpp>

namespace fei {

  /** specialization for FiniteElementData */
  template<>
  struct MatrixTraits<FiniteElementData> {

    /** Return a string type-name for the underlying matrix */
    static const char* typeName()
      { return("FiniteElementData"); }

    static double* getBeginPointer(FiniteElementData* fed)
      {
        return NULL;
      }

    static int getOffset(FiniteElementData* /*fed*/, int row, int col)
      {
        return -1;
      }

    /** Set a specified scalar value throughout the matrix.
     */
    static int setValues(FiniteElementData* fed, double scalar)
      {
	return( -1 );
      }

    /** Query the number of rows. This is expected to be the number of rows
        on the local processor.
    */
    static int getNumLocalRows(FiniteElementData* fed, int& numRows)
    {
      numRows = -1;
      return(-1);
    }

    /** Given a global (zero-based) row number, query the length of that row.
     */
    static int getRowLength(FiniteElementData* fed, int row, int& length)
      {
	return( -1 );
      }

    /** Given a global (zero-based) row number, pass out a copy of the contents
        of that row.
        @param fed
        @param row
        @param len Length of the user-allocated arrays coefs and indices.
        @param coefs User-allocated array which will hold matrix coefficients
        on output.
        @param indices User-allocated array which will hold column-indices on
        output.
        @return error-code 0 if successful. Non-zero return-value may indicate
        that the specified row is not locally owned.
    */
    static int copyOutRow(FiniteElementData* fed,
		      int row, int len, double* coefs, int* indices)
      {
	return( -1 );
      }

    /** Sum a C-style table of coefficient data into the underlying matrix.
     */
    static int putValuesIn(FiniteElementData* fed,
		     int numRows, const int* rows,
		     int numCols, const int* cols,
		     const double* const* values,
                          bool sum_into)
      {
	return( -1 );
      }

    /** Perform any necessary internal communications/synchronizations or other
        operations appropriate at end of data input. For some implementations this
        will be a no-op, so this "default implementation" will return 0.
    */
    static int globalAssemble(FiniteElementData* fed)
    {
      return( fed->loadComplete() );
    }

    /** Compute the matrix-vector product y = A*x */
    static int matvec(FiniteElementData* fed,
		      fei::Vector* x,
		      fei::Vector* y)
    {
      return(-1);
    }

  };//struct MatrixTraits
}//namespace fei


#endif // _fei_MatrixTraits_FEData_hpp_
