/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_MatrixTraits_SSMat_hpp_
#define _fei_MatrixTraits_SSMat_hpp_

//This file defines matrix traits for SSMat matrices
//(well, "matrix-views" to be more precise).
//

#include "fei_SSVec.hpp"
#include "fei_SSMat.hpp"
#include "fei_Vector_Impl.hpp"

namespace fei {

  /** Specialization for SSMat. */
  template<>
  struct MatrixTraits<SSMat> {

    /** Return a string type-name for the underlying matrix */
    static const char* typeName()
      { return("SSMat"); }

    /** Set a specified scalar value throughout the matrix.
     */
    static int setValues(SSMat* mat, double scalar)
      {
	feiArray<SSVec*>& rows = mat->getRows();
	for(int i=0; i<rows.length(); ++i) {
	  if (rows[i] != NULL) {
	    rows[i]->coefs() = scalar;
	  }
	}
	return(0);
      }

    /** Query the number of rows. This is expected to be the number of rows
        on the local processor.
    */
    static int getNumLocalRows(SSMat* mat, int& numRows)
    {
      numRows = mat->getRowNumbers().length();
      return(0);
    }

    /** Given a global (zero-based) row number, query the length of that row.
     */
    static int getRowLength(SSMat* mat, int row, int& length)
      {
	SSVec* ssrow = mat->getRow(row);
	if (ssrow == NULL) {
	  length = 0; return(0);
	}

	length = ssrow->length();
	return( 0 );
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
    static int copyOutRow(SSMat* mat,
		      int row, int len, double* coefs, int* indices)
      {
	SSVec* ssrow = mat->getRow(row);
	if (ssrow == NULL) {
	  return(-1);
	}

	double* rcoefs = ssrow->coefs().dataPtr();
	int* rindices  = ssrow->indices().dataPtr();

	for(int i=0; i<len; ++i) {
	  coefs[i] = rcoefs[i];
	  indices[i] = rindices[i];
	}

	return( 0 );
      }

    /** Sum a C-style table of coefficient data into the underlying matrix.
     */
    static int putValuesIn(SSMat* mat,
                           int numRows, const int* rows,
                           int numCols, const int* cols,
                           const double* const* values,
                           bool sum_into)
      {
	if (numCols < 1 || numRows < 1) return(0);
	int err = 0;
        if (sum_into) {
          for(int i=0; i<numRows; ++i) {
            err += mat->sumInRow(rows[i], cols, values[i], numCols);
          }
        }
        else {
          for(int i=0; i<numRows; ++i) {
            err += mat->putRow(rows[i], cols, values[i], numCols);
          }
        }

	return( err );
      }

    /** Perform any necessary internal communications/synchronizations or other
        operations appropriate at end of data input. For some implementations this
        will be a no-op, so this "default implementation" will return 0.
    */
    static int globalAssemble(SSMat* mat)
    {
      return(0);
    }

    /** Compute the matrix-vector product y = A*x */
    static int matvec(SSMat* mat,
		      fei::Vector* x,
		      fei::Vector* y)
    {
      fei::Vector_Impl<SSVec>* ssx =
	dynamic_cast<fei::Vector_Impl<SSVec>* >(x);
      fei::Vector_Impl<SSVec>* ssy =
	dynamic_cast<fei::Vector_Impl<SSVec>* >(y);

      if (ssx == NULL || ssy == NULL) {
	return(-1);
      }

      SSVec* sx = ssx->getUnderlyingVector();
      SSVec* sy = ssy->getUnderlyingVector();

      return( mat->matVec(*sx, *sy) );
    }

  };//struct MatrixTraits
}//namespace fei

#endif // _fei_MatrixTraits_SSMat_hpp_

