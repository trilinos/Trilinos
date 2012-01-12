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
