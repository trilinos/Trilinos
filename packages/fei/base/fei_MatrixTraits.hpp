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


#ifndef _fei_MatrixTraits_hpp_
#define _fei_MatrixTraits_hpp_

#include <fei_macros.hpp>

namespace fei {
  class Vector;

  /** Define a struct of matrix access traits. The fei matrix implementation
      class fei::Matrix_Impl is essentially a filter which passes data to
      library-specific matrix objects (such as Trilinos/Epetra's
      Epetra_CrsMatrix).  fei::Matrix_Impl is a template, and the template
      parameter is the matrix object. In order to use an arbitrary matrix
      object with fei::Matrix_Impl, it is only necessary to define a
      specialization of this MatrixTraits struct for the matrix object.

      For an example specialization, see
        support-Trilinos/fei_MatrixTraits_Epetra.hpp.

      This "base" MatrixTraits struct provides function stubs for default
      type "T", which will catch the use of any matrix type for which
      specialized traits have not been defined.

      Note:
      Several functions accept row-numbers and/or column-numbers. These are
      always global indices. In a parallel setting, the fei defines a global
      set of row-indices which are partitioned across processors so that
      each processor owns a unique subset of the global row space.
      (Ownership of column-indices is not unique, a given column-index may
      appear in rows on multiple processors.)

      Note2:
      Implementers can safely assume that these functions will only be
      called with locally-owned row-numbers.
  */
  template<typename T>
  struct MatrixTraits {

    /** Return a string type-name for the underlying matrix. May appear in
	debug-output logs, etc. Does not need to exactly correspond to the type,
	but should be descriptive.
    */
    static const char* typeName()
      { return("unsupported"); }

    static double* getBeginPointer(T* /*mat*/)
      {
        return NULL;
      }

    static int getOffset(T* /*mat*/, int /*row*/, int /*col*/)
      {
        return -1;
      }

    /** Set a specified scalar value throughout the matrix.
     */
    static int setValues(T* mat, double scalar)
      { return(-1); }

    /** Query the number of local rows. This is expected to be the number of
        point-entry rows on the local processor.
    */
    static int getNumLocalRows(T* mat, int& numRows)
      { return(-1); }

    /** Given a locally-owned global row number, query the length (number of
        nonzeros) of that row.
     */
    static int getRowLength(T* mat, int row, int& length)
      { return(-1); }

    /** Given a locally-owned global row number, pass out a copy of the
        contents of that row.
	@param mat
	@param row Global row number
	@param len Length of the user-allocated arrays coefs and indices.
	@param coefs User-allocated array which will hold matrix coefficients
	on output.
	@param indices User-allocated array which will hold column-indices on
	output.
	@return error-code 0 if successful. Non-zero return-value may indicate
	that the specified row is not locally owned.
    */
    static int copyOutRow(T* mat,
		      int row, int len, double* coefs, int* indices)
      { return(-1); }

    /** Sum a C-style table of coefficient data into the underlying matrix.
	This is a rectangular array of coefficients for rows/columns defined by
	the 'rows' and 'cols' lists.
     */
    static int putValuesIn(T* mat,
                           int numRows, const int* rows,
                           int numCols, const int* cols,
                           const double* const* values,
                           bool sum_into)
      { return(-1); }

    /** Perform any necessary internal communications/synchronizations or other
        operations appropriate at end of data input. For some implementations
        this will be a no-op, so this "default implementation" returns 0. (The
        Trilinos/Epetra traits specialization calls A->FillComplete() at this
        point.)
    */
    static int globalAssemble(T* A)
    { return(0); }

    /** Compute the matrix-vector product y = A*x. It is expected that the
     underlying matrix object will form the local portion of the result 'y'
     on the local processor. */
    static int matvec(T* A, fei::Vector* x, fei::Vector* y)
    { return(-1); }
  };//struct MatrixTraits

}//namespace fei

#endif // _fei_MatrixTraits_hpp_

