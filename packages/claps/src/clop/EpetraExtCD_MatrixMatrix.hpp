/* Copyright (2001) Sandia Corportation. Under the terms of Contract
 *  * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 *   * work by or on behalf of the U.S. Government.  Export of this program
 *    * may require a license from the United States Government. */


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

#ifndef EPETRAEXTCD_MATRIXMATRIX_H
#define EPETRAEXTCD_MATRIXMATRIX_H

#include <Epetra_RowMatrix.h>
#include <Epetra_CrsMatrix.h>

namespace EpetraExtCD {

  /** Collection of matrix-matrix operations. This class basically
      functions as a namespace, containing only static methods.
      See the program epetraext/test/MatrixMatrix/cxx_main.cpp for
      a usage example.
   */
  class MatrixMatrix {
  public:
    /** destructor */
    virtual ~MatrixMatrix(){}

    /** Given Epetra_CrsMatrix objects A, B and C, form the product C = A*B.
	In a parallel setting, A and B need not have matching distributions,
	but C needs to have the same row-map as A.

    @param A Input, must already have had 'FillComplete()' called.
    @param B Input, must already have had 'FillComplete()' called.
    @param C Result. On entry to this method, it doesn't matter whether
             FillComplete() has already been called on C or not. If it has,
	     then C's graph must already contain all nonzero locations that
	     will be produced when forming the product A*B. On exit,
	     C.FillComplete() will have been called.

    @return error-code, 0 if successful. non-zero returns may result if A or
             B are not already Filled, or if errors occur in putting values
             into C, etc.
     */
    static int Multiply(const Epetra_CrsMatrix& A,
			const Epetra_CrsMatrix& B,
			Epetra_CrsMatrix* & C);

  };//class MatrixMatrix
}//namespace EpetraExtCD

#endif

