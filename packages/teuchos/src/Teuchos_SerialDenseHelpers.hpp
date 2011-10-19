// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef _TEUCHOS_SERIALDENSEHELPERS_HPP_
#define _TEUCHOS_SERIALDENSEHELPERS_HPP_

/*! \file Teuchos_SerialDenseHelpers.hpp
  \brief Non-member helper functions on the templated serial, dense matrix/vector classes.
*/

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_DataAccess.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"	
#include "Teuchos_SerialSymDenseMatrix.hpp"	
#include "Teuchos_SerialDenseVector.hpp"	

namespace Teuchos {

/*! \relates SerialSymDenseMatrix 
  \brief A templated, non-member, helper function for computing the matrix triple-product:  B = alpha*W^T*A*W or B = alpha*W*A*W^T.
  
  \param transw - [in] Compute B = alpha*W^T*A*W if transw = Teuchos::TRANS, else compute B = alpha*W*A*W^T if transw = Teuchos::NO_TRANS.
  \param alpha - [in] The scaling factor.
  \param A - [in] SerialSymDenseMatrix
  \param W - [in] SerialDenseMatrix
  \param B - [out] SerialSymDenseMatrix
  
  \note The syntax for calling this function is:  <tt> Teuchos::symMatTripleProduct<int,double>( Teuchos::TRANS, alpha, A, W, B ) </tt>
*/
template<typename OrdinalType, typename ScalarType>
void symMatTripleProduct( ETransp transw, const ScalarType alpha, const SerialSymDenseMatrix<OrdinalType, ScalarType>& A, 
			  const SerialDenseMatrix<OrdinalType, ScalarType>& W, SerialSymDenseMatrix<OrdinalType, ScalarType>& B )
{
  // Local variables.
  // Note: dimemensions of W are obtained so we can compute W^T*A*W for either cases.
  OrdinalType A_nrowcols = A.numRows();  // A is a symmetric matrix and is assumed square. 
  OrdinalType B_nrowcols = (ETranspChar[transw]!='N') ? W.numCols() : W.numRows();
  OrdinalType W_nrows = (ETranspChar[transw]!='N') ? W.numRows() : W.numCols();
  OrdinalType W_ncols = (ETranspChar[transw]!='N') ? W.numCols() : W.numRows();
  
  bool isBUpper = B.upper();

  // Check for consistent dimensions.
  TEUCHOS_TEST_FOR_EXCEPTION( B_nrowcols != B.numRows(), std::out_of_range,
    "Teuchos::symMatTripleProduct<>() : "
    "Num Rows/Cols B (" << B.numRows() << ") inconsistent with W ("<< B_nrowcols << ")");
  TEUCHOS_TEST_FOR_EXCEPTION( A_nrowcols != W_nrows, std::out_of_range,
    "Teuchos::symMatTripleProduct<>() : "
    "Num Rows/Cols A (" << A_nrowcols << ") inconsistent with W ("<< W_nrows << ")");

  // Scale by zero, initialized B to zeros and return.
  if ( alpha == ScalarTraits<ScalarType>::zero() )
  {
    B.putScalar();
    return;
  }

  // Workspace.
  SerialDenseMatrix<OrdinalType, ScalarType> AW;

  // BLAS class.
  BLAS<OrdinalType, ScalarType> blas;
  ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
  ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();

  // Separate two cases because BLAS only supports symmetric matrix-matrix multply w/o transposes.
  if (ETranspChar[transw]!='N') {
    // Size AW to compute A*W
    AW.shapeUninitialized(A_nrowcols,W_ncols);
  
    // A*W
    AW.multiply( Teuchos::LEFT_SIDE, alpha, A, W, ScalarTraits<ScalarType>::zero() );
    
    // B = W^T*A*W
    if (isBUpper) {
      for (int j=0; j<B_nrowcols; ++j)
	blas.GEMV( transw, W_nrows, j+1, one, W.values(), W.stride(), AW[j], 1, zero, &B(0,j), 1 );
    }
    else {
      for (int j=0; j<B_nrowcols; ++j)
	blas.GEMV( transw, W_nrows, B_nrowcols-j, one, W[j], W.stride(), AW[j], 1, zero, &B(j,j), 1 );
    }
  }
  else {
    // Size AW to compute W*A
    AW.shapeUninitialized(W_ncols, A_nrowcols);

    // W*A
    AW.multiply( Teuchos::RIGHT_SIDE, alpha, A, W, ScalarTraits<ScalarType>::zero() );

    // B = W*A*W^T
    if (isBUpper) {
      for (int j=0; j<B_nrowcols; ++j)
	for (int i=0; i<=j; ++i) 
	  blas.GEMV( transw, 1, A_nrowcols, one, &AW(i,0), AW.stride(), &W(j,0), W.stride(), zero, &B(i,j), 1 );
    }
    else {
      for (int j=0; j<B_nrowcols; ++j)
	for (int i=j; i<B_nrowcols; ++i) 
	  blas.GEMV( transw, 1, A_nrowcols, one, &AW(i,0), AW.stride(), &W(j,0), W.stride(), zero, &B(i,j), 1 );
    }
  }

  return;
}
	
/*! \relates SerialDenseMatrix 
  \brief A templated, non-member, helper function for viewing or copying a column of a SerialDenseMatrix as a SerialDenseVector.
  
  \param CV - [in] Enumerated type set to Teuchos::Copy or Teuchos::View
  \param A - [in] SerialDenseMatrix
  \param col - [in] Integer indicating which column of A to return
  
  \note The syntax for calling this function is:  <tt>Teuchos::SerialDenseVector<int,double> col_j = Teuchos::getCol<int,double>( Teuchos::View, A, j )</tt>
*/
template<typename OrdinalType, typename ScalarType>
SerialDenseVector<OrdinalType,ScalarType>
getCol( DataAccess CV, SerialDenseMatrix<OrdinalType, ScalarType>& A, const OrdinalType col )
{
  return SerialDenseVector<OrdinalType, ScalarType>(CV, A[col], A.numRows());
}

/*! \relates SerialDenseMatrix 
  \brief A templated, non-member, helper function for setting a SerialDenseMatrix column using a SerialDenseVector.
  
  \param v - [in] SerialDenseVector
  \param col - [in] Integer indicating which column of A to replace with v 
  \param A - [out] SerialDenseMatrix
  
  \note The syntax for calling this function is:  bool err = Teuchos::setCol<int,double>( v, j, A )</tt>
*/
template<typename OrdinalType, typename ScalarType>
bool setCol( const SerialDenseVector<OrdinalType, ScalarType>& v,
             const OrdinalType col,
             SerialDenseMatrix<OrdinalType, ScalarType>& A )
{
  if (v.length() != A.numRows()) return false;

  std::copy(v.values(),v.values()+v.length(),A[col]);

  return true;
}

} // namespace Teuchos

#endif /* _TEUCHOS_SERIALDENSEHELPERS_HPP_ */
