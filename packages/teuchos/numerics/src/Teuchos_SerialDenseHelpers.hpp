// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "Teuchos_SerialBandDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"

#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_DefaultSerialComm.hpp"

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
      for (OrdinalType j=0; j<B_nrowcols; ++j)
	blas.GEMV( transw, W_nrows, j+1, one, W.values(), W.stride(), AW[j], 1, zero, &B(0,j), 1 );
    }
    else {
      for (OrdinalType j=0; j<B_nrowcols; ++j)
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
      for (OrdinalType j=0; j<B_nrowcols; ++j)
	for (OrdinalType i=0; i<=j; ++i)
	  blas.GEMV( transw, 1, A_nrowcols, one, &AW(i,0), AW.stride(), &W(j,0), W.stride(), zero, &B(i,j), 1 );
    }
    else {
      for (OrdinalType j=0; j<B_nrowcols; ++j)
	for (OrdinalType i=j; i<B_nrowcols; ++i)
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

/*! \related SerialDenseMatrix
  \brief A templated, non-member, helper function for generating a random SerialDenseMatrix that is synchronized in parallel.

  \param A - [in/out] SerialDenseMatrix to be filled with random numbers that are the same across processors
*/
template <typename OrdinalType, typename ScalarType>
void randomSyncedMatrix( Teuchos::SerialDenseMatrix<OrdinalType, ScalarType>& A )
{
  Teuchos::RCP<const Teuchos::Comm<OrdinalType> > comm;

#ifdef HAVE_MPI
  int mpiStarted = 0;
  MPI_Initialized(&mpiStarted);
  if (mpiStarted)
    comm = Teuchos::DefaultComm<OrdinalType>::getComm();
  else 
    comm = rcp(new Teuchos::SerialComm<OrdinalType>);
#else
  comm = Teuchos::DefaultComm<OrdinalType>::getComm();
#endif

  const OrdinalType procRank = rank(*comm);

  // Construct a separate serial dense matrix and synchronize it to get around
  // input matrices that are subviews of a larger matrix.
  Teuchos::SerialDenseMatrix<OrdinalType, ScalarType> newMatrix( A.numRows(), A.numCols() );
  if (procRank == 0)
    newMatrix.random();
  else
    newMatrix.putScalar( Teuchos::ScalarTraits<ScalarType>::zero() );

  broadcast(*comm, 0, A.numRows()*A.numCols(), newMatrix.values());

  // Assign the synchronized matrix to the input.
  A.assign( newMatrix );
}


/*! \relates SerialBandDenseMatrix
  \brief A templated, non-member, helper function for converting a SerialDenseMatrix to a SerialBandDenseMatrix.

  \param A - [in] SerialDenseMatrix to be converted
  \param kl - [in] Integer indicating desired lower bandwidth of band matrix.
  \param ku - [in] Integer indicating desired upper bandwidth of band matrix.
  \param factorFormat - [in] Bool indicating whether kl extra superdiagonals should be stored to be used by factorization.

  \note The syntax for calling this function is:  <tt>Teuchos::SerialBandDenseMatrix<int,double> AB = Teuchos::generalToBanded<int,double>( A, kl, ku, true )</tt>
*/
template<typename OrdinalType, typename ScalarType>
Teuchos::RCP<SerialBandDenseMatrix<OrdinalType, ScalarType> >
generalToBanded(const RCP<SerialDenseMatrix<OrdinalType,ScalarType> >& A,
		const OrdinalType kl, const OrdinalType ku,
		const bool factorFormat)
{
  OrdinalType m = A->numRows();
  OrdinalType n = A->numCols();

  // Check that the new matrix is consistent.
  TEUCHOS_TEST_FOR_EXCEPTION(A->values()==0, std::invalid_argument,
		     "SerialBandDenseSolver<T>::generalToBanded: A is an empty SerialDenseMatrix<T>!");
  TEUCHOS_TEST_FOR_EXCEPTION(kl<0 || kl>m, std::invalid_argument,
		     "SerialBandDenseSolver<T>::generalToBanded: The lower bandwidth kl is invalid!");
  TEUCHOS_TEST_FOR_EXCEPTION(ku<0 || ku>n, std::invalid_argument,
		     "SerialBandDenseSolver<T>::generalToBanded: The upper bandwidth ku is invalid!");

  OrdinalType extraBands = (factorFormat ? kl : 0);
  Teuchos::RCP<SerialBandDenseMatrix<OrdinalType, ScalarType> > AB =
    rcp( new SerialBandDenseMatrix<OrdinalType,ScalarType>(m,n,kl,extraBands+ku,true));

  for (OrdinalType j = 0; j < n; j++) {
    for (OrdinalType i=TEUCHOS_MAX(0,j-ku); i<=TEUCHOS_MIN(m-1,j+kl); i++) {
      (*AB)(i,j) = (*A)(i,j);
    }
  }
  return(AB);
}

/*! \relates SerialBandDenseMatrix
  \brief A templated, non-member, helper function for converting a SerialBandDenseMatrix to a SerialDenseMatrix.

  \param A - [in] SerialBandDenseMatrix to be converted

  \note The syntax for calling this function is:  <tt>Teuchos::SerialDenseMatrix<int,double> A = Teuchos::bandedToGeneral<int,double>( AB )</tt>
*/
template<typename OrdinalType, typename ScalarType>
Teuchos::RCP<SerialDenseMatrix<OrdinalType, ScalarType> >
bandedToGeneral(const RCP<SerialBandDenseMatrix<OrdinalType,ScalarType> >& AB)
{

  OrdinalType m = AB->numRows();
  OrdinalType n = AB->numCols();
  OrdinalType kl = AB->lowerBandwidth();
  OrdinalType ku = AB->upperBandwidth();

  // Check that the new matrix is consistent.
  TEUCHOS_TEST_FOR_EXCEPTION(AB->values()==0, std::invalid_argument,
		     "SerialBandDenseSolver<T>::bandedToGeneral: AB is an empty SerialBandDenseMatrix<T>!");

  Teuchos::RCP<SerialDenseMatrix<OrdinalType, ScalarType> > A = rcp( new SerialDenseMatrix<OrdinalType,ScalarType>(m,n) );
  for (OrdinalType j = 0; j < n; j++) {
    for (OrdinalType i=TEUCHOS_MAX(0,j-ku); i<=TEUCHOS_MIN(m-1,j+kl); i++) {
      (*A)(i,j) = (*AB)(i,j);
    }
  }
  return(A);
}

} // namespace Teuchos

#endif /* _TEUCHOS_SERIALDENSEHELPERS_HPP_ */
