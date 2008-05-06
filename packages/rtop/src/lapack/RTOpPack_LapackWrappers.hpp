// @HEADER
// ***********************************************************************
// 
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//                Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef RTOPPACK_LAPACK_WRAPPERS_HPP
#define RTOPPACK_LAPACK_WRAPPERS_HPP


#include "RTOpPack_Types.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_as.hpp"


namespace RTOpPack {


/** \brief . */
enum ETransp {
  NOTRANS, TRANS, CONJTRANS
};
const int NUM_ETRANS_ARGS = 3;

/** \brief. */
extern const Teuchos::Tuple<char,NUM_ETRANS_ARGS> transpMap;


/** \brief Peform an in-place factorization of a square or rectangular matrix.
 *
 * \param A [in/out] On input, contains the entries of the square matrix.  On
 * output, contains the L and U factors.
 *
 * \param ipiv [in] On output, contains the pivots used in the factorization.
 * Note: This will be a 1-based valued array since this is a Fortran routine!
 *
 * \param rank [out] On output, gives the rank of the factorization.
 */
template<class Scalar>
void getrf(
  const SubMultiVectorView<Scalar> &A,
  const ArrayView<Teuchos_Index> &ipiv,
  const Ptr<Teuchos_Index> &rank
  );


/** \brief . */
template<class Scalar>
void getrs(
  const ConstSubMultiVectorView<Scalar> &A,
  const ArrayView<const Teuchos_Index> &ipiv,
  const ETransp transp,
  const Ptr<const SubMultiVectorView<Scalar> > &BX
  );


} // namespace RTOpPack


//
// Implementations
//


template<class Scalar>
void RTOpPack::getrf(
  const SubMultiVectorView<Scalar> &A,
  const ArrayView<Teuchos_Index> &ipiv,
  const Ptr<Teuchos_Index> &rank
  )
{
  using Teuchos::as;
  const Teuchos_Index maxRank = TEUCHOS_MIN( A.subDim(), A.numSubCols() );
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( A.subDim() == 0  );
  TEST_FOR_EXCEPT( A.numSubCols() == 0  );
  TEST_FOR_EXCEPT( is_null(A.values()) );
  TEUCHOS_ASSERT_EQUALITY( as<Teuchos_Index>(ipiv.size()), maxRank );
#endif

  Teuchos::LAPACK<Teuchos_Index, Scalar> lapack;
  Teuchos_Index info = -1;
  lapack.GETRF( A.subDim(), A.numSubCols(), A.values().get(), A.leadingDim(),
    &ipiv[0], &info );
  *rank = maxRank;
  TEST_FOR_EXCEPTION(
    info < 0, std::invalid_argument
    ,"getrf(...): Error, Invalid argument "
    << -info << " sent to LAPACK function xGETRF(...)" );
  if (info > 0)
    *rank = info - 1;
}


template<class Scalar>
void RTOpPack::getrs(
  const ConstSubMultiVectorView<Scalar> &A,
  const ArrayView<const Teuchos_Index> &ipiv,
  const ETransp transp,
  const Ptr<const SubMultiVectorView<Scalar> > &BX
  )
{
  using Teuchos::as;
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT( !is_null(BX) );
  TEUCHOS_ASSERT_EQUALITY( A.subDim(), BX->subDim() );
  TEUCHOS_ASSERT_EQUALITY( A.subDim(), A.numSubCols() );
  TEST_FOR_EXCEPT( A.subDim() == 0  );
  TEST_FOR_EXCEPT( A.numSubCols() == 0  );
  TEST_FOR_EXCEPT( is_null(A.values()) );
  TEUCHOS_ASSERT_EQUALITY( A.subDim(), ipiv.size() );
#endif
  Teuchos::LAPACK<Teuchos_Index, Scalar> lapack;
  Teuchos_Index info = -1;
  lapack.GETRS(
    transpMap[transp],
    A.subDim(), BX->numSubCols(), A.values().get(), A.leadingDim(),
    &ipiv[0], BX->values().get(), BX->leadingDim(), &info
    );
  TEST_FOR_EXCEPTION(
    info < 0, std::invalid_argument
    ,"getrs(...): Error, Invalid argument "
    << -info << " sent to LAPACK function xGETRS(...)" );
  // If we get here B is the solution to the linear system.
}


#endif // RTOPPACK_LAPACK_WRAPPERS_HPP
