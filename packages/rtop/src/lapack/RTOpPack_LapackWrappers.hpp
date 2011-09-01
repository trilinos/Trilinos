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
  const ArrayView<int> &ipiv,
  const Ptr<int> &rank
  );


/** \brief . */
template<class Scalar>
void getrs(
  const ConstSubMultiVectorView<Scalar> &A,
  const ArrayView<const int> &ipiv,
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
  const ArrayView<int> &ipiv,
  const Ptr<int> &rank
  )
{
  using Teuchos::as;
  const int maxRank = TEUCHOS_MIN( A.subDim(), A.numSubCols() );
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( A.subDim() == 0  );
  TEST_FOR_EXCEPT( A.numSubCols() == 0  );
  TEST_FOR_EXCEPT( is_null(A.values()) );
  TEUCHOS_ASSERT_EQUALITY( as<int>(ipiv.size()), maxRank );
#endif

  Teuchos::LAPACK<int, Scalar> lapack;
  int info = -1;
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
  const ArrayView<const int> &ipiv,
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
  Teuchos::LAPACK<int, Scalar> lapack;
  int info = -1;
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
