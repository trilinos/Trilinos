/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
//@HEADER
*/

#ifndef IFPACK2_CONDEST_HPP
#define IFPACK2_CONDEST_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_CondestType.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include <Teuchos_Ptr.hpp>

namespace Ifpack2 {

/// \fn Condest
/// \brief Estimate the condition number of the matrix.
///
/// The template parameters of this function are the same and occur in
/// the same order as the template parameters for the Preconditioner
/// class.
///
/// \warning This method is DEPRECATED.  It was inherited from Ifpack,
///   and Ifpack never clearly stated what this method computes.
///   Furthermore, Ifpack's method just estimates the condition number
///   of the matrix A, and ignores the preconditioner -- which is
///   probably not what users thought it did.  If there is sufficient
///   interest, we might reintroduce this method with a different
///   meaning and a better algorithm.
///
/// \param TIFP [in] The Ifpack2 preconditioner.  We need this if
///   <tt>matrix_in</tt> is null.
///
/// \param CT [in] The method to use for computing the condition
///   number estimate.  Currently, the only supported option is
///   <tt>Cheap</tt>.  Unsupported options will throw
///   <tt>std::logic_error</tt>.
///
/// \param MaxIters [in] The number of iterations used to compute the
///   condition number estimate.  Currently, this parameter is ignored.
///
/// \param Tol [in] The convergence tolerance used to compute the
///   condition number estimate.  Currently, this parameter is
///   ignored.
///
/// \param matrix_in [in] Pointer to a Tpetra::RowMatrix.  If nonnull,
///   estimate the condition number of this matrix.  If null, estimate
///   the condition number of TIFP's matrix (as returned by its
///   getMatrix() method).
///
/// The "Cheap" condition number estimate is just \f$\max_i |y_i|\f$,
/// where \f$y = A*[1, \dots, 1]^T\f$.  That is, if the input matrix
/// is \f$A\f$, we multiply it on the right by a vector of ones, and
/// return the infinity norm (maximum absolute value) of the result.
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
Condest (const Ifpack2::Preconditioner<Scalar, LocalOrdinal, GlobalOrdinal, Node>& TIFP,
         const Ifpack2::CondestType CT,
         const int MaxIters = 1550,
         const typename Teuchos::ScalarTraits<Scalar>::magnitudeType& Tol = Teuchos::as<Scalar> (1e-9),
         const Teuchos::Ptr<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& matrix_in = Teuchos::null)
{
  using Teuchos::Ptr;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType MT;
  typedef Teuchos::ScalarTraits<MT> STM;
  typedef Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> row_matrix_type;
  typedef Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> vec_type;

  MT condNumEst = -STS::magnitude( STS::one() );

  // Users may either provide a matrix for which to estimate the
  // condition number, or use the Preconditioner's built-in matrix.
  Ptr<const row_matrix_type> matrix = matrix_in;
  if (matrix_in == Teuchos::null) {
    matrix = TIFP.getMatrix ().ptr ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      matrix == Teuchos::null,
      std::logic_error,
      "Ifpack2::Condest: Both the input matrix (matrix_in) and the Ifpack2 "
      "preconditioner's matrix are null, so we have no matrix with which to "
      "compute a condition number estimate.  This probably indicates a bug "
      "in Ifpack2, since no Ifpack2::Preconditioner subclass should accept a"
      "null matrix.");
  }

  if (CT == Ifpack2::Cheap) {
    vec_type ones (TIFP.getDomainMap ()); // Vector of ones
    ones.putScalar (STS::one ());
    vec_type onesResult (TIFP.getRangeMap ()); // A*ones
    onesResult.putScalar (STS::zero ());
    TIFP.apply (ones, onesResult);
    condNumEst = onesResult.normInf (); // max (abs (A*ones))
    TEUCHOS_TEST_FOR_EXCEPTION(
      STM::isnaninf (condNumEst),
      std::runtime_error,
      "Ifpack2::Condest: $\\|A*[1, ..., 1]^T\\|_{\\infty}$ = " << condNumEst << " is NaN or Inf.");
  } else if (CT == Ifpack2::CG) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "Ifpack2::Condest: Condition number estimation using CG is not currently supported.");
  } else if (CT == Ifpack2::GMRES) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "Ifpack2::Condest: Condition number estimation using GMRES is not currently supported.");
  }
  return condNumEst;
}

}//namespace Ifpack2

#endif // IFPACK2_CONDEST_HPP

