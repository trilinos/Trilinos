/*
//@HEADER
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

#ifndef IFPACK2_POWERMETHOD_HPP
#define IFPACK2_POWERMETHOD_HPP

/// \file Ifpack2_PowerMethod.hpp
/// \brief Definition of power methods
/// \author Graham Harper
///
/// This file describes power methods for use
/// throughout Ifpack2

#include "Kokkos_ArithTraits.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Tpetra_Details_residual.hpp"
#include <cmath>
#include <iostream>

namespace Ifpack2 {
namespace PowerMethod {

namespace { // (anonymous)

// Functor for making sure the real parts of all entries of a vector
// are nonnegative.  We use this in computeInitialGuessForPowerMethod
// below.
template<class OneDViewType,
         class LocalOrdinal = typename OneDViewType::size_type>
class PositivizeVector {
  static_assert (Kokkos::is_view<OneDViewType>::value,
                 "OneDViewType must be a 1-D Kokkos::View.");
  static_assert (static_cast<int> (OneDViewType::rank) == 1,
                 "This functor only works with a 1-D View.");
  static_assert (std::is_integral<LocalOrdinal>::value,
                 "The second template parameter, LocalOrdinal, "
                 "must be an integer type.");
public:
  PositivizeVector (const OneDViewType& x) : x_ (x) {}

  KOKKOS_INLINE_FUNCTION void
  operator () (const LocalOrdinal& i) const
  {
    typedef typename OneDViewType::non_const_value_type IST;
    typedef Kokkos::Details::ArithTraits<IST> STS;
    typedef Kokkos::Details::ArithTraits<typename STS::mag_type> STM;

    if (STS::real (x_(i)) < STM::zero ()) {
      x_(i) = -x_(i);
    }
  }

private:
  OneDViewType x_;
};

} // namespace (anonymous)



/// \brief Use the power method to estimate the maximum eigenvalue
///   of A*D_inv, given an initial guess vector x.
///
/// \param A [in] The Operator to use.
/// \param D_inv [in] Vector to use as implicit right scaling of A.
/// \param numIters [in] Maximum number of iterations of the power
///   method.
/// \param x [in/out] On input: Initial guess Vector for the power
///   method.  Its Map must be the same as that of the domain Map of
///   A.  This method may use this Vector as scratch space.
/// \param y [out] The resulting eigenvector.
/// \param tolerance [in] The relative eigenvalue tolerance. (default: 1e-7)
/// \param eigNormalizationFreq [in] The frequency of normalization. (default: 1)
/// \param out [in] The stream to send verbose output to. (default: null)
///
/// \return Estimate of the maximum eigenvalue of A*D_inv.
template<class OperatorType, class V>
typename V::scalar_type
powerMethodWithInitGuess (const OperatorType& A,
                          const V& D_inv,
                          const int numIters,
                          Teuchos::RCP<V> x,
                          Teuchos::RCP<V> y,
                          const typename Teuchos::ScalarTraits<typename V::scalar_type>::magnitudeType tolerance = 1e-7,
                          const int eigNormalizationFreq = 1,
                          Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::null)
{
  typedef typename V::scalar_type ST;
  typedef Teuchos::ScalarTraits<typename V::scalar_type> STS;
  typedef typename Teuchos::ScalarTraits<typename V::scalar_type>::magnitudeType MT;

  bool verbose = (out != Teuchos::null);

  using std::endl;
  if (verbose) {
    *out << " powerMethodWithInitGuess:" << endl;
  }

  const ST zero = static_cast<ST> (0.0);
  const ST one = static_cast<ST> (1.0);
  ST lambdaMax = zero;
  ST lambdaMaxOld = one;
  ST norm;

  norm = x->norm2 ();
  TEUCHOS_TEST_FOR_EXCEPTION
    (norm == zero, std::runtime_error,
     "Ifpack2::PowerMethod::powerMethodWithInitGuess: The initial guess "
     "has zero norm.  This could be either because Tpetra::Vector::"
     "randomize filled the vector with zeros (if that was used to "
     "compute the initial guess), or because the norm2 method has a "
     "bug.  The first is not impossible, but unlikely.");

  if (verbose) {
    *out << "  Original norm1(x): " << x->norm1 ()
          << ", norm2(x): " << norm << endl;
  }

  x->scale (one / norm);

  if (verbose) {
    *out << "  norm1(x.scale(one/norm)): " << x->norm1 () << endl;
  }

  if(y.is_null())
    y = Teuchos::rcp(new V(A.getRangeMap ()));

  for (int iter = 0; iter < numIters-1; ++iter) {
    if (verbose) {
      *out << "  Iteration " << iter << endl;
    }
    A.apply (*x, *y);
    x->elementWiseMultiply (STS::one(), D_inv, *y, STS::zero());

    if (((iter+1) % eigNormalizationFreq == 0) && (iter < numIters-2)) {
      norm = x->norm2 ();
      if (norm == zero) { // Return something reasonable.
        if (verbose) {
          *out << "   norm is zero; returning zero" << endl;
          *out << "   Power method terminated after "<< iter << " iterations." << endl;
        }
        return zero;
      } else {
        lambdaMaxOld = lambdaMax;
        lambdaMax = pow(norm, Teuchos::ScalarTraits<MT>::one() / eigNormalizationFreq);
        if (Teuchos::ScalarTraits<ST>::magnitude(lambdaMax-lambdaMaxOld) < tolerance * Teuchos::ScalarTraits<ST>::magnitude(lambdaMax)) {
          if (verbose) {
            *out << "  lambdaMax: " << lambdaMax << endl;
            *out << "  Power method terminated after "<< iter << " iterations." << endl;
          }
          return lambdaMax;
        } else if (verbose) {
          *out << "  lambdaMaxOld: " << lambdaMaxOld << endl;
          *out << "  lambdaMax: " << lambdaMax << endl;
          *out << "  |lambdaMax-lambdaMaxOld|/|lambdaMax|: " << Teuchos::ScalarTraits<ST>::magnitude(lambdaMax-lambdaMaxOld)/Teuchos::ScalarTraits<ST>::magnitude(lambdaMax) << endl;
        }
      }
      x->scale (one / norm);
    }
  }
  if (verbose) {
    *out << "  lambdaMax: " << lambdaMax << endl;
  }

  norm = x->norm2 ();
  if (norm == zero) { // Return something reasonable.
    if (verbose) {
      *out << "   norm is zero; returning zero" << endl;
      *out << "   Power method terminated after "<< numIters << " iterations." << endl;
    }
    return zero;
  }
  x->scale (one / norm);
  A.apply (*x, *y);
  y->elementWiseMultiply (STS::one(), D_inv, *y, STS::zero());
  lambdaMax = y->dot (*x);
  if (verbose) {
    *out << "  lambdaMax: " << lambdaMax << endl;
    *out << "  Power method terminated after "<< numIters << " iterations." << endl;
  }

  return lambdaMax;
}


/// \brief Fill x with random initial guess for power method
///
/// \param x [out] Initial guess vector; a domain Map vector of the
///   matrix.
/// \param nonnegativeRealParts [in] Whether to force all entries of
///   x (on output) to have nonnegative real parts.  Defaults to
///   false (don't force).
///
/// This is an implementation detail of powerMethod() below.  For a
/// justification of the second parameter, see Github Issues #64 and
/// #567.
template<class V>
void
computeInitialGuessForPowerMethod (V& x, const bool nonnegativeRealParts)
{
  typedef typename V::device_type::execution_space dev_execution_space;
  typedef typename V::local_ordinal_type LO;

  x.randomize ();

  if (nonnegativeRealParts) {
    auto x_lcl = x.getLocalViewDevice (Tpetra::Access::ReadWrite);
    auto x_lcl_1d = Kokkos::subview (x_lcl, Kokkos::ALL (), 0);

    const LO lclNumRows = static_cast<LO> (x.getLocalLength ());
    Kokkos::RangePolicy<dev_execution_space, LO> range (0, lclNumRows);
    PositivizeVector<decltype (x_lcl_1d), LO> functor (x_lcl_1d);
    Kokkos::parallel_for (range, functor);
  }
}


/// \brief Use the power method to estimate the maximum eigenvalue
///   of A*D_inv.
///
/// \param A [in] The Operator to use.
/// \param D_inv [in] Vector to use as implicit right scaling of A.
/// \param numIters [in] Maximum number of iterations of the power
///   method.
/// \param y [out] The resulting eigenvector.
/// \param tolerance [in] The relative eigenvalue tolerance. (default: 1e-7)
/// \param eigNormalizationFreq [in] The frequency of normalization. (default: 1)
/// \param out [in] The stream to send verbose output to. (default: null)
///
/// \return Estimate of the maximum eigenvalue of A*D_inv.
template<class OperatorType, class V>
typename V::scalar_type
powerMethod (const OperatorType& A,
             const V& D_inv, 
             const int numIters,
             Teuchos::RCP<V> y,
             const typename Teuchos::ScalarTraits<typename V::scalar_type>::magnitudeType tolerance = 1e-7,
             const int eigNormalizationFreq = 1,
             Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::null)
{
  typedef typename V::scalar_type ST;
  typedef Teuchos::ScalarTraits<typename V::scalar_type> STS;

  bool verbose = (out != Teuchos::null);

  if (verbose) {
    *out << "powerMethod:" << std::endl;
  }

  const ST zero = static_cast<ST> (0.0);

  Teuchos::RCP<V> x = Teuchos::rcp(new V(A.getDomainMap ()));
  // For the first pass, just let the pseudorandom number generator
  // fill x with whatever values it wants; don't try to make its
  // entries nonnegative.
  computeInitialGuessForPowerMethod (*x, false);

  ST lambdaMax = powerMethodWithInitGuess (A, D_inv, numIters, x, y, tolerance, eigNormalizationFreq, out);

  // mfh 07 Jan 2015: Taking the real part here is only a concession
  // to the compiler, so that this can build with ScalarType =
  // std::complex<T>.  Our Chebyshev implementation only works with
  // real, symmetric positive definite matrices.  The right thing to
  // do would be what Belos does, which is provide a partial
  // specialization for ScalarType = std::complex<T> with a stub
  // implementation (that builds, but whose constructor throws).
  if (STS::real (lambdaMax) < STS::real (zero)) {
    if (verbose) {
      *out << "real(lambdaMax) = " << STS::real (lambdaMax) << " < 0; "
        "try again with a different random initial guess" << std::endl;
    }
    // Max eigenvalue estimate was negative.  Perhaps we got unlucky
    // with the random initial guess.  Try again with a different (but
    // still random) initial guess.  Only try again once, so that the
    // run time is bounded.

    // For the second pass, make all the entries of the initial guess
    // vector have nonnegative real parts.
    computeInitialGuessForPowerMethod (*x, true);
    lambdaMax = powerMethodWithInitGuess (A, D_inv, numIters, x, y, tolerance, eigNormalizationFreq, out);
  }
  return lambdaMax;
}

} // namespace PowerMethod
} // namespace Ifpack2

#endif // IFPACK2_POWERMETHOD_HPP