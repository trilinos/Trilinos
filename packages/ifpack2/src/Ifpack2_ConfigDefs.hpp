/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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

#ifndef _IFPACK2_CONFIGDEFS_HPP_
#define _IFPACK2_CONFIGDEFS_HPP_

#include <Ifpack2_config.h>
#include <Teuchos_ScalarTraits.hpp>
#include <Tpetra_ConfigDefs.hpp>

#if defined(HAVE_TPETRA_DEBUG) && ! defined(HAVE_IFPACK2_DEBUG)
#  define HAVE_IFPACK2_DEBUG 1
#endif

//The sgn function isn't well defined for complex.
//Is it correct to operate on the real part of x as is done below?
template<class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
IFPACK2_SGN(const Scalar& x)
{
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType magnitudeType;
  static const magnitudeType one = STS::magnitude(STS::one());
  return STS::real(x) > 0.0 ? magnitudeType(-one) : one;
}

/// \namespace Ifpack2
/// \brief Preconditioners and smoothers for Tpetra sparse matrices.
namespace Ifpack2 {
  typedef Tpetra::global_size_t global_size_t;

  /// \namespace Details
  /// \brief Ifpack2 implementation details
  ///
  /// This namespace contains implementation details of Ifpack2.
  /// It is <i>not</i> meant for users.  Users should not rely on
  /// anything in this namespace.
  namespace Details {
    /// \brief Type of relaxation to use.
    /// \warning This is an implementation detail of Ifpack2.  Users
    ///   should not depend on this enum or its integer values.  The
    ///   enum may disappear or change name, or the values may change
    ///   names or values or disappear, at any time.
    enum RelaxationType {
      JACOBI, //!< Jacobi
      GS,     //!< Gauss-Seidel
      SGS     //!< Symmetric Gauss-Seidel
    };
  } // namespace Details
} // namespace Ifpack2

#endif /*_IFPACK2_CONFIGDEFS_HPP_*/
