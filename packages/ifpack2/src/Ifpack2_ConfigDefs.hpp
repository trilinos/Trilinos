// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _IFPACK2_CONFIGDEFS_HPP_
#define _IFPACK2_CONFIGDEFS_HPP_

#include <Ifpack2_config.h>
#include <Teuchos_ScalarTraits.hpp>
#include <Tpetra_ConfigDefs.hpp>

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
      JACOBI,  //!< Jacobi
      GS,      //!< Gauss-Seidel
      SGS,     //!< Symmetric Gauss-Seidel
      MTGS,    //!< Multicore Gauss-Seidel
      MTSGS,   //!< Multicore Symmetric Gauss-Seidel
      MTSPLITJACOBI, //!< Multicore split Jacobi; "split" refers to splitting A = D + R
      RICHARDSON, //!< Richardson
      GS2,     //!< Two-stage Gauss-Seidel
      SGS2     //!< Two-stage Symmetric Gauss-Seidel
    };
  } // namespace Details

  /*!
    @namespace Experimental
    @brief Ifpack2 features that are experimental.  Use at your own risk.
  */
  namespace Experimental {
  }

  /// \namespace DeprecatedAndMayDisappearAtAnyTime
  /// \brief Ifpack2 features that have been DEPRECATED and may
  ///   DISAPPEAR AT ANY TIME.  USE AT YOUR OWN RISK.
  namespace DeprecatedAndMayDisappearAtAnyTime {}

} // namespace Ifpack2

#endif /*_IFPACK2_CONFIGDEFS_HPP_*/
