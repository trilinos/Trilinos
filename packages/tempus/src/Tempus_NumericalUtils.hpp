//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_NumericalUtils_hpp
#define Tempus_NumericalUtils_hpp

#include "Teuchos_ScalarTraits.hpp"

#include "Tempus_config.hpp"

namespace Tempus {

/// Numerical Tolerance (approx. max. significant digits minus two)
template <typename Scalar>
const Scalar numericalTol()
{
  const Scalar numericalTol =
      Teuchos::ScalarTraits<Scalar>::eps() *
      typename Teuchos::ScalarTraits<Scalar>::magnitudeType(100.0);
  return numericalTol;
}

/// Test if value is approximately zero within tolerance.
template <typename Scalar>
bool approxZero(Scalar value,
                Scalar tol = Teuchos::ScalarTraits<Scalar>::sfmin())
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  return ST::magnitude(value) <= ST::magnitude(tol);
}

/// Test if values are approximately equal within the absolute tolerance.
template <typename Scalar>
bool approxEqualAbsTol(Scalar a, Scalar b, Scalar absTol)
{
  return Teuchos::ScalarTraits<Scalar>::magnitude(a - b) < std::abs(absTol);
}

/// Test if values are approximately equal within the relative tolerance.
template <typename Scalar>
bool approxEqual(Scalar a, Scalar b, Scalar relTol = numericalTol<Scalar>())
{
  const Scalar min = std::min(std::abs(a), std::abs(b));
  Scalar absTol    = min * std::abs(relTol);
  if (absTol == Scalar(0.0)) absTol = Teuchos::ScalarTraits<Scalar>::sfmin();
  return approxEqualAbsTol(a, b, absTol);
}

/// Test if values are approximately equal within the relative tolerance given a
/// scale.
template <typename Scalar>
bool approxEqualScale(Scalar a, Scalar b, Scalar scale,
                      Scalar relTol = numericalTol<Scalar>())
{
  return approxEqualAbsTol(a, b, scale * relTol);
}

}  // namespace Tempus

#endif  // Tempus_NumericalUtils_hpp
