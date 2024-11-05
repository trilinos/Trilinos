// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_SET_SCIENTIFIC_HPP
#define TEUCHOS_SET_SCIENTIFIC_HPP

#include <Teuchos_as.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <string>
#include <ios>

namespace Teuchos {

/// \class SetScientific
/// \brief Temporarily make an output stream use scientific
///   notation with sufficient precision.
///
/// On construction, apply the necessary flags to the given
/// output stream so that floating-point numbers are written in
/// scientific notation with precision (dependent on the Scalar
/// type) sufficient to ensure that they can be read in with the
/// same value.  On destruction, restore the original
/// (pre-construction) flags to the output stream.
///
/// This makes SetScientific good for scope-protected alteration
/// of the output stream's flags; no matter how the scope exits
/// (normally or by a thrown exception), the original flags will
/// be restored.  Hence, "temporarily" (or even "politely"): we
/// restore the original flags on scope exit.
///
/// \tparam Scalar A floating-point type, either real or
///   complex, for which Teuchos::ScalarTraits<Scalar> has a
///   specialization.
template<typename Scalar, const bool isFloatingPoint = ! Teuchos::ScalarTraits<Scalar>::isOrdinal>
class SetScientific;


// Partial specialization of SetScientific for floating-point types.
//
// This class currently requires that std::log10() take
// arguments of type Scalar.  This may be relaxed in the future
// if Teuchos::ScalarTraits gets its own log10() method.
template<typename Scalar>
class SetScientific<Scalar, true> {
 public:
  typedef Scalar scalar_type;

  SetScientific(std::ostream& out, int prec = -1):
    out_(out),
    originalFlags_(out.flags()),
    originalPrecision_(out.precision())
  {
    // Print floating-point values in scientific notation.
    out << std::scientific;

    if (prec == -1) prec = Teuchos::SetScientific<Scalar, true>::getDefaultPrecision();

    // Set the number of (decimal) digits after the decimal
    // point to print.
    out.precision(static_cast<std::streamsize>(prec));
  }

  static inline int getDefaultPrecision() {
    typedef Teuchos::ScalarTraits<scalar_type> STS;
    typedef typename STS::magnitudeType magnitude_type;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;

    // We're writing decimal digits, so compute the number of
    // digits we need to get reasonable accuracy when reading
    // values back in.
    //
    // There is actually an algorithm, due to Guy Steele (yes,
    // Java's Guy Steele) et al., for idempotent printing of
    // finite-length floating-point values.  We should actually
    // implement that algorithm, but I don't have time for that
    // now.  Currently, I just print no more than (one decimal
    // digit more than (the number of decimal digits justified
    // by the precision of magnitude_type)).
    //
    // We need to use STM's log10() rather than (say) std::log10
    // here, because STM::base() returns a magnitude_type, not
    // one of C++'s standard integer types.
    const magnitude_type numDecDigits = STM::t() * STM::log10 (STM::base());

    // Round and add one.  The cast to int should not overflow
    // unless STM::t() is _extremely_ large, so we don't need to
    // check for that case here.
    const magnitude_type one = STM::one();
    const magnitude_type two = one + one;
    // Cast from magnitude_type to int, since std::ostream's
    // precision() method expects an int input.
    const int prec = 1 +
      Teuchos::as<int>(magnitude_type((two*numDecDigits + one) / two));
    return prec;
  }

  ~SetScientific () {
    out_.flags (originalFlags_);
  }

 private:
   //! The output stream to which to apply flags.
   std::ostream& out_;

   //! The output stream's original flags.
   std::ios_base::fmtflags originalFlags_;

   //! The output stream's original precision.
   std::streamsize originalPrecision_;
};

//! Partial specialization of SetScientific for non-floating-point types.
template<class Scalar>
class SetScientific<Scalar, false> {
 public:
  typedef Scalar scalar_type;
  SetScientific(std::ostream&) {}
  ~SetScientific() {}
};

} // namespace Teuchos

#endif // TEUCHOS_SET_SCIENTIFIC_HPP
