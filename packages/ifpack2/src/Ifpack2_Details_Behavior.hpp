// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef IFPACK2_DETAILS_BEHAVIOR_HPP
#define IFPACK2_DETAILS_BEHAVIOR_HPP

/// \file Ifpack2_Details_Behavior.hpp
/// \brief Declaration of Ifpack2::Details::Behavior, a class that
///   describes Ifpack2's run-time behavior.

namespace Ifpack2 {
namespace Details {

/// \brief Description of Ifpack2's behavior.
///
/// "Behavior" means things like whether to do extra debug checks or
/// print debug output.  These depend both on build options and on
/// environment variables.  Build options generally control the
/// default behavior.
///
/// This class' methods have the following properties:
///
/// <ul>
/// <li>They may read the environment.</li>
/// <li>They read the environment at most once, on first call.</li>
/// <li>They are idempotent.</li>
/// <li>They are reentrant.</li>
/// </ul>
///
/// We intend for it to be inexpensive to call this class' methods
/// repeatedly.  The idea is that callers need not cache the results;
/// they may simply call the functions freely.
///
/// Ifpack2 does <i>not</i> promise to see changes to environment
/// variables made after first use of this class.  Best practice is to
/// set any desired environment variables before starting the executable.
///
/// Currently, this class supports the following environment variable:
///
/// <tt>IFPACK2_DEBUG</tt>: flags Ifpack2 to turn on debug checking and
/// debug output.
///
/// The environment variable is understood to be "on" or "off", for example:
///
/// <tt>IFPACK2_DEBUG=ON</tt>
/// <tt>IFPACK2_DEBUG=OFF</tt>
///
/// The default value of <tt>IFPACK2_DEBUG</tt> is ON if Ifpack2 is
/// configured with <tt>HAVE_IFPACK2_DEBUG</tt>, otherwise it is OFF.
class Behavior {
 public:
  /// \brief Whether Ifpack2 is in debug mode.
  ///
  /// "Debug mode" means that Ifpack2 may do extra error checks and may
  /// produce more detailed error messages or debug output.
  static bool debug();
};

}  // namespace Details
}  // namespace Ifpack2

#endif  // IFPACK2_DETAILS_BEHAVIOR_HPP