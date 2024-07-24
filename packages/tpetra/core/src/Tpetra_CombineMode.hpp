// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_COMBINEMODE_HPP
#define TPETRA_COMBINEMODE_HPP

/// \file Tpetra_CombineMode.hpp
/// \brief Declaration of Tpetra::CombineMode enum, and a function for
///   setting a Tpetra::CombineMode parameter in a
///   Teuchos::ParameterList.

#include <string>

// Forward declaration of Teuchos::ParameterList.
namespace Teuchos {
  class ParameterList;
} // namespace Teuchos

namespace Tpetra {

  /// \enum CombineMode
  /// \brief Rule for combining data in an Import or Export
  ///
  /// Import or Export (data redistribution) operations might need to
  /// combine data received from other processes with existing data on
  /// the calling process.  This enum tells Tpetra how to do that for
  /// a specific Import or Export operation.  Each Tpetra object may
  /// interpret the CombineMode in a different way, so you should
  /// check the Tpetra object's documentation for details.
  ///
  /// Here is the list of supported combine modes:
  ///   - ADD: Sum new values 
  ///   - INSERT: Insert new values that don't currently exist
  ///   - REPLACE: Replace existing values with new values
  ///   - ABSMAX: If \f$x_{old}\f$ is the old value and \f$x_{new}\f$
  ///     the incoming new value, replace \f$x_{old}\f$ with
  ///     \f$\max\{ x_{old}, x_{new} \}\f$.
  ///   - ZERO: Replace old values with zero
  ///   - ADD_ASSIGN: Do addition assignment (+=) of values into existing values
  ///     \f$x_{target} = x_{target} + \sum(x_{source}\f$  
  ///     May not be supported in all classes
  ///
  /// ADD, REPLACE and ADD_ASSIGN are intended for modifying values that already
  /// exist.  Tpetra objects will generally work correctly if those
  /// values don't already exist.  (For example, ADD will behave like
  /// INSERT if the entry does not yet exist on the calling process.)
  /// However, performance may suffer.
  ///
  /// The ZERO combine mode is a special case that bypasses
  /// communication.  It may seem odd to include a "combine mode" that
  /// doesn't actually combine.  However, this is useful for
  /// computations like domain decomposition with overlap.  A ZERO
  /// combine mode with overlap is different than an ADD combine mode
  /// without overlap.  (See Ifpack2::AdditiveSchwarz, which inspired
  /// inclusion of this combine mode.)  Furthermore, Import and Export
  /// also encapsulate a local permutation; if you want only to
  /// execute the local permutation without communication, you may use
  /// the ZERO combine mode.
  enum CombineMode {
    ADD,     //!< Sum new values 
    INSERT,  //!< Insert new values that don't currently exist
    REPLACE, //!< Replace existing values with new values
    ABSMAX,  //!< Replace old value with maximum of magnitudes of old and new values
    ZERO,    //!< Replace old values with zero
    ADD_ASSIGN //!< Accumulate new values into existing values (may not be supported in all classes)
  };

  /// \brief Set CombineMode parameter in a Teuchos::ParameterList.
  ///
  /// If you are constructing a Teuchos::ParameterList with a
  /// CombineMode parameter, set the parameter by using this function.
  /// This will use a special feature of Teuchos -- custom parameter
  /// list validation -- so that users can specify CombineMode values
  /// by string, rather than enum value.  The strings are the same as
  /// the enum names: "ADD", "INSERT", "REPLACE", "ABSMAX", and
  /// "ZERO".  They are <i>not</i> case sensitive.
  ///
  /// Using this function to set a CombineMode parameter will ensure
  /// that the XML serialization of the resulting
  /// Teuchos::ParameterList will refer to the CombineMode enum values
  /// using human-readable string names, rather than raw integers.
  ///
  /// \param plist [out] Teuchos::ParameterList to which you want to
  ///   add the Tpetra::CombineMode parameter.
  ///
  /// \param paramName [in] String name to use for the parameter.  For
  ///   example, you might wish to call the parameter "Combine Mode",
  ///   "Tpetra::CombineMode", or "combine mode".  The parameter's
  ///   name <i>is</i> case sensitive, even though the string values
  ///   are <i>not</i>.
  void
  setCombineModeParameter (Teuchos::ParameterList& plist,
                           const std::string& paramName);

  //! Human-readable string representation of the given CombineMode.
  std::string combineModeToString (const CombineMode combineMode);

} // namespace Tpetra

#endif // TPETRA_COMBINEMODE_HPP
