// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_EXCEPTIONS_DECL_HPP
#define TPETRA_EXCEPTIONS_DECL_HPP

#include <Tpetra_ConfigDefs.hpp>
#include <stdexcept>

namespace Tpetra {

/// \file Tpetra_Exceptions_decl.hpp
/// \brief Declarations of Tpetra-specific exceptions.
///
/// \warning Users should not depend on the existence or interface of
///   anything in the Details namespace.
///
/// This file includes declarations of exceptions specific to Tpetra.
/// In particular, Details::InvalidGlobalIndex indicates an invalid
/// global index (as the name suggests), and
/// Details::InvalidGlobalRowIndex indicates an invalid global row
/// index (e.g., of a CrsGraph or CrsMatrix).  "Invalid" generally
/// means "not owned by the calling process."
///
/// We do not include Details::InvalidGlobalColumnIndex because no
/// Tpetra class or function currently throws this exception.  It
/// would be natural to add such a class (derived from
/// Details::InvalidGlobalIndex) if this were to change in the future.

/// \namespace Details
/// \brief Implementation details of Tpetra.
///
/// \warning Users must not rely on anything in this namespace.
namespace Details {

/// \class InvalidGlobalIndex
/// \brief Exception thrown by CrsMatrix on invalid global index.
///
/// \tparam GlobalOrdinal Same as the GlobalOrdinal template
///   parameter of Map, CrsGraph, CrsMatrix, MultiVector, etc.
template<class GlobalOrdinal>
class InvalidGlobalIndex : public std::domain_error {
public:
  /// \brief Constructor.
  ///
  /// \param msg [in] The exception message.
  /// \param globalIndex [in] The offending global index.
  InvalidGlobalIndex (const std::string& msg, const GlobalOrdinal globalIndex)
    : std::domain_error (msg), glInd_ (globalIndex) 
    {}
  
  //! The offending global index.
  GlobalOrdinal offendingIndex () const { return glInd_; }
  
private:
  //! The offending global index.
  const GlobalOrdinal glInd_; 
};

/// \class InvalidGlobalRowIndex
/// \brief Exception thrown by CrsMatrix on invalid global row index.
///
/// \tparam GlobalOrdinal Same as the GlobalOrdinal template
///   parameter of Map, CrsGraph, CrsMatrix, MultiVector, etc.
template<class GlobalOrdinal>
class InvalidGlobalRowIndex : public InvalidGlobalIndex<GlobalOrdinal> {
public:
  /// \brief Constructor.
  ///
  /// \param msg [in] The exception message.
  /// \param globalIndex [in] The offending global index.
  InvalidGlobalRowIndex (const std::string& msg, const GlobalOrdinal globalIndex)
    : InvalidGlobalIndex<GlobalOrdinal> (msg, globalIndex)
    {}
};

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_EXCEPTIONS_DECL_HPP
