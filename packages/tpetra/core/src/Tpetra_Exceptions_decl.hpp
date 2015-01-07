// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
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
