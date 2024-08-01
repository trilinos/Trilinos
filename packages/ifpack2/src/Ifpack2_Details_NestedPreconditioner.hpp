// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_DETAILS_NESTEDPRECONDITIONER_HPP
#define IFPACK2_DETAILS_NESTEDPRECONDITIONER_HPP

/// \file Ifpack2_Details_NestedPreconditioner.hpp
/// \brief Declaration of interface for nested preconditioners
///
/// \warning This file is an implementation detail of Ifpack2.  Users
///   must not depend on this file's name or contents.

#include <Ifpack2_Preconditioner.hpp>

namespace Ifpack2 {
namespace Details {

  /// \class NestedPreconditioner
  /// \brief Mix-in interface for nested preconditioners
  /// \tparam PrecType Specialization of Ifpack2::Preconditioner
  ///
  /// \warning This class is an implementation detail of Ifpack2.
  ///   Users must not depend on this class.
  ///
  /// An Ifpack2 preconditioner (a subclass of
  /// Ifpack2::Preconditioner) may inherit from this interface in
  /// order to express the following
  /// <ol>
  /// <li> It has a nested structure (it is an "outer preconditioner"
  ///      that uses an "inner preconditioner") </li>
  /// <li> It can accept an arbitrary "inner preconditioner" at any
  ///      time. </li>
  /// </ol>
  /// For example, AdditiveSchwarz is an "outer preconditioner" that
  /// uses an "inner preconditioner" as the subdomain solver.
  ///
  /// When an outer preconditioner gets a new inner preconditioner, it
  /// must pass the "inner matrix" (for example, a local filter of the
  /// overlap matrix, in the case of AdditiveSchwarz) to the inner
  /// preconditioner.  This necessitates calling the inner
  /// preconditioner's setMatrix() method.  That puts the inner
  /// preconditioner in a "pre-initialized" state (since the structure
  /// of the matrix may have changed).  It is up to the outer
  /// preconditioner whether to call initialize() and compute() on the
  /// inner preconditioner.  As a result, users should consider
  /// setInnerPreconditioner() a collective.  Furthermore,
  /// implementations are responsible for making sure (by attempting a
  /// dynamic cast to CanChangeMatrix) that the inner preconditioner's
  /// matrix can be changed.
  ///
  /// \warning CIRCULAR DEPENDENCIES ARE FORBIDDEN.  You may NOT give
  ///   this object (<tt>*this</tt>) to itself as an inner solver.
  ///   You MAY use an inner solver of the same TYPE as this object
  ///   (as long as this makes sense mathematically), but it must be a
  ///   different instance.  Implementations are strongly encouraged
  ///   to check for this case.
  template<class PrecType>
  class NestedPreconditioner {
  public:
    virtual ~NestedPreconditioner() { }

    /// \brief Set the inner preconditioner.
    ///
    /// \param innerPrec [in/out] The inner preconditioner.  Its
    ///   matrix (if it has one) may be replaced by a matrix specified
    ///   by the outer (this) preconditioner.
    ///
    /// \warning CIRCULAR DEPENDENCIES ARE FORBIDDEN.  You may NOT
    ///   give this object (<tt>*this</tt>) to itself as an inner
    ///   solver.  You MAY use an inner solver of the same TYPE as
    ///   this object (as long as this makes sense mathematically),
    ///   but it must be a different instance.
    ///
    /// \pre <tt>&*innerPrec != this</tt>.
    virtual void
    setInnerPreconditioner (const Teuchos::RCP<PrecType>& innerPrec) = 0;
  };

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_NESTEDPRECONDITIONER_HPP
