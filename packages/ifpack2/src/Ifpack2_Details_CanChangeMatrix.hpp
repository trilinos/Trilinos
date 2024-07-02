// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_DETAILS_CANCHANGEMATRIX_HPP
#define IFPACK2_DETAILS_CANCHANGEMATRIX_HPP

/// \file Ifpack2_Details_CanChangeMatrix.hpp
/// \brief Declaration of interface for preconditioners that can
///   change their matrix after construction
///
/// \warning This file is an implementation detail of Ifpack2.  Users
///   must not depend on this file's name or contents.

#include <Ifpack2_ConfigDefs.hpp>
#include <Tpetra_RowMatrix_decl.hpp>

namespace Ifpack2 {
namespace Details {

  /// \class CanChangeMatrix
  /// \brief Mix-in interface for preconditioners that can change
  ///   their matrix after construction
  /// \tparam RowMatrixType Specialization of Tpetra::RowMatrix
  ///
  /// \warning This class is an implementation detail of Ifpack2.
  ///   Users must not depend on this class.
  ///
  /// An Ifpack2 preconditioner (a subclass of
  /// Ifpack2::Preconditioner) may inherit from this interface in
  /// order to express that users can change its matrix (the matrix
  /// that it preconditions) after construction.
  ///
  /// We introduced this class to facilitate gradual porting of
  /// Ifpack2 preconditioners to support this feature.  Originally,
  /// Ifpack2 preconditioners did not provide an interface to change
  /// their matrix after construction.  If a preconditioner (a
  /// subclass of Ifpack2::Preconditioner) inherits from this
  /// interface, that tells other Ifpack2 classes that the
  /// preconditioner accepts matrix changes.  This is particularly
  /// useful for preconditioners that implement NestedPreconditioner,
  /// since changing the inner preconditioner requires the
  /// NestedPreconditioner to give its new inner preconditioner the
  /// "inner matrix."  The inner matrix depends on the particular
  /// nested preconditioner.  For example, with AdditiveSchwarz, it is
  /// the result of LocalFilter on an overlap matrix.
  ///
  /// Changing the matrix puts the preconditioner back in an
  /// "pre-initialized" state.  You must first call initialize(), then
  /// compute(), before you may call apply() on this preconditioner.
  /// Depending on the implementation, it may be legal to set the
  /// matrix to null.  In that case, you may not call initialize() or
  /// compute() until you have subsequently set a nonnull matrix.
  template<class RowMatrixType>
  class CanChangeMatrix {
  public:
    /// \brief Set the new matrix.
    ///
    /// \param A [in] The new matrix.  This is const by Ifpack2
    ///   convention, for the same reason that the input matrix to any
    ///   Ifpack2 preconditioner's constructor is always const.
    ///
    /// Calling this method with a matrix different than the
    /// preconditioner's current matrix resets the preconditioner's
    /// state.  After calling this method with a nonnull input, you
    /// must first call initialize() and compute() (in that order)
    /// before you may call apply().
    ///
    /// You may call this method with a null input.  If A is null, then
    /// you may not call initialize() or compute() until you first call
    /// this method again with a nonnull input.  This method invalidates
    /// any previous factorization whether or not A is null, so calling
    /// setMatrix() with a null input is one way to clear the
    /// preconditioner's state (and free any memory that it may be
    /// using).
    ///
    /// The new matrix A need not necessarily have the same Maps or even
    /// the same communicator as the original matrix.
    virtual void
    setMatrix (const Teuchos::RCP<const RowMatrixType>& A) = 0;

    //! Destructor
    virtual ~CanChangeMatrix () {}
  };

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_CANCHANGEMATRIX_HPP
