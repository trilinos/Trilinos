// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_TPETRA_ADAPTER_HPP
#define BELOS_TPETRA_ADAPTER_HPP

/// \file BelosTpetraAdapter.hpp
/// \brief Partial specialization of Belos::MultiVecTraits and
///   Belos::OperatorTraits for Tpetra objects.
///
/// \section Belos_TpetraAdapter_sum Summary
///
/// If you want to use Belos solvers with Tpetra objects, include this
/// header file, along with the header file(s) for the solver(s) you
/// want to use.  "Tpetra objects" means the following:
///   - Tpetra::MultiVector for the multivector type (MV)
///   - Tpetra::Operator for the operator type (OP)
///
/// You may use any subclass of Tpetra::Operator here, as long as its
/// template parameters match those of the Tpetra::MultiVector type.
/// Many different Trilinos packages implement Tpetra::Operator
/// subclasses.  For example, when solving a linear system Ax=b, you
/// could use a Tpetra::CrsMatrix or Tpetra::RowMatrix for the matrix
/// A, and a preconditioner from Ifpack2, Amesos2, or MueLu.
///
/// \section Belos_TpetraAdapter_dev Note to Belos developers
///
/// This partial specialization assumes that the first (Scalar)
/// template parameter of Belos::MultiVecTraits and
/// Belos::OperatorTraits matches the first template parameters of
/// Tpetra::MultiVector and Tpetra::Operator.  In terms of Belos
/// solvers, this means that the specialization assumes that the
/// result of an inner product has the same type as any entry of the
/// multivector or matrix.  This is true for most Scalar types of
/// interest, but may not necessarily be true for certain Scalar types
/// implemented in the Stokhos package, or when implementing
/// mixed-precision solvers in certain ways.  If you don't know what
/// this means, don't worry about it.  If you <i>do</i> know what this
/// means, you might need to write your own partial specialization of
/// Belos::MultiVecTraits and Belos::OperatorTraits, for a Scalar type
/// different than that of the Tpetra::MultiVector or
/// Tpetra::Operator.

// #include "BelosMultiVecTraits_Tpetra.hpp"
// #include "BelosOperatorTraits_Tpetra.hpp"
#include "BelosSolverFactory_Tpetra.hpp"

#endif // BELOS_TPETRA_ADAPTER_HPP
