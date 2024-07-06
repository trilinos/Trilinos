// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file   Ifpack2_Details_LinearSolverFactory_decl.hpp
/// \author Mark Hoemmen
/// \brief  Declaration of Ifpack2::Details::LinearSolverFactory.

#ifndef IFPACK2_DETAILS_LINEARSOLVERFACTORY_DECL_HPP
#define IFPACK2_DETAILS_LINEARSOLVERFACTORY_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Trilinos_Details_LinearSolverFactory.hpp"
#include "Tpetra_Operator.hpp"

namespace Ifpack2 {
namespace Details {

/// \class LinearSolverFactory
/// \brief Interface for a "factory" that creates Ifpack2 solvers.
///
/// We use Tpetra's template parameters here, instead of MV, OP, and
/// NormType, because this is not a public-facing class.  We also want
/// to avoid mix-ups between MV and OP.
template<class SC, class LO, class GO, class NT>
class LinearSolverFactory :
  public Trilinos::Details::LinearSolverFactory<Tpetra::MultiVector<SC, LO, GO, NT>,
                                                Tpetra::Operator<SC, LO, GO, NT>,
                                                typename Tpetra::MultiVector<SC, LO, GO, NT>::mag_type>
{
public:
  typedef Trilinos::Details::LinearSolver<Tpetra::MultiVector<SC, LO, GO, NT>,
                                          Tpetra::Operator<SC, LO, GO, NT>,
                                          typename Tpetra::MultiVector<SC, LO, GO, NT>::mag_type> solver_type;

  /// \brief Get an instance of a Ifpack2 solver.
  ///
  /// The solver is wrapped in a Trilinos::Details::LinearSolver
  /// interface.
  ///
  /// \param solverName [in] The solver's name.  Names are case
  ///   sensitive.
  /// \return A pointer to the solver, if the name was valid; else,
  ///   a null pointer (Teuchos::null).
  virtual Teuchos::RCP<solver_type>
  getLinearSolver (const std::string& solverName);

  /// \brief Register this LinearSolverFactory with the central registry.
  ///
  /// Register this LinearSolverFactory with the central registry, for
  /// the given SC, LO, GO, NT template parameters.  This will let any
  /// clients of Trilinos::Details::getLinearSolver create Ifpack2
  /// solvers with those template parameters.
  ///
  /// You may call this function multiple times; it will only have an
  /// effect the first time (it is idempotent).
  ///
  /// Users do not normally have to call this function.  Ifpack2
  /// automatically registers its LinearSolverFactory with the central
  /// repository, for all enabled template parameter combinations.
  static void registerLinearSolverFactory ();
};

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_LINEARSOLVERFACTORY_DECL_HPP
