// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_DETAILS_LINEARSOLVERFACTORY_HPP
#define BELOS_DETAILS_LINEARSOLVERFACTORY_HPP

/// \file Belos_Details_LinearSolverFactory.hpp
/// \brief Implementation of Trilinos::Details::LinearSolverFactory.

#include "BelosSolverFactory.hpp"
#include "Belos_Details_LinearSolver.hpp"
#include "Trilinos_Details_LinearSolverFactory.hpp"

namespace Belos {
namespace Details {

/// \class LinearSolver
/// \brief Belos' implementation of Trilinos::Details::LinearSolverFactory.
template<class MV, class OP, class ScalarType, class NormType>
class LinearSolverFactory :
    public Trilinos::Details::LinearSolverFactory<MV, OP, NormType>
{
public:
  /// \brief Get an instance of a Belos solver.
  ///
  /// The solver is wrapped in a Trilinos::Details::LinearSolver
  /// interface.
  ///
  /// \param solverName [in] The solver's name.  Names are case
  ///   sensitive.
  /// \return A pointer to the solver, if the name was valid; else,
  ///   a null pointer (Teuchos::null).
  virtual Teuchos::RCP<Trilinos::Details::LinearSolver<MV, OP, NormType> >
  getLinearSolver (const std::string& solverName)
  {
    using Teuchos::rcp;
    return rcp (new Belos::Details::LinearSolver<MV, OP, ScalarType, NormType> (solverName));
  }

  /// \brief Register this LinearSolverFactory with the central registry.
  ///
  /// Register this LinearSolverFactory with the central registry, for
  /// the given SC, LO, GO, NT template parameters.  This will let any
  /// clients of Trilinos::Details::getLinearSolver create Belos
  /// solvers with those template parameters.
  ///
  /// You may call this function multiple times; it will only have an
  /// effect the first time (it is idempotent).
  ///
  /// Users do not normally have to call this function.  Belos
  /// automatically registers its LinearSolverFactory with the central
  /// repository, for all enabled template parameter combinations.
  static void registerLinearSolverFactory ()
  {
    typedef Belos::Details::LinearSolverFactory<MV, OP, ScalarType, NormType> this_type;

#ifdef HAVE_TEUCHOSCORE_CXX11
    typedef std::shared_ptr<this_type> ptr_type;
    //typedef std::shared_ptr<Trilinos::Details::LinearSolverFactory<MV, OP> > base_ptr_type;
#else
    typedef Teuchos::RCP<this_type> ptr_type;
    //typedef Teuchos::RCP<Trilinos::Details::LinearSolverFactory<MV, OP> > base_ptr_type;
#endif // HAVE_TEUCHOSCORE_CXX11

    ptr_type factory (new this_type ());
    Trilinos::Details::registerLinearSolverFactory<MV, OP, NormType> ("Belos", factory);
  }
};

} // namespace Details
} // namespace Belos

#endif /* BELOS_DETAILS_LINEARSOLVERFACTORY_HPP */
