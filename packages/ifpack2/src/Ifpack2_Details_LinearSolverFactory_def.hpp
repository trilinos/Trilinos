// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file   Ifpack2_Details_LinearSolverFactory_def.hpp
/// \author Mark Hoemmen
/// \brief  Definition of Ifpack2::Details::LinearSolverFactory.

#ifndef IFPACK2_DETAILS_LINEARSOLVERFACTORY_DEF_HPP
#define IFPACK2_DETAILS_LINEARSOLVERFACTORY_DEF_HPP

#include "Trilinos_Details_LinearSolverFactory.hpp"
#include "Ifpack2_Details_CanChangeMatrix.hpp"
#include "Ifpack2_Details_LinearSolver.hpp"
#include "Ifpack2_Factory.hpp"
#include "Tpetra_RowMatrix.hpp"
#include <type_traits> // std::is_same

namespace Ifpack2 {
namespace Details {

template<class SC, class LO, class GO, class NT>
Teuchos::RCP<typename LinearSolverFactory<SC, LO, GO, NT>::solver_type>
LinearSolverFactory<SC, LO, GO, NT>::
getLinearSolver (const std::string& solverName)
{
  using Teuchos::null;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::TypeNameTraits;
  typedef Ifpack2::Preconditioner<SC, LO, GO, NT> prec_type;
  typedef Tpetra::RowMatrix<SC, LO, GO, NT> ROW;
  const char prefix[] = "Ifpack2::Details::LinearSolverFactory::getLinearSolver: ";

  RCP<prec_type> solver;
  try {
    // The solver to create must be a subclass of
    // Ifpack2::Details::CanChangeMatrix (see documentation of
    // Ifpack2::Details::LinearSolver).  As a result, it should be
    // possible to create the solver with a null matrix.
    solver = Ifpack2::Factory::template create<ROW> (solverName, null);
  }
  catch (std::exception& e) {
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::invalid_argument, prefix << "Failed to create Ifpack2 "
       "preconditioner named \"" << solverName << "\", for the following "
       "template parameters: "
       << "SC = " << TypeNameTraits<SC>::name ()
       << ", LO = " << TypeNameTraits<LO>::name ()
       << ", GO = " << TypeNameTraits<GO>::name ()
       << ", NT = " << TypeNameTraits<NT>::name ()
       << ".  Ifpack2::Factory::create threw an exception: " << e.what ());
  }
  TEUCHOS_TEST_FOR_EXCEPTION
    (solver.is_null (), std::invalid_argument, prefix << "Failed to create "
     "Ifpack2 preconditioner named \"" << solverName << "\", for the "
     "following template parameters: "
     << "SC = " << TypeNameTraits<SC>::name ()
     << ", LO = " << TypeNameTraits<LO>::name ()
     << ", GO = " << TypeNameTraits<GO>::name ()
     << ", NT = " << TypeNameTraits<NT>::name ()
     << ".  Ifpack2::Factory::create returned null.");

  typedef Ifpack2::Details::LinearSolver<SC, LO, GO, NT> impl_type;
  return Teuchos::rcp (new impl_type (solver, solverName));
}

template<class SC, class LO, class GO, class NT>
void
LinearSolverFactory<SC, LO, GO, NT>::
registerLinearSolverFactory ()
{
  typedef Tpetra::MultiVector<SC, LO, GO, NT> MV;
  typedef Tpetra::Operator<SC, LO, GO, NT> OP;
  typedef typename MV::mag_type mag_type;
  typedef Trilinos::Details::LinearSolverFactory<MV, OP, mag_type> factory_base_type;
  typedef Ifpack2::Details::LinearSolverFactory<SC, LO, GO, NT> factory_impl_type;

#ifdef HAVE_TEUCHOSCORE_CXX11
  typedef std::shared_ptr<factory_base_type> base_ptr_type;
  typedef std::shared_ptr<factory_impl_type> impl_ptr_type;
#else
  typedef Teuchos::RCP<factory_base_type> base_ptr_type;
  typedef Teuchos::RCP<factory_impl_type> impl_ptr_type;
#endif // HAVE_TEUCHOSCORE_CXX11

  impl_ptr_type factory (new factory_impl_type ());
  base_ptr_type factoryBase = factory; // implicit cast to base class

  TEUCHOS_TEST_FOR_EXCEPTION
    (factoryBase.get () == NULL, std::logic_error, "Factory is null!  This "
     "should never happen!  Please report this bug to the Ifpack2 developers.");

// #ifdef HAVE_IFPACK2_DEBUG
//   {
//     using std::cerr;
//     using std::endl;
//     using Teuchos::TypeNameTraits;
//     cerr << "Registering Ifpack2 LinearSolverFactory for"
//          << " SC = " << TypeNameTraits<SC>::name ()
//          << ", LO = " << TypeNameTraits<LO>::name ()
//          << ", GO = " << TypeNameTraits<GO>::name ()
//          << ", NT = " << TypeNameTraits<NT>::name ()
//          << ", and mag_type = " << TypeNameTraits<mag_type>::name ()
//          << endl;
//   }
// #endif // HAVE_IFPACK2_DEBUG
  Trilinos::Details::registerLinearSolverFactory<MV, OP, mag_type> ("Ifpack2", factoryBase);
}

} // namespace Details
} // namespace Ifpack2

// Do explicit instantiation of Ifpack2::Details::LinearSolverFactory,
// for Tpetra objects, with the given Tpetra template parameters.
#define IFPACK2_DETAILS_LINEARSOLVERFACTORY_INSTANT( SC, LO, GO, NT ) \
  template class Ifpack2::Details::LinearSolverFactory<SC, LO, GO, NT>;

#endif // IFPACK2_DETAILS_LINEARSOLVERFACTORY_DEF_HPP
