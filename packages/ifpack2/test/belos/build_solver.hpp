// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _build_solver_hpp_
#define _build_solver_hpp_

#include "Teuchos_RefCountPtr.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosSolverFactory.hpp"

template<class Scalar,class MV, class OP>
Teuchos::RCP<Belos::SolverManager<Scalar,MV,OP> >
build_solver(Teuchos::ParameterList& test_params,
             Teuchos::RCP<Belos::LinearProblem<Scalar,MV,OP> > problem)
{
  Teuchos::ParameterList bparams;
  if (test_params.isSublist("Belos")) {
    bparams = test_params.sublist("Belos");
  }
  Teuchos::RCP<Teuchos::ParameterList> rcpparams = Teuchos::rcp(&bparams,false);

  std::string solver_type("not specified");

  Ifpack2::getParameter(test_params, "solver_type", solver_type);

  if (solver_type == "not specified") {
    throw std::runtime_error("Error in build_solver: solver_type not specified.");
  }

  Teuchos::RCP<Belos::SolverManager<Scalar,MV,OP> > solver;
  Belos::SolverFactory<Scalar, MV, OP> factory;
  try {
    solver = factory.create (solver_type, rcpparams);
  } catch (std::exception& e) {
    std::cout << "*** FAILED: Belos::SolverFactory::create threw an exception: "
              << e.what () << std::endl;
    return Teuchos::null;
  }

  solver->setProblem( problem );

  return solver;
}
#endif

