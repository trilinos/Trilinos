// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include "NOX_Observer_ReusePreconditionerFactory.hpp"
#include "NOX_Observer_ReusePreconditioner.hpp"
#include "Teuchos_ParameterList.hpp"

Teuchos::RCP<NOX::Observer>
NOX::createReusePreconditionerObserver(Teuchos::ParameterList& pl)
{
  Teuchos::ParameterList validParams;
  validParams.set<bool>("Update prec at start of nonlinear solve",true);
  validParams.set<int>("Update prec after this many nonlinear solves",0);
  validParams.set<bool>("Reset nonlinear solve count on failed solve",true);
  validParams.set<int>("Update prec after this many nonlinear iterations",0);
  validParams.set<int>("Update prec after this many stalled linear solves",0);
  validParams.set<int>("Max linear iterations for stall",50);

  pl.validateParametersAndSetDefaults(validParams);

  bool update_at_start = pl.get<bool>("Update prec at start of nonlinear solve");
  int update_after_n_nonlinear_solves = pl.get<int>("Update prec after this many nonlinear solves");
  bool reset_nonlinear_solve_count_on_failed_solve = pl.get<bool>("Reset nonlinear solve count on failed solve");
  int update_after_n_iters = pl.get<int>("Update prec after this many nonlinear iterations");
  int update_after_n_stalls = pl.get<int>("Update prec after this many stalled linear solves");
  int max_iters_for_stall = pl.get<int>("Max linear iterations for stall");

  auto observer = Teuchos::rcp(new NOX::ObserverReusePreconditioner);

  if (update_at_start)
    observer->updateAtStartOfSolve();

  // Zero or negative value disables
  if (update_after_n_nonlinear_solves > 0)
    observer->updateAfterNNonlinearSolves(update_after_n_nonlinear_solves,
                                          reset_nonlinear_solve_count_on_failed_solve);

  // Zero or negative value disables
  if (update_after_n_iters > 0)
    observer->updateAfterNIterations(update_after_n_iters);

  // Zero or negative value disables
  if (update_after_n_stalls > 0)
    observer->updateOnLinearSolverStall(max_iters_for_stall,update_after_n_stalls);

  return observer;
}
