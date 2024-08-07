// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_NOXSOLVER_DEF_HPP
#define PIRO_NOXSOLVER_DEF_HPP

#include "Piro_NOXSolver.hpp"
#include "Piro_MatrixFreeDecorator.hpp"

#include "Thyra_ModelEvaluatorHelpers.hpp"

#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultInverseLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include <stdexcept>
#include <cstddef>
#include <ostream>

template <typename Scalar>
Piro::NOXSolver<Scalar>::
NOXSolver(const Teuchos::RCP<Teuchos::ParameterList> &appParams_,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &in_model_,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &in_adjointModel_,
    const Teuchos::RCP<ObserverBase<Scalar> > &observer_) :
    SteadyStateSolver<Scalar>(in_model_, in_adjointModel_),
    appParams(appParams_),
    observer(observer_),
    solver(new Thyra::NOXNonlinearSolver),
    out(Teuchos::VerboseObjectBase::getDefaultOStream()),
    model(in_model_),
    adjointModel(in_adjointModel_),
    writeOnlyConvergedSol(appParams_->get("Write Only Converged Solution", true)),
    solveState(true),
    current_iteration(-1)
    {
  using Teuchos::RCP;

  std::string sensitivity_method_str = appParams->get("Sensitivity Method", "Forward");
  this->setSensitivityMethod(sensitivity_method_str);

  const RCP<Teuchos::ParameterList> noxParams =
      Teuchos::sublist(appParams, "NOX", /*mustAlreadyExist =*/ false);

  solver->setParameterList(noxParams);

  std::string jacobianSource = appParams->get("Jacobian Operator", "Have Jacobian");

  exitUponFailedNOXSolve = appParams->get("Exit on Failed NOX Solve", false); 

  reComputeWithZeroInitialGuess = appParams->get("On Failure Solve With Zero Initial Guess", false);

  if (jacobianSource == "Matrix-Free") {
    if (appParams->isParameter("Matrix-Free Perturbation")) {
      model = Teuchos::rcp(new Piro::MatrixFreeDecorator<Scalar>(in_model_,
          appParams->get<double>("Matrix-Free Perturbation")));
    }
    else 
      model = Teuchos::rcp(new Piro::MatrixFreeDecorator<Scalar>(in_model_));
  }
  solver->setModel(model);
    }

template <typename Scalar>
void Piro::NOXSolver<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
    {
  using Teuchos::RCP;
  bool observeFinalSolution = true;
  const int num_p = this->num_p();

  //For Analysis problems we typically do not want to write the solution at each NOX solver.
  //Instead we write at every "write interval" iterations of optimization solver.
  //If write_interval == -1, we print after every successful NOX solver
  //If write_inteval == 0 we never print.
  //This relies on the fact that sensitivities are always called by ROL at each iteration to asses whether the solver is converged
  //TODO: when write_interval>1, at the moment there is no guarantee that the final iteration of the optimization (i.e. the converged solution) gets printed

  if(appParams->isSublist("Optimization Status")){
    auto optimizationParams = appParams->sublist("Optimization Status");
    if(optimizationParams.isParameter("Optimizer Iteration Number"))
      observeFinalSolution = false;

    solveState = optimizationParams.isParameter("Compute State") ? optimizationParams.template get<bool>("Compute State") : true;
  }

  // Forward all parameters to underlying model
  Thyra::ModelEvaluatorBase::InArgs<Scalar> modelInArgs = this->getModel().createInArgs();
  for (int l = 0; l < num_p; ++l) {
    if (Teuchos::nonnull(inArgs.get_p(l)))
      modelInArgs.set_p(l, inArgs.get_p(l));
    else
      modelInArgs.set_p(l, this->getModel().getNominalValues().get_p(l));

    modelInArgs.set_p_direction(l, inArgs.get_p_direction(l));
  }

  // Find the solution of the implicit underlying model
  Thyra::SolveStatus<Scalar> solve_status;
  const Thyra::SolveCriteria<Scalar> solve_criteria;

  if(solveState)
  {
    const auto timer = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Piro::NOXSolver::evalModelImpl::solve")));
    solver->setBasePoint(modelInArgs);

    const RCP<const Thyra::VectorBase<Scalar> > modelNominalState =
        this->getModel().getNominalValues().get_x();


    RCP<Thyra::VectorBase<Scalar> > initial_guess;
    
    //use nominal initial guess only if inArgs.get_x() is not defined
    bool nominalInitialGuess = Teuchos::is_null(inArgs.get_x());
    if(nominalInitialGuess) {
      initial_guess = modelNominalState->clone_v();
    } else {
      initial_guess = inArgs.get_x()->clone_v();
    }

    solve_status = solver->solve(initial_guess.get(), &solve_criteria, /*delta =*/ NULL);
    if((solve_status.solveStatus != ::Thyra::SOLVE_STATUS_CONVERGED) && reComputeWithZeroInitialGuess) {
      auto norm_prev_init_guess = nominalInitialGuess ? initial_guess->norm_1() : inArgs.get_x()->norm_1();
      if(norm_prev_init_guess != 0.0) { //recompute only if the previous initial guess was not already zero
        *out << "Piro::NOXSolver, Trying to solve the nonlinear problem again with a zero initial guess" << std::endl;
        initial_guess->assign(0.0);
        solve_status = solver->solve(initial_guess.get(), &solve_criteria, /*delta =*/ NULL);
      }
    }

    //  MPerego: I think it is better not to throw an error when the solver does not converge.
    //  One can look at the solver status to check whether the solution is converged.
    if (exitUponFailedNOXSolve == true) {
      TEUCHOS_TEST_FOR_EXCEPTION(
          solve_status.solveStatus != ::Thyra::SOLVE_STATUS_CONVERGED,
          std::runtime_error,
          "Nonlinear solver failed to converge");
    }

    if(appParams->isSublist("Optimization Status")){
      appParams->sublist("Optimization Status").set("State Solve Converged", solve_status.solveStatus==Thyra::SOLVE_STATUS_CONVERGED);
      appParams->sublist("Optimization Status").set("Compute State", false);
    }

    auto final_point = model->createInArgs();
    if (final_point.supports(Thyra::ModelEvaluatorBase::IN_ARG_x)) {
      final_point.set_x(solver->get_current_x());
    }
    model->reportFinalPoint(final_point,solve_status.solveStatus==Thyra::SOLVE_STATUS_CONVERGED);

    observeFinalSolution = observeFinalSolution && (solve_status.solveStatus == ::Thyra::SOLVE_STATUS_CONVERGED || !writeOnlyConvergedSol);
  }

  // Retrieve final solution to evaluate underlying model
  const RCP<const Thyra::VectorBase<Scalar> > finalSolution = solver->get_current_x();
  modelInArgs.set_x(finalSolution);

  this->evalConvergedModelResponsesAndSensitivities(modelInArgs, outArgs, *appParams);
  
  bool computeReducedHessian = false;
#ifdef Thyra_BUILD_HESSIAN_SUPPORT
  for (int g_index=0; g_index<this->num_g(); ++g_index) {
    for (int p_index=0; p_index<this->num_p(); ++p_index)
      if (outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_pp, g_index, p_index, p_index))
        if(Teuchos::nonnull(outArgs.get_hess_vec_prod_g_pp(g_index, p_index, p_index))) {
          computeReducedHessian = true;
          break;
        }
  }
#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT

  if(computeReducedHessian == true)   
    this->evalReducedHessian(modelInArgs, outArgs, *appParams);

  if (Teuchos::nonnull(this->observer) && observeFinalSolution) {
    this->observer->observeSolution(*finalSolution);
  }

    }

#endif /*PIRO_NOXSOLVER_DEF_HPP*/
