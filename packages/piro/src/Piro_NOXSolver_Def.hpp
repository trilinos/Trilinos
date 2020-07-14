// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
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
    const Teuchos::RCP<ObserverBase<Scalar> > &observer_) :
    SteadyStateSolver<Scalar>(in_model_),
    appParams(appParams_),
    observer(observer_),
    solver(new Thyra::NOXNonlinearSolver),
    out(Teuchos::VerboseObjectBase::getDefaultOStream()),
    model(in_model_),
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
  const int num_g = this->num_g();

  //For Analysis problems we typically do not want to write the solution at each NOX solver.
  //Instead we write at every "write interval" iterations of optimization solver.
  //If write_interval == -1, we print after every successful NOX solver
  //If write_inteval == 0 we ever print.
  //This relies on the fact that sensitivities are always called by ROL at each iteration to asses whether the solver is converged
  //TODO: when write_interval>1, at the moment there is no guarantee that the final iteration of the optimization (i.e. the converged solution) gets printed

  //When write_interval>0 we print only after computing the sensitivities,
  //to make sure to print them updated if required.
  bool solving_sensitivities = false;
  for (int i=0; i<num_p; i++) {
    for (int j=0; j<=num_g; j++) {
      if (!outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, i).none() && !outArgs.get_DgDp(j,i).isEmpty()) {
        solving_sensitivities = true;
        break;
      }
    }
  }
  if(appParams->isSublist("Analysis")){
    auto analysisParams = appParams->sublist("Analysis");
    if(analysisParams.isSublist("Optimization Status")){
      auto optimizationParams = analysisParams.sublist("Optimization Status");
      if(optimizationParams.isParameter("Optimizer Iteration Number")) {
        int iteration = optimizationParams.template get<int>("Optimizer Iteration Number");
        int write_interval = analysisParams.get("Write Interval",1);
        Teuchos::RCP<Teuchos::ParameterList> writeParams = Teuchos::rcp(new Teuchos::ParameterList());
        if(analysisParams.isSublist("ROL") && analysisParams.sublist("ROL").get("Use Old Reduced Space Interface", false)) {
          if (write_interval > 0) {
            observeFinalSolution = (iteration >= 0 && solving_sensitivities && iteration != current_iteration) ? (iteration%write_interval == 0) : false;
            if(observeFinalSolution)
              current_iteration = iteration;
          }
          else if (write_interval == 0)
            observeFinalSolution = false;
          else if (write_interval == -1)
            observeFinalSolution = true;
        }
        else
          observeFinalSolution = false;
      }

      // Call observer method for each parameter if optimization variables have changed
      if (optimizationParams.isParameter("Optimization Variables Changed")) {
        if (optimizationParams.template get<bool>("Optimization Variables Changed")) {
          if (analysisParams.isSublist("ROL")) {
            auto rolParams = analysisParams.sublist("ROL");
            if(rolParams.get("Use Old Reduced Space Interface", false)) {
              int num_parameters = rolParams.get("Number of Parameters", 1);
              for(int i=0; i<num_parameters; ++i) {
                std::ostringstream ss; ss << "Parameter Vector Index " << i;
                const int p_ind = rolParams.get(ss.str(), i);
                const auto paramNames = *this->getModel().get_p_names(p_ind);
                if(Teuchos::nonnull(this->observer))
                  this->observer->parameterChanged(paramNames[0]);
              }
              optimizationParams.set("Optimization Variables Changed", false);
            }
          }
        }
      }
      solveState = optimizationParams.isParameter("Compute State") ? optimizationParams.template get<bool>("Compute State") : true;
    }
  }

  // Forward all parameters to underlying model
  Thyra::ModelEvaluatorBase::InArgs<Scalar> modelInArgs = this->getModel().createInArgs();
  for (int l = 0; l < this->num_p(); ++l) {
    modelInArgs.set_p(l, inArgs.get_p(l));
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
    if(Teuchos::nonnull(inArgs.get_x())) { //used in optimization
      initial_guess = inArgs.get_x()->clone_v();
    } else {
      initial_guess = modelNominalState->clone_v();
    }

    solve_status = solver->solve(initial_guess.get(), &solve_criteria, /*delta =*/ NULL);
    if((solve_status.solveStatus != ::Thyra::SOLVE_STATUS_CONVERGED) && reComputeWithZeroInitialGuess && (initial_guess->norm_1()!=0.0)) {
      *out << "Piro::NOXSolver, Trying to solve the nonlinear problem again with a zero initial guess" << std::endl;
      initial_guess->assign(0.0);
      solve_status = solver->solve(initial_guess.get(), &solve_criteria, /*delta =*/ NULL);
    }


    //  MPerego: I think it is better not to throw an error when the solver does not converge.
    //  One can look at the solver status to check whether the solution is converged.
    if (exitUponFailedNOXSolve == true) {
      TEUCHOS_TEST_FOR_EXCEPTION(
          solve_status.solveStatus != ::Thyra::SOLVE_STATUS_CONVERGED,
          std::runtime_error,
          "Nonlinear solver failed to converge");
    }

    if(appParams->isSublist("Analysis")){
      Teuchos::ParameterList& analysisParams = appParams->sublist("Analysis");
      if(analysisParams.isSublist("Optimization Status")) {
        analysisParams.sublist("Optimization Status").set("State Solve Converged", solve_status.solveStatus==Thyra::SOLVE_STATUS_CONVERGED);
        analysisParams.sublist("Optimization Status").set("Compute State", false);
      }
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

  this->evalConvergedModelResponsesAndSensitivities(modelInArgs, outArgs);

  if (Teuchos::nonnull(this->observer) && observeFinalSolution) {
    this->observer->observeSolution(*finalSolution);
  }

    }

#endif /*PIRO_NOXSOLVER_DEF_HPP*/
