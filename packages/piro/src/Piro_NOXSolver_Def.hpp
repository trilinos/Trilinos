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

#include <stdexcept>
#include <cstddef>
#include <ostream>

template <typename Scalar>
Piro::NOXSolver<Scalar>::
NOXSolver(const Teuchos::RCP<Teuchos::ParameterList> &appParams_,
	  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model_,
          const Teuchos::RCP<ObserverBase<Scalar> > &observer_) :
  SteadyStateSolver<Scalar>(model_),
  appParams(appParams_),
  observer(observer_),
  solver(new Thyra::NOXNonlinearSolver),
  out(Teuchos::VerboseObjectBase::getDefaultOStream()),
  model(model_),
  writeOnlyConvergedSol(appParams_->get("Write Only Converged Solution", true))
{
  using Teuchos::RCP;

  const RCP<Teuchos::ParameterList> noxParams =
    Teuchos::sublist(appParams, "NOX", /*mustAlreadyExist =*/ false);

  solver->setParameterList(noxParams);

  std::string jacobianSource = appParams->get("Jacobian Operator", "Have Jacobian");

  if (jacobianSource == "Matrix-Free") {
    if (appParams->isParameter("Matrix-Free Perturbation")) {
      model = Teuchos::rcp(new Piro::MatrixFreeDecorator<Scalar>(model_,
                           appParams->get<double>("Matrix-Free Perturbation")));
    }
    else 
      model = Teuchos::rcp(new Piro::MatrixFreeDecorator<Scalar>(model_));
  }
  solver->setModel(model);
}

template <typename Scalar>
void Piro::NOXSolver<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
{
  using Teuchos::RCP;

  // Forward all parameters to underlying model
  Thyra::ModelEvaluatorBase::InArgs<Scalar> modelInArgs = this->getModel().createInArgs();
  for (int l = 0; l < this->num_p(); ++l) {
    modelInArgs.set_p(l, inArgs.get_p(l));
  }

  // Find the solution of the implicit underlying model
  Thyra::SolveStatus<Scalar> solve_status;
  {
    solver->setBasePoint(modelInArgs);

    const RCP<const Thyra::VectorBase<Scalar> > modelNominalState =
      this->getModel().getNominalValues().get_x();

    const RCP<Thyra::VectorBase<Scalar> > initial_guess = modelNominalState->clone_v();

    const Thyra::SolveCriteria<Scalar> solve_criteria;
    solve_status = solver->solve(initial_guess.get(), &solve_criteria, /*delta =*/ NULL);

//  MPerego: I think it is better not to throw an error when the solver does not converge. 
//  One can look at the solver status to check whether the solution is converged.
//    TEUCHOS_TEST_FOR_EXCEPTION(
//        solve_status.solveStatus != ::Thyra::SOLVE_STATUS_CONVERGED,
//        std::runtime_error,
//       "Nonlinear solver failed to converge");
  }

  // Retrieve final solution to evaluate underlying model
  const RCP<const Thyra::VectorBase<Scalar> > finalSolution = solver->get_current_x();
  modelInArgs.set_x(finalSolution);

  if (Teuchos::nonnull(this->observer) && (solve_status.solveStatus == ::Thyra::SOLVE_STATUS_CONVERGED || !writeOnlyConvergedSol)) {
    this->observer->observeSolution(*finalSolution);
  }

  this->evalConvergedModel(modelInArgs, outArgs);
}

#endif /*PIRO_NOXSOLVER_DEF_HPP*/
