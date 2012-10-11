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

#include "Piro_NOXSolver.hpp"

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TestForException.hpp"

#include <stdexcept>
#include <cstddef>
#include <ostream>

template <typename Scalar>
Piro::NOXSolver<Scalar>::
NOXSolver(Teuchos::RCP<Teuchos::ParameterList> appParams_,
	  Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<Scalar> > model_) :
  appParams(appParams_),
  model(model_),
  num_p(model->Np()),
  num_g(model->Ng()),
  solver(new Thyra::NOXNonlinearSolver),
  out(Teuchos::VerboseObjectBase::getDefaultOStream())
{
  using Teuchos::RCP;

  const RCP<Teuchos::ParameterList> noxParams =
    Teuchos::sublist(appParams, "NOX", /*mustAlreadyExist =*/ false);

  solver->setParameterList(noxParams);
  solver->setModel(model);
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::NOXSolver<Scalar>::get_p_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l >= num_p || l < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in Piro::NOXSolver::get_p_map():  " <<
                     "Invalid parameter index l = " <<
                     l << std::endl);
  return model->get_p_space(l);
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::NOXSolver<Scalar>::get_g_space(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(j > num_g || j < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in Piro::NOXSolver::get_g_map():  " <<
                     "Invalid response index j = " <<
                     j << std::endl);

  if (j < num_g) {
    return model->get_g_space(j);
  } else {
    // j == num_g, corresponding to the state by convention
    return model->get_x_space();
  }
}

template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::NOXSolver<Scalar>::getNominalValues() const
{
  return model->getNominalValues();
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> Piro::NOXSolver<Scalar>::createInArgs() const
{
  Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p);
  return inArgs;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar> Piro::NOXSolver<Scalar>::createOutArgsImpl() const
{
  Thyra::ModelEvaluatorBase::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());

  // One additional response slot for the solution vector
  outArgs.set_Np_Ng(num_p, num_g+1);

  Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs = model->createOutArgs();
  for (int i=0; i<num_g; i++)
    for (int j=0; j<num_p; j++)
      if (!modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, i, j).none())
        outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, i, j,
             Thyra::ModelEvaluatorBase::DerivativeSupport(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM));

  return outArgs;
}

template <typename Scalar>
void Piro::NOXSolver<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
{
  using Teuchos::RCP;

  // Solve the underlying model
  // TODO: pass the parameters from inArgs down to the underlying model used by the solver.
  // Currently, only the default values of p set in createInArgs are used in the nonlinear solver !
  {
    const RCP<Thyra::VectorBase<double> > initial_guess =
      model->getNominalValues().get_x()->clone_v();

    const Thyra::SolveCriteria<double> solve_criteria;
    const Thyra::SolveStatus<double> solve_status =
      solver->solve(initial_guess.get(), &solve_criteria, /*delta =*/ NULL);

    TEUCHOS_TEST_FOR_EXCEPTION(
        solve_status.solveStatus != ::Thyra::SOLVE_STATUS_CONVERGED,
        std::runtime_error,
        "Nonlinear solver failed to converge");
  }

  // Retrieve the solution at convergence to perform underlying model evaluation
  const RCP<const Thyra::VectorBase<Scalar> > finalSolution = solver->get_current_x();
  // Solution at convergence is the response at index num_g
  const RCP<Thyra::VectorBase<Scalar> > gx_out = outArgs.get_g(num_g);
  if (Teuchos::nonnull(gx_out)) {
    Thyra::copy(*finalSolution, gx_out.ptr());
  }

  // Setup input for underlying model
  Thyra::ModelEvaluatorBase::InArgs<Scalar> modelInArgs = model->createInArgs();
  {
    // Evaluation to be done at final solution
    modelInArgs.set_x(finalSolution);

    // Forward all parameters
    for (int l = 0; l < num_p; ++l) {
      modelInArgs.set_p(l, inArgs.get_p(l));
    }
  }

  // Setup output for underlying model
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs = model->createOutArgs();
  {
    // Forward all responses
    for (int j = 0; j < num_g; ++j) {
      const RCP<Thyra::VectorBase<Scalar> > g_out = outArgs.get_g(j);
      modelOutArgs.set_g(j, g_out);
    }

    // Request derivatives required to compute forward sensitivities
    for (int l = 0; l < num_p; ++l) {
      for (int j = 0; j < num_g; ++j) {
	if (!outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l).none()) {
	  const Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdp_deriv =
	    outArgs.get_DgDp(j, l);
	  const RCP<Thyra::MultiVectorBase<Scalar> > dgdp_mv =
	    dgdp_deriv.getMultiVector();
	  
	  if (Teuchos::nonnull(dgdp_mv)) {
	    const Thyra::ModelEvaluatorBase::DerivativeSupport dpdg_support =
	      modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l);
	    TEUCHOS_TEST_FOR_EXCEPT(!dpdg_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM));
	    
	    modelOutArgs.set_DgDp(j, l, dgdp_deriv);
	    
	    const RCP<Thyra::MultiVectorBase<Scalar> > dfdp_mv =
	      Thyra::createMembers(model->get_f_space(), model->get_p_space(l));
	    modelOutArgs.set_DfDp(l,
				  Thyra::ModelEvaluatorBase::Derivative<Scalar>(dfdp_mv, Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM));
	    
	    const RCP<Thyra::MultiVectorBase<Scalar> > dgdx_mv =
	      Thyra::createMembers(model->get_x_space(), model->get_g_space(j));
	    modelOutArgs.set_DgDx(j,
				  Thyra::ModelEvaluatorBase::Derivative<Scalar>(dgdx_mv, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM));
	    
	    const RCP<Thyra::LinearOpWithSolveBase<Scalar> > jacobian =
	      model->create_W();
	    modelOutArgs.set_W(jacobian);
	  }
        }
      }
    }
  }

  model->evalModel(modelInArgs, modelOutArgs);

  // Compute forward sensitivities
  for (int l = 0; l < num_p; ++l) {
    for (int j = 0; j < num_g; ++j) {
      if (!outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l).none()) {
	const Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdp_deriv =
	  outArgs.get_DgDp(j, l);
	const RCP<Thyra::MultiVectorBase<Scalar> > dgdp_mv =
	  dgdp_deriv.getMultiVector();
	
	if (Teuchos::nonnull(dgdp_mv)) {
	  const RCP<Thyra::MultiVectorBase<Scalar> > dxdp_mv =
	    Thyra::createMembers(model->get_x_space(), model->get_p_space(l));
	  assign(dxdp_mv.ptr(), Teuchos::ScalarTraits<Scalar>::zero());
	  
	  const RCP<Thyra::MultiVectorBase<Scalar> > dfdp_mv =
	    modelOutArgs.get_DfDp(l).getMultiVector();
	  const RCP<Thyra::LinearOpWithSolveBase<Scalar> > jacobian =
	    modelOutArgs.get_W();
	  const RCP<Thyra::MultiVectorBase<Scalar> > dgdx_mv =
	    modelOutArgs.get_DgDx(j).getMultiVector();
	  
	  const Thyra::SolveCriteria<Scalar> defaultSolveCriteria;
	  const Thyra::SolveStatus<Scalar> solveStatus =
	    Thyra::solve(*jacobian, Thyra::NOTRANS, *dfdp_mv, dxdp_mv.ptr(), Teuchos::ptr(&defaultSolveCriteria));
	  TEUCHOS_TEST_FOR_EXCEPTION(
				     solveStatus.solveStatus == Thyra::SOLVE_STATUS_UNCONVERGED,
				     std::runtime_error,
				     "Jacobian solver failed to converge");
	  
	  Thyra::apply(
		       *dgdx_mv,
		       Thyra::TRANS,
		       *dxdp_mv,
		       dgdp_mv.ptr(),
		       -Teuchos::ScalarTraits<Scalar>::one(),
		       Teuchos::ScalarTraits<Scalar>::one());
	}
      }
    }
  }
}
