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

// TODO: Remove with debugging statements
#include <iostream>

template <typename Scalar>
Piro::NOXSolver<Scalar>::
NOXSolver(Teuchos::RCP<Teuchos::ParameterList> appParams_,
	  Teuchos::RCP< Thyra::ModelEvaluatorDefaultBase<Scalar> > model_) :
  appParams(appParams_),
  model(model_)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out = Teuchos::VerboseObjectBase::getDefaultOStream();

  RCP<Teuchos::ParameterList> noxParams =
	rcp(&(appParams->sublist("NOX")),false);

//   string jacobianSource = appParams->get("Jacobian Operator", "Have Jacobian");

//   if (jacobianSource == "Matrix-Free") {
//     TEUCHOS_TEST_FOR_EXCEPTION(jacobianSource == "Matrix-Free", std::logic_error,
//        "MATRIX_free not yet implemented for Piro Thyra");
//     model = rcp(new Piro::Thyra::MatrixFreeDecorator(model));
//   }

  // Grab some modelEval stuff from underlying model
  num_p = model->createInArgs().Np();
  num_g = model->createOutArgs().Ng();

  solver = rcp(new Thyra::NOXNonlinearSolver);
  solver->setParameterList(noxParams);
  solver->setModel(model);
}

template <typename Scalar>
Piro::NOXSolver<Scalar>::~NOXSolver()
{
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

  if      (j < num_g) return model->get_g_space(j);
  else return model->get_x_space(); // j == num_g
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

  // Ng is 1 bigger then model-Ng so that the solution vector can be an outarg
  outArgs.set_Np_Ng(num_p, num_g+1);

  Thyra::ModelEvaluatorBase::OutArgs<Scalar> model_outargs = model->createOutArgs();
  for (int i=0; i<num_g; i++)
    for (int j=0; j<num_p; j++)
      if (!model_outargs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, i, j).none())
        outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, i, j,
             Thyra::ModelEvaluatorBase::DerivativeSupport(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM));

  return outArgs;
}

template <typename Scalar>
void Piro::NOXSolver<Scalar>::evalModelImpl(
       const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
       const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs ) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  TEUCHOS_TEST_FOR_EXCEPTION(
      num_p > 1,
      std::logic_error,
      "More than 1 parameter not yet supported");

  TEUCHOS_TEST_FOR_EXCEPTION(
      num_g > 1,
      std::logic_error,
      "More than 1 response not yet supported");

  // Parse InArgs
  // TODO: Handle case num_p > 1
  RCP<const Thyra::VectorBase<Scalar> > p_in;
  if (num_p > 0) {
    p_in = inArgs.get_p(0);
  }

  // Parse OutArgs
  RCP< Thyra::VectorBase<Scalar> > g_out;
  // TODO: Handle case num_g > 1
  if (num_g > 0) {
    g_out = outArgs.get_g(0);
  }
  // By convention, response at index num_g is the solution vector
  RCP< Thyra::VectorBase<Scalar> > gx_out = outArgs.get_g(num_g);

  // TODO: Parse out-args for sensitivity calculation

  // Solve the underlying model
  {
    RCP< ::Thyra::VectorBase<double> >
      initial_guess = model->getNominalValues().get_x()->clone_v();

    const ::Thyra::SolveCriteria<double> solve_criteria;
    const ::Thyra::SolveStatus<double> solve_status =
      solver->solve(initial_guess.get(), &solve_criteria, NULL);

    TEUCHOS_TEST_FOR_EXCEPTION(
        solve_status.solveStatus != ::Thyra::SOLVE_STATUS_CONVERGED,
        std::runtime_error,
        "Nonlinear solver failed to converge");
  }

  // Return the final solution as the last response
  const RCP<const Thyra::VectorBase<Scalar> > finalSolution = solver->get_current_x();
  if (Teuchos::nonnull(gx_out)) {
    Thyra::copy(*finalSolution, gx_out.ptr());
  }

  // Compute responses at final solution
  // TODO: Handle case num_g > 1
  if (Teuchos::nonnull(g_out)) {
     // Setup input
     Thyra::ModelEvaluatorBase::InArgs<Scalar> model_inargs = model->createInArgs();
     {
       model_inargs.set_x(finalSolution);
       if (num_p > 0) {
         model_inargs.set_p(0, p_in);
       }
     }

     // Setup output
     Thyra::ModelEvaluatorBase::OutArgs<Scalar> model_outargs = model->createOutArgs();
     {
       model_outargs.set_g(0, g_out);
     }

     model->evalModel(model_inargs, model_outargs);
  }

  // TODO: Handle case num_g > 1 or num_p > 1
  if (num_p > 0 && num_g > 0) {
    const Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdp_deriv = outArgs.get_DgDp(0, 0);
    const RCP<Thyra::MultiVectorBase<Scalar> > dgdp_mv = dgdp_deriv.getMultiVector();

    if (Teuchos::nonnull(dgdp_mv)) {
      Thyra::ModelEvaluatorBase::InArgs<Scalar> model_inargs = model->createInArgs();
      {
        model_inargs.set_x(finalSolution);
        if (num_p > 0) {
          model_inargs.set_p(0, p_in);
        }
      }

      Thyra::ModelEvaluatorBase::OutArgs<Scalar> model_outargs = model->createOutArgs();
      const Thyra::ModelEvaluatorBase::DerivativeSupport dpdg_support =
        model_outargs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, 0, 0);
      TEUCHOS_TEST_FOR_EXCEPT(!dpdg_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM));

      model_outargs.set_DgDp(0, 0, dgdp_deriv);

      const RCP<Thyra::MultiVectorBase<Scalar> > dfdp_mv =
        Thyra::createMembers(model->get_f_space(), model->get_p_space(0));
      model_outargs.set_DfDp(0,
          Thyra::ModelEvaluatorBase::Derivative<Scalar>(dfdp_mv, Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM));

      const RCP<Thyra::MultiVectorBase<Scalar> > dgdx_mv =
        Thyra::createMembers(model->get_x_space(), model->get_g_space(0));
      model_outargs.set_DgDx(0,
          Thyra::ModelEvaluatorBase::Derivative<Scalar>(dgdx_mv, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM));

      const RCP<Thyra::LinearOpWithSolveBase<Scalar> > jacobian =
        model->create_W();
      model_outargs.set_W(jacobian);

      model->evalModel(model_inargs, model_outargs);

      const RCP<Thyra::MultiVectorBase<Scalar> > dxdp_mv =
        Thyra::createMembers(model->get_x_space(), model->get_p_space(0));
      assign(dxdp_mv.ptr(), Teuchos::ScalarTraits<Scalar>::zero());

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
