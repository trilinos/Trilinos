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

#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultInverseLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Array.hpp"

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
                     "Error in Piro::NOXSolver::get_p_space():  " <<
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
                     "Error in Piro::NOXSolver::get_g_space():  " <<
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
  Thyra::ModelEvaluatorBase::InArgs<Scalar> result = this->createInArgs();
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> modelNominalValues = model->getNominalValues();
  for (int l = 0; l < num_p; ++l) {
    result.set_p(l, modelNominalValues.get_p(l));
  }
  return result;
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
  outArgs.set_Np_Ng(num_p, num_g + 1);

  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs = model->createOutArgs();

  // Sensitivity support
  // 1) Jacobian operator required
  if (modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_W)) {
    for (int j = 0; j < num_g; ++j) {
      // 2) DgDx required
      const Thyra::ModelEvaluatorBase::DerivativeSupport dgdx_support =
        modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, j);
        const bool dgdx_forwardOpSupport =
          dgdx_support.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
        const bool dgdx_mvGradSupport =
          dgdx_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
      if (dgdx_forwardOpSupport || dgdx_mvGradSupport) {
        for (int l = 0; l < num_p; ++l) {
          // 3) DfDp required
          const Thyra::ModelEvaluatorBase::DerivativeSupport dfdp_support =
            modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp, l);
          const bool dfdp_forwardOpSupport =
            dfdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
          const bool dfdp_mvJacSupport =
            dfdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
          if (dfdp_forwardOpSupport || dfdp_mvJacSupport) {
            // 4) Dgdp required
            const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support =
              modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l);
            Thyra::ModelEvaluatorBase::DerivativeSupport sensitivitySupport;
            if (dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP) &&
                dgdx_forwardOpSupport && dfdp_forwardOpSupport) {
              sensitivitySupport.plus(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
            }
            if (dfdp_mvJacSupport) {
              if (dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
                sensitivitySupport.plus(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
              }
              if (dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM) && dgdx_mvGradSupport) {
                sensitivitySupport.plus(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
              }
            }
            outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l, sensitivitySupport);
          }
        }
      }
    }
  }

  return outArgs;
}

template <typename Scalar>
void Piro::NOXSolver<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

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

  // Compute responses
  {
    Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs = model->createOutArgs();
    // Forward to underlying model
    for (int j = 0; j < num_g; ++j) {
      const RCP<Thyra::VectorBase<Scalar> > g_out = outArgs.get_g(j);
      modelOutArgs.set_g(j, g_out);
    }
    model->evalModel(modelInArgs, modelOutArgs);
  }

  // Sensitivities
  for (int j = 0; j < num_g; ++j) {
    for (int l = 0; l < num_p; ++l) {
      const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support =
        outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l);
      if (!dgdp_support.none()) {
        Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs = model->createOutArgs();
        {
          const Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdp_deriv =
            outArgs.get_DgDp(j, l);
          if (!dgdp_deriv.isEmpty()) {
            // Solver sensitivity DpDg(j, l) requested by user
            const RCP<Thyra::MultiVectorBase<Scalar> > dgdp_mv =
              dgdp_deriv.getMultiVector();
            if (Teuchos::nonnull(dgdp_mv)) {
              // Model DgDp(j, l) required in multivector form
              modelOutArgs.set_DgDp(j, l, dgdp_deriv);
              // Model DfDp(l) required in multivector jacobian form
              const RCP<Thyra::MultiVectorBase<Scalar> > dfdp_mv =
                Thyra::createMembers(model->get_f_space(), model->get_p_space(l));
              modelOutArgs.set_DfDp(l,
                  Thyra::ModelEvaluatorBase::Derivative<Scalar>(dfdp_mv, Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM));

              const bool dgdp_mvGrad_required =
                dgdp_deriv.getMultiVectorOrientation() == Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM;
              Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdx_deriv;
              if (dgdp_mvGrad_required) {
                // Model DgDx(j) required in multivector gradient form
                dgdx_deriv =
                  Thyra::create_DgDx_mv(*model, j, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
              } else {
                // Model DgDx(j) required in linear operator form
                const RCP<Thyra::LinearOpBase<Scalar> > dgdx_op = model->create_DgDx_op(j);
                dgdx_deriv = Thyra::ModelEvaluatorBase::Derivative<Scalar>(dgdx_op);
              }
              modelOutArgs.set_DgDx(j, dgdx_deriv);
            } else {
              // Model DgDp(j, l) required in linear operator form
              const RCP<Thyra::DefaultAddedLinearOp<Scalar> > dgdp_op =
                Teuchos::rcp_dynamic_cast<Thyra::DefaultAddedLinearOp<Scalar> >(dgdp_deriv.getLinearOp());
              TEUCHOS_TEST_FOR_EXCEPTION(
                  Teuchos::is_null(dgdp_op),
                  std::invalid_argument,
                  "Illegal operator for DgDp(" <<
                  "j = " << j << ", " <<
                  "index l = " << l << ")\n");

              // Model DgDp(j, l) required in linear operator form
              const RCP<Thyra::LinearOpBase<Scalar> > model_dgdp_op = model->create_DgDp_op(j, l);
              modelOutArgs.set_DgDp(j, l, model_dgdp_op);

              // Model DfDp(l) required in linear operator form
              const RCP<Thyra::LinearOpBase<Scalar> > model_dfdp_op = model->create_DfDp_op(l);
              modelOutArgs.set_DfDp(l, model_dfdp_op);
              // Model DgDx(j) required in linear operator form
              const RCP<Thyra::LinearOpBase<Scalar> > model_dgdx_op = model->create_DgDx_op(j);
              modelOutArgs.set_DgDx(l, model_dgdx_op);
            }

            // Model Jacobian with solve required
            {
              const RCP<Thyra::LinearOpWithSolveBase<Scalar> > jacobian =
                model->create_W();
              modelOutArgs.set_W(jacobian);
            }

            // Call underlying model
            model->evalModel(modelInArgs, modelOutArgs);

            // Assemble solver sensitivity
            {
              const Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdp_deriv =
                outArgs.get_DgDp(j, l);
              const RCP<Thyra::MultiVectorBase<Scalar> > dgdp_mv =
                dgdp_deriv.getMultiVector();
              if (Teuchos::nonnull(dgdp_mv)) {
                // DxDp = J^-1 * DfDp
                const RCP<Thyra::MultiVectorBase<Scalar> > dxdp_mv =
                  Thyra::createMembers(model->get_x_space(), model->get_p_space(l));
                {
                  assign(dxdp_mv.ptr(), Teuchos::ScalarTraits<Scalar>::zero());

                  const RCP<Thyra::LinearOpWithSolveBase<Scalar> > jacobian =
                    modelOutArgs.get_W();
                  const RCP<const Thyra::MultiVectorBase<Scalar> > dfdp_mv =
                    modelOutArgs.get_DfDp(l).getMultiVector();

                  const Thyra::SolveCriteria<Scalar> defaultSolveCriteria;
                  const Thyra::SolveStatus<Scalar> solveStatus =
                    Thyra::solve(
                        *jacobian,
                        Thyra::NOTRANS,
                        *dfdp_mv,
                        dxdp_mv.ptr(),
                        Teuchos::ptr(&defaultSolveCriteria));
                  TEUCHOS_TEST_FOR_EXCEPTION(
                      solveStatus.solveStatus == Thyra::SOLVE_STATUS_UNCONVERGED,
                      std::runtime_error,
                      "Jacobian solver failed to converge");
                }
                // DgDp -= < DgDx, J^-1 * DfDp >
                if (dgdp_deriv.getMultiVectorOrientation() == Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM) {
                  const RCP<const Thyra::MultiVectorBase<Scalar> > dgdx_mv =
                    modelOutArgs.get_DgDx(j).getMultiVector();
                  Thyra::apply(
                      *dxdp_mv,
                      Thyra::TRANS,
                      *dgdx_mv,
                      dgdp_mv.ptr(),
                      -Teuchos::ScalarTraits<Scalar>::one(),
                      Teuchos::ScalarTraits<Scalar>::one());
                } else {
                  const RCP<const Thyra::LinearOpBase<Scalar> > dgdx_op =
                    modelOutArgs.get_DgDx(j).getLinearOp();
                  Thyra::apply(
                      *dgdx_op,
                      Thyra::NOTRANS,
                      *dxdp_mv,
                      dgdp_mv.ptr(),
                      -Teuchos::ScalarTraits<Scalar>::one(),
                      Teuchos::ScalarTraits<Scalar>::one());
                }
              } else {
                //Linop
                const RCP<Thyra::DefaultAddedLinearOp<Scalar> > dgdp_op =
                  Teuchos::rcp_dynamic_cast<Thyra::DefaultAddedLinearOp<Scalar> >(dgdp_deriv.getLinearOp());
                // Redundant check
                TEUCHOS_ASSERT(Teuchos::nonnull(dgdp_op));

                // Build DgDx * (DfDx)^-1 * DfDp
                const RCP<const Thyra::LinearOpBase<Scalar> > dfdp_op =
                  modelOutArgs.get_DfDp(l).getLinearOp();
                const RCP<const Thyra::LinearOpBase<Scalar> > dfdx_inv_op =
                  Thyra::inverse<Scalar>(modelOutArgs.get_W());
                const RCP<const Thyra::LinearOpBase<Scalar> > dgdx_t_op =
                  modelOutArgs.get_DgDx(j).getLinearOp();
                const RCP<const Thyra::LinearOpBase<Scalar> > mult_op =
                  Thyra::multiply<Scalar>(dgdx_t_op, dfdx_inv_op, dfdp_op);

                // Build DpDg - (DgDx * (DfDx)^-1 * DfDp)
                const RCP<const Thyra::LinearOpBase<Scalar> > minus_mult_op =
                  Thyra::scale<Scalar>(-Teuchos::ScalarTraits<Scalar>::one(), mult_op);
                const RCP<const Thyra::LinearOpBase<Scalar> > model_dgdp_op =
                  modelOutArgs.get_DgDp(j, l).getLinearOp();
                Teuchos::Array<RCP<const Thyra::LinearOpBase<Scalar> > > op_args(2);
                op_args[0] = model_dgdp_op;
                op_args[1] = minus_mult_op;
                dgdp_op->initialize(op_args);
              }
            }
          }
        }
      }
    }
  }
}

template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Piro::NOXSolver<Scalar>::create_DgDp_op_impl(int j, int l) const
{
  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::rcp;

  const Array<RCP<const Thyra::LinearOpBase<Scalar> > > dummy(1, zero(this->get_g_space(j), this->get_p_space(l)));
  return rcp(new Thyra::DefaultAddedLinearOp<Scalar>(dummy));
}
