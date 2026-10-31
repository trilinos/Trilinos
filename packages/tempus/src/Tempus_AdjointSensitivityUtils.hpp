//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_AdjointSensitivityUtils_hpp
#define Tempus_AdjointSensitivityUtils_hpp

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TestForException.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"

#include "Tempus_TimeStepControl.hpp"
#include "Tempus_TimeEventList.hpp"

namespace Tempus {

/** \brief Type of response for transient adjoint sensitivities.
 *
 *  Selected with the "Response Type" parameter of the "Sensitivities"
 *  sublist. With g(x,p,t) the model response:
 *  - FinalTime:    J(p) = g(x(T),p,T)
 *  - TimeIntegral: J(p) = int_{t0}^{T} g(x(t),p,t) dt
 *  - SumOverTimes: J(p) = sum_k g(x(t_k),p,t_k), with the times t_k given by
 *                  the "Response Time List" parameter.
 */
enum class AdjointResponseType { FinalTime,
                                 TimeIntegral,
                                 SumOverTimes };

namespace AdjointSensitivityUtils {

/// Convert the "Response Type" parameter value to the enum.
inline AdjointResponseType parseResponseType(const std::string& name)
{
  if (name == "Final Time") return AdjointResponseType::FinalTime;
  if (name == "Time Integral") return AdjointResponseType::TimeIntegral;
  if (name == "Sum Over Times") return AdjointResponseType::SumOverTimes;
  TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::invalid_argument,
      "Error - invalid \"Response Type\" = \""
          << name
          << "\". Valid values are \"Final Time\", \"Time Integral\" and "
             "\"Sum Over Times\".\n");
  return AdjointResponseType::FinalTime;  // not reached
}

/// Convert a comma-separated list of times (as in "Output Time List") to a
/// sorted vector.
template <class Scalar>
std::vector<Scalar> parseTimeList(const std::string& str)
{
  std::vector<Scalar> times;
  std::string::size_type lastPos = str.find_first_not_of(", ", 0);
  std::string::size_type pos     = str.find_first_of(",", lastPos);
  while (lastPos != std::string::npos) {
    std::string token = str.substr(lastPos, pos - lastPos);
    token.erase(token.find_last_not_of(" \t\n\r") + 1);
    if (!token.empty()) {
      std::size_t numParsed = 0;
      double value          = 0.0;
      try {
        value = std::stod(token, &numParsed);
      }
      catch (const std::exception&) {
        numParsed = 0;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(
          numParsed != token.size(), std::invalid_argument,
          "Error - could not parse \"" << token << "\" in the time list \""
                                       << str << "\".\n");
      times.push_back(Scalar(value));
    }
    if (pos == std::string::npos) break;
    lastPos = str.find_first_not_of(", ", pos);
    pos     = str.find_first_of(",", lastPos);
  }
  std::sort(times.begin(), times.end());
  return times;
}

/** \brief Evaluate dg/dx^T (gradient form, x_space x num_g) of response
 *  g_index at the given inArgs.
 *
 *  If dgdx_dot_check is true and the model supports dg/dx_dot in gradient
 *  form, dg/dx_dot is evaluated in the same call and an exception is thrown
 *  if it is nonzero, since the adjoint formulation requires responses that do
 *  not depend on x_dot.
 */
template <class Scalar>
void evalResponseGradientX(const Thyra::ModelEvaluator<Scalar>& model,
                           const int g_index,
                           const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
                           Thyra::MultiVectorBase<Scalar>& dgdx,
                           const bool dgdx_dot_check, const std::string& caller)
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::OutArgs<Scalar> outArgs = model.createOutArgs();
  TEUCHOS_TEST_FOR_EXCEPTION(
      !outArgs.supports(MEB::OUT_ARG_DgDx, g_index)
           .supports(MEB::DERIV_MV_GRADIENT_FORM),
      std::logic_error,
      "Error - " << caller << ": the model must support dg/dx in gradient "
                 << "(DERIV_MV_GRADIENT_FORM) form for response " << g_index
                 << ".\n");
  outArgs.set_DgDx(g_index, MEB::Derivative<Scalar>(Teuchos::rcpFromRef(dgdx),
                                                    MEB::DERIV_MV_GRADIENT_FORM));

  Teuchos::RCP<Thyra::MultiVectorBase<Scalar>> dgdx_dot;
  if (dgdx_dot_check &&
      outArgs.supports(MEB::OUT_ARG_DgDx_dot, g_index).supports(MEB::DERIV_MV_GRADIENT_FORM)) {
    dgdx_dot = Thyra::createMembers(model.get_x_space(), dgdx.domain()->dim());
    outArgs.set_DgDx_dot(g_index, MEB::Derivative<Scalar>(dgdx_dot, MEB::DERIV_MV_GRADIENT_FORM));
  }
  model.evalModel(inArgs, outArgs);

  if (dgdx_dot != Teuchos::null) {
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;
    Teuchos::Array<Mag> norms(dgdx_dot->domain()->dim());
    Thyra::norms_inf(*dgdx_dot, norms());
    const Mag maxNorm = *std::max_element(norms.begin(), norms.end());
    TEUCHOS_TEST_FOR_EXCEPTION(
        maxNorm != Teuchos::ScalarTraits<Mag>::zero(), std::logic_error,
        "Error - " << caller << ": response " << g_index
                   << " depends on x_dot (|dg/dx_dot| = " << maxNorm
                   << "). The \"Time Integral\" and \"Sum Over Times\" "
                   << "adjoint responses require g = g(x,p,t).\n");
  }
}

/** \brief Evaluate dg/dp^T (gradient form, p_space x num_g) of response
 *  g_index with respect to parameter p_index at the given inArgs.
 *
 *  Models that only provide the Jacobian form are supported by transposing,
 *  as done for the final-time response.
 */
template <class Scalar>
void evalResponseGradientP(const Thyra::ModelEvaluator<Scalar>& model,
                           const int g_index, const int p_index,
                           const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
                           Thyra::MultiVectorBase<Scalar>& dgdp)
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::OutArgs<Scalar> outArgs = model.createOutArgs();
  MEB::DerivativeSupport dgdp_support =
      outArgs.supports(MEB::OUT_ARG_DgDp, g_index, p_index);
  if (dgdp_support.supports(MEB::DERIV_MV_GRADIENT_FORM)) {
    outArgs.set_DgDp(g_index, p_index,
                     MEB::Derivative<Scalar>(Teuchos::rcpFromRef(dgdp),
                                             MEB::DERIV_MV_GRADIENT_FORM));
    model.evalModel(inArgs, outArgs);
  }
  else if (dgdp_support.supports(MEB::DERIV_MV_JACOBIAN_FORM)) {
    const int num_g = dgdp.domain()->dim();
    const int num_p = dgdp.range()->dim();
    Teuchos::RCP<Thyra::MultiVectorBase<Scalar>> dgdp_trans =
        Thyra::createMembers(model.get_g_space(g_index), num_p);
    outArgs.set_DgDp(g_index, p_index, MEB::Derivative<Scalar>(dgdp_trans, MEB::DERIV_MV_JACOBIAN_FORM));
    model.evalModel(inArgs, outArgs);
    Thyra::DetachedMultiVectorView<Scalar> dgdp_view(dgdp);
    Thyra::ConstDetachedMultiVectorView<Scalar> dgdp_trans_view(*dgdp_trans);
    for (int i = 0; i < num_p; ++i)
      for (int j = 0; j < num_g; ++j) dgdp_view(i, j) = dgdp_trans_view(j, i);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Invalid dg/dp support");
}

/** \brief Compute y = alpha * df/dp^T * v + beta * y.
 *
 *  v is a multivector in the residual (f) space and y a multivector in the
 *  space of parameter p_index, with the same number of columns.  df/dp is
 *  evaluated at inArgs in the form supported by the model (linear operator,
 *  Jacobian or gradient multivector).  dfdp_op and dfdp_mv cache the
 *  derivative objects between calls (pass null RCPs on the first call).
 */
template <class Scalar>
void applyDfDpTranspose(const Thyra::ModelEvaluator<Scalar>& model,
                        const int p_index,
                        const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
                        const Thyra::MultiVectorBase<Scalar>& v,
                        Thyra::MultiVectorBase<Scalar>& y, const Scalar alpha,
                        const Scalar beta,
                        Teuchos::RCP<Thyra::LinearOpBase<Scalar>>& dfdp_op,
                        Teuchos::RCP<Thyra::MultiVectorBase<Scalar>>& dfdp_mv)
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::OutArgs<Scalar> outArgs        = model.createOutArgs();
  MEB::DerivativeSupport dfdp_support = outArgs.supports(MEB::OUT_ARG_DfDp, p_index);
  Thyra::EOpTransp trans              = Thyra::CONJTRANS;
  if (dfdp_support.supports(MEB::DERIV_LINEAR_OP)) {
    if (dfdp_op == Teuchos::null) dfdp_op = model.create_DfDp_op(p_index);
    outArgs.set_DfDp(p_index, MEB::Derivative<Scalar>(dfdp_op));
  }
  else if (dfdp_support.supports(MEB::DERIV_MV_JACOBIAN_FORM)) {
    if (dfdp_mv == Teuchos::null)
      dfdp_mv = Thyra::createMembers(model.get_f_space(), model.get_p_space(p_index)->dim());
    outArgs.set_DfDp(p_index, MEB::Derivative<Scalar>(dfdp_mv, MEB::DERIV_MV_JACOBIAN_FORM));
    dfdp_op = dfdp_mv;
  }
  else if (dfdp_support.supports(MEB::DERIV_MV_GRADIENT_FORM)) {
    if (dfdp_mv == Teuchos::null)
      dfdp_mv = Thyra::createMembers(model.get_p_space(p_index), model.get_f_space()->dim());
    outArgs.set_DfDp(p_index, MEB::Derivative<Scalar>(dfdp_mv, MEB::DERIV_MV_GRADIENT_FORM));
    dfdp_op = dfdp_mv;
    trans   = Thyra::CONJ;
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Invalid df/dp support");
  model.evalModel(inArgs, outArgs);
  dfdp_op->apply(trans, v, Teuchos::ptrFromRef(y), alpha, beta);
}

/** \brief Create a solver for the adjoint model's W = (alpha*df/dx_dot + beta*df/dx)^T at inArgs.
 *
 *  Uses the W out-arg of the adjoint model if supported, otherwise its W_op,
 *  its W_factory, and the preconditioner of the W_factory (or of the model).
 */
template <class Scalar>
Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar>>
createAdjointLinearOpWithSolve(
    const Thyra::ModelEvaluator<Scalar>& adjoint_model,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs, const Scalar alpha,
    const Scalar beta)
{
  using Teuchos::RCP;
  using Thyra::LinearOpBase;
  using Thyra::LinearOpWithSolveBase;
  using Thyra::LinearOpWithSolveFactoryBase;
  using Thyra::PreconditionerBase;
  using Thyra::PreconditionerFactoryBase;
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgs<Scalar> inargs = inArgs;
  inargs.set_alpha(alpha);
  inargs.set_beta(beta);
  MEB::OutArgs<Scalar> adj_outargs = adjoint_model.createOutArgs();
  RCP<LinearOpWithSolveBase<Scalar>> W;
  if (adj_outargs.supports(MEB::OUT_ARG_W)) {
    // Model supports W
    W = adjoint_model.create_W();
    adj_outargs.set_W(W);
    adjoint_model.evalModel(inargs, adj_outargs);
  }
  else {
    // Otherwise model must support a W_op and W factory
    RCP<const LinearOpWithSolveFactoryBase<Scalar>> lowsfb =
        adjoint_model.get_W_factory();
    TEUCHOS_TEST_FOR_EXCEPTION(lowsfb == Teuchos::null, std::logic_error,
                               "Adjoint ME must support W out-arg or provide "
                               "a W_factory for non-identity mass matrix");

    // Compute W_op (and W_prec if supported)
    RCP<LinearOpBase<Scalar>> W_op = adjoint_model.create_W_op();
    adj_outargs.set_W_op(W_op);
    RCP<PreconditionerFactoryBase<Scalar>> prec_factory = lowsfb->getPreconditionerFactory();
    RCP<PreconditionerBase<Scalar>> W_prec;
    if (prec_factory != Teuchos::null)
      W_prec = prec_factory->createPrec();
    else if (adj_outargs.supports(MEB::OUT_ARG_W_prec)) {
      W_prec = adjoint_model.create_W_prec();
      adj_outargs.set_W_prec(W_prec);
    }
    adjoint_model.evalModel(inargs, adj_outargs);

    // Create and initialize W
    W = lowsfb->createOp();
    if (W_prec != Teuchos::null) {
      if (prec_factory != Teuchos::null)
        prec_factory->initializePrec(Thyra::defaultLinearOpSource<Scalar>(W_op), W_prec.get());
      Thyra::initializePreconditionedOp<Scalar>(*lowsfb, W_op, W_prec, W.ptr());
    }
    else
      Thyra::initializeOp<Scalar>(*lowsfb, W_op, W.ptr());
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
      W == Teuchos::null, std::logic_error,
      "A null W has been encountered in "
      "Tempus::AdjointSensitivityUtils::createAdjointLinearOpWithSolve!\n");
  return W;
}

/** \brief Add the response times to the time step control as time events on
 *  which the time steps land exactly.
 *
 *  This only adjusts the time steps for "Step Type" = "Variable"; with
 *  "Constant" time steps the response times must be on the time grid.
 *  The time steps are not flagged for output.
 */
template <class Scalar>
void addResponseTimeEvents(TimeStepControl<Scalar>& tsc,
                           const std::vector<Scalar>& times)
{
  if (times.empty()) return;
  tsc.getTimeEvents()->add(Teuchos::rcp(
      new TimeEventList<Scalar>(times, "Response Time List", true)));
}

}  // namespace AdjointSensitivityUtils
}  // namespace Tempus

#endif  // Tempus_AdjointSensitivityUtils_hpp
