// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_ModelEvaluatorHelpers.hpp"

namespace Piro {

namespace Detail {

Teuchos::Array<bool> parseResponseParameters(Teuchos::ParameterList &params, int responseCount);
int parseResponseIndex(Teuchos::ParameterList &params);
bool parseSensitivityParameters(Teuchos::ParameterList &params);
Teuchos::Array<bool> createResponseTableFromIndex(int index, int responseCount);

template <typename Scalar, typename VectorType, typename MultiVectorType>
void PerformSolveImpl(
    const Thyra::ModelEvaluator<Scalar> &model,
    const Teuchos::Array<bool> &computeResponses,
    bool computeSensitivities,
    Teuchos::Array<Teuchos::RCP<VectorType> > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<MultiVectorType> > > &sensitivities,
    Teuchos::Array<Teuchos::RCP<MultiVectorType> > &directions,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<MultiVectorType> > > &reducedHessian,
    Teuchos::RCP<SolutionObserverBase<Scalar, VectorType> > observer)
{
  const int responseCount = model.Ng();
  const int parameterCount = model.Np();

  // Empty function in/out arguments
  responses.resize(responseCount);
  sensitivities.resize(responseCount, Teuchos::Array<Teuchos::RCP<MultiVectorType> >(parameterCount));
  reducedHessian.resize(responseCount, Teuchos::Array<Teuchos::RCP<MultiVectorType> >(parameterCount));

  // Setup input arguments
  Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs = model.createInArgs();
  for (int l = 0; l < parameterCount; ++l) {
    inArgs.set_p_direction(l, directions[l]);
  }

  // Setup output arguments
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs = model.createOutArgs();
  {
    for (int j = 0; j < responseCount; ++j) {
      if (computeResponses[j]) {
        const Teuchos::RCP<Thyra::VectorBase<Scalar> > g = Thyra::createMember(*model.get_g_space(j));
        outArgs.set_g(j, g);
        responses[j] = g;

        if(Teuchos::nonnull(observer))
           observer->observeResponse(j, Teuchos::rcpFromRef(outArgs),
                Teuchos::rcpFromRef(responses), g);

        if (computeSensitivities) {
          for (int l = 0; l < parameterCount; ++l) {
            const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support =
              outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l);
            Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dgdp_orient;
            if (dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM)||dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
              //Somehow, when DERIV_MV_JACOBIAN_FORM is supported, also  DERIV_MV_GRADIENT_FORM is.
              //When p is a distributed parameter, only DERIV_MV_GRADIENT_FORM is supported.
              dgdp_orient =  (dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM)) ?
                Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM:
                Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM;

              const Thyra::ModelEvaluatorBase::DerivativeMultiVector<Scalar> dgdp =
                Thyra::create_DgDp_mv(model, j, l, dgdp_orient);
              outArgs.set_DgDp(j, l, dgdp);
              sensitivities[j][l] = dgdp.getMultiVector();
            }
          }
        }
#ifdef Thyra_BUILD_HESSIAN_SUPPORT
        for (int l = 0; l < parameterCount; ++l) {
          if (Teuchos::nonnull(inArgs.get_p_direction(l))) {
            if (outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_pp, j, l, l)) {
              const Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > hess_g_pp =
                Thyra::createMembers(inArgs.get_p_direction(l)->col(0)->space(),inArgs.get_p_direction(l)->domain()->dim());
              outArgs.set_hess_vec_prod_g_pp(j, l, l, hess_g_pp);
              reducedHessian[j][l] = hess_g_pp;
            }
          }
        }
#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT
      }
    }
  }

/*
   TODO: This improper use of "reportFinalPoint" needs to be redone, but will
   require redesign of this section of code.
*/

  Thyra::ModelEvaluator<Scalar>* mutable_modelPtr =
         const_cast<Thyra::ModelEvaluator<Scalar>* >(&model);

  Thyra::ModelEvaluatorBase::InArgs<Scalar> finalPoint;

  mutable_modelPtr->reportFinalPoint(finalPoint, /*unused*/ true);

  // Solve the problem using the default values for the parameters
  model.evalModel(inArgs, outArgs);

  // If the model has set that it provides IN_ARG_x, get the final solution
  if(finalPoint.supports(Thyra::ModelEvaluatorBase::IN_ARG_x)){
        Teuchos::RCP<const Thyra::VectorBase<Scalar> > x_init = finalPoint.get_x();
        responses[responseCount-1] = Teuchos::rcp_const_cast< Thyra::VectorBase<Scalar> >(x_init);
  }

}

template <typename Scalar>
void PerformSolveImpl(
    const Thyra::ModelEvaluator<Scalar> &piroModel,
    int responseIndex,
    Teuchos::RCP<Thyra::VectorBase<Scalar> > &response)
{
  const Teuchos::Array<bool> computeResponses =
    createResponseTableFromIndex(responseIndex, piroModel.Ng());
  const bool computeSensitivities = false;

  Teuchos::Array<Teuchos::RCP<Thyra::VectorBase<Scalar> > > responses;
  Teuchos::Array<Teuchos::Array<Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > > > sensitivities;
  Teuchos::Array<Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > > directions(piroModel.Np());
  Teuchos::Array<Teuchos::Array<Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > > > reducedHessian;
  Teuchos::RCP<SolutionObserverBase<Scalar, Thyra::VectorBase<Scalar> > > observer;
  PerformSolveImpl(piroModel, computeResponses, computeSensitivities, responses, sensitivities, directions, reducedHessian, observer);

  response = responses[responseIndex];
}

template <typename Scalar, typename VectorType, typename MultiVectorType>
void PerformSolveImpl(
    const Thyra::ModelEvaluator<Scalar> &model,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<VectorType> > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<MultiVectorType> > > &sensitivities)
{
  const int responseCount = model.Ng();

  // Parse parameters
  const Teuchos::Array<bool> computeResponses = parseResponseParameters(solveParams, responseCount);
  const bool computeSensitivities = parseSensitivityParameters(solveParams);

  Teuchos::Array<Teuchos::RCP<MultiVectorType> > directions(model.Np());
  Teuchos::Array<Teuchos::Array<Teuchos::RCP<MultiVectorType> > > reducedHessian;
  Teuchos::RCP<SolutionObserverBase<Scalar, VectorType> > observer;
  PerformSolveImpl(model, computeResponses, computeSensitivities, responses, sensitivities, directions, reducedHessian, observer);
}

template <typename Scalar, typename VectorType, typename MultiVectorType>
void PerformSolveImpl(
    const Thyra::ModelEvaluator<Scalar> &model,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<VectorType> > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<MultiVectorType> > > &sensitivities,
    Teuchos::RCP<SolutionObserverBase<Scalar, VectorType> > observer)
{
  const int responseCount = model.Ng();

  // Parse parameters
  const Teuchos::Array<bool> computeResponses = parseResponseParameters(solveParams, responseCount);
  const bool computeSensitivities = parseSensitivityParameters(solveParams);

  Teuchos::Array<Teuchos::RCP<MultiVectorType> > directions(model.Np());
  Teuchos::Array<Teuchos::Array<Teuchos::RCP<MultiVectorType> > > reducedHessian;
  PerformSolveImpl(model, computeResponses, computeSensitivities, responses, sensitivities, directions, reducedHessian, observer);
}

template <typename Scalar, typename VectorType, typename MultiVectorType>
void PerformSolveImpl(
    const Thyra::ModelEvaluator<Scalar> &model,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<VectorType> > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<MultiVectorType> > > &sensitivities,
    Teuchos::Array<Teuchos::RCP<MultiVectorType> > &directions,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<MultiVectorType> > > &reducedHessian)
{
  const int responseCount = model.Ng();

  // Parse parameters
  const Teuchos::Array<bool> computeResponses = parseResponseParameters(solveParams, responseCount);
  const bool computeSensitivities = parseSensitivityParameters(solveParams);

  Teuchos::RCP<SolutionObserverBase<Scalar, VectorType> > observer;
  PerformSolveImpl(model, computeResponses, computeSensitivities, responses, sensitivities, directions, reducedHessian, observer);
}

template <typename Scalar, typename VectorType, typename MultiVectorType>
void PerformSolveImpl(
    const Thyra::ModelEvaluator<Scalar> &model,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<VectorType> > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<MultiVectorType> > > &sensitivities,
    Teuchos::Array<Teuchos::RCP<MultiVectorType> > &directions,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<MultiVectorType> > > &reducedHessian,
    Teuchos::RCP<SolutionObserverBase<Scalar, VectorType> > observer)
{
  const int responseCount = model.Ng();

  // Parse parameters
  const Teuchos::Array<bool> computeResponses = parseResponseParameters(solveParams, responseCount);
  const bool computeSensitivities = parseSensitivityParameters(solveParams);

  PerformSolveImpl(model, computeResponses, computeSensitivities, responses, sensitivities, directions, reducedHessian, observer);
}

} // namespace Detail

// PerformSolve overlaods (Simple wrappers around PerformSolveBase)

template <typename Scalar>
void PerformSolve(
    const Thyra::ResponseOnlyModelEvaluatorBase<Scalar> &piroModel,
    Teuchos::RCP<Thyra::VectorBase<Scalar> > &response)
{
  PerformSolveBase(piroModel, response);
}

template <typename Scalar>
void PerformSolve(
    const Thyra::ResponseOnlyModelEvaluatorBase<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::RCP<Thyra::VectorBase<Scalar> > &response)
{
  PerformSolveBase(piroModel, solveParams, response);
}

template <typename Scalar>
void PerformSolve(
    const Thyra::ResponseOnlyModelEvaluatorBase<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::RCP<Thyra::VectorBase<Scalar> > &response,
    Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > &sensitivity)
{
  PerformSolveBase(piroModel, solveParams, response, sensitivity);
}

template <typename Scalar>
void PerformSolve(
    const Thyra::ResponseOnlyModelEvaluatorBase<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<Thyra::VectorBase<Scalar> > > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > > > &sensitivities)
{
  PerformSolveBase(piroModel, solveParams, responses, sensitivities);
}

template <typename Scalar>
void PerformSolve(
    const Thyra::ResponseOnlyModelEvaluatorBase<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > > &sensitivities)
{
  PerformSolveBase(piroModel, solveParams, responses, sensitivities);
}

template <typename Scalar>
void PerformSolve(
    const Thyra::ResponseOnlyModelEvaluatorBase<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<Thyra::VectorBase<Scalar> > > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > > > &sensitivities,
    Teuchos::Array<Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > > &directions,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > > > &reducedHessian)
{
  PerformSolveBase(piroModel, solveParams, responses, sensitivities, directions, reducedHessian);
}

template <typename Scalar>
void PerformSolve(
    const Thyra::ResponseOnlyModelEvaluatorBase<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > > &sensitivities,
    Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > &directions,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > > &reducedHessian)
{
  PerformSolveBase(piroModel, solveParams, responses, sensitivities, directions, reducedHessian);
}

// PerformSolveBase overloads (less statically safe versions)

template <typename Scalar>
void PerformSolveBase(
    const Thyra::ModelEvaluator<Scalar> &piroModel,
    Teuchos::RCP<Thyra::VectorBase<Scalar> > &response)
{
  Detail::PerformSolveImpl(piroModel, 0, response);
}

template <typename Scalar>
void PerformSolveBase(
    const Thyra::ModelEvaluator<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::RCP<Thyra::VectorBase<Scalar> > &response)
{
  const int index = Detail::parseResponseIndex(solveParams);
  Detail::PerformSolveImpl(piroModel, index, response);
}

template <typename Scalar>
void PerformSolveBase(
    const Thyra::ModelEvaluator<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::RCP<Thyra::VectorBase<Scalar> > &response,
    Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > &sensitivity)
{
  const int responseCount = piroModel.Ng();
  const int index = Detail::parseResponseIndex(solveParams);
  const Teuchos::Array<bool> computeResponses = Detail::createResponseTableFromIndex(index, responseCount);
  const bool computeSensitivities = Detail::parseSensitivityParameters(solveParams);

  Teuchos::Array<Teuchos::RCP<Thyra::VectorBase<Scalar> > > responses;
  Teuchos::Array<Teuchos::Array<Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > > > sensitivities;
  Teuchos::RCP<SolutionObserverBase<Scalar, Thyra::VectorBase<Scalar> > > observer;
  Detail::PerformSolveImpl(piroModel, computeResponses, computeSensitivities, responses, sensitivities, observer);

  response = responses[index];
  sensitivity = (piroModel.Np() > 0) ? sensitivities[index].front() : Teuchos::null;
}

template <typename Scalar>
void PerformSolveBase(
    const Thyra::ModelEvaluator<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<Thyra::VectorBase<Scalar> > > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > > > &sensitivities)
{
  Detail::PerformSolveImpl(piroModel, solveParams, responses, sensitivities);
}

template <typename Scalar>
void PerformSolveBase(
    const Thyra::ModelEvaluator<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > > &sensitivities)
{
  Detail::PerformSolveImpl(piroModel, solveParams, responses, sensitivities);
}

template <typename Scalar>
void PerformSolveBase(
    const Thyra::ModelEvaluator<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > > &sensitivities,
    Teuchos::RCP<SolutionObserverBase<Scalar, const Thyra::VectorBase<Scalar> > > observer)
{
  Detail::PerformSolveImpl(piroModel, solveParams, responses, sensitivities, observer);
}

template <typename Scalar>
void PerformSolveBase(
    const Thyra::ModelEvaluator<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<Thyra::VectorBase<Scalar> > > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > > > &sensitivities,
    Teuchos::Array<Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > > &directions,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > > > &reducedHessian)
{
  Detail::PerformSolveImpl(piroModel, solveParams, responses, sensitivities, directions, reducedHessian);
}

template <typename Scalar>
void PerformSolveBase(
    const Thyra::ModelEvaluator<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > > &sensitivities,
    Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > &directions,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > > &reducedHessian)
{
  Detail::PerformSolveImpl(piroModel, solveParams, responses, sensitivities, directions, reducedHessian);
}

template <typename Scalar>
void PerformSolveBase(
    const Thyra::ModelEvaluator<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > > &sensitivities,
    Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > &directions,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > > &reducedHessian,
    Teuchos::RCP<SolutionObserverBase<Scalar, const Thyra::VectorBase<Scalar> > > observer)
{
  Detail::PerformSolveImpl(piroModel, solveParams, responses, sensitivities, directions, reducedHessian, observer);
}

} // namespace Piro
