//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_WrapCombinedFSAModelEvaluator_hpp
#define Tempus_WrapCombinedFSAModelEvaluator_hpp

#include "Tempus_config.hpp"
#include "Tempus_WrapperModelEvaluatorPairIMEX_CombinedFSA.hpp"
#include "Tempus_WrapperModelEvaluatorPairPartIMEX_CombinedFSA.hpp"

namespace Tempus {

/*! Helper function for creating a CombinedForwardSensitivityModelEvaluator
 * from a given application model evaluator.  It handles the complexity
 * introducted by IMEX steppers where the sensitivity model evaluator needs
 * to be put inside the IMEX pair model evaluators.
 */
template <typename Scalar>
Teuchos::RCP<SensitivityModelEvaluatorBase<Scalar> >
wrapCombinedFSAModelEvaluator(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >&
        sens_residual_model,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& sens_solve_model,
    const Teuchos::RCP<const Teuchos::ParameterList>& pList = Teuchos::null)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  RCP<SensitivityModelEvaluatorBase<Scalar> > wrapped_model;

  // Test if model is an IMEX pair
  RCP<const WrapperModelEvaluatorPairIMEX_Basic<Scalar> > modelPairIMEX =
      rcp_dynamic_cast<const WrapperModelEvaluatorPairIMEX_Basic<Scalar> >(
          model);
  RCP<const WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >
      modelPairPartIMEX = rcp_dynamic_cast<
          const WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >(model);

  // It isn't clear how to handle separate sensitivity and solve model
  // evaluators in the IMEX case, so punt for now.  Would they also be IMEX MEs?
  if ((modelPairIMEX != Teuchos::null || modelPairPartIMEX != Teuchos::null) &&
      (model.ptr() != sens_residual_model.ptr()) &&
      (model.ptr() != sens_solve_model.ptr()))
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Unique model evaluators for state, sensitivity residual, and "
        "sensitivity solve is not supported for IMEX");

  if (modelPairIMEX != Teuchos::null) {
    wrapped_model = rcp(new WrapperModelEvaluatorPairIMEX_CombinedFSA<Scalar>(
        modelPairIMEX, pList));
  }
  else if (modelPairPartIMEX != Teuchos::null) {
    wrapped_model =
        rcp(new WrapperModelEvaluatorPairPartIMEX_CombinedFSA<Scalar>(
            modelPairPartIMEX, pList));
  }
  else {
    wrapped_model = rcp(new CombinedForwardSensitivityModelEvaluator<Scalar>(
        model, sens_residual_model, sens_solve_model, pList));
  }

  return wrapped_model;
}

template <typename Scalar>
Teuchos::RCP<SensitivityModelEvaluatorBase<Scalar> >
wrapCombinedFSAModelEvaluator(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& sens_residual_model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& sens_solve_model,
    const Teuchos::RCP<const Teuchos::ParameterList>& pList = Teuchos::null)
{
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > cmodel = model;
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > csens_residual_model =
      sens_residual_model;
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > csens_solve_model =
      sens_solve_model;
  return wrapCombinedFSAModelEvaluator(cmodel, csens_residual_model,
                                       csens_solve_model, pList);
}

}  // namespace Tempus

#endif
