// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_WrapCombinedFSAModelEvaluator_hpp
#define Tempus_WrapCombinedFSAModelEvaluator_hpp

#include "Tempus_WrapperModelEvaluatorPairIMEX_CombinedFSA.hpp"
#include "Tempus_WrapperModelEvaluatorPairPartIMEX_CombinedFSA.hpp"

namespace Tempus {

/*! Helper function for creating a CombinedForwardSensitivityModelEvaluator
 * from a given application model evaluator.  It handles the complexity
 * introducted by IMEX steppers where the sensitivity model evaluator needs
 * to be put inside the IMEX pair model evaluators.
 */
template <typename Scalar>
Teuchos::RCP< SensitivityModelEvaluatorBase<Scalar> >
wrapCombinedFSAModelEvaluator(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > & model,
  const Teuchos::RCP<const Teuchos::ParameterList>& pList = Teuchos::null)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  RCP<SensitivityModelEvaluatorBase<Scalar> > wrapped_model;

  // Test if model is an IMEX pair
  RCP<const WrapperModelEvaluatorPairIMEX_Basic<Scalar> > modelPairIMEX =
    rcp_dynamic_cast<const WrapperModelEvaluatorPairIMEX_Basic<Scalar> >(model);
  RCP<const WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> > modelPairPartIMEX =
    rcp_dynamic_cast<const WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >(model);
  if (modelPairIMEX != Teuchos::null) {
    wrapped_model =
      rcp(new WrapperModelEvaluatorPairIMEX_CombinedFSA<Scalar>(
            modelPairIMEX, pList));
  }
  else if (modelPairPartIMEX != Teuchos::null) {
    wrapped_model =
      rcp(new WrapperModelEvaluatorPairPartIMEX_CombinedFSA<Scalar>(
            modelPairPartIMEX, pList));
  }
  else {
    wrapped_model =
      rcp(new CombinedForwardSensitivityModelEvaluator<Scalar>(model, pList));
  }

  return wrapped_model;
}

template <typename Scalar>
Teuchos::RCP< SensitivityModelEvaluatorBase<Scalar> >
wrapCombinedFSAModelEvaluator(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > & model,
  const Teuchos::RCP<const Teuchos::ParameterList>& pList = Teuchos::null)
{
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > cmodel = model;
  return wrapCombinedFSAModelEvaluator(cmodel, pList);
}

} // namespace Tempus

#endif
