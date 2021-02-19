// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Tempus_ExplicitTemplateInstantiation.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION
#include "Tempus_StepperSubcycling.hpp"
#include "Tempus_StepperSubcycling_impl.hpp"

namespace Tempus {

  TEMPUS_INSTANTIATE_TEMPLATE_CLASS(StepperSubcycling)

  // Nonmember constructor
  template Teuchos::RCP<StepperSubcycling<double> >
  createStepperSubcycling(
    const Teuchos::RCP<const Thyra::ModelEvaluator<double> >& model,
    Teuchos::RCP<Teuchos::ParameterList> pl);

}

#endif
