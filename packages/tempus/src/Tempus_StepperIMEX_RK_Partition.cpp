// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Tempus_ExplicitTemplateInstantiation.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION
#include "Tempus_StepperIMEX_RK_Partition.hpp"
#include "Tempus_StepperIMEX_RK_Partition_impl.hpp"

namespace Tempus {

  TEMPUS_INSTANTIATE_TEMPLATE_CLASS(StepperIMEX_RK_Partition)

  // Nonmember constructor
  template Teuchos::RCP<StepperIMEX_RK_Partition<double> >
  createStepperIMEX_RK_Partition(
    const Teuchos::RCP<const Thyra::ModelEvaluator<double> >& model,
    std::string stepperType,
    Teuchos::RCP<Teuchos::ParameterList> pl);

}

#endif
