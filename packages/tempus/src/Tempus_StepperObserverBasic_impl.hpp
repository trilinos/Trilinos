// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperObserverBasic_impl_hpp
#define Tempus_StepperObserverBasic_impl_hpp

#include "Tempus_RKButcherTableauBuilder.hpp"
#include "Tempus_config.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Tempus_WrapperModelEvaluatorBasic.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "NOX_Thyra.H"


namespace Tempus {

template<class Scalar> StepperObserverBasic<Scalar>::StepperObserverBasic(){}

template<class Scalar> StepperObserverBasic<Scalar>::~StepperObserverBasic(){}

template<class Scalar>
void StepperObserverBasic<Scalar>::observeBeginTakeStep(
  Teuchos::RCP<SolutionHistory<Scalar> > sh,
  Stepper<Scalar> & stepper)
{}

template<class Scalar>
void StepperObserverBasic<Scalar>::observeEndTakeStep(
  Teuchos::RCP<SolutionHistory<Scalar> > sh,
  Stepper<Scalar> & stepper)
{}

} // namespace Tempus
#endif // Tempus_StepperObserverBasic_impl_hpp
