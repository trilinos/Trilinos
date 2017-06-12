// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperFactory_hpp
#define Tempus_StepperFactory_hpp

#include "Teuchos_ParameterList.hpp"
#include "Tempus_StepperForwardEuler.hpp"
#include "Tempus_StepperBackwardEuler.hpp"
#include "Tempus_StepperNewmark.hpp"
#include "Tempus_StepperNewmarkExplicit.hpp"
#include "Tempus_StepperHHTAlpha.hpp"
#include "Tempus_StepperExplicitRK.hpp"
#include "Tempus_StepperDIRK.hpp"


namespace Tempus {

/** \brief Stepper factory.
 *
 * <b>Adding Steppers</b>
 *    -#
 */
template<class Scalar>
class StepperFactory
{
public:

  /// Constructor
  StepperFactory(){}

  /// Destructor
  virtual ~StepperFactory() {}

  /// Create default stepper from stepper type (e.g., "Forward Euler").
  Teuchos::RCP<Stepper<Scalar> > createStepper(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    std::string stepperType = "")
  {
    if (stepperType == "") stepperType = "Forward Euler";
    return this->createStepper(model, stepperType, Teuchos::null);
  }

  /// Create stepper from ParameterList with its details.
  Teuchos::RCP<Stepper<Scalar> > createStepper(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    std::string stepperType = stepperPL->get<std::string>("Stepper Type");
    return this->createStepper(model, stepperType, stepperPL);
  }

private:
  /// Very simple factory method
  Teuchos::RCP<Stepper<Scalar> > createStepper(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    std::string stepperType,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    using Teuchos::rcp;
    if (stepperType == "Forward Euler")
      return rcp(new StepperForwardEuler<Scalar>(model, stepperPL));
    else if (stepperType == "Backward Euler")
      return rcp(new StepperBackwardEuler<Scalar>(model, stepperPL));
    else if (stepperType == "Newmark Beta")
      return rcp(new StepperNewmark<Scalar>(model, stepperPL));
    else if (stepperType == "Newmark Beta Explicit")
      return rcp(new StepperNewmarkExplicit<Scalar>(model, stepperPL));
    else if (stepperType == "HHT-Alpha")
      return rcp(new StepperHHTAlpha<Scalar>(model, stepperPL));
    else if (
      stepperType == "RK Forward Euler" ||
      stepperType == "RK Explicit 4 Stage" ||
      stepperType == "RK Explicit 3/8 Rule" ||
      stepperType == "RK Explicit 4 Stage 3rd order by Runge" ||
      stepperType == "RK Explicit 5 Stage 3rd order by Kinnmark and Gray"||
      stepperType == "RK Explicit 3 Stage 3rd order" ||
      stepperType == "RK Explicit 3 Stage 3rd order TVD" ||
      stepperType == "RK Explicit 3 Stage 3rd order by Heun" ||
      stepperType == "RK Explicit 2 Stage 2nd order by Runge" ||
      stepperType == "RK Explicit Trapezoidal" ||
      stepperType == "General ERK" )
      return rcp(new StepperExplicitRK<Scalar>(model, stepperType, stepperPL));
    else if (
      stepperType == "SDIRK 1 Stage 1st order" ||
      stepperType == "SDIRK 2 Stage 2nd order" ||
      stepperType == "SDIRK 2 Stage 3rd order" ||
      stepperType == "SDIRK 3 Stage 4th order" ||
      stepperType == "SDIRK 5 Stage 4th order" ||
      stepperType == "SDIRK 5 Stage 5th order"
      //stepperType == "DIRK 2 Stage 3rd order"
      )
      return rcp(new StepperDIRK<Scalar>(model, stepperType, stepperPL));
    else if (
      stepperType == "RK Backward Euler" ||
      stepperType == "RK Implicit 3 Stage 6th Order Kuntzmann & Butcher" ||
      stepperType == "RK Implicit 4 Stage 8th Order Kuntzmann & Butcher" ||
      stepperType == "RK Implicit 2 Stage 4th Order Hammer & Hollingsworth" ||
      stepperType == "IRK 1 Stage Theta Method" ||
      stepperType == "IRK 2 Stage Theta Method" ||
      stepperType == "RK Implicit 1 Stage 2nd order Gauss" ||
      stepperType == "RK Implicit 2 Stage 4th order Gauss" ||
      stepperType == "RK Implicit 3 Stage 6th order Gauss" ||
      stepperType == "RK Implicit 1 Stage 1st order Radau left" ||
      stepperType == "RK Implicit 2 Stage 3rd order Radau left" ||
      stepperType == "RK Implicit 3 Stage 5th order Radau left" ||
      stepperType == "RK Implicit 1 Stage 1st order Radau right" ||
      stepperType == "RK Implicit 2 Stage 3rd order Radau right" ||
      stepperType == "RK Implicit 3 Stage 5th order Radau right" ||
      stepperType == "RK Implicit 2 Stage 2nd order Lobatto A" ||
      stepperType == "RK Implicit 3 Stage 4th order Lobatto A" ||
      stepperType == "RK Implicit 4 Stage 6th order Lobatto A" ||
      stepperType == "RK Implicit 2 Stage 2nd order Lobatto B" ||
      stepperType == "RK Implicit 3 Stage 4th order Lobatto B" ||
      stepperType == "RK Implicit 4 Stage 6th order Lobatto B" ||
      stepperType == "RK Implicit 2 Stage 2nd order Lobatto C" ||
      stepperType == "RK Implicit 3 Stage 4th order Lobatto C" ||
      stepperType == "RK Implicit 4 Stage 6th order Lobatto C" ) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - Implicit RK not implemented yet!.\n");
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Unknown 'Stepper Type' = " << stepperType);
    }
  }

};


} // namespace Tempus
#endif // Tempus_StepperFactory_hpp
