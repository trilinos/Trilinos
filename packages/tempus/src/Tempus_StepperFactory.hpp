#ifndef Tempus_StepperFactory_hpp
#define Tempus_StepperFactory_hpp

#include "Teuchos_ParameterList.hpp"
#include "Tempus_StepperForwardEuler.hpp"
#include "Tempus_StepperBackwardEuler.hpp"
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

  /// Very simple factory method
  Teuchos::RCP<Stepper<Scalar> > createStepper(
    Teuchos::RCP<Teuchos::ParameterList>                      stepperPL,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model     )
  {
    using Teuchos::rcp;
    std::string stepper = stepperPL->get<std::string>("Stepper Type");
    if (stepper == "Forward Euler")
      return rcp(new StepperForwardEuler<Scalar>(stepperPL, model));
    else if (stepper == "Backward Euler")
      return rcp(new StepperBackwardEuler<Scalar>(stepperPL, model));
    else if (
      stepper == "RK Forward Euler" ||
      stepper == "RK Explicit 4 Stage" ||
      stepper == "RK Explicit 3/8 Rule" ||
      stepper == "RK Explicit 4 Stage 3rd order by Runge" ||
      stepper == "RK Explicit 5 Stage 3rd order by Kinnmark and Gray"||
      stepper == "RK Explicit 3 Stage 3rd order" ||
      stepper == "RK Explicit 3 Stage 3rd order TVD" ||
      stepper == "RK Explicit 3 Stage 3rd order by Heun" ||
      stepper == "RK Explicit 2 Stage 2nd order by Runge" ||
      stepper == "RK Explicit Trapezoidal" ||
      stepper == "General ERK" )
      return rcp(new StepperExplicitRK<Scalar>(stepperPL, model));
    else if (
      stepper == "SDIRK 1 Stage 1st order" ||
      stepper == "SDIRK 2 Stage 2nd order" ||
      stepper == "SDIRK 2 Stage 3rd order" ||
      stepper == "SDIRK 3 Stage 4th order" ||
      stepper == "SDIRK 5 Stage 4th order" ||
      stepper == "SDIRK 5 Stage 5th order"
      //stepper == "DIRK 2 Stage 3rd order"
      )
      return rcp(new StepperDIRK<Scalar>(stepperPL, model));
    else if (
      stepper == "RK Backward Euler" ||
      stepper == "RK Implicit 3 Stage 6th Order Kuntzmann & Butcher" ||
      stepper == "RK Implicit 4 Stage 8th Order Kuntzmann & Butcher" ||
      stepper == "RK Implicit 2 Stage 4th Order Hammer & Hollingsworth" ||
      stepper == "IRK 1 Stage Theta Method" ||
      stepper == "IRK 2 Stage Theta Method" ||
      stepper == "RK Implicit 1 Stage 2nd order Gauss" ||
      stepper == "RK Implicit 2 Stage 4th order Gauss" ||
      stepper == "RK Implicit 3 Stage 6th order Gauss" ||
      stepper == "RK Implicit 1 Stage 1st order Radau left" ||
      stepper == "RK Implicit 2 Stage 3rd order Radau left" ||
      stepper == "RK Implicit 3 Stage 5th order Radau left" ||
      stepper == "RK Implicit 1 Stage 1st order Radau right" ||
      stepper == "RK Implicit 2 Stage 3rd order Radau right" ||
      stepper == "RK Implicit 3 Stage 5th order Radau right" ||
      stepper == "RK Implicit 2 Stage 2nd order Lobatto A" ||
      stepper == "RK Implicit 3 Stage 4th order Lobatto A" ||
      stepper == "RK Implicit 4 Stage 6th order Lobatto A" ||
      stepper == "RK Implicit 2 Stage 2nd order Lobatto B" ||
      stepper == "RK Implicit 3 Stage 4th order Lobatto B" ||
      stepper == "RK Implicit 4 Stage 6th order Lobatto B" ||
      stepper == "RK Implicit 2 Stage 2nd order Lobatto C" ||
      stepper == "RK Implicit 3 Stage 4th order Lobatto C" ||
      stepper == "RK Implicit 4 Stage 6th order Lobatto C" ) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - Implicit RK not implemented yet!.\n");
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Unknown 'Stepper Type' = " << stepper);
    }
    return Teuchos::null; // Should never get here!
  }

};


} // namespace Tempus
#endif // Tempus_StepperFactory_hpp
