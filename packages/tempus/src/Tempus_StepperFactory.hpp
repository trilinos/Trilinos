#ifndef TEMPUS_STEPPERFACTORY_HPP
#define TEMPUS_STEPPERFACTORY_HPP

//Teuchos
//#include "Teuchos_VerboseObject.hpp"
//#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterList.hpp"
//#include "Teuchos_StandardParameterEntryValidators.hpp"
//// Thyra
//#include "Thyra_ModelEvaluator.hpp"
//// Tempus
#include "Tempus_Stepper.hpp"
#include "Tempus_StepperForwardEuler.hpp"
#include "Tempus_StepperBackwardEuler.hpp"


namespace Tempus {

enum StepperType {
  FORWARD_EULER,
  BACKWARD_EULER
};

static std::string ForwardEuler_name  = "Forward Euler";
static std::string BackwardEuler_name = "Backward Euler";
static std::string Stepper_name       = "Stepper";
static std::string Stepper_default    = ForwardEuler_name;

Teuchos::Array<std::string> Stepper_names = Teuchos::tuple<std::string>(
  ForwardEuler_name,
  BackwardEuler_name);

const Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<StepperType> >
StepperValidator = Teuchos::rcp(
  new Teuchos::StringToIntegralParameterEntryValidator<StepperType>(
    Stepper_names,
    Teuchos::tuple<Tempus::StepperType>(
      FORWARD_EULER,
      BACKWARD_EULER),
    Stepper_name));


const std::string toString(const StepperType s)
{
  switch(s) {
    case FORWARD_EULER:
      return ForwardEuler_name;
    case BACKWARD_EULER:
      return BackwardEuler_name;
    default:
      TEUCHOS_TEST_FOR_EXCEPT("Invalid StepperType!");
  }
  return 0; // Should never get here!
}


const StepperType fromString(const std::string ss)
{
  if (ss == ForwardEuler_name)
    return FORWARD_EULER;
  else if (ss == BackwardEuler_name)
    return BACKWARD_EULER;
  else
    TEUCHOS_TEST_FOR_EXCEPT("Invalid String for StepperType!");

  return FORWARD_EULER; // Should never get here!
}


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

  /// Factory Method
  Teuchos::RCP<Stepper<Scalar> > createStepper(
    std::string                                         stepperName,
    Teuchos::RCP<Teuchos::ParameterList>                stepper_pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model         )
  {
    StepperType stepperType = fromString(stepperName);
    switch(stepperType) {
      case FORWARD_EULER:
        return Teuchos::rcp(
          new StepperForwardEuler<Scalar>(stepper_pList, model));
      case BACKWARD_EULER:
        return Teuchos::rcp(
          new StepperBackwardEuler<Scalar>(stepper_pList, model));
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Unknown StepperType = " << stepperName);
    }
    return Teuchos::null; // Should never get here!
  }

};


/// Return Stepper_name used for selecting the Stepper from the ParameterList
const std::string getStepperName() {return Stepper_name; }

/// Return the default Stepper
const std::string getStepperDefault() {return ForwardEuler_name; }


} // namespace Tempus
#endif // TEMPUS_STEPPERFACTORY_HPP
