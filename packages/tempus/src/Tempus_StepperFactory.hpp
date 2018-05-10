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
#include "Tempus_StepperBDF2.hpp"
#include "Tempus_StepperNewmarkImplicitAForm.hpp"
#include "Tempus_StepperNewmarkImplicitDForm.hpp"
#include "Tempus_StepperNewmarkExplicitAForm.hpp"
#include "Tempus_StepperHHTAlpha.hpp"
#include "Tempus_StepperExplicitRK.hpp"
#include "Tempus_StepperDIRK.hpp"
#include "Tempus_StepperIMEX_RK.hpp"
#include "Tempus_StepperIMEX_RK_Partition.hpp"
#include "Tempus_StepperLeapfrog.hpp"
#include "Tempus_StepperOperatorSplit.hpp"
#include "Tempus_StepperTrapezoidal.hpp"


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
    std::string stepperType = "Forward Euler")
  {
    if (stepperType == "") stepperType = "Forward Euler";
    return this->createStepper(model, stepperType, Teuchos::null);
  }

  /// Create stepper from ParameterList with its details.
  Teuchos::RCP<Stepper<Scalar> > createStepper(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    std::string stepperType = "Forward Euler";
    if (stepperPL != Teuchos::null)
      stepperType = stepperPL->get<std::string>("Stepper Type","Forward Euler");
    return this->createStepper(model, stepperType, stepperPL);
  }

  /// Create stepper from ParameterList with its details.
  Teuchos::RCP<Stepper<Scalar> > createStepper(
    std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > models,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    std::string stepperType = stepperPL->get<std::string>("Stepper Type");
    return this->createStepper(models, stepperType, stepperPL);
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
    else if (stepperType == "Trapezoidal Method")
      return rcp(new StepperTrapezoidal<Scalar>(model, stepperPL));
    else if (stepperType == "BDF2")
      return rcp(new StepperBDF2<Scalar>(model, stepperPL));
    else if (stepperType == "Newmark Implicit a-Form")
      return rcp(new StepperNewmarkImplicitAForm<Scalar>(model, stepperPL));
    else if (stepperType == "Newmark Implicit d-Form")
      return rcp(new StepperNewmarkImplicitDForm<Scalar>(model, stepperPL));
    else if (stepperType == "Newmark Explicit a-Form")
      return rcp(new StepperNewmarkExplicitAForm<Scalar>(model, stepperPL));
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
      stepperType == "Bogacki-Shampine 3(2) Pair" ||
      stepperType == "General ERK" )
      return rcp(new StepperExplicitRK<Scalar>(model, stepperType, stepperPL));
    else if (
      stepperType == "RK Backward Euler" ||
      stepperType == "IRK 1 Stage Theta Method" ||
      stepperType == "SDIRK 1 Stage 1st order" ||
      stepperType == "SDIRK 2 Stage 2nd order" ||
      stepperType == "SDIRK 2 Stage 3rd order" ||
      stepperType == "EDIRK 2 Stage 3rd order" ||
      stepperType == "EDIRK 2 Stage Theta Method" ||
      stepperType == "SDIRK 3 Stage 4th order" ||
      stepperType == "SDIRK 5 Stage 4th order" ||
      stepperType == "SDIRK 5 Stage 5th order" ||
      stepperType == "General DIRK"
      )
      return rcp(new StepperDIRK<Scalar>(model, stepperType, stepperPL));
    else if (
      stepperType == "RK Implicit 3 Stage 6th Order Kuntzmann & Butcher" ||
      stepperType == "RK Implicit 4 Stage 8th Order Kuntzmann & Butcher" ||
      stepperType == "RK Implicit 2 Stage 4th Order Hammer & Hollingsworth" ||
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
    else if (
      stepperType == "IMEX RK 1st order" ||
      stepperType == "IMEX RK SSP2"      ||
      stepperType == "IMEX RK ARS 233"   ||
      stepperType == "General IMEX RK" )
      return rcp(new StepperIMEX_RK<Scalar>(model, stepperType, stepperPL));
    else if (
      stepperType == "Partitioned IMEX RK 1st order" ||
      stepperType == "Partitioned IMEX RK SSP2"      ||
      stepperType == "Partitioned IMEX RK ARS 233"   ||
      stepperType == "General Partitioned IMEX RK" )
      return rcp(new StepperIMEX_RK_Partition<Scalar>(
                        model, stepperType, stepperPL));
    else if (stepperType == "Leapfrog")
      return rcp(new StepperLeapfrog<Scalar>(model, stepperPL));
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Unknown 'Stepper Type' = " << stepperType);
    }
  }

  Teuchos::RCP<Stepper<Scalar> > createStepper(
    std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > models,
    std::string stepperType,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    if (stepperType == "Operator Split")
      return rcp(new StepperOperatorSplit<Scalar>(models, stepperPL));
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Unknown 'Stepper Type' = " << stepperType);
    }
  }

};


} // namespace Tempus
#endif // Tempus_StepperFactory_hpp
