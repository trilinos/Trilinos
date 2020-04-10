// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperForwardEulerModifierDefault_hpp
#define Tempus_StepperForwardEulerModifierDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperForwardEulerModifierBase.hpp"


namespace Tempus {

/** \brief Default modifier for StepperForwardEuler.
 *
 *  The default modifier provides no-op functionality for the modifier.
 *  See StepperForwardEulerModifierBase for details on the algorithm.
 */
template<class Scalar>
class StepperForwardEulerModifierDefault
  : virtual public Tempus::StepperForwardEulerModifierBase<Scalar>
{
public:

  /// Constructor
  StepperForwardEulerModifierDefault(){}

  /// Destructor
  virtual ~StepperForwardEulerModifierDefault(){}

  /// Modify ForwardEuler Stepper.
  virtual void modify(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Teuchos::RCP<StepperForwardEuler<Scalar> > /* stepper */,
    const typename StepperForwardEulerAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    switch(actLoc) {
      case StepperForwardEulerAppAction<Scalar>::BEGIN_STEP:
      case StepperForwardEulerAppAction<Scalar>::BEFORE_EXPLICIT_EVAL:
      case StepperForwardEulerAppAction<Scalar>::END_STEP:
      {
        // No-op.
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - unknown action location.\n");
    }
  }

};

} // namespace Tempus

#endif // Tempus_StepperForwardEulerModifierDefault_hpp
