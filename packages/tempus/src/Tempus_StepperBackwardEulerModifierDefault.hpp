// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperBackwardEulerModifierDefault_hpp
#define Tempus_StepperBackwardEulerModifierDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperBackwardEulerModifierBase.hpp"


namespace Tempus {

/** \brief Default modifier for StepperBackwardEuler.
 *
 *  The default modifier provides no-op functionality for the modifier.
 *  See StepperBackwardEulerModifierBase for details on the algorithm.
 */
template<class Scalar>
class StepperBackwardEulerModifierDefault
  : virtual public Tempus::StepperBackwardEulerModifierBase<Scalar>
{
public:

  /// Constructor
  StepperBackwardEulerModifierDefault(){}

  /// Destructor
  virtual ~StepperBackwardEulerModifierDefault(){}

  /// Modify BackwardEuler Stepper.
  virtual void modify(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Teuchos::RCP<StepperBackwardEuler<Scalar> > /* stepper */,
    const typename StepperBackwardEulerAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    switch(actLoc) {
      case StepperBackwardEulerAppAction<Scalar>::BEGIN_STEP:
      case StepperBackwardEulerAppAction<Scalar>::BEFORE_SOLVE:
      case StepperBackwardEulerAppAction<Scalar>::AFTER_SOLVE:
      case StepperBackwardEulerAppAction<Scalar>::END_STEP:
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

#endif // Tempus_StepperBackwardEulerModifierDefault_hpp
