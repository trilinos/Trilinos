// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperBackwardEulerObserverDefault_hpp
#define Tempus_StepperBackwardEulerObserverDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperBackwardEulerObserverBase.hpp"


namespace Tempus {

/** \brief Default observer for StepperBackwardEuler.
 *
 *  The default observer provides no-op functionality for the observer.
 *  See StepperBackwardEulerObserverBase for details on the algorithm.
 */
template<class Scalar>
class StepperBackwardEulerObserverDefault
  : virtual public Tempus::StepperBackwardEulerObserverBase<Scalar>
{
public:

  /// Constructor
  StepperBackwardEulerObserverDefault(){}

  /// Destructor
  virtual ~StepperBackwardEulerObserverDefault(){}

  /// Observe BackwardEuler Stepper at end of takeStep.
  virtual void observe(
    Teuchos::RCP<const SolutionHistory<Scalar> > /* sh */,
    Teuchos::RCP<const StepperBackwardEuler<Scalar> > /* stepper */,
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

#endif // Tempus_StepperBackwardEulerObserverDefault_hpp
