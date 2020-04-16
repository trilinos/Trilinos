// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperRKObserverDefault_hpp
#define Tempus_StepperRKObserverDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperRKObserverBase.hpp"


namespace Tempus {

/** \brief Default observer for StepperRK.
 *
 *  The default observer provides no-op functionality for the observer.
 *  See StepperRKObserverBase for details on the algorithm.
 */
template<class Scalar>
class StepperRKObserverDefault
  : virtual public Tempus::StepperRKObserverBase<Scalar>
{
public:

  /// Constructor
  StepperRKObserverDefault(){}

  /// Destructor
  virtual ~StepperRKObserverDefault(){}

  /// Observe RK Stepper at end of takeStep.
  virtual void observe(
    Teuchos::RCP<const SolutionHistory<Scalar> > /* sh */,
    Teuchos::RCP<const StepperERK<Scalar> > /* stepper */,
    const typename StepperRKAppAction<Scalar>::ACTION_LOCATION actLoc) const
  {
    switch(actLoc) {
      case StepperRKAppAction<Scalar>::BEGIN_STEP:
      case StepperRKAppAction<Scalar>::BEFORE_SOLVE:
      case StepperRKAppAction<Scalar>::AFTER_SOLVE:
      case StepperRKAppAction<Scalar>::END_STEP:
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

#endif // Tempus_StepperRKObserverDefault_hpp
