// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperOperatorSplitObserverDefault_hpp
#define Tempus_StepperOperatorSplitObserverDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperOperatorSplitObserverBase.hpp"

namespace Tempus {

/** \brief Default observer for StepperOperatorSplit.
 *
 *  The default observer provides no-op functionality for the observer.
 *  See StepperOperatorSplitObserverBase for details on the algorithm.
 */
template<class Scalar>
class StepperOperatorSplitObserverDefault
  : virtual public Tempus::StepperOperatorSplitObserverBase<Scalar>
{
public:

  /// Constructor
  StepperOperatorSplitObserverDefault(){}

  /// Destructor
  virtual ~StepperOperatorSplitObserverDefault(){}

  /// Observe OperatorSplit Stepper at end of takeStep.
  virtual void observe(
    Teuchos::RCP<const SolutionHistory<Scalar> > /* sh */,
    Teuchos::RCP<const StepperOperatorSplit<Scalar> > /* stepper */,
    const typename StepperOperatorSplitAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    switch(actLoc) {
      case StepperOperatorSplitAppAction<Scalar>::BEGIN_STEP:
      case StepperOperatorSplitAppAction<Scalar>::BEFORE_STEPPER:
      case StepperOperatorSplitAppAction<Scalar>::AFTER_STEPPER:
      case StepperOperatorSplitAppAction<Scalar>::END_STEP:
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

#endif // Tempus_StepperOperatorSplitObserverDefault_hpp
