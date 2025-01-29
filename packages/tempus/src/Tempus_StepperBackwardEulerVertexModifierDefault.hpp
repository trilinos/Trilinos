// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperBackwardEulerVertexModifierDefault_hpp
#define Tempus_StepperBackwardEulerVertexModifierDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperBackwardEulerVertexModifierBase.hpp"

// Applications can uncomment this include in their implementation,
// if they need access to the stepper methods.
//#include "Tempus_StepperBackwardEulerVertex.hpp"


namespace Tempus {

/** \brief Default modifier for StepperBackwardEulerVertex.
 *
 *  The default modifier provides no-op functionality for the modifier.
 *  See StepperBackwardEulerVertexModifierBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template<class Scalar>
class StepperBackwardEulerVertexModifierDefault
  : virtual public Tempus::StepperBackwardEulerVertexModifierBase<Scalar>
{
public:

  /// Constructor
  StepperBackwardEulerVertexModifierDefault(){}

  /// Destructor
  virtual ~StepperBackwardEulerVertexModifierDefault(){}

  /// Modify BackwardEulerVertex Stepper.
  virtual void modify(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Teuchos::RCP<StepperBackwardEulerVertex<Scalar> > /* stepper */,
    const typename StepperBackwardEulerVertexAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    switch(actLoc) {
      case StepperBackwardEulerVertexAppAction<Scalar>::BEGIN_STEP:
      case StepperBackwardEulerVertexAppAction<Scalar>::BEFORE_SOLVE:
      case StepperBackwardEulerVertexAppAction<Scalar>::AFTER_SOLVE:
      case StepperBackwardEulerVertexAppAction<Scalar>::END_STEP:
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

#endif // Tempus_StepperBackwardEulerVertexModifierDefault_hpp
