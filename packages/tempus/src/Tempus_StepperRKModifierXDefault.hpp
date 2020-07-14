// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperRKModifierX_hpp
#define Tempus_StepperRKModifierX_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperRKModifierXBase.hpp"


namespace Tempus {

/** \brief Default ModifierX for StepperRK.
 *
 *  The default provides no-op functionality for ModifierX.
 *  See StepperRKModifierXBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template<class Scalar>
class StepperRKModifierXDefault
  : virtual public Tempus::StepperRKModifierXBase<Scalar>
{
public:

  /// Constructor
  StepperRKModifierXDefault(){}

  /// Destructor
  virtual ~StepperRKModifierXDefault(){}

  /// Modify solution based on the MODIFIER_TYPE.
  virtual void modify(
    Teuchos::RCP<Thyra::VectorBase<Scalar> > /* x */,
    const Scalar /* time */, const Scalar /* dt */,
    const int /* stageNumber */,
    const typename StepperRKModifierXBase<Scalar>::MODIFIER_TYPE modType)
  {
    switch(modType) {
      case StepperRKModifierXBase<Scalar>::X_BEGIN_STEP:
      case StepperRKModifierXBase<Scalar>::X_BEGIN_STAGE:
      case StepperRKModifierXBase<Scalar>::X_BEFORE_SOLVE:
      case StepperRKModifierXBase<Scalar>::X_AFTER_SOLVE:
      case StepperRKModifierXBase<Scalar>::X_BEFORE_EXPLICIT_EVAL:
      case StepperRKModifierXBase<Scalar>::X_END_STAGE:
      case StepperRKModifierXBase<Scalar>::X_END_STEP:
      {
        // No-op.
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - unknown modifier type.\n");
    }
  }

};

} // namespace Tempus

#endif // Tempus_StepperRKModifierX_hpp
