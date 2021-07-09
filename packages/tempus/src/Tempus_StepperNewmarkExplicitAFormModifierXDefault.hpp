// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperNewmarkExplicitAFormModifierX_hpp
#define Tempus_StepperNewmarkExplicitAFormModifierX_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperNewmarkExplicitAFormModifierXBase.hpp"


namespace Tempus {

/** \brief Default ModifierX for StepperNewmarkExplicitAForm.
 *
 *  The default provides no-op functionality for ModifierX.
 *  See StepperNewmarkExplicitAFormModifierXBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template<class Scalar>
class StepperNewmarkExplicitAFormModifierXDefault
  : virtual public Tempus::StepperNewmarkExplicitAFormModifierXBase<Scalar>
{
public:

  /// Constructor
  StepperNewmarkExplicitAFormModifierXDefault(){}

  /// Destructor
  virtual ~StepperNewmarkExplicitAFormModifierXDefault(){}

  /// Modify solution based on the MODIFIER_TYPE.
  virtual void modify(
    Teuchos::RCP<Thyra::VectorBase<Scalar> > /* x */,
    const Scalar /* time */, const Scalar /* dt */,
    const typename StepperNewmarkExplicitAFormModifierXBase<Scalar>::MODIFIER_TYPE modType)
  {
    switch(modType) {
      case StepperNewmarkExplicitAFormModifierXBase<Scalar>::X_BEGIN_STEP:
      case StepperNewmarkExplicitAFormModifierXBase<Scalar>::X_BEFORE_EXPLICIT_EVAL:
      case StepperNewmarkExplicitAFormModifierXBase<Scalar>::X_AFTER_EXPLICIT_EVAL:
      case StepperNewmarkExplicitAFormModifierXBase<Scalar>::X_END_STEP:
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

#endif // Tempus_StepperNewmarkExplicitAFormModifierX_hpp
