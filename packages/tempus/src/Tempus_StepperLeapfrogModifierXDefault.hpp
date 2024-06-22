//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperLeapfrogModifierX_hpp
#define Tempus_StepperLeapfrogModifierX_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperLeapfrogModifierXBase.hpp"

namespace Tempus {

/** \brief Default ModifierX for StepperLeapfrog.
 *
 *  The default provides no-op functionality for ModifierX.
 *  See StepperLeapfrogModifierXBase for details on the algorithm.
 */
template <class Scalar>
class StepperLeapfrogModifierXDefault
  : virtual public Tempus::StepperLeapfrogModifierXBase<Scalar> {
 public:
  /// Constructor
  StepperLeapfrogModifierXDefault() {}

  /// Destructor
  virtual ~StepperLeapfrogModifierXDefault() {}

  /// Modify solution based on the MODIFIER_TYPE.
  virtual void modify(
      Teuchos::RCP<Thyra::VectorBase<Scalar> > /* x */, const Scalar /* time */,
      const Scalar /* dt */,
      const typename StepperLeapfrogModifierXBase<Scalar>::MODIFIER_TYPE
          modType)
  {
    switch (modType) {
      case StepperLeapfrogModifierXBase<Scalar>::X_BEGIN_STEP:
      case StepperLeapfrogModifierXBase<Scalar>::X_BEFORE_X_UPDATE:
      case StepperLeapfrogModifierXBase<Scalar>::X_BEFORE_EXPLICIT_EVAL:
      case StepperLeapfrogModifierXBase<Scalar>::X_BEFORE_XDOT_UPDATE:
      case StepperLeapfrogModifierXBase<Scalar>::X_END_STEP: {
        // No-op.
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                   "Error - unknown modifier type.\n");
    }
  }
};

}  // namespace Tempus

#endif  // Tempus_StepperLeapfrogModifierX_hpp
