//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperBDF2ModifierX_hpp
#define Tempus_StepperBDF2ModifierX_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperBDF2ModifierXBase.hpp"

namespace Tempus {

/** \brief Default ModifierX for StepperBDF2.
 *
 *  The default provides no-op functionality for ModifierX.
 *  See StepperBDF2ModifierXBase for details on the algorithm.
 */
template <class Scalar>
class StepperBDF2ModifierXDefault
  : virtual public Tempus::StepperBDF2ModifierXBase<Scalar> {
 public:
  /// Constructor
  StepperBDF2ModifierXDefault() {}

  /// Destructor
  virtual ~StepperBDF2ModifierXDefault() {}

  /// Modify solution based on the MODIFIER_TYPE.
  virtual void modify(
      Teuchos::RCP<Thyra::VectorBase<Scalar> > /* x */, const Scalar /* time */,
      const Scalar /* dt */,
      const typename StepperBDF2ModifierXBase<Scalar>::MODIFIER_TYPE modType)
  {
    switch (modType) {
      case StepperBDF2ModifierXBase<Scalar>::X_BEGIN_STEP:
      case StepperBDF2ModifierXBase<Scalar>::X_BEFORE_SOLVE:
      case StepperBDF2ModifierXBase<Scalar>::X_AFTER_SOLVE:
      case StepperBDF2ModifierXBase<Scalar>::X_END_STEP: {
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

#endif  // Tempus_StepperBDF2ModifierX_hpp
