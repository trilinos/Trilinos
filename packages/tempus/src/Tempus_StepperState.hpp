//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperState_hpp
#define Tempus_StepperState_hpp

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include <string>

#include "Tempus_config.hpp"

namespace Tempus {

template <class Scalar>
/** \brief StepperState is a simple class to hold state information about the
 * stepper.
 *
 * <b>Design Considerations</b>
 *   - The purpose of StepperState is to provide a means to store any stepper
 *     information that is required to restart from a checkpoint.
 *   - The StepperState should be held by the SolutionState so that it has
 *     all the information needed to restart the solution.
 *   - Many Steppers will simply need this base class, because
 *     they do not have any other additional state information.
 *   - StepperState can be inherited to expand the state information.
 *   - Examples of other information that could be included in derived
 *     StepperStates:
 *     - Metadata and SolutionHistory for a BDF method
 *   - The base class currently has the Stepper name so the Stepper can
 *     check if the StepperState is usable.
 */
class StepperState
  : public Teuchos::Describable,
    public Teuchos::VerboseObject<Tempus::StepperState<Scalar> >

{
 public:
  /// Constructor
  StepperState(std::string name = "Default") : stepperName_(name) {}

  /// Clone copy constructor
  virtual Teuchos::RCP<StepperState<Scalar> > clone() const
  {
    Teuchos::RCP<StepperState<Scalar> > ss_out =
        Teuchos::rcp(new StepperState<Scalar>(this->stepperName_));
    return ss_out;
  }

  /// This is a deep copy
  virtual void copy(const Teuchos::RCP<const StepperState<Scalar> >& ss)
  {
    stepperName_ = ss->stepperName_;
  }

  /// \name Overridden from Teuchos::Describable
  //@{
  virtual std::string description() const
  {
    return "Tempus::StepperState - '" + stepperName_ + "'";
  }

  virtual void describe(Teuchos::FancyOStream& out,
                        const Teuchos::EVerbosityLevel /* verbLevel */) const
  {
    auto l_out = Teuchos::fancyOStream(out.getOStream());
    Teuchos::OSTab ostab(*l_out, 2, this->description());
    l_out->setOutputToRootOnly(0);

    *l_out << "\n--- " << this->description() << " ---" << std::endl;
  }
  //@}

  std::string stepperName_;  ///< Name of the creating Stepper.
};
}  // namespace Tempus
#endif  // Tempus_StepperState_hpp
