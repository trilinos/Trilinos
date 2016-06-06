#ifndef TEMPUS_STEPPERSTATE_HPP
#define TEMPUS_STEPPERSTATE_HPP

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include <string>


namespace Tempus {

template<class Scalar>
/** \brief StepperState is a simple class to hold state information about the stepper.
 *
 * <b>Design Considerations</b>
 *   - The purpose of StepperState is to provide a means to store any stepper
 *     information that is required to restart from a checkpoint.
 *   - The StepperState should be held by the SolutionState so that it has
 *     all the information needed to restart the solution.
 *   - Many time integrators will simply need this base class, because
 *     they do not have any other additional state information.
 *   - StepperState can be inherited to expand the state information.
 *   - Examples of other information that could be included in derived
 *     StepperStates:
 *     - Theta value for a theta method
 *     - Metadata and SolutionHistory for a BDF method
 *   - The base class currently has the Stepper name so the Stepper can
 *     check if the StepperState is usable.
 *   - The StepperState also keeps track of the the Stepper status, i.e.,
 *     PASSED or FAILED.
 */
class StepperState :
  public Teuchos::Describable,
  public Teuchos::VerboseObject<Tempus::StepperState<Scalar> >

{
public:
  /// Constructor
  StepperState(std::string name_, Status stepperStatus_ = WORKING)
    : stepperName(name_), stepperStatus(stepperStatus_){}

  /// This is a deep copy
  virtual Teuchos::RCP<StepperState<Scalar> > clone() const
  {
     Teuchos::RCP<StepperState<Scalar> > ss_out =
       Teuchos::rcp(new StepperState<Scalar> (this->stepperName,
                                              this->stepperStatus));
     return ss_out;
  }

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual std::string description() const
    {
      std::string name = "Tempus::StepperState";
      return(name);
    }
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const
    {
      out << description() << "::describe" << std::endl
          << "  stepperName   = " << stepperName << std::endl
          << "  stepperStatus = " << toString(stepperStatus) << std::endl;
    }
  //@}

  std::string stepperName;    ///< Name of the creating Stepper.

  /** The stepperStatus is used to indicate whether the Stepper has PASSED or
   FAILED.  WORKING is used for prior and during the Stepper.
   */
  Status      stepperStatus;

};
} // namespace Tempus
#endif // TEMPUS_STEPPERSTATE_HPP
