//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_Integrator_hpp
#define Tempus_Integrator_hpp

#include "Tempus_config.hpp"
#include "Tempus_Types.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"

#include <string>

namespace Teuchos {
class Time;
}

namespace Tempus {
template <typename Scalar>
class Stepper;
template <typename Scalar>
class SolutionHistory;
template <typename Scalar>
class TimeStepControl;
}  // namespace Tempus

namespace Tempus {

/** \brief Thyra Base interface for time integrators.
 *  Time integrators are designed to advance the solution from an initial
 *  time, \f$t_0\f$, to a final time, \f$t_f\f$.
 *
 * <b>Design Considerations</b>
 *   - Integrators manage multiple time steps
 *   - Integrators have a single Stepper
 *   - Time-step ramping and startup are handled by TimeStepControl.
 *   - Solution output, e.g., solution plotting and check pointing,
 *     is coordinated in the Integrator.
 *   - Solution stability is handled in the timeStepControl, e.g., CFL
 *     constraint.
 *   - Error control over multiple time steps is handled in the Integrators,
 *     while error control over a single time step is handled in the Steppers.
 *   - Integrators will collect error control information from the Stepper
 *     and determine the next time step size and order.
 *   - Integrator maintains its own copy of the time history in the
 *     SolutionHistory, which may be just a single time step up to
 *     the entire solution history.
 *   - Integrators should compute the next time, and after accepting
 *     the time step advance the solution.  This allows a simple undo
 *     capability, if a solution is not acceptable.
 *
 * <b> CS Design Considerations</b>
 *   - Integrators will be fully constructed with input or default parameters.
 *   - All input parameters (i.e., ParameterList) can be set by public methods.
 *   - The Integrator ParameterList must be consistent.
 *     - The "set" methods which update parameters in the ParameterList
 *       must update the Integrator ParameterList.
 */
template <class Scalar>
class Integrator
  : virtual public Teuchos::Describable,
    virtual public Teuchos::VerboseObject<Tempus::Integrator<Scalar> > {
 public:
  /// \name Basic integrator methods
  //@{
  /// Advance the solution to time, and return true if successful.
  virtual bool advanceTime(const Scalar time_final) = 0;
  /// Get current time
  virtual Scalar getTime() const = 0;
  /// Get current index
  virtual int getIndex() const = 0;
  /// Get the Status
  virtual Tempus::Status getStatus() const = 0;
  /// Set the Status
  virtual void setStatus(const Tempus::Status st) = 0;
  /// Get the stepper
  virtual Teuchos::RCP<Stepper<Scalar> > getStepper() const = 0;
  /// Returns the SolutionHistory for this Integrator
  virtual Teuchos::RCP<const SolutionHistory<Scalar> > getSolutionHistory()
      const = 0;
  /// Returns the SolutionHistory for this Integrator
  virtual Teuchos::RCP<SolutionHistory<Scalar> >
  getNonConstSolutionHistory() = 0;
  /// Returns the TimeStepControl for this Integrator
  virtual Teuchos::RCP<const TimeStepControl<Scalar> > getTimeStepControl()
      const = 0;
  virtual Teuchos::RCP<TimeStepControl<Scalar> >
  getNonConstTimeStepControl() = 0;
  /// Returns the IntegratorTimer_ for this Integrator
  virtual Teuchos::RCP<Teuchos::Time> getIntegratorTimer() const = 0;
  virtual Teuchos::RCP<Teuchos::Time> getStepperTimer() const    = 0;
  //@}
};

}  // namespace Tempus
#endif  // Tempus_Integrator_hpp
