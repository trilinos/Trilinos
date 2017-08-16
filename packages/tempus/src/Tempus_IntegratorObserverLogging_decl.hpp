// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_IntegratorObserverLogging_decl_hpp
#define Tempus_IntegratorObserverLogging_decl_hpp

#include "Tempus_IntegratorObserver.hpp"
#include "Tempus_TimeStepControl.hpp"

namespace Tempus {

/** \brief This observer logs calls to observer functions.
 *  This observer simply logs and counts the calls to each of the
 *  observer functions.  This is useful in monirtoring and debugging
 *  the time integration.
 */
template<class Scalar>
class IntegratorObserverLogging
  : virtual public Tempus::IntegratorObserver<Scalar>
{
public:

  /// Constructor
  IntegratorObserverLogging(
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory,
    const Teuchos::RCP<TimeStepControl<Scalar> >& timeStepControl);

  /// Destructor
  virtual ~IntegratorObserverLogging();

  /// \name Override IntegratorObserver basic methods
  //@{
    /// Observe the beginning of the time integrator.
    virtual void observeStartIntegrator();

    /// Observe the beginning of the time step loop.
    virtual void observeStartTimeStep();

    /// Observe after the next time step size is selected.
    virtual void observeNextTimeStep(Status & integratorStatus);

    /// Observe before Stepper takes step.
    virtual void observeBeforeTakeStep();

    /// Observe after Stepper takes step.
    virtual void observeAfterTakeStep();

    /// Observe after accepting time step.
    virtual void observeAcceptedTimeStep(Status & integratorStatus);

    /// Observe the end of the time integrator.
    virtual void observeEndIntegrator(const Status integratorStatus);

    virtual void setSolutionHistory(Teuchos::RCP<SolutionHistory<Scalar> > sh);

    virtual void setTimeStepControl(Teuchos::RCP<TimeStepControl<Scalar> > tsc);
  //@}

  void resetLogCounters();

  Teuchos::RCP<const std::map<std::string,int> > getCounters();

  Teuchos::RCP<const std::list<std::string> > getOrder();

  /** \name String names logged in map
      Use these strings to validate a call stack with this observer
  */
  //@{
    const std::string nameObserveStartIntegrator_;
    const std::string nameObserveStartTimeStep_;
    const std::string nameObserveNextTimeStep_;
    const std::string nameObserveBeforeTakeStep_;
    const std::string nameObserveAfterTakeStep_;
    const std::string nameObserveAcceptedTimeStep_;
    const std::string nameObserveEndIntegrator_;
  //@}

protected:

  Teuchos::RCP<SolutionHistory<Scalar> > solutionHistory_;
  Teuchos::RCP<TimeStepControl<Scalar> > timeStepControl_;

private:

  /** \brief Asserts next call on the stack is correct and removes from stack

      This is a const method so that it can be called from the
      derived IntegratorObserver methods that are const.
  */
  void logCall(const std::string call) const;

  Teuchos::RCP< std::map<std::string,int> > counters_;
  Teuchos::RCP< std::list<std::string> > order_;

};

} // namespace Tempus
#endif // Tempus_IntegratorObserverLogging_decl_hpp
