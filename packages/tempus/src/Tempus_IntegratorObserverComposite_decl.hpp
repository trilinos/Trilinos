// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_IntegratorObserverComposite_decl_hpp
#define Tempus_IntegratorObserverComposite_decl_hpp

#include "Tempus_IntegratorObserver.hpp"
#include "Tempus_TimeStepControl.hpp"
#include <vector>

namespace Tempus {

/** \brief This observer 
 */
template<class Scalar>
class IntegratorObserverComposite
  : virtual public Tempus::IntegratorObserver<Scalar>
{
public:

  /// Constructor
  IntegratorObserverComposite(
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory,
    const Teuchos::RCP<TimeStepControl<Scalar> >& timeStepControl);

  /// Destructor
  virtual ~IntegratorObserverComposite();

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

    // add observer to the composite observer list
    void addObserver(const Teuchos::RCP<IntegratorObserver<Scalar> > &observer);
  //@}

protected:

  Teuchos::RCP<SolutionHistory<Scalar> > solutionHistory_;
  Teuchos::RCP<TimeStepControl<Scalar> > timeStepControl_;

private:

  std::vector<Teuchos::RCP<IntegratorObserver<Scalar > > > observers_;

};

} // namespace Tempus
#endif // Tempus_IntegratorObserverComposite_decl_hpp
