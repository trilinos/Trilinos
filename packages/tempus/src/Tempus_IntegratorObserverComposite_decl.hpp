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

  /// Default constructor
  IntegratorObserverComposite();

  /// Destructor
  virtual ~IntegratorObserverComposite();

  /// \name Override IntegratorObserver basic methods
  //@{
    /// Observe the beginning of the time integrator.
    virtual void observeStartIntegrator(const Integrator<Scalar>& integrator) override;

    /// Observe the beginning of the time step loop.
    virtual void observeStartTimeStep(const Integrator<Scalar>& integrator) override;

    /// Observe after the next time step size is selected.
    virtual void observeNextTimeStep(const Integrator<Scalar>& integrator) override;

    /// Observe before Stepper takes step.
    virtual void observeBeforeTakeStep(const Integrator<Scalar>& integrator) override;

    /// Observe after Stepper takes step.
    virtual void observeAfterTakeStep(const Integrator<Scalar>& integrator) override;

    /// Observe after accepting time step.
    virtual void observeAcceptedTimeStep(const Integrator<Scalar>& integrator) override;

    /// Observe the end of the time integrator.
    virtual void observeEndIntegrator(const Integrator<Scalar>& integrator) override;

    // add observer to the composite observer list
    void addObserver(const Teuchos::RCP<IntegratorObserver<Scalar> > &observer);
  //@}

private:

  std::vector<Teuchos::RCP<IntegratorObserver<Scalar > > > observers_;

};

} // namespace Tempus
#endif // Tempus_IntegratorObserverComposite_decl_hpp
