// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_IntegratorObserverBasic_decl_hpp
#define Tempus_IntegratorObserverBasic_decl_hpp

#include "Tempus_IntegratorObserver.hpp"

namespace Tempus {

/** \brief IntegratorObserverBasic class for time integrators.
 *  This basic class has simple no-op functions, as all basic
 *  functionality should be handled through other methods.
 */
template<class Scalar>
class IntegratorObserverBasic
  : virtual public Tempus::IntegratorObserver<Scalar>
{
public:

  /// Constructor
  IntegratorObserverBasic();

  /// Destructor
  virtual ~IntegratorObserverBasic();

  /// \name Basic IntegratorObserver methods
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
  //@}

};
} // namespace Tempus
#endif // Tempus_IntegratorObserverBasic_decl_hpp
