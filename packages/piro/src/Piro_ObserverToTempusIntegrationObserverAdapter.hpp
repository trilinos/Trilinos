// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_OBSERVERTOTEMPUSINTEGRATIONOBSERVERADAPTER_HPP
#define PIRO_OBSERVERTOTEMPUSINTEGRATIONOBSERVERADAPTER_HPP

#include "Tempus_IntegratorObserverBasic.hpp"

#include "Piro_Helpers.hpp" 
#include "Piro_ObserverBase.hpp"

#include "Teuchos_RCP.hpp"

namespace Piro {

template <typename Scalar>
class ObserverToTempusIntegrationObserverAdapter : public Tempus::IntegratorObserverBasic<Scalar> {


public:


  // Constructor
  ObserverToTempusIntegrationObserverAdapter(
    const Teuchos::RCP<const Tempus::SolutionHistory<Scalar> >& solutionHistory,
    const Teuchos::RCP<const Tempus::TimeStepControl<Scalar> >& timeStepControl,
    const Teuchos::RCP<Piro::ObserverBase<Scalar> > &wrappedObserver,
    const bool supports_x_dotdot = false, 
    const bool abort_on_fail_at_min_dt = false, 
    const SENS_METHOD sens_method = NONE); 

  // Overridden from Tempus::IntegratorObserver

  //@{
  /// Destructor

  virtual ~ObserverToTempusIntegrationObserverAdapter() = default;

  /// Observe the beginning of the time integrator.
  virtual void observeStartIntegrator(const Tempus::Integrator<Scalar>& integrator) override;

  /// Observe the end of the time step loop.
  virtual void observeEndTimeStep(const Tempus::Integrator<Scalar>& integrator) override;
  //@}

private:

  void observeTimeStep();
  Teuchos::RCP<const Tempus::SolutionHistory<Scalar> > solutionHistory_;
  Teuchos::RCP<const Tempus::TimeStepControl<Scalar> > timeStepControl_;
  Teuchos::RCP<Teuchos::FancyOStream> out_;
  bool hasSensitivities_;
  Teuchos::RCP<ObserverBase<Scalar> > wrappedObserver_;
  bool supports_x_dotdot_;
  Scalar previous_dt_; 
  bool abort_on_fail_at_min_dt_;
  
  SENS_METHOD sens_method_;
};

} // namespace Piro

#include "Piro_ObserverToTempusIntegrationObserverAdapter_Def.hpp"

#endif /* PIRO_OBSERVERTOTEMPUSINTEGRATIONOBSERVERADAPTER_HPP */
