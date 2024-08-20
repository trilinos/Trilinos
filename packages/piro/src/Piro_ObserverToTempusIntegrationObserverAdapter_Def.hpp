// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_ObserverToTempusIntegrationObserverAdapter.hpp"

#include "Teuchos_Ptr.hpp"

//#define DEBUG_OUTPUT

// Constructor
template <typename Scalar>
Piro::ObserverToTempusIntegrationObserverAdapter<Scalar>::ObserverToTempusIntegrationObserverAdapter(
    const Teuchos::RCP<const Tempus::SolutionHistory<Scalar> >& solutionHistory,
    const Teuchos::RCP<const Tempus::TimeStepControl<Scalar> >& timeStepControl,
    const Teuchos::RCP<Piro::ObserverBase<Scalar> > &wrappedObserver,
    const bool supports_x_dotdot, const bool abort_on_fail_at_min_dt, 
    const SENS_METHOD sens_method)
    : solutionHistory_(solutionHistory),
      timeStepControl_(timeStepControl),
      out_(Teuchos::VerboseObjectBase::getDefaultOStream()),
      wrappedObserver_(wrappedObserver),
      supports_x_dotdot_(supports_x_dotdot),
      abort_on_fail_at_min_dt_(abort_on_fail_at_min_dt), 
      sens_method_(sens_method)  
{
  previous_dt_ = 0.0;
}

template <typename Scalar>
Piro::ObserverToTempusIntegrationObserverAdapter<Scalar>::~ObserverToTempusIntegrationObserverAdapter()
{
  //Nothing to do
}

template <typename Scalar>
void
Piro::ObserverToTempusIntegrationObserverAdapter<Scalar>::
observeStartIntegrator(const Tempus::Integrator<Scalar>& integrator)
{
  // store off the solution history and time step control
  solutionHistory_ = integrator.getSolutionHistory();
  timeStepControl_ = integrator.getTimeStepControl();

  std::time_t begin = std::time(nullptr);
  const Teuchos::RCP<Teuchos::FancyOStream> out = integrator.getOStream();
  Teuchos::OSTab ostab(out,0,"ScreenOutput");
  *out << "\nTempus - IntegratorBasic\n"
       << std::asctime(std::localtime(&begin)) << "\n"
       << "  Stepper = " << integrator.getStepper()->description() << "\n"
       << "  Simulation Time Range  [" << integrator.getTimeStepControl()->getInitTime()
       << ", " << integrator.getTimeStepControl()->getFinalTime() << "]\n"
       << "  Simulation Index Range [" << integrator.getTimeStepControl()->getInitIndex()
       << ", " << integrator.getTimeStepControl()->getFinalIndex() << "]\n"
       << "============================================================================\n"
       << "  Step       Time         dt  Abs Error  Rel Error  Order  nFail  dCompTime"
       << std::endl;

  this->observeTimeStep();
}

template <typename Scalar>
void
Piro::ObserverToTempusIntegrationObserverAdapter<Scalar>::
observeStartTimeStep(const Tempus::Integrator<Scalar>& )
{
  //Nothing to do
}

template <typename Scalar>
void
Piro::ObserverToTempusIntegrationObserverAdapter<Scalar>::
observeNextTimeStep(const Tempus::Integrator<Scalar>& )
{
  //Nothing to do
}

template <typename Scalar>
void
Piro::ObserverToTempusIntegrationObserverAdapter<Scalar>::
observeBeforeTakeStep(const Tempus::Integrator<Scalar>& )
{
  //Nothing to do
}


template <typename Scalar>
void
Piro::ObserverToTempusIntegrationObserverAdapter<Scalar>::
observeAfterTakeStep(const Tempus::Integrator<Scalar>& )
{
  //Nothing to do
}


template<class Scalar>
void
Piro::ObserverToTempusIntegrationObserverAdapter<Scalar>::
observeAfterCheckTimeStep(const Tempus::Integrator<Scalar>& integrator)
{
  //Nothing to do
}


template<class Scalar>
void
Piro::ObserverToTempusIntegrationObserverAdapter<Scalar>::
observeEndTimeStep(const Tempus::Integrator<Scalar>& integrator)
{
  using Teuchos::RCP;
  RCP<Tempus::SolutionStateMetaData<Scalar> > csmd =
    integrator.getSolutionHistory()->getCurrentState()->getMetaData();

  if ((csmd->getOutputScreen() == true) or
      (csmd->getOutput() == true) or
      (csmd->getTime() == integrator.getTimeStepControl()->getFinalTime())) {

     const Scalar steppertime = integrator.getStepperTimer()->totalElapsedTime();
     // reset the stepper timer
     integrator.getStepperTimer()->reset();

     const Teuchos::RCP<Teuchos::FancyOStream> out = integrator.getOStream();
     Teuchos::OSTab ostab(out,0,"ScreenOutput");
     *out<<std::scientific<<std::setw( 6)<<std::setprecision(3)<<csmd->getIStep()
        <<std::setw(11)<<std::setprecision(3)<<csmd->getTime()
        <<std::setw(11)<<std::setprecision(3)<<csmd->getDt()
        <<std::setw(11)<<std::setprecision(3)<<csmd->getErrorAbs()
        <<std::setw(11)<<std::setprecision(3)<<csmd->getErrorRel()
        <<std::fixed     <<std::setw( 7)<<std::setprecision(1)<<csmd->getOrder()
        <<std::scientific<<std::setw( 7)<<std::setprecision(3)<<csmd->getNFailures()
        <<std::setw(11)<<std::setprecision(3)<<steppertime
        <<std::endl;
  }
  this->observeTimeStep();
}


template <typename Scalar>
void
Piro::ObserverToTempusIntegrationObserverAdapter<Scalar>::
observeEndIntegrator(const Tempus::Integrator<Scalar>& integrator)
{
  //this->observeTimeStep();

  std::string exitStatus;
  //const Scalar runtime = integrator.getIntegratorTimer()->totalElapsedTime();
  if (integrator.getSolutionHistory()->getCurrentState()->getSolutionStatus() ==
      Tempus::Status::FAILED or integrator.getStatus() == Tempus::Status::FAILED) {
    exitStatus = "Time integration FAILURE!";
  } else {
    exitStatus = "Time integration complete.";
  }
  std::time_t end = std::time(nullptr);
  const Scalar runtime = integrator.getIntegratorTimer()->totalElapsedTime();
  const Teuchos::RCP<Teuchos::FancyOStream> out = integrator.getOStream();
  Teuchos::OSTab ostab(out,0,"ScreenOutput");
  *out << "============================================================================\n"
       << "  Total runtime = " << runtime << " sec = "
       << runtime/60.0 << " min\n"
       << std::asctime(std::localtime(&end))
       << exitStatus << "\n"
       << std::endl;

}

template <typename Scalar>
void
Piro::ObserverToTempusIntegrationObserverAdapter<Scalar>::observeTimeStep()
{
  Scalar current_dt; 
  if (Teuchos::nonnull(solutionHistory_->getCurrentState())) {
    current_dt = solutionHistory_->getCurrentState()->getTimeStep();
  }
  else {
    current_dt = solutionHistory_->getCurrentState()->getTimeStep();
  }
  
  //Don't observe solution if step failed to converge
  if ((solutionHistory_->getCurrentState() != Teuchos::null) &&
     (solutionHistory_->getCurrentState()->getSolutionStatus() == Tempus::Status::FAILED)) {
    Scalar min_dt = timeStepControl_->getMinTimeStep(); 
    if ((previous_dt_ == current_dt) && (previous_dt_ == min_dt)) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, 
        "\n Error!  Time integration FAILURE!  Time integrator has failed using Minimum Time Step = " << min_dt << "\n" 
        << "and is unable to reduce time step further.\n");  
    }
    previous_dt_ = current_dt; 
    return;
  }

  //Get current state and its derivatives
  Teuchos::RCP<Thyra::VectorBase<Scalar>> x = solutionHistory_->getCurrentState()->getX(); 
  Teuchos::RCP<Thyra::VectorBase<Scalar>> xdot = solutionHistory_->getCurrentState()->getXDot(); 
  Teuchos::RCP<Thyra::VectorBase<Scalar>> xdotdot = solutionHistory_->getCurrentState()->getXDotDot(); 
  
  //Declare solution and its time derivatives
  Teuchos::RCP<Thyra::VectorBase<Scalar>> solution = Teuchos::null; 
  Teuchos::RCP<Thyra::VectorBase<Scalar>> solution_dot = Teuchos::null; 
  Teuchos::RCP<Thyra::VectorBase<Scalar>> solution_dotdot = Teuchos::null; 
  
  //Declare solution derivative dx/dp (for forward sens) 
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > solution_dxdp_mv = Teuchos::null; 
  
  if ((sens_method_ == NONE) || (sens_method_ == ADJOINT)) { //this is a hack for now for 
	                                                     //the adjoint case
  //if (sens_method_ == NONE) { 
    solution = x; 
    solution_dot = xdot; 
    solution_dotdot = xdotdot; 
  } 
  else if (sens_method_ == FORWARD) {
    //Get solution and its derivatives
    typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
    Teuchos::RCP<DMVPV> x_fwd = Teuchos::rcp_dynamic_cast<DMVPV>(x);
    Teuchos::RCP<DMVPV> xdot_fwd = Teuchos::rcp_dynamic_cast<DMVPV>(xdot);
    Teuchos::RCP<DMVPV> xdotdot_fwd = Teuchos::rcp_dynamic_cast<DMVPV>(xdotdot);
    //First vector in x_fwd is solution, 2nd vector is dx/dp
    solution = x_fwd->getNonconstMultiVector()->col(0);
    solution_dot = x_fwd->getNonconstMultiVector()->col(0);
    solution_dotdot = x_fwd->getNonconstMultiVector()->col(0);
    const int num_param = x_fwd->getMultiVector()->domain()->dim()-1;
    const Teuchos::Range1D rng(1,num_param);
    solution_dxdp_mv = x_fwd->getNonconstMultiVector()->subView(rng);
#ifdef DEBUG_OUTPUT
    for (int np = 0; np < num_param; np++) {
      *out_ << "\n*** Piro::ObserverToTempusIntegrationObserverAdapter dxdp" << np << " ***\n";
      Teuchos::RCP<const Thyra::VectorBase<Scalar>> solution_dxdp = solution_dxdp_mv->col(np);
      Teuchos::Range1D range;
      RTOpPack::ConstSubVectorView<Scalar> dxdpv;
      solution_dxdp->acquireDetachedView(range, &dxdpv);
      auto dxdpa = dxdpv.values();
      for (auto i = 0; i < dxdpa.size(); ++i) *out_ << dxdpa[i] << " ";
      *out_ << "\n*** Piro::ObserverToTempusIntegrationObserverAdapter dxdp" << np << " ***\n";
    }
#endif
  }
  //IKT FIXME: Commenting out for now since cast doesn't work...
 /* else if (sens_method_ == ADJOINT) {
    //Get solution and its derivatives
    //First vector in x_fwd is solution, 2nd vector is dg/dp
    //TODO?  put in observation of dg/dp
    typedef Thyra::DefaultProductVector<Scalar> DPV; 
    std::cout << "IKT x = " << x << "\n"; 
    Teuchos::RCP<DPV> x_adj = Teuchos::rcp_dynamic_cast<DPV>(x);
    std::cout << "IKT x = " << x_adj << "\n"; 
    Teuchos::RCP<DPV> xdot_adj = Teuchos::rcp_dynamic_cast<DPV>(xdot);
    Teuchos::RCP<DPV> xdotdot_fwd = Teuchos::rcp_dynamic_cast<DPV>(xdotdot);
    solution = x_adj->getNonconstVectorBlock(0);
    solution_dot = x_adj->getNonconstVectorBlock(0);
    solution_dotdot = x_adj->getNonconstVectorBlock(0);
  }*/

  const Scalar scalar_time = solutionHistory_->getCurrentState()->getTime();
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType StampScalar;
  const StampScalar time = Teuchos::ScalarTraits<Scalar>::real(scalar_time);

  //Observe solution
  //IKT, 7/21/2021: this may need to be extended to observe DgDp for 
  //transient adjoint sensitivities
  if (solution_dot != Teuchos::null) { //have solution_dot 
    if (supports_x_dotdot_) { //Transient problem suppors second derivs
      if (solution_dxdp_mv != Teuchos::null) {
        wrappedObserver_->observeSolution(*solution, *solution_dxdp_mv, *solution_dot, *solution_dotdot, time);
      }
      else {
        wrappedObserver_->observeSolution(*solution, *solution_dot, *solution_dotdot, time);
      }
    }
    else { //Trasient problem that does not support second derivs in ttime
      if (solution_dxdp_mv != Teuchos::null) {
        wrappedObserver_->observeSolution(*solution, *solution_dxdp_mv, *solution_dot, time);
      }
      else {
        wrappedObserver_->observeSolution(*solution, *solution_dot, time);
      }
    }
  }
  else { //no solution_dot 
    if (solution_dxdp_mv != Teuchos::null) {
      wrappedObserver_->observeSolution(*solution, *solution_dxdp_mv, time);
    }
    else {
      wrappedObserver_->observeSolution(*solution, time);
    }
  }
  previous_dt_ = current_dt; 
}

