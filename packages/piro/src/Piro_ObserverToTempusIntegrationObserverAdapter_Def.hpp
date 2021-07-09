// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
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

  //Get solution
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
  Teuchos::RCP<Thyra::VectorBase<Scalar>> x = solutionHistory_->getCurrentState()->getX(); 
  Teuchos::RCP<DMVPV> X = Teuchos::rcp_dynamic_cast<DMVPV>(x);
  //IKT, 5/12/2020: getX() returns a vector containing the sensitivities for the 
  //case of sensitivity calculations for sensitivity integrator. 
  //Hence, we extract the first column, which is the solution.
  Teuchos::RCP<Thyra::VectorBase<Scalar>> solution = (sens_method_ == NONE) ? x : X->getNonconstMultiVector()->col(0);
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > solution_dxdp_mv = Teuchos::null; 
  if (sens_method_ == FORWARD) {
    const int num_param = X->getMultiVector()->domain()->dim()-1;
    const Teuchos::Range1D rng(1,num_param);
    solution_dxdp_mv = X->getNonconstMultiVector()->subView(rng);
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
  //Get solution_dot 
  Teuchos::RCP<const Thyra::VectorBase<Scalar>> xdot = solutionHistory_->getCurrentState()->getXDot(); 
  Teuchos::RCP<const DMVPV> XDot = Teuchos::rcp_dynamic_cast<const DMVPV>(xdot);

  const Scalar scalar_time = solutionHistory_->getCurrentState()->getTime();
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType StampScalar;
  const StampScalar time = Teuchos::ScalarTraits<Scalar>::real(scalar_time);

  //Get solution_dotdot 
  Teuchos::RCP<const Thyra::VectorBase<Scalar>> xdotdot = solutionHistory_->getCurrentState()->getXDotDot(); 
  Teuchos::RCP<const DMVPV> XDotDot = Teuchos::rcp_dynamic_cast<const DMVPV>(xdotdot);
  if (Teuchos::nonnull(xdot)) {
    //IKT, 5/11/2020: getXDot() returns a vector containing the sensitivities for 
    //the case of sensitivity calculations for sentivity integrator 
    //Hence, we extract the first column, which is the solution_dot.
    Teuchos::RCP<const Thyra::VectorBase<Scalar>> solution_dot 
          = (sens_method_ == NONE) ? xdot : XDot->getMultiVector()->col(0);
    if (supports_x_dotdot_) {
      //IKT, 5/11/2020: getXDotDot() returns a vector containing the sensitivities for 
      //the case of sensitivity calculations. 
      //Hence, we extract the first column, which is the solution_dotdot.
      Teuchos::RCP<const Thyra::VectorBase<Scalar>> solution_dotdot 
          = (sens_method_ == NONE) ? xdot : XDotDot->getMultiVector()->col(0);
      if (solution_dxdp_mv != Teuchos::null) {
        wrappedObserver_->observeSolution(*solution, *solution_dxdp_mv, *solution_dot, *solution_dotdot, time);
      }
      else {
        wrappedObserver_->observeSolution(*solution, *solution_dot, *solution_dotdot, time);
      }
    }
    else {
      if (solution_dxdp_mv != Teuchos::null) {
        wrappedObserver_->observeSolution(*solution, *solution_dxdp_mv, *solution_dot, time);
      }
      else {
        wrappedObserver_->observeSolution(*solution, *solution_dot, time);
      }
    }
  }
  else {
    if (solution_dxdp_mv != Teuchos::null) {
      wrappedObserver_->observeSolution(*solution, *solution_dxdp_mv, time);
    }
    else {
      wrappedObserver_->observeSolution(*solution, time);
    }
  }
  previous_dt_ = current_dt; 
}

