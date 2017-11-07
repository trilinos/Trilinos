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


// Constructor
template <typename Scalar>
Piro::ObserverToTempusIntegrationObserverAdapter<Scalar>::ObserverToTempusIntegrationObserverAdapter(
    const Teuchos::RCP<const Tempus::SolutionHistory<Scalar> >& solutionHistory,
    const Teuchos::RCP<const Tempus::TimeStepControl<Scalar> >& timeStepControl,
    const Teuchos::RCP<Piro::ObserverBase<Scalar> > &wrappedObserver, 
    const bool supports_x_dotdot)
    : solutionHistory_(solutionHistory),
      timeStepControl_(timeStepControl),
      out_(Teuchos::VerboseObjectBase::getDefaultOStream()),
      wrappedObserver_(wrappedObserver),
      supports_x_dotdot_(supports_x_dotdot) 
{
  //Currently, sensitivities are not supported in Tempus.
  hasSensitivities_ = false;
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
  this->observeTimeStep();
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


template <typename Scalar>
void
Piro::ObserverToTempusIntegrationObserverAdapter<Scalar>::
observeAcceptedTimeStep(const Tempus::Integrator<Scalar>& )
{
  //Nothing to do
}


template <typename Scalar>
void
Piro::ObserverToTempusIntegrationObserverAdapter<Scalar>::
observeEndIntegrator(const Tempus::Integrator<Scalar>& )
{
  this->observeTimeStep();
}

template <typename Scalar>
void
Piro::ObserverToTempusIntegrationObserverAdapter<Scalar>::observeTimeStep()
{
  //Get solution
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > solution = solutionHistory_->getCurrentState()->getX();
  solution.assert_not_null();
  //Get solution_dot
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > solution_dot = solutionHistory_->getCurrentState()->getXDot();

  const Scalar scalar_time = solutionHistory_->getCurrentState()->getTime();
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType StampScalar;
  const StampScalar time = Teuchos::ScalarTraits<Scalar>::real(scalar_time);

  Teuchos::RCP<const Thyra::VectorBase<Scalar> > solution_dotdot = solutionHistory_->getCurrentState()->getXDotDot();
  if (Teuchos::nonnull(solution_dot))
  {
    if (supports_x_dotdot_) {

      wrappedObserver_->observeSolution(*solution, *solution_dot, *solution_dotdot, time);
    }
   else {
      wrappedObserver_->observeSolution(*solution, *solution_dot, time);
   }
  }
  else {
    wrappedObserver_->observeSolution(*solution, time);
  }
}

