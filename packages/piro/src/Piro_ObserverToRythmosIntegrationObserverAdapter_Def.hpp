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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#include "Piro_ObserverToRythmosIntegrationObserverAdapter.hpp"

#include "Rythmos_ForwardSensitivityStepper.hpp"
#include "Rythmos_extractStateAndSens.hpp"

#include "Teuchos_Ptr.hpp"

template <typename Scalar>
Piro::ObserverToRythmosIntegrationObserverAdapter<Scalar>::ObserverToRythmosIntegrationObserverAdapter(
    const Teuchos::RCP<ObserverBase<Scalar> > &wrappedObserver) :
  wrappedObserver_(wrappedObserver)
{}

template <typename Scalar>
Teuchos::RCP<Rythmos::IntegrationObserverBase<Scalar> >
Piro::ObserverToRythmosIntegrationObserverAdapter<Scalar>::cloneIntegrationObserver() const
{
  return Teuchos::rcp(new ObserverToRythmosIntegrationObserverAdapter(*this));
}

template <typename Scalar>
void
Piro::ObserverToRythmosIntegrationObserverAdapter<Scalar>::resetIntegrationObserver(
    const Rythmos::TimeRange<Scalar> &/*integrationTimeDomain*/)
{
  // Nothing to do
}

template <typename Scalar>
void
Piro::ObserverToRythmosIntegrationObserverAdapter<Scalar>::observeStartTimeIntegration(
    const Rythmos::StepperBase<Scalar> &stepper)
{
  this->observeTimeStep(stepper);
}

template <typename Scalar>
void
Piro::ObserverToRythmosIntegrationObserverAdapter<Scalar>::observeEndTimeIntegration(
    const Rythmos::StepperBase<Scalar> &/*stepper*/)
{
  // Nothing to do
}

template <typename Scalar>
void
Piro::ObserverToRythmosIntegrationObserverAdapter<Scalar>::observeStartTimeStep(
    const Rythmos::StepperBase<Scalar> &/*stepper*/,
    const Rythmos::StepControlInfo<Scalar> &/*stepCtrlInfo*/,
    const int /*timeStepIter*/)
{
  // Nothing to do
}

template <typename Scalar>
void
Piro::ObserverToRythmosIntegrationObserverAdapter<Scalar>::observeCompletedTimeStep(
    const Rythmos::StepperBase<Scalar> &stepper,
    const Rythmos::StepControlInfo<Scalar> &/*stepCtrlInfo*/,
    const int /*timeStepIter*/)
{
  this->observeTimeStep(stepper);
}

template <typename Scalar>
void
Piro::ObserverToRythmosIntegrationObserverAdapter<Scalar>::observeFailedTimeStep(
    const Rythmos::StepperBase<Scalar> &/*stepper*/,
    const Rythmos::StepControlInfo<Scalar> &/*stepCtrlInfo*/,
    const int /*timeStepIter*/)
{
  // Nothing to do
}

template <typename Scalar>
void
Piro::ObserverToRythmosIntegrationObserverAdapter<Scalar>::observeTimeStep(
    const Rythmos::StepperBase<Scalar> &stepper)
{
  typedef Rythmos::ForwardSensitivityStepper<Scalar> FSS;
  const bool hasSensitivities = Teuchos::nonnull(
      Teuchos::ptr_dynamic_cast<const FSS>(Teuchos::ptr(&stepper)));

  this->observeTimeStepStatus(stepper.getStepStatus(), hasSensitivities);
}

template <typename Scalar>
void
Piro::ObserverToRythmosIntegrationObserverAdapter<Scalar>::observeTimeStepStatus(
    const Rythmos::StepStatus<Scalar> &status,
    bool hasSensitivities)
{
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > solution = status.solution;
  solution.assert_not_null();

  Teuchos::RCP<const Thyra::VectorBase<Scalar> > solution_dot = status.solutionDot;

  if (hasSensitivities) {
    Rythmos::extractState(solution, solution_dot, &solution, &solution_dot);
  }

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType StampScalar;
  const StampScalar time = Teuchos::ScalarTraits<Scalar>::real(status.time);

  if (Teuchos::nonnull(solution_dot)) {
    wrappedObserver_->observeSolution(*solution, *solution_dot, time);
  } else {
    wrappedObserver_->observeSolution(*solution, time);
  }
}
