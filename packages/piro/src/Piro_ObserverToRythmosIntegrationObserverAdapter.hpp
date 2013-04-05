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

#ifndef PIRO_OBSERVATORTOINTEGRATIONOBSERVERSDAPTER_HPP
#define PIRO_OBSERVATORTOINTEGRATIONOBSERVERSDAPTER_HPP

#include "Rythmos_IntegrationObserverBase.hpp"

#include "Piro_ObserverBase.hpp"

#include "Teuchos_RCP.hpp"

namespace Piro {

template <typename Scalar>
class ObserverToRythmosIntegrationObserverAdapter : public Rythmos::IntegrationObserverBase<Scalar> {
public:
  explicit ObserverToRythmosIntegrationObserverAdapter(
      const Teuchos::RCP<ObserverBase<Scalar> > &wrappedObserver);

  // Overridden from Rythmos::IntegrationObserverBase

  virtual Teuchos::RCP<Rythmos::IntegrationObserverBase<Scalar> > cloneIntegrationObserver() const;

  virtual void resetIntegrationObserver(
      const Rythmos::TimeRange<Scalar> &integrationTimeDomain);

  virtual void observeStartTimeIntegration(
      const Rythmos::StepperBase<Scalar> &stepper);

  virtual void observeEndTimeIntegration(
      const Rythmos::StepperBase<Scalar> &stepper);

  virtual void observeStartTimeStep(
    const Rythmos::StepperBase<Scalar> &stepper,
    const Rythmos::StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter);

  virtual void observeCompletedTimeStep(
    const Rythmos::StepperBase<Scalar> &stepper,
    const Rythmos::StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter);

  virtual void observeFailedTimeStep(
    const Rythmos::StepperBase<Scalar> &stepper,
    const Rythmos::StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter);

private:
  void observeTimeStep(const Rythmos::StepperBase<Scalar> &stepper);

  void observeTimeStepStatus(
      const Rythmos::StepStatus<Scalar> &status,
      bool hasSensitivities);

  Teuchos::RCP<ObserverBase<Scalar> > wrappedObserver_;
};

} // namespace Piro

#include "Piro_ObserverToRythmosIntegrationObserverAdapter_Def.hpp"

#endif /* PIRO_OBSERVATORTOINTEGRATIONOBSERVERSDAPTER_HPP */
