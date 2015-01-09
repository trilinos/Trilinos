//@HEADER
// ***********************************************************************
//
//                     Rythmos Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef RYTHMOS_COMPOSITE_INTEGRATION_OBSERVER_HPP
#define RYTHMOS_COMPOSITE_INTEGRATION_OBSERVER_HPP


#include "Rythmos_IntegrationObserverBase.hpp"
#include "Teuchos_as.hpp"


namespace Rythmos {


/** \brief Standard composite observer subclass.
 *
 * ToDo: Finish Documentation
 */
template<class Scalar>
class CompositeIntegrationObserver
  : public IntegrationObserverBase<Scalar>
{
public:

  /** \name Constructors/Initializers/Accessors */
  //@{

  /** \brief . */
  CompositeIntegrationObserver();

  /** \brief . */
  void addObserver(
    const RCP<IntegrationObserverBase<Scalar> > &observer
    );

  // ToDo: Add functions to add observers

  //@}

  /** \name Overridden from IntegrationObserverBase */
  //@{

  /** \brief . */
  virtual RCP<IntegrationObserverBase<Scalar> >
  cloneIntegrationObserver() const;

  /** \brief . */
  virtual void resetIntegrationObserver(
    const TimeRange<Scalar> &integrationTimeDomain
    );

  void observeStartTimeIntegration(const StepperBase<Scalar> &stepper);

  void observeEndTimeIntegration(const StepperBase<Scalar> &stepper);

  void observeStartTimeStep(
    const StepperBase<Scalar> &stepper,
    const StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter
    );

  /** \brief . */
  virtual void observeCompletedTimeStep(
    const StepperBase<Scalar> &stepper,
    const StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter
    );

  /** \brief . */
  virtual void observeFailedTimeStep(
    const StepperBase<Scalar> &stepper,
    const StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter
    );

  //@}

private:

  Array<RCP<IntegrationObserverBase<Scalar> > > observers_;

};


/** \brief Non-member constructor.
 *
 * \relates CompositeIntegrationObserver
 */
template<class Scalar>
RCP<CompositeIntegrationObserver<Scalar> > createCompositeIntegrationObserver()
{
  RCP<CompositeIntegrationObserver<Scalar> >
    frsco(new CompositeIntegrationObserver<Scalar>());
  return frsco;
}


//
// Implementations
//


// Constructors/Initializers/Accessors


template<class Scalar>
CompositeIntegrationObserver<Scalar>::CompositeIntegrationObserver()
{}


template<class Scalar>
void CompositeIntegrationObserver<Scalar>::addObserver(
  const RCP<IntegrationObserverBase<Scalar> > &observer
  )
{
  observers_.push_back(observer);
}


// Overridden from IntegrationObserverBase


template<class Scalar>
RCP<IntegrationObserverBase<Scalar> >
CompositeIntegrationObserver<Scalar>::cloneIntegrationObserver() const
{
  using Teuchos::as;
  RCP<CompositeIntegrationObserver<Scalar> >
    compositeObserver = createCompositeIntegrationObserver<Scalar>();
  for (int i = 0; i < as<int>(observers_.size()); ++i ) {
    compositeObserver->addObserver(observers_[i]->cloneIntegrationObserver());
  }
  return compositeObserver;
}


template<class Scalar>
void CompositeIntegrationObserver<Scalar>::resetIntegrationObserver(
  const TimeRange<Scalar> &integrationTimeDomain
  )
{
  using Teuchos::as;
  for (int i = 0; i < as<int>(observers_.size()); ++i ) {
    observers_[i]->resetIntegrationObserver(integrationTimeDomain);
  }
}

template<class Scalar>
void CompositeIntegrationObserver<Scalar>::observeStartTimeIntegration(
  const StepperBase<Scalar> &stepper
  )
{
  using Teuchos::as;

  const RCP<FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();

  for (int i = 0; i < as<int>(observers_.size()); ++i ) {
    RCP<IntegrationObserverBase<Scalar> > observer = observers_[i];
    observer->setOStream(out);
    observer->setVerbLevel(verbLevel);
    observer->observeStartTimeIntegration(stepper);
  }
}

template<class Scalar>
void CompositeIntegrationObserver<Scalar>::observeEndTimeIntegration(
  const StepperBase<Scalar> &stepper
  )
{
  using Teuchos::as;

  const RCP<FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();

  for (int i = 0; i < as<int>(observers_.size()); ++i ) {
    RCP<IntegrationObserverBase<Scalar> > observer = observers_[i];
    observer->setOStream(out);
    observer->setVerbLevel(verbLevel);
    observer->observeEndTimeIntegration(stepper);
  }
}

template<class Scalar>
void CompositeIntegrationObserver<Scalar>::observeStartTimeStep(
  const StepperBase<Scalar> &stepper,
  const StepControlInfo<Scalar> &stepCtrlInfo,
  const int timeStepIter
  )
{
  using Teuchos::as;

  const RCP<FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();

  for (int i = 0; i < as<int>(observers_.size()); ++i ) {
    RCP<IntegrationObserverBase<Scalar> > observer = observers_[i];
    observer->setOStream(out);
    observer->setVerbLevel(verbLevel);
    observer->observeStartTimeStep(stepper,stepCtrlInfo,timeStepIter);
  }
}

template<class Scalar>
void CompositeIntegrationObserver<Scalar>::observeCompletedTimeStep(
  const StepperBase<Scalar> &stepper,
  const StepControlInfo<Scalar> &stepCtrlInfo,
  const int timeStepIter
  )
{
  using Teuchos::as;

  const RCP<FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();

  for (int i = 0; i < as<int>(observers_.size()); ++i ) {
    RCP<IntegrationObserverBase<Scalar> > observer = observers_[i];
    observer->setOStream(out);
    observer->setVerbLevel(verbLevel);
    observer->observeCompletedTimeStep(stepper,stepCtrlInfo,timeStepIter);
  }
}


template<class Scalar>
void CompositeIntegrationObserver<Scalar>::observeFailedTimeStep(
  const StepperBase<Scalar> &stepper,
  const StepControlInfo<Scalar> &stepCtrlInfo,
  const int timeStepIter
  )
{
  using Teuchos::as;

  const RCP<FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();

  for (int i = 0; i < as<int>(observers_.size()); ++i ) {
    RCP<IntegrationObserverBase<Scalar> > observer = observers_[i];
    observer->setOStream(out);
    observer->setVerbLevel(verbLevel);
    observer->observeFailedTimeStep(stepper,stepCtrlInfo,timeStepIter);
  }
}


} // namespace Rythmos


#endif //RYTHMOS_COMPOSITE_INTEGRATOR_OBSERVER_HPP
