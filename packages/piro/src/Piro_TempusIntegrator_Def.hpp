// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_TempusIntegrator.hpp"

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Assert.hpp"

#include <string>
#include <stdexcept>
#include <iostream>

template <typename Scalar>
Piro::TempusIntegrator<Scalar>::TempusIntegrator(Teuchos::RCP< Teuchos::ParameterList > pList,
   const Teuchos::RCP< Thyra::ModelEvaluator< Scalar > > &model,
   const SENS_METHOD sens_method) :
   out_(Teuchos::VerboseObjectBase::getDefaultOStream())
{
  if (sens_method == NONE) {
    //no sensitivities
    basicIntegrator_ = Tempus::createIntegratorBasic<Scalar>(pList, model);
    fwdSensIntegrator_ = Teuchos::null;
    adjSensIntegrator_ = Teuchos::null;
  }
  else if (sens_method == FORWARD) {
    //forward sensitivities
    basicIntegrator_ = Teuchos::null;
    //Remove 'Tempus->Sensitivities->Response Function Index' parameter, if it appears
    //in the input file, as this does not make sense for forward sensitivities
    if (pList->isSublist("Sensitivities")){
      Teuchos::ParameterList& tempusSensPL = pList->sublist("Sensitivities", true);
      if (tempusSensPL.isParameter("Response Function Index")) {
        tempusSensPL.remove("Response Function Index");
      }
    }
    fwdSensIntegrator_ = Tempus::createIntegratorForwardSensitivity<Scalar>(pList, model);
    adjSensIntegrator_ = Teuchos::null;
  }
  else if (sens_method == ADJOINT) {
    //adjoint sensitivities
    basicIntegrator_ = Teuchos::null;
    fwdSensIntegrator_ = Teuchos::null;
    adjSensIntegrator_ = Tempus::createIntegratorAdjointSensitivity<Scalar>(pList, model);
  }
}

template <typename Scalar>
Piro::TempusIntegrator<Scalar>::TempusIntegrator(Teuchos::RCP< Teuchos::ParameterList > pList,
   const Teuchos::RCP< Thyra::ModelEvaluator< Scalar > > &model,
   const Teuchos::RCP< Thyra::ModelEvaluator< Scalar > > &adjoint_model,
   const SENS_METHOD sens_method) :
   out_(Teuchos::VerboseObjectBase::getDefaultOStream())
{
  if (sens_method == ADJOINT) {
    //adjoint sensitivities
    basicIntegrator_ = Teuchos::null;
    fwdSensIntegrator_ = Teuchos::null;
    adjSensIntegrator_ = Tempus::createIntegratorAdjointSensitivity<Scalar>(pList, model, adjoint_model);
  }
  else { //throw
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
                 "Error in Piro::TempusIntegrator: constructor taking adjoint model evaluator " <<
		 "is only valid for Adjoint sensitivities\n");
  }
}

template <typename Scalar>
Teuchos::RCP<Tempus::Stepper<Scalar>>
Piro::TempusIntegrator<Scalar>::getStepper() const
{
  Teuchos::RCP<Tempus::Stepper<Scalar>> stepper;
  if (basicIntegrator_ != Teuchos::null) {
    stepper = basicIntegrator_->getStepper();
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    stepper = fwdSensIntegrator_->getStepper();
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    stepper = adjSensIntegrator_->getStepper();
  }
  return stepper;
}

template <typename Scalar>
bool
Piro::TempusIntegrator<Scalar>::advanceTime(const Scalar time_final)
{
  bool out;
  if (basicIntegrator_ != Teuchos::null) {
    out = basicIntegrator_->advanceTime(time_final);
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    out = fwdSensIntegrator_->advanceTime(time_final);
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    out = adjSensIntegrator_->advanceTime(time_final);
  }
  return out;
}

template <typename Scalar>
Scalar
Piro::TempusIntegrator<Scalar>::getTime() const
{
  Scalar time;
  if (basicIntegrator_ != Teuchos::null) {
    time = basicIntegrator_->getTime();
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    time = fwdSensIntegrator_->getTime();
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    time = adjSensIntegrator_->getTime();
  }
  return time;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar>>
Piro::TempusIntegrator<Scalar>::getX() const
{
  Teuchos::RCP<const Thyra::VectorBase<Scalar>> x;
  if (basicIntegrator_ != Teuchos::null) {
    x = basicIntegrator_->getX();
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    x = fwdSensIntegrator_->getX();
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    x = adjSensIntegrator_->getX();
  }
  return x;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar>>
Piro::TempusIntegrator<Scalar>::getXDot() const
{
  Teuchos::RCP<const Thyra::VectorBase<Scalar>> xdot;
  if (basicIntegrator_ != Teuchos::null) {
    xdot = basicIntegrator_->getXDot();
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    xdot = fwdSensIntegrator_->getXDot();
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    xdot = adjSensIntegrator_->getXDot();
  }
  return xdot;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar>>
Piro::TempusIntegrator<Scalar>::getXDotDot() const
{
  Teuchos::RCP<const Thyra::VectorBase<Scalar>> xdotdot;
  if (basicIntegrator_ != Teuchos::null) {
    xdotdot = basicIntegrator_->getXDotDot();
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    xdotdot = fwdSensIntegrator_->getXDotDot();
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    xdotdot = adjSensIntegrator_->getXDotDot();
  }
  return xdotdot;
}

template <typename Scalar>
Teuchos::RCP<const Tempus::SolutionHistory<Scalar>>
Piro::TempusIntegrator<Scalar>::getSolutionHistory() const
{
  Teuchos::RCP<const Tempus::SolutionHistory<Scalar>> soln_history;
  if (basicIntegrator_ != Teuchos::null) {
    soln_history = basicIntegrator_->getSolutionHistory();
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    soln_history = fwdSensIntegrator_->getSolutionHistory();
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    soln_history = adjSensIntegrator_->getSolutionHistory();
  }
  return soln_history;
}

template <typename Scalar>
Teuchos::RCP<const Tempus::TimeStepControl<Scalar>>
Piro::TempusIntegrator<Scalar>::getTimeStepControl() const
{
  Teuchos::RCP<const Tempus::TimeStepControl<Scalar>> ts_control;
  if (basicIntegrator_ != Teuchos::null) {
    ts_control = basicIntegrator_->getTimeStepControl();
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    ts_control = fwdSensIntegrator_->getTimeStepControl();
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    ts_control = adjSensIntegrator_->getTimeStepControl();
  }
  return ts_control;
}

template <typename Scalar>
void
Piro::TempusIntegrator<Scalar>::clearObservers()
{
  if (basicIntegrator_ != Teuchos::null) {
    basicIntegrator_->setObserver();
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    fwdSensIntegrator_->setObserver();
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    adjSensIntegrator_->setObserver();
  }
}

template <typename Scalar>
void
Piro::TempusIntegrator<Scalar>::setObserver(Teuchos::RCP<Tempus::IntegratorObserver<Scalar>> obs)
{
  if (basicIntegrator_ != Teuchos::null) {
    basicIntegrator_->setObserver(obs);
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    fwdSensIntegrator_->setObserver(obs);
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    adjSensIntegrator_->setObserver(obs);
  }
}

template <typename Scalar>
void
Piro::TempusIntegrator<Scalar>::clearSolutionHistory()
{
  if (basicIntegrator_ != Teuchos::null) {
    basicIntegrator_->getNonConstSolutionHistory()->clear();
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    fwdSensIntegrator_->getNonConstSolutionHistory()->clear();
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    adjSensIntegrator_->getNonConstSolutionHistory()->clear();
  }
}

template <typename Scalar>
void
Piro::TempusIntegrator<Scalar>::initialize()
{
  if (basicIntegrator_ != Teuchos::null) {
    basicIntegrator_->initialize();
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    fwdSensIntegrator_->initialize();
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    adjSensIntegrator_->initialize();
  }
}

template <typename Scalar>
void
Piro::TempusIntegrator<Scalar>::initializeSolutionHistory(Scalar t0,
    Teuchos::RCP< const Thyra::VectorBase< Scalar > > x0,
    Teuchos::RCP< const Thyra::VectorBase< Scalar > > xdot0,
    Teuchos::RCP< const Thyra::VectorBase< Scalar > > xdotdot0,
    Teuchos::RCP< const Thyra::MultiVectorBase< Scalar > > DxDp0,
    Teuchos::RCP< const Thyra::MultiVectorBase< Scalar > > DxdotDp0,
    Teuchos::RCP< const Thyra::MultiVectorBase< Scalar > > DxdotdotDp0)
{
  if (basicIntegrator_ != Teuchos::null) {
    basicIntegrator_->initializeSolutionHistory(t0, x0, xdot0, xdotdot0);
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    fwdSensIntegrator_->initializeSolutionHistory(t0, x0, xdot0, xdotdot0, DxDp0, DxdotDp0, DxdotdotDp0);
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    adjSensIntegrator_->initializeSolutionHistory(t0, x0, xdot0, xdotdot0, DxDp0, DxdotDp0, DxdotdotDp0);
  }
}

template <typename Scalar>
Tempus::Status
Piro::TempusIntegrator<Scalar>::getStatus() const
{
  Tempus::Status status;
  if (basicIntegrator_ != Teuchos::null) {
    status = basicIntegrator_->getStatus();
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    status = fwdSensIntegrator_->getStatus();
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    status = adjSensIntegrator_->getStatus();
  }
  return status;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
Piro::TempusIntegrator<Scalar>::getDxDp() const
{
  if (fwdSensIntegrator_ != Teuchos::null) {
    return fwdSensIntegrator_->getDxDp();
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
                 "Error in Piro::TempusIntegrator: getDxDp() is not valid for the requested integrator type, " <<
		 "which is not of type Tempus::IntegratorForwardSensitivity!\n");
  }
}

template <typename Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
Piro::TempusIntegrator<Scalar>::getDxDotDp() const
{
  if (fwdSensIntegrator_ != Teuchos::null) {
    return fwdSensIntegrator_->getDxDotDp();
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
                 "Error in Piro::TempusIntegrator: getDxDotDp() is not valid for the requested integrator type, " <<
		 "which is not of type Tempus::IntegratorForwardSensitivity!\n");
  }
}

template <typename Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
Piro::TempusIntegrator<Scalar>::getDxDotDotDp() const
{
  if (fwdSensIntegrator_ != Teuchos::null) {
    return fwdSensIntegrator_->getDxDotDotDp();
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
                 "Error in Piro::TempusIntegrator: getDxDotDotDp() is not valid for the requested integrator type, " <<
		 "which is not of type Tempus::IntegratorForwardSensitivity!\n");
  }
}


template <typename Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
Piro::TempusIntegrator<Scalar>::getDgDp() const
{
  if (adjSensIntegrator_ != Teuchos::null) {
    return adjSensIntegrator_->getDgDp();
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
                 "Error in Piro::TempusIntegrator: getDgDp() is not valid for the requested integrator type, " <<
		 "which is not of type Tempus::IntegratorAdjointSensitivity!\n");
  }
}
