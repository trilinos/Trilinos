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
  else { //throw
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
                 "Error in Piro::TempusIntegrator: invalid sensitivity method '" <<
                 sens_method << "'!\n");
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
  if (basicIntegrator_ != Teuchos::null) {
    return basicIntegrator_->getStepper();
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    return fwdSensIntegrator_->getStepper();
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    return adjSensIntegrator_->getStepper();
  }
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
               "Error in Piro::TempusIntegrator: no integrator stored!\n");
}

template <typename Scalar>
bool
Piro::TempusIntegrator<Scalar>::advanceTime(const Scalar time_final)
{
  if (basicIntegrator_ != Teuchos::null) {
    return basicIntegrator_->advanceTime(time_final);
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    return fwdSensIntegrator_->advanceTime(time_final);
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    return adjSensIntegrator_->advanceTime(time_final);
  }
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
               "Error in Piro::TempusIntegrator: no integrator stored!\n");
}

template <typename Scalar>
Scalar
Piro::TempusIntegrator<Scalar>::getTime() const
{
  if (basicIntegrator_ != Teuchos::null) {
    return basicIntegrator_->getTime();
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    return fwdSensIntegrator_->getTime();
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    return adjSensIntegrator_->getTime();
  }
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
               "Error in Piro::TempusIntegrator: no integrator stored!\n");
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar>>
Piro::TempusIntegrator<Scalar>::getX() const
{
  if (basicIntegrator_ != Teuchos::null) {
    return basicIntegrator_->getX();
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    return fwdSensIntegrator_->getX();
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    return adjSensIntegrator_->getX();
  }
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
               "Error in Piro::TempusIntegrator: no integrator stored!\n");
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar>>
Piro::TempusIntegrator<Scalar>::getXDot() const
{
  if (basicIntegrator_ != Teuchos::null) {
    return basicIntegrator_->getXDot();
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    return fwdSensIntegrator_->getXDot();
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    return adjSensIntegrator_->getXDot();
  }
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
               "Error in Piro::TempusIntegrator: no integrator stored!\n");
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar>>
Piro::TempusIntegrator<Scalar>::getXDotDot() const
{
  if (basicIntegrator_ != Teuchos::null) {
    return basicIntegrator_->getXDotDot();
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    return fwdSensIntegrator_->getXDotDot();
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    return adjSensIntegrator_->getXDotDot();
  }
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
               "Error in Piro::TempusIntegrator: no integrator stored!\n");
}

template <typename Scalar>
Teuchos::RCP<const Tempus::SolutionHistory<Scalar>>
Piro::TempusIntegrator<Scalar>::getSolutionHistory() const
{
  if (basicIntegrator_ != Teuchos::null) {
    return basicIntegrator_->getSolutionHistory();
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    return fwdSensIntegrator_->getSolutionHistory();
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    return adjSensIntegrator_->getSolutionHistory();
  }
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
               "Error in Piro::TempusIntegrator: no integrator stored!\n");
}

template <typename Scalar>
Teuchos::RCP<const Tempus::TimeStepControl<Scalar>>
Piro::TempusIntegrator<Scalar>::getTimeStepControl() const
{
  if (basicIntegrator_ != Teuchos::null) {
    return basicIntegrator_->getTimeStepControl();
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    return fwdSensIntegrator_->getTimeStepControl();
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    return adjSensIntegrator_->getTimeStepControl();
  }
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
               "Error in Piro::TempusIntegrator: no integrator stored!\n");
}

template <typename Scalar>
void
Piro::TempusIntegrator<Scalar>::clearObservers()
{
  if (basicIntegrator_ != Teuchos::null) {
    basicIntegrator_->setObserver();
    return;
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    fwdSensIntegrator_->setObserver();
    return;
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    adjSensIntegrator_->setObserver();
    return;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
               "Error in Piro::TempusIntegrator: no integrator stored!\n");
}

template <typename Scalar>
void
Piro::TempusIntegrator<Scalar>::setObserver(Teuchos::RCP<Tempus::IntegratorObserver<Scalar>> obs)
{
  if (basicIntegrator_ != Teuchos::null) {
    basicIntegrator_->setObserver(obs);
    return;
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    fwdSensIntegrator_->setObserver(obs);
    return;
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    adjSensIntegrator_->setObserver(obs);
    return;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
               "Error in Piro::TempusIntegrator: no integrator stored!\n");
}

template <typename Scalar>
void
Piro::TempusIntegrator<Scalar>::clearSolutionHistory()
{
  if (basicIntegrator_ != Teuchos::null) {
    basicIntegrator_->getNonConstSolutionHistory()->clear();
    return;
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    fwdSensIntegrator_->getNonConstSolutionHistory()->clear();
    return;
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    adjSensIntegrator_->getNonConstSolutionHistory()->clear();
    return;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
               "Error in Piro::TempusIntegrator: no integrator stored!\n");
}

template <typename Scalar>
void
Piro::TempusIntegrator<Scalar>::initialize()
{
  if (basicIntegrator_ != Teuchos::null) {
    basicIntegrator_->initialize();
    return;
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    fwdSensIntegrator_->initialize();
    return;
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    adjSensIntegrator_->initialize();
    return;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
               "Error in Piro::TempusIntegrator: no integrator stored!\n");
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
  if (basicIntegrator_ != Teuchos::null) {
    return basicIntegrator_->getStatus();
  }
  if (fwdSensIntegrator_ != Teuchos::null) {
    return fwdSensIntegrator_->getStatus();
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    return adjSensIntegrator_->getStatus();
  }
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
               "Error in Piro::TempusIntegrator: no integrator stored!\n");
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


template <typename Scalar>
Teuchos::RCP<const Tempus::Integrator<Scalar> >
Piro::TempusIntegrator<Scalar>:: getIntegrator () const
{
  Teuchos::RCP<const Tempus::Integrator<Scalar> > out;
  if (basicIntegrator_ != Teuchos::null) {
    out = basicIntegrator_;
  } else if (fwdSensIntegrator_ != Teuchos::null) {
    out = fwdSensIntegrator_;
  } else if (adjSensIntegrator_ != Teuchos::null) {
    out = adjSensIntegrator_;
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                 "Error in Piro::TempusIntegrator: no integrator stored!\n");
  }
  return out;
}
