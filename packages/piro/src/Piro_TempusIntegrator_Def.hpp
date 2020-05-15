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
    basicIntegrator_ = Tempus::integratorBasic<Scalar>(pList, model);
    fwdSensIntegrator_ = Teuchos::null; 
    adjSensIntegrator_ = Teuchos::null; 
  }
  else if (sens_method == FORWARD) {
    //forward sensitivities
    basicIntegrator_ = Teuchos::null;
    fwdSensIntegrator_ = Tempus::integratorForwardSensitivity<Scalar>(pList, model);
    adjSensIntegrator_ = Teuchos::null; 
  }
  else if (sens_method == ADJOINT) {
    //adjoint sensitivities
    basicIntegrator_ = Teuchos::null;
    fwdSensIntegrator_ = Teuchos::null; 
    adjSensIntegrator_ = Tempus::integratorAdjointSensitivity<Scalar>(pList, model);
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
    basicIntegrator_->getObserver()->clearObservers();
  }
  //IKT: fwdSensIntegrator and adjSensIntegrator do not support composite 
  //observer so clearObservers() routine is not relevant.  If a 
  //composite observer is ever added to these integrators, we'd need 
  //to uncomment the code below 
  /*if (fwdSensIntegrator_ != Teuchos::null) {
    fwdSensIntegrator_->getObserver()->clearObservers(); 
  }
  if (adjSensIntegrator_ != Teuchos::null) {
    adjSensIntegrator_->getObserver()->clearObservers(); 
  }*/
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
  //IKT, FIXME: adjSensIntegrator_ has no routine setObsever(obs)
  //Look into.
  /*if (adjSensIntegrator_ != Teuchos::null) {
    adjSensIntegrator_->setObserver(obs); 
  }*/
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
  //IKT FIXME: adjSensIntegrator_ has no initialize() routine.
  //Look into.
  /*if (adjSensIntegrator_ != Teuchos::null) {
    adjSensIntegrator_->initialize(); 
  }*/
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
    return Teuchos::null; 
  }
}

template <typename Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> 
Piro::TempusIntegrator<Scalar>::getDxdotDp() const
{
  if (fwdSensIntegrator_ != Teuchos::null) {
    return fwdSensIntegrator_->getDxdotDp(); 
  }
  else {
    return Teuchos::null; 
  }
}

template <typename Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> 
Piro::TempusIntegrator<Scalar>::getDxdotdotDp() const
{
  if (fwdSensIntegrator_ != Teuchos::null) {
    return fwdSensIntegrator_->getDxdotdotDp(); 
  }
  else {
    return Teuchos::null; 
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
    return Teuchos::null; 
  }
}
