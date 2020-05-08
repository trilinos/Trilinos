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

#define DEBUG_OUTPUT

#include <string>
#include <stdexcept>
#include <iostream>

template <typename Scalar>
Piro::TempusIntegrator<Scalar>::TempusIntegrator(Teuchos::RCP< Teuchos::ParameterList > pList, 
   const Teuchos::RCP< Thyra::ModelEvaluator< Scalar > > &model,
   int sensitivities_requested) 
{
#ifdef DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  //IKT FIXME: make sensitivities_requested an enum type
  //IKT FIXME: remove duplicate integrators - should be able to just have 1 integrator rather than 3, and do 
  //casts to figure out which one we need. 
  if (sensitivities_requested == 0) {
    //no sensitivities
    basicIntegrator_ = Tempus::integratorBasic<Scalar>(pList, model);
    fwdSensIntegrator_ = Teuchos::null; 
    adjSensIntegrator_ = Teuchos::null; 
  }
  else if (sensitivities_requested == 1) {
    //forward sensitivities
    basicIntegrator_ = Teuchos::null;
    fwdSensIntegrator_ = Tempus::integratorForwardSensitivity<Scalar>(pList, model);
    adjSensIntegrator_ = Teuchos::null; 
  }
  else if (sensitivities_requested == 2) {
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
  //IKT FIXME - fill in! 
}

template <typename Scalar>
bool 
Piro::TempusIntegrator<Scalar>::advanceTime(const Scalar time_final)
{
  //IKT FIXME - fill in! 
}

template <typename Scalar>
Scalar 
Piro::TempusIntegrator<Scalar>::getTime() const
{
//IKT FIXME - fill in!
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar>> 
Piro::TempusIntegrator<Scalar>::getX() const
{
  //IKT FIXME - fill in! 
}
 
template <typename Scalar>
Teuchos::RCP<const Tempus::SolutionHistory<Scalar>> 
Piro::TempusIntegrator<Scalar>::getSolutionHistory() const
{
  //IKT FIXME - fill in!
}

template <typename Scalar>
Teuchos::RCP<const Tempus::TimeStepControl<Scalar>> 
Piro::TempusIntegrator<Scalar>::getTimeStepControl() const
{
  //IKT FIXME - fill in!
}

template <typename Scalar>
void 
Piro::TempusIntegrator<Scalar>::clearObservers()
{
  //IKT FIXME - fill in! 
}

template <typename Scalar>
void 
Piro::TempusIntegrator<Scalar>::setObserver(Teuchos::RCP<Tempus::IntegratorObserver<Scalar>> obs)
{
  //IKT FIXME - fill in!
}

template <typename Scalar>
void 
Piro::TempusIntegrator<Scalar>::initialize()
{
  //IKT FIXME - fill in!
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
  //IKT FIXME - fill in!
}

template <typename Scalar>
Tempus::Status 
Piro::TempusIntegrator<Scalar>::getStatus() const
{
  //IKT FIXME - fill in!
}
