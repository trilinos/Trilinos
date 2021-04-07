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

#ifndef PIRO_TRANSIENTINTEGRATOR_H
#define PIRO_TRANSIENTINTEGRATOR_H

#include "Piro_ConfigDefs.hpp"

#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_IntegratorForwardSensitivity.hpp"
#include "Tempus_IntegratorAdjointSensitivity.hpp"
#include "Piro_Helpers.hpp"

#include <map>
#include <string>

namespace Piro {

/** \brief Thyra-based Model Evaluator for Tempus solves using Tempus
 * */
template <typename Scalar>
class TempusIntegrator
{
public:
  /** \name Constructors/initializers */
  //@{

  /** \brief . */
  TempusIntegrator(Teuchos::RCP< Teuchos::ParameterList > pList, const Teuchos::RCP< Thyra::ModelEvaluator< Scalar > > &model,
                   const SENS_METHOD sens_method = NONE);

  Teuchos::RCP<Tempus::Stepper<Scalar>> getStepper() const;

  bool advanceTime(const Scalar time_final);

  Scalar getTime() const;

  Teuchos::RCP<const Thyra::VectorBase<Scalar>> getX() const;

  Teuchos::RCP<const Tempus::SolutionHistory<Scalar>> getSolutionHistory() const;

  Teuchos::RCP<const Tempus::TimeStepControl<Scalar>> getTimeStepControl() const;

  void clearObservers();

  void setObserver(Teuchos::RCP<Tempus::IntegratorObserver<Scalar>> obs = Teuchos::null);

  void initialize();

  void initializeSolutionHistory(Scalar t0,
                                 Teuchos::RCP< const Thyra::VectorBase< Scalar > > x0,
                                 Teuchos::RCP< const Thyra::VectorBase< Scalar > > xdot0 = Teuchos::null,
                                 Teuchos::RCP< const Thyra::VectorBase< Scalar > > xdotdot0 = Teuchos::null,
                                 Teuchos::RCP< const Thyra::MultiVectorBase< Scalar > > DxDp0 = Teuchos::null,
                                 Teuchos::RCP< const Thyra::MultiVectorBase< Scalar > > DxdotDp0 = Teuchos::null,
                                 Teuchos::RCP< const Thyra::MultiVectorBase< Scalar > > DxdotdotDp0 = Teuchos::null);

  Tempus::Status getStatus() const;

  // The following 3 routines are only for forward sensitivities
  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> getDxDp() const;
  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> getDxdotDp() const;
  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> getDxdotdotDp() const;

  //The following routine is only for adjoint sensitivities
  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> getDgDp() const;

private:

  Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > basicIntegrator_;
  Teuchos::RCP<Tempus::IntegratorForwardSensitivity<Scalar> > fwdSensIntegrator_;
  Teuchos::RCP<Tempus::IntegratorAdjointSensitivity<Scalar> > adjSensIntegrator_;
  Teuchos::RCP<Teuchos::FancyOStream> out_;

};

}

#include "Piro_TempusIntegrator_Def.hpp"
#endif
