// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

  //Second constructor which takes in forward and adjoint ME - needed/valid for adjoint transient sensitivities only
  TempusIntegrator(Teuchos::RCP< Teuchos::ParameterList > pList, const Teuchos::RCP< Thyra::ModelEvaluator< Scalar > > &model,
		   const Teuchos::RCP< Thyra::ModelEvaluator< Scalar > > &adjoint_model,
                   const SENS_METHOD sens_method = NONE);

  Teuchos::RCP<Tempus::Stepper<Scalar>> getStepper() const;

  bool advanceTime(const Scalar time_final);

  Scalar getTime() const;

  Teuchos::RCP<const Thyra::VectorBase<Scalar>> getX() const;

  Teuchos::RCP<const Thyra::VectorBase<Scalar>> getXDot() const;

  Teuchos::RCP<const Thyra::VectorBase<Scalar>> getXDotDot() const;

  Teuchos::RCP<const Tempus::SolutionHistory<Scalar>> getSolutionHistory() const;

  Teuchos::RCP<const Tempus::TimeStepControl<Scalar>> getTimeStepControl() const;

  void clearObservers();

  void setObserver(Teuchos::RCP<Tempus::IntegratorObserver<Scalar>> obs = Teuchos::null);

  void clearSolutionHistory();

  void initialize();

  void initializeSolutionHistory(Scalar t0,
                                 Teuchos::RCP< const Thyra::VectorBase< Scalar > > x0,
                                 Teuchos::RCP< const Thyra::VectorBase< Scalar > > xdot0 = Teuchos::null,
                                 Teuchos::RCP< const Thyra::VectorBase< Scalar > > xdotdot0 = Teuchos::null,
                                 Teuchos::RCP< const Thyra::MultiVectorBase< Scalar > > DxDp0 = Teuchos::null,
                                 Teuchos::RCP< const Thyra::MultiVectorBase< Scalar > > DxDotDp0 = Teuchos::null,
                                 Teuchos::RCP< const Thyra::MultiVectorBase< Scalar > > DxdotDotDp0 = Teuchos::null);

  Tempus::Status getStatus() const;

  // The following 3 routines are only for forward sensitivities
  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> getDxDp() const;
  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> getDxDotDp() const;
  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> getDxDotDotDp() const;

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
