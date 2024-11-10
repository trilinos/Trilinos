// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Piro_TempusStepperFactory_hpp__
#define __Piro_TempusStepperFactory_hpp__

#include "Teuchos_RCP.hpp"

#include "Tempus_Stepper.hpp"

#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_NonlinearSolverBase.hpp"

namespace Piro {

template <typename Scalar>
class TempusStepperFactory {
public:

  virtual ~TempusStepperFactory() {}

  virtual Teuchos::RCP<Tempus::Stepper<Scalar> > buildStepper(
                        const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > & model,
                        const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > & solver,
                        const Teuchos::RCP<Teuchos::ParameterList> & paramList) = 0;
  
};

}

#endif
