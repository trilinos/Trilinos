// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_ConfigDefs.hpp"

#ifdef HAVE_PIRO_NOX
#include "Piro_NOXSolver.hpp"
#include "Thyra_AdaptiveSolutionManager.hpp"
#include "Piro_LOCASolver.hpp"
#include "Piro_LOCAAdaptiveSolver.hpp"
#include "Piro_VelocityVerletSolver.hpp"
#include "Piro_TrapezoidRuleSolver.hpp"
#endif /* HAVE_PIRO_NOX */

#ifdef HAVE_PIRO_TEMPUS
#include "Piro_TempusSolver.hpp"  
#endif /*HAVE_PIRO_TEMPUS */

#include "Teuchos_TestForException.hpp"

#include <string>
#include <stdexcept>

namespace Piro {

template <typename Scalar>
Teuchos::RCP<Thyra::ResponseOnlyModelEvaluatorBase<Scalar> > SolverFactory::createSolverAdaptive(
    const Teuchos::RCP<Teuchos::ParameterList> &piroParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &adjointModel,
    const Teuchos::RCP<Thyra::AdaptiveSolutionManager> &solMgr,
    const Teuchos::RCP<Piro::ObserverBase<Scalar> > &observer)
{
  Teuchos::RCP<Thyra::ResponseOnlyModelEvaluatorBase<Scalar> > result;

  const std::string &solverType = piroParams->get("Solver Type", "NOX");

#ifdef HAVE_PIRO_NOX
  if (solverType == "NOX") {
    result = Teuchos::rcp(new NOXSolver<Scalar>(piroParams, model, adjointModel, observer));
  } else
  if (solverType == "Velocity Verlet") {
    result = Teuchos::rcp(new VelocityVerletSolver<Scalar>(
         piroParams, model, solMgr, observer));
  } else
  if (solverType == "Trapezoid Rule") {
    result = Teuchos::rcp(new TrapezoidRuleSolver<Scalar>(
        piroParams, model, solMgr, observer));
  } else
  if (solverType == "LOCA") {
    if(Teuchos::nonnull(solMgr))
      result = observedLocaSolver(piroParams, model, adjointModel, solMgr, observer);
    else
      result = observedLocaSolver(piroParams, model, adjointModel, observer);
  } else
#endif /* HAVE_PIRO_NOX */
#ifdef HAVE_PIRO_TEMPUS
  if (solverType == "Tempus") {
    result = tempusSolver<Scalar>(piroParams, model, adjointModel, observer);
  } else
#endif /* HAVE_PIRO_TEMPUS */
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true,
        std::runtime_error,
        "Error: Unknown Piro Solver Type: " << solverType);
  }

  return result;
}

/*
The below is DEPRECATED!
Please do not use
*/

template <typename Scalar>
Teuchos::RCP<Thyra::ResponseOnlyModelEvaluatorBase<Scalar> > SolverFactory::createSolver(
    const Teuchos::RCP<Teuchos::ParameterList> &piroParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &adjointModel,
    const Teuchos::RCP<Piro::ObserverBase<Scalar> > &observer)
{
  Teuchos::RCP<Thyra::ResponseOnlyModelEvaluatorBase<Scalar> > result;

  const std::string &solverType = piroParams->get("Solver Type", "NOX");

#ifdef HAVE_PIRO_NOX
  if (solverType == "NOX") {
    result = Teuchos::rcp(new NOXSolver<Scalar>(piroParams, model, adjointModel, observer));
  } else
  if (solverType == "Velocity Verlet") {
    result = Teuchos::rcp(new VelocityVerletSolver<Scalar>( 
         piroParams, model, Teuchos::null, observer));
  } else
  if (solverType == "Trapezoid Rule") {
    result = Teuchos::rcp(new TrapezoidRuleSolver<Scalar>(
        piroParams, model, Teuchos::null, observer));
  } else
  if (solverType == "LOCA") {
    result = observedLocaSolver(piroParams, model, adjointModel, observer);
  } else
#endif /* HAVE_PIRO_NOX */
#ifdef HAVE_PIRO_TEMPUS
  if (solverType == "Tempus") {
    result = tempusSolver<Scalar>(piroParams, model, adjointModel, observer);
  } else
#endif /* HAVE_PIRO_TEMPUS */
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true,
        std::runtime_error,
        "Error: Unknown Piro Solver Type: " << solverType);
  }

  return result;
}

} // namespace Piro
