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

#include "Piro_ConfigDefs.hpp"

#ifdef HAVE_PIRO_NOX
#include "Piro_NOXSolver.hpp"
#include "Thyra_AdaptiveSolutionManager.hpp"
#include "Piro_LOCASolver.hpp"
#include "Piro_LOCAAdaptiveSolver.hpp"
#include "Piro_VelocityVerletSolver.hpp"
#include "Piro_TrapezoidRuleSolver.hpp"
#endif /* HAVE_PIRO_NOX */

#ifdef HAVE_PIRO_RYTHMOS
// point.
#include "Piro_RythmosSolver.hpp"
#endif /* HAVE_PIRO_RYTHMOS */

#ifdef HAVE_PIRO_TEMPUS
#include "Piro_TempusSolver.hpp"  
#endif /*HAVE_PIRO_TEMPUS */

#include "Teuchos_TestForException.hpp"

#include <string>
#include <stdexcept>

namespace Piro {

template <typename Scalar>
Teuchos::RCP<Thyra::ResponseOnlyModelEvaluatorBase<Scalar> > SolverFactory::createSolver(
    const Teuchos::RCP<Teuchos::ParameterList> &piroParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
    const Teuchos::RCP<Thyra::AdaptiveSolutionManager> &solMgr,
    const Teuchos::RCP<Piro::ObserverBase<Scalar> > &observer)
{
  Teuchos::RCP<Thyra::ResponseOnlyModelEvaluatorBase<Scalar> > result;

  const std::string &solverType = piroParams->get("Solver Type", "NOX");

#ifdef HAVE_PIRO_NOX
  if (solverType == "NOX") {
    result = Teuchos::rcp(new NOXSolver<Scalar>(piroParams, model, observer));
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
      result = observedLocaSolver(piroParams, model, solMgr, observer);
    else
      result = observedLocaSolver(piroParams, model, observer);
  } else
#endif /* HAVE_PIRO_NOX */
#ifdef HAVE_PIRO_RYTHMOS
  if (solverType == "Rythmos") {
    result = rythmosSolver<Scalar>(piroParams, model, observer);
  } else
#endif /* HAVE_PIRO_RYTHMOS */
#ifdef HAVE_PIRO_TEMPUS
  if (solverType == "Tempus") {
    result = tempusSolver<Scalar>(piroParams, model, observer);
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
    const Teuchos::RCP<Piro::ObserverBase<Scalar> > &observer)
{
  Teuchos::RCP<Thyra::ResponseOnlyModelEvaluatorBase<Scalar> > result;

  const std::string &solverType = piroParams->get("Solver Type", "NOX");

#ifdef HAVE_PIRO_NOX
  if (solverType == "NOX") {
    result = Teuchos::rcp(new NOXSolver<Scalar>(piroParams, model, observer));
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
    result = observedLocaSolver(piroParams, model, observer);
  } else
#endif /* HAVE_PIRO_NOX */
#ifdef HAVE_PIRO_RYTHMOS
  if (solverType == "Rythmos") {
    result = rythmosSolver<Scalar>(piroParams, model, observer);
  } else
#endif /* HAVE_PIRO_RYTHMOS */
#ifdef HAVE_PIRO_TEMPUS
  if (solverType == "Tempus") {
    result = tempusSolver<Scalar>(piroParams, model, observer);
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
