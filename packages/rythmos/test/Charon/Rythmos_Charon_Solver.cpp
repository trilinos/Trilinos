//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2009) Sandia Corporation
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

// Charon includes:
#include "Rythmos_Charon_Solver.hpp"
#include "Rythmos_Charon_IntegrationControlAndObserver.hpp"
#include "Rythmos_Charon_ImplicitBDFStepperStepControl.hpp"
#include "Rythmos_Charon_ImplicitBDFStepperErrWtVecCalc.hpp"

// General Trilinos includes:
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_TimeMonitor.hpp"

// Rythmos includes:
#include "Rythmos_ImplicitBDFStepper.hpp"
#include "Rythmos_BackwardEulerStepper.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "Rythmos_DefaultIntegrator.hpp"
#include "Rythmos_StepperHelpers.hpp"


Teuchos::RCP<const Teuchos::ParameterList> RythmosCharonSolver::getStaticValidParameters()
{
  return Teuchos::null;
}
Teuchos::RCP<Thyra::NonlinearSolverBase<double> >
RythmosCharonSolver::buildTimeStepNonlinearSolver(
  const Teuchos::RCP<Teuchos::ParameterList> &timeStepNonlinearSolverSublist,
  Teuchos::FancyOStream &out
  )
{
  return Teuchos::null;
}
Teuchos::RCP<Rythmos::SolverAcceptingStepperBase<double> >
RythmosCharonSolver::buildRythmosTimeStepper(
  const Teuchos::RCP<Teuchos::ParameterList> &rythmosStepperSelectionPL,
  Teuchos::FancyOStream &out
  )
{
  return Teuchos::null;
}
Teuchos::RCP<RythmosCharon::RythmosCharonIntegrationControlAndObserver>
RythmosCharonSolver::buildIntegrationControlAndObserverStrategy(
  const Teuchos::RCP<Teuchos::ParameterList> &charonStepControlAndObservationSettingsPL,
  const Teuchos::RCP<const Thyra::EpetraModelEvaluator> &epetraThyraModel,
  Teuchos::FancyOStream &out
  )
{
  return Teuchos::null;
}
Thyra::ModelEvaluatorBase::InArgs<double>
RythmosCharonSolver::getStateInitialCondition()
{
  Thyra::ModelEvaluatorBase::InArgs<double> inArgs;
  return inArgs;
}
void RythmosCharonSolver::updateCharonState(
  const Thyra::VectorBase<double> &x_dot,
  const Thyra::VectorBase<double> &x,
  const double currentTime,
  const double timeStep,
  const bool guaranteeUpToDateAuxiliaryData,
  Teuchos::FancyOStream &out
  )
{ }
void RythmosCharonSolver::writeOutput(
  const double currentTime,
  const double timeStep,
  const int outIter,
  Teuchos::FancyOStream &out
  )
{ }
RythmosCharonSolver::RythmosCharonSolver() { }
void RythmosCharonSolver::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& paramList
  )
{ }
Teuchos::RCP<Teuchos::ParameterList> RythmosCharonSolver::getNonconstParameterList()
{
  return Teuchos::null;
}
Teuchos::RCP<Teuchos::ParameterList> RythmosCharonSolver::unsetParameterList()
{
  return Teuchos::null;
}
Teuchos::RCP<const Teuchos::ParameterList> RythmosCharonSolver::getParameterList() const
{
  return Teuchos::null;
}
Teuchos::RCP<const Teuchos::ParameterList> RythmosCharonSolver::getValidParameters() const
{
  return Teuchos::null;
}
void RythmosCharonSolver::setup() { }
bool RythmosCharonSolver::solve() 
{ 
  return false;
}


