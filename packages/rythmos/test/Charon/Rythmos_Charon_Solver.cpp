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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
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
#include "../SinCos/SinCosModel.hpp"
#include "Rythmos_SimpleIntegrationControlStrategy.hpp"


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
  return Rythmos::timeStepNonlinearSolver<double>();
}
Teuchos::RCP<Rythmos::SolverAcceptingStepperBase<double> >
RythmosCharonSolver::buildRythmosTimeStepper(
  const Teuchos::RCP<Teuchos::ParameterList> &rythmosStepperSelectionPL,
  Teuchos::FancyOStream &out
  )
{
  return Rythmos::backwardEulerStepper<double>();
}
Teuchos::RCP<RythmosCharon::CharonIntegrationControlAndObserver>
RythmosCharonSolver::buildIntegrationControlAndObserverStrategy(
  const Teuchos::RCP<Teuchos::ParameterList> &charonStepControlAndObservationSettingsPL,
  const Teuchos::RCP<const Thyra::ModelEvaluator<double> > &epetraThyraModel,
  Teuchos::FancyOStream &out
  )
{
  return Teuchos::rcp(new RythmosCharon::CharonIntegrationControlAndObserver());
}
Thyra::ModelEvaluatorBase::InArgs<double>
RythmosCharonSolver::getStateInitialCondition()
{
  return epetraThyraModel_->getNominalValues();
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
void RythmosCharonSolver::setup() 
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using Teuchos::sublist;

  epetraThyraModel_ = Rythmos::sinCosModel(true);

  const std::string TimeStepNonlinearSolverSelection_name = "Time Step Nonlinear Solver Selection";
  const std::string RythmosStepperSelection_name = "Rythmos Stepper Selection";
  const std::string CharonStepControlAndObservationSettings_name = "Charon Step Control and Observation Settings";

  const RCP<Teuchos::FancyOStream> out = this->getOStream();

  RCP<ParameterList> rythmosSettingsPL = Teuchos::parameterList();
  const RCP<Thyra::NonlinearSolverBase<double> > nonlinearSolver =
    buildTimeStepNonlinearSolver(
      sublist(rythmosSettingsPL, TimeStepNonlinearSolverSelection_name),
      *out
      );

  rythmosStepper_ = buildRythmosTimeStepper(
    sublist(rythmosSettingsPL, RythmosStepperSelection_name),
    *out
    );
  rythmosStepper_->setModel(epetraThyraModel_);
  rythmosStepper_->setSolver(nonlinearSolver);

  ryhmosCharonIntegrationControlAndObserver_ =
    buildIntegrationControlAndObserverStrategy(
      sublist(rythmosSettingsPL, CharonStepControlAndObservationSettings_name),
      epetraThyraModel_, 
      *out
      );

  RCP<Rythmos::DefaultIntegrator<double> > integrator = 
    Rythmos::defaultIntegrator<double>(
      ryhmosCharonIntegrationControlAndObserver_, // control
      ryhmosCharonIntegrationControlAndObserver_ // observer
      );

  {
    Teuchos::RCP<Rythmos::SimpleIntegrationControlStrategy<double> > intCont = 
      Rythmos::simpleIntegrationControlStrategy<double>();
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->set("Take Variable Steps",false);
    pl->set("Fixed dt",0.1);
    intCont->setParameterList(pl);
    integrator->setIntegrationControlStrategy(intCont);
  }
  rythmosStateIntegrator_ = integrator;

  rythmosStepper_->setInitialCondition(
    getStateInitialCondition()
    );

}
bool RythmosCharonSolver::solve() 
{ 
  const double finalTime = 1.0;

  rythmosStateIntegrator_->setStepper(rythmosStepper_, finalTime);

  rythmosStateIntegrator_->getFwdPoints(
      Teuchos::tuple<double>(finalTime),
      0, // x_vec
      0, // xdot_vec
      0  // accuracy_vec
      );

  return true;
}


