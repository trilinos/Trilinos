//@HEADER
// ***********************************************************************
//
//                     Rythmos Package
//                 Copyright (2006) Sandia Corporation
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

#include "Rythmos_GAASPInterface.hpp"
#include "GAdjointSolve.h"
#include "GErrorEstimate.h"
#include "Rythmos_GAASPGModel_ThyraModelEvaluator.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

namespace Rythmos {

// Static members

const std::string GAASPInterface::sTime_name_ = "sTime";
const double GAASPInterface::sTime_default_ = 0.0;
const std::string GAASPInterface::eTime_name_ = "eTime";
const double GAASPInterface::eTime_default_ = 10.0;
const std::string GAASPInterface::timeStep_name_ = "timeStep";
const double GAASPInterface::timeStep_default_ = 0.001;
const std::string GAASPInterface::method_name_ = "Time Integration Method Pair";
const std::string GAASPInterface::method_default_ = "DG1";
const std::string GAASPInterface::qtyOfInterest_name_ = "Quantity of Interest";
const std::string GAASPInterface::qtyOfInterest_default_ = "Average";
const std::string GAASPInterface::uTOL_name_ = "uTOL";
const double GAASPInterface::uTOL_default_ = 0.00001;
const std::string GAASPInterface::meshRefineMethod_name_ = "Mesh Refinement Method";
const std::string GAASPInterface::meshRefineMethod_default_ = "Weighted Adaptive";

GAASPInterface::GAASPInterface() {
  isInitialized_ = false;
  forwardSolveCompleted_ = false;
  adjointSolveCompleted_ = false;
  errorEstimateCompleted_ = false;
  refineMeshCompleted_ = false;
}

void GAASPInterface::initialize() {
  TEST_FOR_EXCEPTION(is_null(tModel_), std::logic_error,
      "Error, no model has been set, this must be done before anything can be run!"
      );

  // Create GAASP objects for reuse.
  forwardSolver_ = Teuchos::rcp(new GAASP::GForwardSolve);
  adjointSolver_ = Teuchos::rcp(new GAASP::GAdjointSolve);
  errorEstimator_ = Teuchos::rcp(new GAASP::GErrorEstimate);
  meshRefiner_ = Teuchos::rcp(new GAASP::GMeshRefine);
  isInitialized_ = true;
}

void GAASPInterface::forwardSolve() {
  if (!isInitialized_) {
    initialize();
  }
  if (!forwardSolveCompleted_) {
    GAASP::TSimControl fsim = {sTime_, eTime_, timeStep_, 0, method_};
    forwardSolver_->Setup(gModel_,fsim,initialCondition_);
  }
  if (refineMeshCompleted_) {
    forwardSolver_->updateMesh(newMesh_,nNodes_);
  }
  forwardSolver_->solve();
	forwardSolver_->printSolution();
  // Get data out of the solver to save for the adjoint and the error estimation
	mesh_ = forwardSolver_->getMesh();
	nSteps_ = forwardSolver_->getNSteps();
	fSoln_ = forwardSolver_->getSolution();
	fDer_ = forwardSolver_->getDerivatives();
  forwardSolveCompleted_ = true;
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  if (Teuchos::as<int>(verbLevel) != Teuchos::VERB_NONE) {
    *out << "GAASPInterface::forwardSolve nSteps_ = " << nSteps_+1 << std::endl;
  }
}

void GAASPInterface::adjointSolve() {
  TEST_FOR_EXCEPTION(!forwardSolveCompleted_, std::logic_error,
      "Error, a forward solve must be completed before an adjoint solve can be run!"
      );
  if (!adjointSolveCompleted_) {
    GAASP::TAdjointControl asim = {
      sTime_, 
      eTime_, 
      timeStep_, 
      0, 
      qtyOfInterest_, 
      method_, 
      1
    };
    adjointSolver_->Setup(gModel_, asim, fSoln_, mesh_);
  } else {
    adjointSolver_->updateData(mesh_, nSteps_, fSoln_);
  }
	adjointSolver_->solve();
	adjointSolver_->printSolution();
	aSoln_ = adjointSolver_->getSolution();
  adjointSolveCompleted_ = true;
}

Teuchos::RCP<const GAASPErrorEstimate> GAASPInterface::computeErrorEstimate() {
  TEST_FOR_EXCEPTION(!adjointSolveCompleted_, std::logic_error, 
      "Error, an adjoint solve must be completed before an error estimate can be calculated!"
      );
  if (!errorEstimateCompleted_) {
    errorEstimator_->Setup(
        fSoln_, 
        aSoln_, 
        fDer_, 
        mesh_, 
        initialCondition_.size(), 
        nSteps_, 
        gModel_, 
        method_
        );
  } else {
    errorEstimator_->updateData( 
        mesh_,
        nSteps_,
        fSoln_,
        fDer_,
        aSoln_
        );
  }
	errorEstimator_->compute();
	errorEstimator_->printErrorEstimates();
  // Copy data into GAASPErrorEstimate object
  errorEstimate_ = errorEstimator_->getErrorEstimate();
  intError_ = errorEstimator_->getIntervalContributions();
  Teuchos::RCP<GAASPErrorEstimate> ee = Teuchos::rcp(new GAASPErrorEstimate);
  ee->setErrorEstimate(errorEstimate_);
  ee->setIntervalErrorContributions(intError_);
  errorEstimateCompleted_ = true;
  return(ee);
}

void GAASPInterface::refineMesh() {
  TEST_FOR_EXCEPTION(!errorEstimateCompleted_, std::logic_error,
      "Error, the mesh cannot be refined until an error estimate has been completed!"
      );
  if (!refineMeshCompleted_) {
    meshRefiner_->Setup(
        mesh_,
        intError_,
        meshRefineMethod_,
        method_,
        uTOL_,
        dim_,
        nSteps_+1
        );
  } else {
    meshRefiner_->updateData(
        intError_,
        mesh_,
        nSteps_+1
        );
  }
  meshRefiner_->refine();
  newMesh_ = meshRefiner_->getNewMesh();
  nNodes_ = meshRefiner_->getNumberofNodes();
  refineMeshCompleted_ = true;
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  if (Teuchos::as<int>(verbLevel) != Teuchos::VERB_NONE) {
    *out << "GAASPInterface::refineMesh nNodes_ = " << nNodes_ << std::endl;
  }
}

void GAASPInterface::setThyraModelEvaluator(Teuchos::RCP<Thyra::ModelEvaluator<double> > tModel) {
  TEST_FOR_EXCEPTION(is_null(tModel), std::logic_error,
      "Error, a null Model Evaluator was passed in!");
  tModel_ = tModel;
  // Get the initial condition to the ModelEvaluator
  Teuchos::RCP<const Thyra::VectorBase<double> > x_init = tModel_->getNominalValues().get_x();

//  Teuchos::RCP<Teuchos::FancyOStream>
//    out = Teuchos::VerboseObjectBase::getDefaultOStream();
//  *out << "GAASPInterface::setThyraModelEvaluator:  x_init = " << std::endl;
//  x_init->describe(*out,Teuchos::VERB_EXTREME);

  dim_ = tModel_->get_x_space()->dim();
  initialCondition_.clear();
  for (int i=0; i<dim_ ; ++i) {
    // TODO:  Change this to Serial double * view
//    *out << "GAASPInterface::setThyraModelEvaluator:             x_init["<<i<<"] = " << Thyra::get_ele<double>(*x_init,i) << std::endl;
    initialCondition_.push_back(Thyra::get_ele<double>(*x_init,i));
//    *out << "GAASPInterface::setThyraModelEvaluator:  initialCondition_["<<i<<"] = " << initialCondition_[i] << std::endl;
  }

  gModel_ = gModel_ThyraModelEvaluator(tModel_);
}

std::string GAASPInterface::description() const
{
  std::string name = "Rythmos::GAASPInterface";
  return(name);
}

void GAASPInterface::describe(
  Teuchos::FancyOStream                &out
  ,const Teuchos::EVerbosityLevel      verbLevel
  ) const
{
  out << description() << "::describe" << std::endl;
}

void GAASPInterface::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList)
{
  TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParametersAndSetDefaults(*this->getValidParameters(),0);
  paramList_ = paramList;
  sTime_ = Teuchos::get<double>(*paramList_,sTime_name_);
  eTime_ = Teuchos::get<double>(*paramList_,eTime_name_);
  timeStep_ = Teuchos::get<double>(*paramList_,timeStep_name_);
  method_ = Teuchos::getIntegralValue<GAASP::TMethodFlag>(*paramList_,method_name_);
  qtyOfInterest_ = Teuchos::getIntegralValue<GAASP::TErrorFlag>(*paramList_,qtyOfInterest_name_);
  uTOL_ = Teuchos::get<double>(*paramList_,uTOL_name_);
  meshRefineMethod_ = 
    Teuchos::getIntegralValue<GAASP::TRefinementFlag>(*paramList_,meshRefineMethod_name_);
  Teuchos::readVerboseObjectSublist(&*paramList_,this);
}

Teuchos::RCP<Teuchos::ParameterList>
GAASPInterface::getNonconstParameterList()
{
  return(paramList_);
}

Teuchos::RCP<Teuchos::ParameterList> GAASPInterface::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_param_list = paramList_;
  paramList_ = Teuchos::null;
  return(temp_param_list);
}

Teuchos::RCP<const Teuchos::ParameterList> GAASPInterface::getValidParameters() const {
  static Teuchos::RCP<Teuchos::ParameterList> validPL;

  if (is_null(validPL)) {

    Teuchos::RCP<Teuchos::ParameterList>
      pl = Teuchos::parameterList();
    
    pl->set(sTime_name_, sTime_default_,
        "Set the start time for the simulation"
        );
    pl->set(eTime_name_, eTime_default_,
        "Set the end time for the simulation"
        );
    pl->set(timeStep_name_, timeStep_default_,
        "Set the fixed time step size for the simulation"
        );
    Teuchos::setStringToIntegralParameter<GAASP::TMethodFlag>(
        method_name_, method_default_, 
        "Set the pair of integration methods for the forward and adjoint solves",
        Teuchos::tuple<std::string>("DG0","DG1"),
        Teuchos::tuple<GAASP::TMethodFlag>(GAASP::Method_DG0,GAASP::Method_DG1),
        &*pl
        );
    Teuchos::setStringToIntegralParameter<GAASP::TErrorFlag>(
        qtyOfInterest_name_, qtyOfInterest_default_, 
        "Set the quantity of interest for the global error estimate",
        Teuchos::tuple<std::string>("End","Average"),
        Teuchos::tuple<GAASP::TErrorFlag>(GAASP::Error_End,GAASP::Error_Avg),
        &*pl
        );
    pl->set(uTOL_name_, uTOL_default_,
        "Set the global error tolerance for global error control"
        );
    Teuchos::setStringToIntegralParameter<GAASP::TRefinementFlag>(
        meshRefineMethod_name_, meshRefineMethod_default_, 
        "Set the mesh refinement method",
        Teuchos::tuple<std::string>(
          "By Contribution",
          "Probabilistic 1",
          "Probabilistic 2", 
          "Weighted Adaptive"
          ),
        Teuchos::tuple<GAASP::TRefinementFlag>(
          GAASP::Refine_ByContribution,
          GAASP::Refine_Probabilistic1,
          GAASP::Refine_Probabilistic2,
          GAASP::Refine_WeightedAdaptive
          ),
        &*pl
        );

    Teuchos::setupVerboseObjectSublist(&*pl);

    validPL = pl;

  }

  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"getValidParameters");
  if (Teuchos::as<int>(verbLevel) == Teuchos::VERB_HIGH) {
    *out << "Setting up valid parameterlist." << std::endl;
    validPL->print(*out);
  }

  return (validPL);
}

} // namespace Rythmos


