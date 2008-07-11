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

namespace Rythmos {

GAASPInterface::GAASPInterface() {
  forwardSolveCompleted_ = false;
  adjointSolveCompleted_ = false;
  // Move all of these to the parameter list
  sTime_ = 0.0;
  eTime_ = 10.0;
  timeStep_ = 0.001;
  method_ = GAASP::Method_DG1;
  qtyOfInterest_ = GAASP::Error_Avg; // Convert from Rythmos::ERROR_QUANTITY_OF_INTEREST to GAASP::TErrorFlag
}

void GAASPInterface::forwardSolve() {
  TEST_FOR_EXCEPTION(is_null(tModel_), std::logic_error, 
      "Error, no model has been set, this must be done before a forward solve can be run!");
  GAASP::TSimControl fsim = {sTime_, eTime_, timeStep_, 0, method_};
  GAASP::GForwardSolve forwardSolver(gModel_,fsim,initialCondition_);
  forwardSolver.solve();
  // Get data out of the solver to save for the adjoint and the error estimation
	mesh_ = forwardSolver.getMesh();
	forwardSolver.printSolution();
	nSteps_ = forwardSolver.getNSteps();
	fSoln_ = forwardSolver.getSolution();
	fDer_ = forwardSolver.getDerivatives();
  forwardSolveCompleted_ = true;
}

void GAASPInterface::adjointSolve() {
  TEST_FOR_EXCEPTION(!forwardSolveCompleted_, std::logic_error,
      "Error, a forward solve must be completed before an adjoint solve can be run!");
  GAASP::TAdjointControl asim = {
    sTime_, 
    eTime_, 
    timeStep_, 
    0, 
    qtyOfInterest_, 
    method_, 
    1
  };
  GAASP::GAdjointSolve adjointSolver(gModel_, asim, fSoln_, mesh_);
	adjointSolver.solve();
	aSoln_ = adjointSolver.getSolution();
	adjointSolver.printSolution();
  adjointSolveCompleted_ = true;
}

void GAASPInterface::computeErrorEstimate() {
  TEST_FOR_EXCEPTION(!adjointSolveCompleted_, std::logic_error, 
      "Error, an adjoint solve must be completed before an error estimate can be calculated!");
  GAASP::GErrorEstimate errorEstimator(
      fSoln_, 
      aSoln_, 
      fDer_, 
      mesh_, 
      initialCondition_.size(), 
      nSteps_, 
      gModel_, 
      method_
      );
	errorEstimator.compute();
	errorEstimator.printErrorEstimates();
  // Copy data into GAASPErrorEstimate object
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

  int dim = tModel_->get_x_space()->dim();
  initialCondition_.clear();
  for (int i=0; i<dim ; ++i) {
    // TODO:  Change this to Serial double * view
//    *out << "GAASPInterface::setThyraModelEvaluator:             x_init["<<i<<"] = " << Thyra::get_ele<double>(*x_init,i) << std::endl;
    initialCondition_.push_back(Thyra::get_ele<double>(*x_init,i));
//    *out << "GAASPInterface::setThyraModelEvaluator:  initialCondition_["<<i<<"] = " << initialCondition_[i] << std::endl;
  }

  gModel_ = gModel_ThyraModelEvaluator(tModel_);
}

} // namespace Rythmos


