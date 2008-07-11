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

#include "Rythmos_GAASPErrorEstimator.hpp"
#include "Rythmos_GAASPHelpers.hpp"
#include "GAdjointSolve.h"
#include "GErrorEstimate.h"

namespace Rythmos {
  
// GAASPErrorEstimator Definitions

GAASPErrorEstimator::GAASPErrorEstimator():
  qtyOfInterest_(INVALID_ERROR_QTY),
  isInitialized_(false)
{}

void GAASPErrorEstimator::setModel( Teuchos::RCP<Thyra::ModelEvaluator<double> > model ) {
  TEST_FOR_EXCEPT(is_null(model));
  model_ = model;
}

void GAASPErrorEstimator::setQuantityOfInterest( 
      ERROR_QUANTITY_OF_INTEREST qtyOfInterest
      ) {
  qtyOfInterest_ = qtyOfInterest;
}

void GAASPErrorEstimator::initialize_() {
  TEST_FOR_EXCEPT(is_null(model_));
  // Create GAASP interface
  gaaspInterfacePtr_ = Teuchos::rcp(new GAASPInterface );
  // Verify we have a valid model.
  gaaspInterfacePtr_->setThyraModelEvaluator(model_);
  // Verify we have a valid quantity of interest
  isInitialized_ = true;
}

Teuchos::RCP<ErrorEstimateBase<double> > GAASPErrorEstimator::getErrorEstimate() {
  if (!isInitialized_) {
    initialize_();
  }
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"GAASPErrorEstimator::");
  *out << "getErrorEstimate:  Calling GAASPInterface::forwardSolve()..." << std::endl;
  gaaspInterfacePtr_->forwardSolve();
  *out << "getErrorEstimate:  Calling GAASPInterface::adjointSolve()..." << std::endl;
  gaaspInterfacePtr_->adjointSolve();
  *out << "getErrorEstimate:  Calling GAASPInterface::computErrorEstimate()..." << std::endl;
  gaaspInterfacePtr_->computeErrorEstimate();
  // Copy Error Estimates into ErrorEstimate object
  Teuchos::RCP<GAASPErrorEstimate> gaaspEE = Teuchos::rcp(new GAASPErrorEstimate);
  // Pass out ErrorEstimate object
  return(gaaspEE);
}

std::string GAASPErrorEstimator::description() const
{
  std::string name = "Rythmos::GAASPErrorEstimator";
  return(name);
}

void GAASPErrorEstimator::describe(
  Teuchos::FancyOStream                &out
  ,const Teuchos::EVerbosityLevel      verbLevel
  ) const
{
  out << description() << "::describe" << std::endl;
}

void GAASPErrorEstimator::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList)
{
  paramList_ = paramList;
}

Teuchos::RCP<Teuchos::ParameterList> GAASPErrorEstimator::getParameterList()
{
  return(paramList_);
}

Teuchos::RCP<Teuchos::ParameterList> GAASPErrorEstimator::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_param_list = paramList_;
  paramList_ = Teuchos::null;
  return(temp_param_list);
}


} // namespace Rythmos

