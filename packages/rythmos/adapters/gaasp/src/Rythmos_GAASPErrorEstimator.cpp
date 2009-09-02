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
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

namespace Rythmos {
  
// Static members

const std::string GAASPErrorEstimator::GAASPInterface_name_ = "GAASP Interface Parameters";
  
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
  // Pass parameter sublist to GAASP Interface
  Teuchos::RCP<Teuchos::ParameterList> pl = sublist(paramList_, GAASPInterface_name_);
  gaaspInterfacePtr_->setParameterList(pl);
  
  // Verify we have a valid model.
  gaaspInterfacePtr_->setThyraModelEvaluator(model_);
  // Verify we have a valid quantity of interest
  isInitialized_ = true;
}

Teuchos::RCP<const ErrorEstimateBase<double> > GAASPErrorEstimator::getErrorEstimate() {
  if (!isInitialized_) {
    initialize_();
  }
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"GAASPErrorEstimator::");
  if (Teuchos::as<int>(verbLevel) != Teuchos::VERB_NONE) {
    *out << "getErrorEstimate:  Calling GAASPInterface::forwardSolve()..." << std::endl;
  }
  gaaspInterfacePtr_->forwardSolve();
  if (Teuchos::as<int>(verbLevel) != Teuchos::VERB_NONE) {
    *out << "getErrorEstimate:  Calling GAASPInterface::adjointSolve()..." << std::endl;
  }
  gaaspInterfacePtr_->adjointSolve();
  if (Teuchos::as<int>(verbLevel) != Teuchos::VERB_NONE) {
    *out << "getErrorEstimate:  Calling GAASPInterface::computeErrorEstimate()..." << std::endl;
  }
  Teuchos::RCP<const GAASPErrorEstimate> gaaspEE = gaaspInterfacePtr_->computeErrorEstimate();
  if (Teuchos::as<int>(verbLevel) != Teuchos::VERB_NONE) {
    *out << "getErrorEstimate:  Global Error Estimate = " << fabs(gaaspEE->getTotalError()) << std::endl;
  }
  return(gaaspEE);
}

Teuchos::RCP<const ErrorEstimateBase<double> > GAASPErrorEstimator::controlGlobalError(double uTOL) {
  if (!isInitialized_) {
    initialize_();
  }
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"GAASPErrorEstimator::");

  paramList_->sublist(GAASPInterface_name_).set("uTOL",uTOL);
  gaaspInterfacePtr_->setParameterList(sublist(paramList_,GAASPInterface_name_));
  if (Teuchos::as<int>(verbLevel) != Teuchos::VERB_NONE) {
    *out << "controlGlobalError:  Global error tolerance = " << uTOL << std::endl;
  }

  Teuchos::RCP<const ErrorEstimateBase<double> > gaaspEE;

  gaaspEE = getErrorEstimate();
  while (fabs(gaaspEE->getTotalError()) > uTOL) {
    if (Teuchos::as<int>(verbLevel) != Teuchos::VERB_NONE) {
      *out << "controlGlobalError:  Calling GAASPInterface::refineMesh()..." << std::endl;
    }
    gaaspInterfacePtr_->refineMesh();
    gaaspEE = getErrorEstimate();
  }

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

void GAASPErrorEstimator::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& paramList)
{
  TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParametersAndSetDefaults(*this->getValidParameters(),0);
  paramList_ = paramList;
  Teuchos::readVerboseObjectSublist(&*paramList_,this);
}

Teuchos::RCP<Teuchos::ParameterList>
GAASPErrorEstimator::getNonconstParameterList()
{
  return(paramList_);
}

Teuchos::RCP<Teuchos::ParameterList>
GAASPErrorEstimator::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_param_list = paramList_;
  paramList_ = Teuchos::null;
  return(temp_param_list);
}

Teuchos::RCP<const Teuchos::ParameterList> GAASPErrorEstimator::getValidParameters() const {
  static Teuchos::RCP<Teuchos::ParameterList> validPL;

  if (is_null(validPL)) {

    Teuchos::RCP<Teuchos::ParameterList>
      pl = Teuchos::parameterList();

    Teuchos::RCP<Teuchos::ParameterList> gaaspPL = sublist(pl,GAASPInterface_name_);
    GAASPInterface gaaspI;
    gaaspPL->setParameters(*(gaaspI.getValidParameters()));
    
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

