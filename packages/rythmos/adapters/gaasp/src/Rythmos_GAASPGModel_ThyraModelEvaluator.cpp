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


#include "Rythmos_GAASPGModel_ThyraModelEvaluator.hpp"
#include "Rythmos_GAASPHelpers.hpp"
#include "Teuchos_RCPBoostSharedPtrConversions.hpp"

namespace Rythmos {

int GModel_ThyraModelEvaluator::calcDerivs(double *yin, double *yout, double t) {
  TEST_FOR_EXCEPT(is_null(thyraModel_));
  // Set up const view of yin as Thyra::VectorBase<double>

  // Convert double * to const Thyra::VectorBase view
  Teuchos::RCP<const Thyra::VectorBase<double> > thyra_yin = 
    constThyraVectorBaseDouble(thyraModel_->get_x_space(), yin, dim_);
  
  // Set x on the inArgs.
  inArgs_.set_x(thyra_yin);
  
  // Set t on the inArgs.
  inArgs_.set_t(t);
  
  // Set up non const view of yout as Thyra::VectorBase<double>

  // Convert double * to nonconst Thyra::VectorBase view
  Teuchos::RCP<Thyra::VectorBase<double> > thyra_yout = 
    thyraVectorBaseDouble( thyraModel_->get_f_space(), yout, dim_ );

  // Set f on the outArgs.
  outArgs_.set_f(thyra_yout);

  // Evaluate model
  thyraModel_->evalModel(inArgs_,outArgs_);
  
  // All the magic happens when we go out of scope and the views are deleted.
  return(0);
}


void GModel_ThyraModelEvaluator::updateParameters(std::vector<double>) {
  TEST_FOR_EXCEPT(true);
}

int GModel_ThyraModelEvaluator::getDim() {
  TEST_FOR_EXCEPT(is_null(thyraModel_));
  return(dim_);
}

void GModel_ThyraModelEvaluator::setThyraModel(Teuchos::RCP<Thyra::ModelEvaluator<double> > model) {
  TEST_FOR_EXCEPT(is_null(model));
  thyraModel_ = model;
  dim_ = thyraModel_->get_x_space()->dim();
  inArgs_ = thyraModel_->createInArgs();
  outArgs_ = thyraModel_->createOutArgs();
}
  
// Helper functions:

// Convert a Thyra::ModelEvaluator into a GModelBase
boost::shared_ptr<GModel_ThyraModelEvaluator> gModel_ThyraModelEvaluator(Teuchos::RCP<Thyra::ModelEvaluator<double> > tModel) {
  TEST_FOR_EXCEPT(is_null(tModel));
  boost::shared_ptr<GModel_ThyraModelEvaluator> gModel(new GModel_ThyraModelEvaluator);
  gModel->setThyraModel(tModel);
  return(gModel);
}

} // namespace Rythmos

