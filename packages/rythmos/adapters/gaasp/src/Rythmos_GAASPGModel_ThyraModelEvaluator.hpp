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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_GAASP_GMODEL_THYRA_MODEL_EVALUATOR_H
#define Rythmos_GAASP_GMODEL_THYRA_MODEL_EVALUATOR_H

#include "Thyra_ModelEvaluator.hpp"
#include "GModelBase.h"
#include "Teuchos_RCPBoostSharedPtrConversions.hpp"

namespace Rythmos {

class GModel_ThyraModelEvaluator : public virtual GAASP::GModelBase {
  public:
    int calcDerivs(double *yin, double *yout, double t);
    void updateParameters(std::vector<double>);
    int getDim();
    void setThyraModel(Teuchos::RCP<Thyra::ModelEvaluator<double> > model);
  private:
    Teuchos::RCP<Thyra::ModelEvaluator<double> > thyraModel_;
    int dim_;
    Thyra::ModelEvaluator<double>::InArgs<double> inArgs_;
    Thyra::ModelEvaluator<double>::OutArgs<double> outArgs_;
};

// Helper Function

// Convert a Thyra::ModelEvaluator into a GModelBase
boost::shared_ptr<GModel_ThyraModelEvaluator> gModel_ThyraModelEvaluator(Teuchos::RCP<Thyra::ModelEvaluator<double> > tModel);

} // namespace Rythmos


#endif // Rythmos_GAASP_GMODEL_THYRA_MODEL_EVALUATOR_H
