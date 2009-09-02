//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
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

#include "Rythmos_UnitTestModels.hpp"
#include "../SinCos/SinCosModel.hpp"
#include "../VanderPol/VanderPolModel.hpp"

namespace Rythmos {

//RCP<TestModelBuilder> testModelBuilder()
//{
//  RCP<TestModelBuilder> tmb = rcp(new TestModelBuilder );
//  return tmb;
//}
//
//
//RCP<Thyra::ModelEvaluator<double> > createTestModel(
//    const std::string& modelName,
//    const RCP<ParameterList>& modelPL
//    )
//{
//  RCP<TestModelBuilder> tmb = testModelBuilder();
//  RCP<Thyra::ModelEvaluator<double> > model;
//  if (!Teuchos::is_null(modelPL)) {
//    RCP<ParameterList> tmbPL = Teuchos::parameterList();
//    tmbPL->set("Test Model Type",modelName);
//    ParameterList& modelSubList = tmbPL->sublist(modelName);
//    modelSubList.setParameters(*modelPL);
//    tmb->setParameterList(tmbPL);
//    model = tmb->create();
//  }
//  else {
//    model = tmb->create(modelName);
//  }
//  return model;
//}
//
//TestModelBuilder::TestModelBuilder()
//{
//  this->initializeDefaults_();
//}
//
//TestModelBuilder::~TestModelBuilder()
//{
//}
//
//void TestModelBuilder::setModelFactory(
//  const RCP<const Teuchos::AbstractFactory<Thyra::ModelEvaluator<double> > > &modelFactory,
//  const std::string &modelName
//  )
//{
//  builder_.setObjectFactory(modelFactory, modelName);
//}
//
//
//std::string
//TestModelBuilder::getModelName() const
//{
//  return builder_.getObjectName();
//}
//
//void TestModelBuilder::setParameterList(
//  RCP<Teuchos::ParameterList> const& paramList
//  )
//{
//  builder_.setParameterList(paramList);
//}
//
//
//RCP<Teuchos::ParameterList>
//TestModelBuilder::getNonconstParameterList()
//{
//  return builder_.getNonconstParameterList();
//}
//
//
//RCP<Teuchos::ParameterList>
//TestModelBuilder::unsetParameterList()
//{
//  return builder_.unsetParameterList();
//}
//
//
//RCP<const Teuchos::ParameterList>
//TestModelBuilder::getParameterList() const
//{
//  return builder_.getParameterList();
//}
//
//
//RCP<const Teuchos::ParameterList>
//TestModelBuilder::getValidParameters() const
//{
//  return builder_.getValidParameters();
//}
//
//
//RCP<Thyra::ModelEvaluator<double> >
//TestModelBuilder::create(
//  const std::string &modelName
//  ) const
//{
//  return builder_.create(modelName);
//}
//
//
//void TestModelBuilder::initializeDefaults_()
//{
//
//  using Teuchos::abstractFactoryStd;
//
//  builder_.setObjectName("Rythmos::TestModel");
//  builder_.setObjectTypeName("Test Model Type");
//
//  //
//  // Test Models
//  //
//  
//  builder_.setObjectFactory(
//      abstractFactoryStd< Thyra::ModelEvaluator<double>, SinCosModel >(),
//      "SinCos"
//      );
//
//  builder_.setObjectFactory(
//      abstractFactoryStd< Thyra::ModelEvaluator<double>, VanderPolModel >(),
//      "VanderPol"
//      );
//
//  //builder_.setObjectFactory(
//  //    abstractFactoryStd< Thyra::ModelEvaluator<double>, DiagonalTransientModel >(),
//  //    "DiagonalTransient"
//  //    );
//
//  builder_.setDefaultObject("SinCos");
//  
//}

} // namespace Rythmos

