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
#ifndef Rythmos_ERK_CONVERGENCETEST_H
#define Rythmos_ERK_CONVERGENCETEST_H

#include "Rythmos_Types.hpp"
#include "Rythmos_ConvergenceTestHelpers.hpp"
#include "Rythmos_ExplicitRKStepper.hpp"

namespace Rythmos {

using Thyra::ModelEvaluator;

template<class Scalar>
Array<std::string> getERKButcherTableauNames()
{
  Array<std::string> explicitRKTableauNames;
  RCP<DefaultRKButcherTableauFactory<Scalar> > rkbtFactory = rKButcherTableauFactory<Scalar>();
  RCP<const ParameterList> validPL = rkbtFactory->getValidParameters();
  Teuchos::ParameterList::ConstIterator plIt = validPL->begin();
  for (;plIt != validPL->end() ; plIt++) {
    std::string rkbt_name = validPL->name(plIt);
    if (rkbt_name == "RKButcherTableau Type") { continue; }
    RCP<RKButcherTableauBase<Scalar> > rkbt = rkbtFactory->create(rkbt_name);
    if (determineRKBTType(*rkbt) == RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_ERK) {
      explicitRKTableauNames.push_back(rkbt_name);
    }
  }
  return explicitRKTableauNames;
}

template<class Scalar>
class ExplicitRKStepperFactory : public virtual StepperFactoryBase<Scalar>
{
  public:
    ExplicitRKStepperFactory(RCP<ModelFactoryBase<Scalar> > modelFactory)
    {
      modelFactory_ = modelFactory;
      index_ = 0;
      explicitRKNames_ = getERKButcherTableauNames<Scalar>();
    }
    virtual ~ExplicitRKStepperFactory() {}
    RCP<StepperBase<Scalar> > getStepper() const
    {
      RCP<ModelEvaluator<Scalar> > model = modelFactory_->getModel();
      RCP<RKButcherTableauBase<Scalar> > rkbt = rkbtFactory_.create(explicitRKNames_[index_]);
      RCP<ExplicitRKStepper<Scalar> > stepper = explicitRKStepper<Scalar>(model,rkbt);
      return(stepper);
    }
    void setIndex(int index)
    {
      index_ = index;
    }
    int maxIndex()
    {
      return Teuchos::as<int>(explicitRKNames_.size());
    }
  private:
    RCP<ModelFactoryBase<Scalar> > modelFactory_;
    int index_;
    Array<std::string> explicitRKNames_;
    DefaultRKButcherTableauFactory<Scalar> rkbtFactory_;
};
// non-member constructor
template<class Scalar>
RCP<ExplicitRKStepperFactory<Scalar> > explicitRKStepperFactory( 
    RCP<ModelFactoryBase<Scalar> > modelFactory)
{
  RCP<ExplicitRKStepperFactory<Scalar> > erkFactory = Teuchos::rcp(
      new ExplicitRKStepperFactory<Scalar>(modelFactory)
      );
  return erkFactory;
}

} // namespace Rythmos 

#endif // Rythmos_ERK_CONVERGENCETEST_H

