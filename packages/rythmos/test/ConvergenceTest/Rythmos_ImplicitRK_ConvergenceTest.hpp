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
#ifndef Rythmos_IRK_CONVERGENCETEST_H
#define Rythmos_IRK_CONVERGENCETEST_H

#include "Rythmos_Types.hpp"
#include "Rythmos_ConvergenceTestHelpers.hpp"
#include "Rythmos_ImplicitRKStepper.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"

namespace Rythmos {

template<class Scalar>
Array<std::string> getIRKButcherTableauNames()
{
  Array<std::string> allRKTableauNames = getS_RKButcherTableauMethodNames();
  Array<std::string> iRKTableauNames;

  DefaultRKButcherTableauFactory<Scalar> rkbtFactory;
  Teuchos::ParameterList pl;
  pl.set("Selection Type", "Method by name");

  int N = Teuchos::as<int>(allRKTableauNames.size());
  for (int i=0 ; i<N ; ++i) {
    pl.set("Method by name", allRKTableauNames[i]);
    RKButcherTableau<Scalar> rkbt = rkbtFactory.create(pl);
    if (determineRKBTType(rkbt) == RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_IRK) {
      iRKTableauNames.push_back(allRKTableauNames[i]);
    }
  }
  return iRKTableauNames;
}

template<class Scalar>
Array<std::string> getDIRKButcherTableauNames()
{
  Array<std::string> allRKTableauNames = getS_RKButcherTableauMethodNames();
  Array<std::string> dIRKTableauNames;

  DefaultRKButcherTableauFactory<Scalar> rkbtFactory;
  Teuchos::ParameterList pl;
  pl.set("Selection Type", "Method by name");

  int N = Teuchos::as<int>(allRKTableauNames.size());
  for (int i=0 ; i<N ; ++i) {
    pl.set("Method by name", allRKTableauNames[i]);
    RKButcherTableau<Scalar> rkbt = rkbtFactory.create(pl);
    if (    (determineRKBTType(rkbt) == RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_DIRK)
         || (determineRKBTType(rkbt) == RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_SDIRK)
        ) {
      dIRKTableauNames.push_back(allRKTableauNames[i]);
    }
  }
  return dIRKTableauNames;
}

template<class Scalar>
class ImplicitRKStepperFactory : public virtual StepperFactoryBase<Scalar>
{
  public:
    ImplicitRKStepperFactory(RCP<ModelFactoryBase<Scalar> > modelFactory) 
    { 
      modelFactory_ = modelFactory;
      index_ = 0;
      implicitRKNames_ = getIRKButcherTableauNames<Scalar>();
    }
    virtual ~ImplicitRKStepperFactory() {}
    void setIndex(int index) { index_ = index; }
    int maxIndex() { return Teuchos::as<int>(implicitRKNames_.size()); }
    RCP<StepperBase<Scalar> > getStepper() const 
    { 
      // Get the model:
      RCP<ModelEvaluator<Scalar> > model = modelFactory_->getModel();
      // Get the RKBT:
      Teuchos::ParameterList paramList;
      paramList.set("Selection Type", "Method by name");
      paramList.set("Method by name", implicitRKNames_[index_]);
      RKButcherTableau<Scalar> rkbt = rkbtFactory_.create(paramList);
      // Create the nonlinear solver
      RCP<Rythmos::TimeStepNonlinearSolver<double> >
        nonlinearSolver = Rythmos::timeStepNonlinearSolver<double>();
      // Create the W_factory
      RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > irk_W_factory = modelFactory_->get_W_factory();
      RCP<ImplicitRKStepper<Scalar> > stepper = implicitRKStepper<Scalar>(model,nonlinearSolver,irk_W_factory, rkbt);
      // Set verbosity on Stepper:
      RCP<Teuchos::ParameterList> stepperPL = Teuchos::parameterList();
      RCP<Teuchos::ParameterList> stepperVOPL = Teuchos::sublist(stepperPL,"VerboseObject");
      stepperVOPL->set("Verbosity Level","none");
      stepper->setParameterList(stepperPL);
      return(stepper);
    }
  private:
    RCP<ModelFactoryBase<Scalar> > modelFactory_;
    int index_;
    Array<std::string> implicitRKNames_;
    DefaultRKButcherTableauFactory<Scalar> rkbtFactory_;
};
template<class Scalar>
RCP<ImplicitRKStepperFactory<Scalar> > implicitRKStepperFactory(RCP<ModelFactoryBase<Scalar> > modelFactory) 
{
  RCP<ImplicitRKStepperFactory<Scalar> > irkFactory = Teuchos::rcp(
      new ImplicitRKStepperFactory<Scalar>(modelFactory)
      );
  return irkFactory;
}

template<class Scalar>
class DiagonalImplicitRKStepperFactory : public virtual StepperFactoryBase<Scalar>
{
  public:
    DiagonalImplicitRKStepperFactory(RCP<ModelFactoryBase<Scalar> > modelFactory) 
    { 
      modelFactory_ = modelFactory;
      index_ = 0;
      implicitRKNames_ = getDIRKButcherTableauNames<Scalar>();
    }
    virtual ~DiagonalImplicitRKStepperFactory() {}
    void setIndex(int index) { index_ = index; }
    int maxIndex() { return Teuchos::as<int>(implicitRKNames_.size()); }
    RCP<StepperBase<Scalar> > getStepper() const 
    { 
      // Get the model:
      RCP<ModelEvaluator<Scalar> > model = modelFactory_->getModel();
      // Get the RKBT:
      Teuchos::ParameterList paramList;
      paramList.set("Selection Type", "Method by name");
      paramList.set("Method by name", implicitRKNames_[index_]);
      RKButcherTableau<Scalar> rkbt = rkbtFactory_.create(paramList);
      // Create the nonlinear solver
      RCP<Rythmos::TimeStepNonlinearSolver<double> >
        nonlinearSolver = Rythmos::timeStepNonlinearSolver<double>();
      // Create the W_factory
      RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > irk_W_factory = modelFactory_->get_W_factory();
      RCP<ImplicitRKStepper<Scalar> > stepper = implicitRKStepper<Scalar>(model,nonlinearSolver,irk_W_factory, rkbt);
      return(stepper);
    }
  private:
    RCP<ModelFactoryBase<Scalar> > modelFactory_;
    int index_;
    Array<std::string> implicitRKNames_;
    DefaultRKButcherTableauFactory<Scalar> rkbtFactory_;
};
// non-member constructor
template<class Scalar>
RCP<DiagonalImplicitRKStepperFactory<Scalar> > diagonalImplicitRKStepperFactory(
    RCP<ModelFactoryBase<Scalar> > modelFactory
    )
{
  RCP<DiagonalImplicitRKStepperFactory<Scalar> > dirkStepper = Teuchos::rcp(
      new DiagonalImplicitRKStepperFactory<Scalar>(modelFactory)
      );
  return dirkStepper;
}

} // namespace Rythmos 

#endif // Rythmos_IRK_CONVERGENCETEST_H

