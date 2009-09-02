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


#include "Teuchos_UnitTestHarness.hpp"

#include "Rythmos_ExplicitRK_ConvergenceTest.hpp"


namespace Rythmos {

using Thyra::VectorBase;
using Thyra::VectorSpaceBase;
using Teuchos::is_null;

TEUCHOS_UNIT_TEST( Rythmos_ExplicitRKStepper, GlobalErrorConvergenceStudy ) {

  RCP<SinCosModelFactory> modelFactory = sinCosModelFactory(false);
  RCP<SinCosModelExactSolutionObject> exactSolution = sinCosModelExactSolutionObject(modelFactory);
  RCP<ExplicitRKStepperFactory<double> > stepperFactory = explicitRKStepperFactory<double>(modelFactory);
  StepperFactoryAndExactSolutionObject<double> stepperFactoryAndExactSolution(stepperFactory,exactSolution);

  int N = stepperFactory->maxIndex();
  for (int index=0; index<N ; ++index) {
    stepperFactory->setIndex(index);

//    RCP<Teuchos::FancyOStream> fancyOut = Teuchos::VerboseObjectBase::getDefaultOStream();
//    RCP<ExplicitRKStepper<double> > erkStepper = rcp_dynamic_cast<ExplicitRKStepper<double> >(stepperFactoryAndExactSolution.getStepper(),true);
//    erkStepper->getRKButcherTableau().describe(*fancyOut,Teuchos::VERB_EXTREME);

    double slope = computeOrderByGlobalErrorConvergenceStudy(stepperFactoryAndExactSolution);

    int order = stepperFactoryAndExactSolution.getStepper()->getOrder();
    double tol = 1.0e-1;
    TEST_FLOATING_EQUALITY( slope, 1.0*order, tol ); // is slope close to order?
  }
}


TEUCHOS_UNIT_TEST( Rythmos_ExplicitRKStepper, LocalErrorConvergenceStudy ) {

  RCP<SinCosModelFactory> modelFactory = sinCosModelFactory(false);
  RCP<SinCosModelExactSolutionObject> exactSolution = sinCosModelExactSolutionObject(modelFactory);
  RCP<ExplicitRKStepperFactory<double> > stepperFactory = explicitRKStepperFactory<double>(modelFactory);
  StepperFactoryAndExactSolutionObject<double> stepperFactoryAndExactSolution(stepperFactory,exactSolution);

  int N = stepperFactory->maxIndex();
  for (int index=0 ; index<N ; ++index) {
    stepperFactory->setIndex(index);

//    RCP<Teuchos::FancyOStream> fancyOut = Teuchos::VerboseObjectBase::getDefaultOStream();
//    RCP<ExplicitRKStepper<double> > erkStepper = rcp_dynamic_cast<ExplicitRKStepper<double> >(stepperFactoryAndExactSolution.getStepper(),true);
//    erkStepper->getRKButcherTableau().describe(*fancyOut,Teuchos::VERB_EXTREME);
//
    double slope = computeOrderByLocalErrorConvergenceStudy(stepperFactoryAndExactSolution);

    int order = stepperFactoryAndExactSolution.getStepper()->getOrder();
    int localOrder = order+1; // I don't know why the order is coming out one higher than it should!?!
    double tol = 1.0e-2;
    TEST_FLOATING_EQUALITY( slope, 1.0*localOrder, tol ); // is slope close to order?
  }
}



} // namespace Rythmos

