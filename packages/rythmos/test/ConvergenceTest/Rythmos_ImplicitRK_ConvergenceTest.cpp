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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER


#include "Teuchos_UnitTestHarness.hpp"

#include "Rythmos_ImplicitRK_ConvergenceTest.hpp"
#include "Rythmos_RKButcherTableau.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"


namespace Rythmos {

using Thyra::VectorBase;
using Thyra::VectorSpaceBase;
using Teuchos::is_null;


TEUCHOS_UNIT_TEST( Rythmos_ImplicitRKStepper, GlobalErrorConvergenceStudy ) {

  RCP<DiagonalModelFactory> modelFactory = diagonalModelFactory();
  RCP<DiagonalModelExactSolutionObject> exactSolution = diagonalModelExactSolutionObject(modelFactory);

  RCP<ImplicitRKStepperFactory<double> > stepperFactory = implicitRKStepperFactory<double>(modelFactory);
  StepperFactoryAndExactSolutionObject<double> stepperFactoryAndExactSolution(stepperFactory,exactSolution);

  int N = stepperFactory->maxIndex();
  for (int index=0 ; index<N ; ++index) {
    stepperFactory->setIndex(index);

    //RCP<Teuchos::FancyOStream> fancyOut = Teuchos::VerboseObjectBase::getDefaultOStream();
    //RCP<ImplicitRKStepper<double> > dirkStepper = rcp_dynamic_cast<ImplicitRKStepper<double> >(stepperFactoryAndExactSolution.getStepper(),true);
    //dirkStepper->getRKButcherTableau().describe(*fancyOut,Teuchos::VERB_EXTREME);

    int order = stepperFactoryAndExactSolution.getStepper()->getOrder();

    int numCuts = 1;
    //int numCuts = 4;
    if (order > 4) { numCuts = 2; }
    double slope = computeOrderByGlobalErrorConvergenceStudy(stepperFactoryAndExactSolution,numCuts);

    double tol = 2.0e-1;
    TEST_FLOATING_EQUALITY( slope, 1.0*order, tol ); // is slope close to order?
  }
}

/*
TEUCHOS_UNIT_TEST( Rythmos_ImplicitRKStepper, LocalErrorConvergenceStudy ) {

  RCP<DiagonalModelFactory> modelFactory = diagonalModelFactory();
  RCP<DiagonalModelExactSolutionObject> exactSolution = diagonalModelExactSolutionObject(modelFactory);

  RCP<ImplicitRKStepperFactory<double> > stepperFactory = implicitRKStepperFactory<double>(modelFactory);
  StepperFactoryAndExactSolutionObject<double> stepperFactoryAndExactSolution(stepperFactory,exactSolution);

  int N = stepperFactory->maxIndex();
  for (int index=0 ; index<N ; ++index) {
    stepperFactory->setIndex(index);

    RCP<Teuchos::FancyOStream> fancyOut = Teuchos::VerboseObjectBase::getDefaultOStream();
    RCP<ImplicitRKStepper<double> > dirkStepper = rcp_dynamic_cast<ImplicitRKStepper<double> >(stepperFactoryAndExactSolution.getStepper(),true);
    dirkStepper->getRKButcherTableau().describe(*fancyOut,Teuchos::VERB_EXTREME);

    int order = stepperFactoryAndExactSolution.getStepper()->getOrder();

    int numCuts = 9;
    if (order > 4) { numCuts = 8; }
    double slope = computeOrderByLocalErrorConvergenceStudy(stepperFactoryAndExactSolution,numCuts);


    int localOrder = order+1; // I don't know why the order is coming out one higher than it should!?!
    double tol = 1.0e-10;
    TEST_FLOATING_EQUALITY( slope, 1.0*localOrder, tol ); // is slope close to order?
  }
}
*/


TEUCHOS_UNIT_TEST( Rythmos_ImplicitRKStepper, DIRKGlobalErrorConvergenceStudy ) {

  RCP<SinCosModelFactory> modelFactory = sinCosModelFactory(true);
  RCP<SinCosModelExactSolutionObject> exactSolution = sinCosModelExactSolutionObject(modelFactory);
  RCP<DiagonalImplicitRKStepperFactory<double> > stepperFactory = diagonalImplicitRKStepperFactory<double>(modelFactory);
  StepperFactoryAndExactSolutionObject<double> stepperFactoryAndExactSolution(stepperFactory,exactSolution);

  int N = stepperFactory->maxIndex();
  for (int index=0 ; index<N ; ++index) {
    stepperFactory->setIndex(index);

    //RCP<Teuchos::FancyOStream> fancyOut = Teuchos::VerboseObjectBase::getDefaultOStream();
    //RCP<ImplicitRKStepper<double> > dirkStepper = rcp_dynamic_cast<ImplicitRKStepper<double> >(stepperFactoryAndExactSolution.getStepper(),true);
    //dirkStepper->getRKButcherTableau().describe(*fancyOut,Teuchos::VERB_EXTREME);

    int order = stepperFactoryAndExactSolution.getStepper()->getOrder();

    int numCuts = 3;
    double slope = computeOrderByGlobalErrorConvergenceStudy(stepperFactoryAndExactSolution,numCuts);


    double tol = 1.0e-1;
    TEST_FLOATING_EQUALITY( slope, 1.0*order, tol ); // is slope close to order?
  }
}

/*
TEUCHOS_UNIT_TEST( Rythmos_ImplicitRKStepper, DIRKLocalErrorConvergenceStudy ) {
  RCP<SinCosModelFactory> modelFactory = sinCosModelFactory(true);
  RCP<SinCosModelExactSolutionObject> exactSolution = sinCosModelExactSolutionObject(modelFactory);
  RCP<DiagonalImplicitRKStepperFactory<double> > stepperFactory = diagonalImplicitRKStepperFactory<double>(modelFactory);
  StepperFactoryAndExactSolutionObject<double> stepperFactoryAndExactSolution(stepperFactory,exactSolution);

  int N = stepperFactory->maxIndex();
  for (int index=0 ; index<N ; ++index) {
    stepperFactory->setIndex(index);

    RCP<Teuchos::FancyOStream> fancyOut = Teuchos::VerboseObjectBase::getDefaultOStream();
    RCP<ImplicitRKStepper<double> > dirkStepper = rcp_dynamic_cast<ImplicitRKStepper<double> >(stepperFactoryAndExactSolution.getStepper(),true);
    dirkStepper->getRKButcherTableau().describe(*fancyOut,Teuchos::VERB_EXTREME);

    double slope = computeOrderByLocalErrorConvergenceStudy(stepperFactoryAndExactSolution);

    int order = stepperFactoryAndExactSolution.getStepper()->getOrder();

    int localOrder = order+1; // I don't know why the order is coming out one higher than it should!?!
    double tol = 1.0e-10;
    TEST_FLOATING_EQUALITY( slope, 1.0*localOrder, tol ); // is slope close to order?
  }
}
*/


} // namespace Rythmos

