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

#include "Rythmos_ImplicitRK_ConvergenceTest.hpp"
#include "Rythmos_RKButcherTableau.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"


namespace Rythmos {

using Thyra::VectorBase;
using Thyra::VectorSpaceBase;
using Teuchos::is_null;

RCP<ImplicitRKStepper<double> > getIRKStepperByOrder(int order) {
  RCP<SinCosModel> model = sinCosModel(true);
  DefaultRKButcherTableauFactory<double> rkbtFactory;
  Teuchos::ParameterList pl;
  pl.set("Selection Type", "Implicit method by order");
  pl.set("Implicit method by order", order);
  RKButcherTableau<double> rkbt = rkbtFactory.create(pl);
  // Create the base point
  Thyra::ModelEvaluatorBase::InArgs<double> basePoint = model->createInArgs();
  RCP<VectorBase<double> > base_x = Thyra::createMember(model->get_x_space());
  V_S(&*base_x,0.0);
  basePoint.set_x(base_x);
  RCP<VectorBase<double> > base_x_dot = Thyra::createMember(model->get_x_space());
  V_S(&*base_x_dot,1.0);
  basePoint.set_x_dot(base_x_dot);
  double base_t = 0.0;
  basePoint.set_t(base_t);
  // Create the nonlinear solver
  RCP<Rythmos::TimeStepNonlinearSolver<double> >
    nonlinearSolver = Rythmos::timeStepNonlinearSolver<double>();
  // Create the IRK W factory
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory = 
    Thyra::defaultSerialDenseLinearOpWithSolveFactory<double>();
  RCP<ImplicitRKStepper<double> > stepper = implicitRKStepper<double>(model,nonlinearSolver, irk_W_factory, rkbt);
  return stepper;
}

RCP<ImplicitRKStepper<double> > getSDIRKStepperByIndex(int index) {
  RCP<SinCosModel> model = sinCosModel(true);
  DefaultRKButcherTableauFactory<double> rkbtFactory;
  Teuchos::ParameterList pl;
  pl.set("Selection Type", "Method by name");
  if (index == 0) {
    pl.set("Method by name", "Backward Euler");
  } else if (index == 1) {
    pl.set("Method by name", "SDIRK 3 Stage 4th order");
  } else if (index == 2) {
    pl.set("Method by name", "SDIRK 5 Stage 4th order");
  } else if (index == 3) {
    pl.set("Method by name", "SDIRK 5 Stage 5th order");
  } else {
    TEST_FOR_EXCEPT_MSG(true, "Invalid index!\n");
  }
  RKButcherTableau<double> rkbt = rkbtFactory.create(pl);
  // Create the base point
  Thyra::ModelEvaluatorBase::InArgs<double> basePoint = model->createInArgs();
  RCP<VectorBase<double> > base_x = Thyra::createMember(model->get_x_space());
  V_S(&*base_x,0.0);
  basePoint.set_x(base_x);
  RCP<VectorBase<double> > base_x_dot = Thyra::createMember(model->get_x_space());
  V_S(&*base_x_dot,1.0);
  basePoint.set_x_dot(base_x_dot);
  double base_t = 0.0;
  basePoint.set_t(base_t);
  // Create the nonlinear solver
  RCP<Rythmos::TimeStepNonlinearSolver<double> >
    nonlinearSolver = Rythmos::timeStepNonlinearSolver<double>();
  // Create the IRK W factory
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory = 
    Thyra::defaultSerialDenseLinearOpWithSolveFactory<double>();
  RCP<ImplicitRKStepper<double> > stepper = implicitRKStepper<double>(model,nonlinearSolver, irk_W_factory, rkbt);
  return stepper;
}

/*
TEUCHOS_UNIT_TEST( Rythmos_ImplicitRKStepper, GlobalErrorConvergenceStudy ) {

  SinCosModelIRKStepperFactory irkFactory;
  for (int order=1 ; order<=1 ; ++order) {
    irkFactory.setOrder(order);
    double slope = computeOrderByGlobalErrorConvergenceStudy(irkFactory);

    double tol = 1.0e-10;
    TEST_FLOATING_EQUALITY( slope, 1.0*order, tol ); // is slope close to order?
  }
}

TEUCHOS_UNIT_TEST( Rythmos_ImplicitRKStepper, LocalErrorConvergenceStudy ) {
  SinCosModelIRKStepperFactory irkFactory;
  for (int order=1 ; order<=1 ; ++order) {
    irkFactory.setOrder(order);
    double slope = computeOrderByLocalErrorConvergenceStudy(irkFactory);

    int localOrder = order+1; // I don't know why the order is coming out one higher than it should!?!
    double tol = 1.0e-10;
    TEST_FLOATING_EQUALITY( slope, 1.0*localOrder, tol ); // is slope close to order?
  }
}
*/

TEUCHOS_UNIT_TEST( Rythmos_ImplicitRKStepper, SDIRKGlobalErrorConvergenceStudy ) {

  SinCosModelSDIRKStepperFactory sdirkFactory;
  for (int index=0 ; index<4 ; ++index) {
    sdirkFactory.setIndex(index);
    double slope = computeOrderByGlobalErrorConvergenceStudy(sdirkFactory);

    int order = sdirkFactory.create()->getOrder();

    double tol = 1.0e-2;
    if (index == 3) { tol = 1.0e-1; }
    TEST_FLOATING_EQUALITY( slope, 1.0*order, tol ); // is slope close to order?
  }
}

TEUCHOS_UNIT_TEST( Rythmos_ImplicitRKStepper, SDIRKLocalErrorConvergenceStudy ) {
  SinCosModelSDIRKStepperFactory sdirkFactory;
  for (int index=0 ; index<4 ; ++index) {
    sdirkFactory.setIndex(index);
    double slope = computeOrderByLocalErrorConvergenceStudy(sdirkFactory);

    int order = sdirkFactory.create()->getOrder();

    int localOrder = order+1; // I don't know why the order is coming out one higher than it should!?!
    double tol = 1.0e-1;
    if (index == 3) { tol = 1.0e-3; }
    TEST_FLOATING_EQUALITY( slope, 1.0*localOrder, tol ); // is slope close to order?
  }
}


} // namespace Rythmos

