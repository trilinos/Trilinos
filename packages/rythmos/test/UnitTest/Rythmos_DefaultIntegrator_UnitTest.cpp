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

#include "Rythmos_Types.hpp"
#include "Rythmos_UnitTestHelpers.hpp"

#include "Rythmos_ExplicitRKStepper.hpp"
#include "Rythmos_RKButcherTableauBuilder.hpp"

#include "../SinCos/SinCosModel.hpp"

#include "Rythmos_DefaultIntegrator.hpp"
#include "Rythmos_InterpolationBuffer.hpp"
#include "Rythmos_SimpleIntegrationControlStrategy.hpp"
#include "Rythmos_IntegratorBuilder.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Thyra_DetachedVectorView.hpp"

namespace Rythmos {

using Thyra::VectorBase;

// Test the ERK stepper through the integrator
TEUCHOS_UNIT_TEST( Rythmos_DefaultIntegrator, ExplicitRKStepper ) {
  // Integrator
  RCP<DefaultIntegrator<double> > integrator = defaultIntegrator<double>();

  // Stepper
  double finalTime = 1.0;
  RCP<SinCosModel> model = sinCosModel(false);
  RCP<RKButcherTableauBase<double> > rkbt = createRKBT<double>("Forward Euler");
  TEST_ASSERT( !is_null(rkbt) );
  RCP<ExplicitRKStepper<double> > stepper = explicitRKStepper<double>(model,rkbt);
  // Set initial condition on stepper
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  stepper->setInitialCondition(ic);
  // Set stepper on integrator:
  integrator->setStepper(stepper, finalTime);

  // Trailing Interpolation Buffer
  // Do not use with ExplicitRK yet as the InterpolationBuffer does not work without x_dot
  
  // IntegrationControlStrategy to specify fixed steps for this stepper
  double dt = 0.1;
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set("Take Variable Steps",false);
  pl->set("Fixed dt", dt);
  RCP<SimpleIntegrationControlStrategy<double> > intCont = simpleIntegrationControlStrategy<double>(pl);
  TimeRange<double> tr(0.0,finalTime);
  intCont->resetIntegrationControlStrategy(tr); // ??? Why do I need to do this?
  integrator->setIntegrationControlStrategy(intCont);

  // Ask integrator for points
  Array<double> time_vec;
  Array<RCP<const VectorBase<double> > > x_vec;
  double N = 10;
  for (int i=0 ; i<=N ; ++i) {
    double t = 0.0 + i*finalTime/N;
    time_vec.push_back(t);
  }
  integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

  // Verify that these points are accurate
  // Since we're using Forward Euler, we can write down the exact numerical solution
  double tol = 1.0e-10;
  double exact_x0 = 0.0; // nominal values on SinCosModel
  double exact_x1 = 1.0;
  for (int i=0 ; i<=N ; ++i) {
    {
      Thyra::ConstDetachedVectorView<double> x_vec_view( *(x_vec[i]) );
      TEST_FLOATING_EQUALITY( exact_x0, x_vec_view[0], tol );
      TEST_FLOATING_EQUALITY( exact_x1, x_vec_view[1], tol );
    }
    double x0 = exact_x0;
    double x1 = exact_x1;
    exact_x0 += dt*x1;
    exact_x1 -= dt*x0;
  }
}


TEUCHOS_UNIT_TEST( Rythmos_DefaultIntegrator, maxNumTimeSteps ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<SinCosModel> model = sinCosModel(true);
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*(ib->getValidParameters()));
  pl->sublist("Integrator Settings").sublist("Integrator Selection").set("Integrator Type","Default Integrator");
  pl->sublist("Integrator Settings").sublist("Integrator Selection").sublist("Default Integrator").set("Max Number Time Steps",10);
  pl->sublist("Stepper Settings").sublist("Stepper Selection").set("Stepper Type","Backward Euler");
  pl->sublist("Integration Control Strategy Selection").set("Integration Control Strategy Type","Simple Integration Control Strategy");
  pl->sublist("Integration Control Strategy Selection").sublist("Simple Integration Control Strategy").set("Take Variable Steps",false);
  pl->sublist("Integration Control Strategy Selection").sublist("Simple Integration Control Strategy").set("Fixed dt",0.01);
  ib->setParameterList(pl);
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver = timeStepNonlinearSolver<double>();
  RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
  Teuchos::Array<double> time_vec;
  time_vec.push_back(pl->sublist("Integrator Settings").get<double>("Final Time"));
  TEST_NOTHROW(
    integrator->getFwdPoints(time_vec,NULL,NULL,NULL)
    );
}


/*
TEUCHOS_UNIT_TEST( Rythmos_DefaultIntegrator, momento ) {
  RCP<Momento<double> > integrator_momento;
  // place to store solution at time = 1.0
  RCP<const VectorBase<double> > x_norestart; 
  RCP<const VectorBase<double> > xdot_norestart;
  // Create an integrator with Backward Euler 
  // Integrate to t = 0.5, pull out the momento, then integrate the rest of the way to t = 1.0.
  {
    RCP<SinCosModel> model = sinCosModel(true);
    Thyra::ModelEvaluator::InArgs<double> model_ic = model->getNominalValues();
    RCP<Thyra::NonLinearSolver<double> > nlSolver = timeStepNonlinearSolver<double>();
    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
    {
      RCP<ParameterList> pl = Teuchos::parameterList();
      pl->setParameters(*(ib->getValidParameters()));
      pl->sublist("Stepper Settings").sublist("Stepper Selection").set("Stepper Type","Backward Euler");
      pl->sublist("Integration Control Strategy Selection").set("Integration Control Strategy Type","Simple Integration Control Strategy");
      pl->sublist("Integration Control Strategy Selection").sublist("Simple Integration Control Strategy").set("Take Variable Steps",false);
      pl->sublist("Integration Control Strategy Selection").sublist("Simple Integration Control Strategy").set("Fixed dt",0.1);
      ib->setParameterList(pl);
    }
    RCP<IntegratorBase<double> > integrator = ib->create(model,model_ic,nlSolver);
    {
      // Integrate to t=0.5
      Array<double> time_vec;
      time_vec.push_back(0.5);
      Array<RCP<const VectorBase<double> > > x_vec;
      Array<RCP<const VectorBase<double> > > xdot_vec;
      Array<double> accuracy_vec;
      integrator->getFwdPoints(time_vec,&x_vec,&xdot_vec,&accuracy_vec);
    }
    // Pull out the momento.
    integrator_momento = integrator->getMomento();
    {
      // Finish integrating to t=1.0
      time_vec.push_back(1.0);
      Array<RCP<const VectorBase<double> > > x_vec;
      Array<RCP<const VectorBase<double> > > xdot_vec;
      Array<double> accuracy_vec;
      integrator->getFwdPoints(time_vec,&x_vec,&xdot_vec,&accuracy_vec);
      x_norestart = x_vec[0];
      xdot_norestart = xdot_vec[0];
    }
  }
  // Create an integrator with Backward Euler, pass the momento back, and verify the data is the same
  {
    RCP<SinCosModel> model = sinCosModel(true);
    Thyra::ModelEvaluator::InArgs<double> model_ic = model->getNominalValues();
    RCP<Thyra::NonLinearSolver<double> > nlSolver = timeStepNonlinearSolver<double>();
    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
    {
      RCP<ParameterList> pl = Teuchos::parameterList();
      pl->setParameters(*(ib->getValidParameters()));
      pl->sublist("Stepper Settings").sublist("Stepper Selection").set("Stepper Type","Backward Euler");
      pl->sublist("Integration Control Strategy Selection").set("Integration Control Strategy Type","Simple Integration Control Strategy");
      pl->sublist("Integration Control Strategy Selection").sublist("Simple Integration Control Strategy").set("Take Variable Steps",false);
      pl->sublist("Integration Control Strategy Selection").sublist("Simple Integration Control Strategy").set("Fixed dt",0.1);
      ib->setParameterList(pl);
    }
    RCP<IntegratorBase<double> > integrator = ib->create(model,model_ic,nlSolver);
    integrator->setMomento(*integrator_momento);
    {
      // Finish integrating to t=1.0
      time_vec.push_back(1.0);
      Array<RCP<const VectorBase<double> > > x_vec;
      Array<RCP<const VectorBase<double> > > xdot_vec;
      Array<double> accuracy_vec;
      integrator->getFwdPoints(time_vec,&x_vec,&xdot_vec,&accuracy_vec);
      
      RCP<const VectorBase<double> > x_restart = x_vec[0];
      RCP<const VectorBase<double> > xdot_restart = xdot_vec[0];
      // Verify that x_restart == x_norestart
      // Verify that xdot_restart == xdot_norestart
      TEST_ASSERT(
        Thyra::testRelNormDiffErr(
          "x_norestart", *x_norestart, "x_restart", *x_restart,
          "epsilon", SMT::eps(),
          "epsilon", SMT::eps(),
          0
          )
        );
      TEST_ASSERT(
        Thyra::testRelNormDiffErr(
          "xdot_norestart", *xdot_norestart, "xdot_restart", *xdot_restart,
          "epsilon", SMT::eps(),
          "epsilon", SMT::eps(),
          0
          )
        );
    }
  }
}
*/

} // namespace Rythmos

