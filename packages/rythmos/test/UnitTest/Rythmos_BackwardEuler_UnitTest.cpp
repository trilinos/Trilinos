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
#include "Rythmos_BackwardEulerStepper.hpp"
#include "../SinCos/SinCosModel.hpp"
#include "Rythmos_UnitTestModels.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "Rythmos_IntegratorBuilder.hpp"

namespace Rythmos {

TEUCHOS_UNIT_TEST( Rythmos_BackwardEulerStepper, clone ) {
  RCP<BackwardEulerStepper<double> > stepper = backwardEulerStepper<double>();
  TEST_ASSERT( stepper->supportsCloning() );
  RCP<StepperBase<double> > newStepper = stepper->cloneStepperAlgorithm();
  TEST_ASSERT( !Teuchos::is_null(newStepper) );
  {
    RCP<BackwardEulerStepper<double> > beStepper = Teuchos::rcp_dynamic_cast<BackwardEulerStepper<double> >(newStepper,false);
    TEST_ASSERT(!is_null(beStepper));
  }
}

TEUCHOS_UNIT_TEST( Rythmos_BackwardEulerStepper, momento_create ) {
  RCP<const MomentoBase<double> > momento;
  {
    RCP<BackwardEulerStepper<double> > stepper = backwardEulerStepper<double>();
    RCP<SinCosModel> model = sinCosModel(true);
    stepper->setModel(model);
    momento = stepper->getMomento();
    TEST_ASSERT( !is_null(momento) );
  }
  {
    RCP<SinCosModel> model = sinCosModel(true);
    RCP<BackwardEulerStepper<double> > stepper = backwardEulerStepper<double>();
    stepper->setMomento(momento.ptr(),model,Teuchos::null);
    momento = Teuchos::null;
    RCP<const Thyra::ModelEvaluator<double> > modelOut = stepper->getModel();
    TEST_ASSERT( !is_null(modelOut) );
    RCP<const SinCosModel> scModel = Teuchos::rcp_dynamic_cast<const SinCosModel>(modelOut,false);
    TEST_ASSERT( !is_null(scModel) );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_BackwardEulerStepper, restart ) {
  RCP<const MomentoBase<double> > stepper_momento;
  // place to store solution at time = 1.0
  RCP<const VectorBase<double> > x_norestart; 
  RCP<const VectorBase<double> > xdot_norestart;
  // Create Backward Euler stepper
  // Step to t = 0.5, pull out the momento, then step the rest of the way to t = 1.0.
  {
    RCP<SinCosModel> model = sinCosModel(true);
    Thyra::ModelEvaluatorBase::InArgs<double> model_ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > neSolver = timeStepNonlinearSolver<double>();
    RCP<BackwardEulerStepper<double> > stepper = backwardEulerStepper<double>(model,neSolver);
    stepper->setInitialCondition(model_ic);
    double dt = 0.1;
    // Step to t=0.5
    for (int i=0 ; i<5 ; ++i) {
      double dt_taken = stepper->takeStep(dt,STEP_TYPE_FIXED);
      TEST_ASSERT( dt_taken == dt );
    }
    // Pull out the momento
    stepper_momento = stepper->getMomento();
    // Step to t=1.0
    for (int i=0 ; i<5 ; ++i) {
      double dt_taken = stepper->takeStep(dt,STEP_TYPE_FIXED);
      TEST_ASSERT( dt_taken == dt );
    }
    {
      // Pull out the final solution
      Array<double> time_vec;
      time_vec.push_back(1.0);
      Array<RCP<const VectorBase<double> > > x_vec;
      Array<RCP<const VectorBase<double> > > xdot_vec;
      Array<double> accuracy_vec;
      stepper->getPoints(time_vec,&x_vec,&xdot_vec,&accuracy_vec);
      x_norestart = x_vec[0];
      xdot_norestart = xdot_vec[0];
    }
  }
  // Create a new Backward Euler Stepper
  {
    RCP<BackwardEulerStepper<double> > stepper = backwardEulerStepper<double>();
    RCP<SinCosModel> model = sinCosModel(true);
    RCP<Thyra::NonlinearSolverBase<double> > neSolver = timeStepNonlinearSolver<double>();
    // Put the momento back into the stepper
    stepper->setMomento(stepper_momento.ptr(),model,neSolver);

    double dt = 0.1;
    // Step to t=1.0
    for (int i=0 ; i<5 ; ++i) {
      double dt_taken = stepper->takeStep(dt,STEP_TYPE_FIXED);
      TEST_ASSERT( dt_taken == dt );
    }
    {
      // Pull out the final solution after restart
      Array<double> time_vec;
      time_vec.push_back(1.0);
      Array<RCP<const VectorBase<double> > > x_vec;
      Array<RCP<const VectorBase<double> > > xdot_vec;
      Array<double> accuracy_vec;
      stepper->getPoints(time_vec,&x_vec,&xdot_vec,&accuracy_vec);
      RCP<const VectorBase<double> > x_restart = x_vec[0];
      RCP<const VectorBase<double> > xdot_restart = xdot_vec[0];

      // Verify that x_restart == x_norestart
      RCP<VectorBase<double> > x_diff = createMember(x_restart->space());
      V_VmV( x_diff.ptr(), *x_norestart, *x_restart );
      double x_normDiff = norm(*x_diff);
      double tol = 1.0e-10;
      TEST_COMPARE( x_normDiff, <, tol );

      // Verify that xdot_restart == xdot_norestart
      RCP<VectorBase<double> > xdot_diff = createMember(x_restart->space());
      V_VmV( xdot_diff.ptr(), *xdot_norestart, *xdot_restart );
      double xdot_normDiff = norm(*xdot_diff);
      TEST_COMPARE( xdot_normDiff, <, tol );

    }
  }
}

TEUCHOS_UNIT_TEST( Rythmos_BackwardEulerStepper, checkConsistentState ) {
  {
    RCP<SinCosModel> model = sinCosModel(true);
    Thyra::ModelEvaluatorBase::InArgs<double> model_ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > neSolver = timeStepNonlinearSolver<double>();
    RCP<BackwardEulerStepper<double> > stepper = backwardEulerStepper<double>(model,neSolver);
    stepper->setInitialCondition(model_ic);
    double dt = 0.1;
    double dt_taken = stepper->takeStep(dt,STEP_TYPE_FIXED);
    TEST_ASSERT( dt == dt_taken );
    RCP<const MomentoBase<double> > momento = stepper->getMomento();
    TEST_THROW(stepper->setMomento(momento.ptr(),Teuchos::null,Teuchos::null), std::logic_error);
    TEST_THROW(stepper->setMomento(momento.ptr(),model,Teuchos::null), std::logic_error);
    TEST_THROW(stepper->setMomento(momento.ptr(),Teuchos::null,neSolver), std::logic_error);
    TEST_NOTHROW(stepper->setMomento(momento.ptr(),model,neSolver));
  }
  {
    // Initialize a valid non-const momento:
    RCP<BackwardEulerStepperMomento<double> > momento;
    {
      RCP<SinCosModel> model = sinCosModel(true);
      Thyra::ModelEvaluatorBase::InArgs<double> model_ic = model->getNominalValues();
      RCP<Thyra::NonlinearSolverBase<double> > neSolver = timeStepNonlinearSolver<double>();
      RCP<BackwardEulerStepper<double> > stepper = backwardEulerStepper<double>(model,neSolver);
      stepper->setInitialCondition(model_ic);
      double dt = 0.1;
      double dt_taken = stepper->takeStep(dt,STEP_TYPE_FIXED);
      TEST_ASSERT( dt == dt_taken );
      RCP<const MomentoBase<double> > m = stepper->getMomento();
      RCP<const BackwardEulerStepperMomento<double> > beMomento = Teuchos::rcp_dynamic_cast<const BackwardEulerStepperMomento<double> >(m,true);
      momento = Teuchos::rcp_dynamic_cast<BackwardEulerStepperMomento<double> >(beMomento->clone(),true);
    }
    {
      // Check if isInitialized_ == true, but 
      // model_ = null or solver = null or haveInitialCondition_ = false or interpolator = null
      RCP<BackwardEulerStepperMomento<double> > m = Teuchos::rcp_dynamic_cast<BackwardEulerStepperMomento<double> >(momento->clone(),true);
      RCP<BackwardEulerStepper<double> > stepper = backwardEulerStepper<double>();
      RCP<SinCosModel> model = sinCosModel(true);
      RCP<Thyra::NonlinearSolverBase<double> > neSolver = timeStepNonlinearSolver<double>();
      TEST_NOTHROW(stepper->setMomento(m.ptr(),model,neSolver));
      TEST_THROW(stepper->setMomento(m.ptr(),Teuchos::null,neSolver),std::logic_error);
      TEST_NOTHROW(stepper->setMomento(m.ptr(),model,neSolver));
      TEST_THROW(stepper->setMomento(m.ptr(),model,Teuchos::null),std::logic_error);
      TEST_NOTHROW(stepper->setMomento(m.ptr(),model,neSolver));
      m->set_haveInitialCondition(false);
      TEST_THROW(stepper->setMomento(m.ptr(),model,neSolver),std::logic_error);
      m->set_haveInitialCondition(true);
      TEST_NOTHROW(stepper->setMomento(m.ptr(),model,neSolver));
      RCP<InterpolatorBase<double> > interp = m->get_interpolator();
      m->set_interpolator(Teuchos::null);
      TEST_THROW(stepper->setMomento(m.ptr(),model,neSolver),std::logic_error);
      m->set_interpolator(interp);
      TEST_NOTHROW(stepper->setMomento(m.ptr(),model,neSolver));
    }
    {
      // Check if haveInitialCondition_ == true, but 
      // t_ == nan or 
      // t_old_ == nan or 
      // scaled_x_old == null or 
      // x_dot_old == null or
      // x == null or 
      // x_dot == null
      typedef Teuchos::ScalarTraits<double> ST;
      RCP<BackwardEulerStepperMomento<double> > m = Teuchos::rcp_dynamic_cast<BackwardEulerStepperMomento<double> >(momento->clone(),true);
      RCP<BackwardEulerStepper<double> > stepper = backwardEulerStepper<double>();
      RCP<SinCosModel> model = sinCosModel(true);
      RCP<Thyra::NonlinearSolverBase<double> > neSolver = timeStepNonlinearSolver<double>();
      TEST_NOTHROW(stepper->setMomento(m.ptr(),model,neSolver));

      double t = m->get_t();
      m->set_t(ST::nan());
      TEST_THROW(stepper->setMomento(m.ptr(),model,neSolver),std::logic_error);
      m->set_t(t);
      TEST_NOTHROW(stepper->setMomento(m.ptr(),model,neSolver));

      double t_old = m->get_t_old();
      m->set_t_old(ST::nan());
      TEST_THROW(stepper->setMomento(m.ptr(),model,neSolver),std::logic_error);
      m->set_t_old(t_old);
      TEST_NOTHROW(stepper->setMomento(m.ptr(),model,neSolver));

      RCP<VectorBase<double> > scaled_x_old = m->get_scaled_x_old();
      m->set_scaled_x_old(Teuchos::null);
      TEST_THROW(stepper->setMomento(m.ptr(),model,neSolver),std::logic_error);
      m->set_scaled_x_old(scaled_x_old);
      TEST_NOTHROW(stepper->setMomento(m.ptr(),model,neSolver));

      RCP<VectorBase<double> > x_dot_old = m->get_x_dot_old();
      m->set_x_dot_old(Teuchos::null);
      TEST_THROW(stepper->setMomento(m.ptr(),model,neSolver),std::logic_error);
      m->set_x_dot_old(x_dot_old);
      TEST_NOTHROW(stepper->setMomento(m.ptr(),model,neSolver));

      RCP<VectorBase<double> > x = m->get_x();
      m->set_x(Teuchos::null);
      TEST_THROW(stepper->setMomento(m.ptr(),model,neSolver),std::logic_error);
      m->set_x(x);
      TEST_NOTHROW(stepper->setMomento(m.ptr(),model,neSolver));

      RCP<VectorBase<double> > x_dot = m->get_x_dot();
      m->set_x_dot(Teuchos::null);
      TEST_THROW(stepper->setMomento(m.ptr(),model,neSolver),std::logic_error);
      m->set_x_dot(x_dot);
      TEST_NOTHROW(stepper->setMomento(m.ptr(),model,neSolver));
    }
    {
      // Check that if numSteps_ > 0 then isInitialized_ == true and haveInitialCondition_ == true
      RCP<BackwardEulerStepperMomento<double> > m = Teuchos::rcp_dynamic_cast<BackwardEulerStepperMomento<double> >(momento->clone(),true);
      TEST_ASSERT( m->get_numSteps() == 1 );
      RCP<BackwardEulerStepper<double> > stepper = backwardEulerStepper<double>();
      RCP<SinCosModel> model = sinCosModel(true);
      RCP<Thyra::NonlinearSolverBase<double> > neSolver = timeStepNonlinearSolver<double>();
      TEST_NOTHROW(stepper->setMomento(m.ptr(),model,neSolver));

      m->set_isInitialized(false);
      TEST_THROW(stepper->setMomento(m.ptr(),model,neSolver),std::logic_error);
      m->set_isInitialized(true);
      TEST_NOTHROW(stepper->setMomento(m.ptr(),model,neSolver));
    }
  }
}

} // namespace Rythmos

