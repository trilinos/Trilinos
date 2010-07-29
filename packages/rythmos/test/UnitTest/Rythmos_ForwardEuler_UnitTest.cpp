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
#include "Rythmos_ForwardEulerStepper.hpp"
#include "../SinCos/SinCosModel.hpp"
#include "Rythmos_UnitTestModels.hpp"
#include "Thyra_TestingTools.hpp"

namespace Rythmos {

TEUCHOS_UNIT_TEST( Rythmos_ForwardEulerStepper, clone ) {
  RCP<ForwardEulerStepper<double> > stepper = forwardEulerStepper<double>();
  TEST_ASSERT( stepper->supportsCloning() );
  RCP<StepperBase<double> > newStepper = stepper->cloneStepperAlgorithm();
  TEST_ASSERT( !Teuchos::is_null(newStepper) );
  {
    RCP<ForwardEulerStepper<double> > feStepper = Teuchos::rcp_dynamic_cast<ForwardEulerStepper<double> >(newStepper,false);
    TEST_ASSERT(!is_null(feStepper));
  }
}

TEUCHOS_UNIT_TEST( Rythmos_ForwardEulerStepper, momento_create ) {
  RCP<const MomentoBase<double> > momento;
  {
    RCP<ForwardEulerStepper<double> > stepper = forwardEulerStepper<double>();
    RCP<SinCosModel> model = sinCosModel(false);
    stepper->setModel(model);
    momento = stepper->getMomento();
    TEST_ASSERT( !is_null(momento) );
  }
  {
    RCP<SinCosModel> model = sinCosModel(false);
    RCP<ForwardEulerStepper<double> > stepper = forwardEulerStepper<double>();
    stepper->setMomento(momento.ptr());
    momento = Teuchos::null;
    RCP<const Thyra::ModelEvaluator<double> > modelOut = stepper->getModel();
    TEST_ASSERT( !is_null(modelOut) );
    RCP<const SinCosModel> scModel = Teuchos::rcp_dynamic_cast<const SinCosModel>(modelOut,false);
    TEST_ASSERT( !is_null(scModel) );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_ForwardEulerStepper, restart ) {
  Array<RCP<StateSerializerStrategy<double> > > sss_array;
  sss_array.push_back(rcp(new XMLStateSerializerStrategy<double>()));
  sss_array.push_back(rcp(new BinaryStateSerializerStrategy<double>()));
  int N = Teuchos::as<int>(sss_array.size());
  for (int n=0 ; n<N; ++n) {
    StateSerializerStrategy<double>& sss = *sss_array[n];
    std::string data;
    // place to store solution at time = 1.0
    RCP<const VectorBase<double> > x_norestart; 
    // Create Forward Euler stepper
    // Step to t = 0.5, pull out the momento, then step the rest of the way to t = 1.0.
    {
      RCP<SinCosModel> model = sinCosModel(false);
      Thyra::ModelEvaluatorBase::InArgs<double> model_ic = model->getNominalValues();
      RCP<ForwardEulerStepper<double> > stepper = forwardEulerStepper<double>(model);
      stepper->setInitialCondition(model_ic);
      double dt = 0.1;
      // Step to t=0.5
      for (int i=0 ; i<5 ; ++i) {
        double dt_taken = stepper->takeStep(dt,STEP_TYPE_FIXED);
        TEST_ASSERT( dt_taken == dt );
      }
      // Pull out the momento
      RCP<const MomentoBase<double> > stepper_momento = stepper->getMomento();
      TEST_ASSERT( !Teuchos::is_null(stepper_momento) );
      std::ostringstream oStream;
      stepper_momento->serialize(sss,oStream);
      data = oStream.str();
      stepper_momento = Teuchos::null;
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
        x_norestart = x_vec[0]->clone_v();
      }
    }
    out << "data = >>" << data << "<<" << std::endl;
    // Create a new Forward Euler Stepper
    {
      RCP<SinCosModel> model = sinCosModel(false);
      Thyra::ModelEvaluatorBase::InArgs<double> model_ic = model->getNominalValues();

      std::istringstream iStream(data);
      RCP<ForwardEulerStepperMomento<double> > stepper_momento = rcp(new ForwardEulerStepperMomento<double>());
      stepper_momento->set_model(model);
      stepper_momento->deSerialize(sss,iStream);
      RCP<ForwardEulerStepper<double> > stepper = forwardEulerStepper<double>();
      RCP<Thyra::ModelEvaluatorBase::InArgs<double> > model_ic_ptr = Teuchos::rcp(new Thyra::ModelEvaluatorBase::InArgs<double>(model_ic));
      stepper_momento->set_basePoint(model_ic_ptr);
      // Put the momento back into the stepper
      stepper->setMomento(stepper_momento.ptr());

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

        // Verify that x_restart == x_norestart
        RCP<VectorBase<double> > x_diff = createMember(x_restart->space());
        V_VmV( x_diff.ptr(), *x_norestart, *x_restart );
        double x_normDiff = norm(*x_diff);
        double tol = 1.0e-10;
        TEST_COMPARE( x_normDiff, <, tol );
      }
    }
  }
}

TEUCHOS_UNIT_TEST( Rythmos_ForwardEulerStepper, checkConsistentState ) {
  {
    RCP<SinCosModel> model = sinCosModel(false);
    Thyra::ModelEvaluatorBase::InArgs<double> model_ic = model->getNominalValues();
    RCP<ForwardEulerStepper<double> > stepper = forwardEulerStepper<double>(model);
    stepper->setInitialCondition(model_ic);
    double dt = 0.1;
    double dt_taken = stepper->takeStep(dt,STEP_TYPE_FIXED);
    TEST_ASSERT( dt == dt_taken );
    RCP<const MomentoBase<double> > momento = stepper->getMomento();
    RCP<ForwardEulerStepperMomento<double> > feMomento = Teuchos::rcp_dynamic_cast<ForwardEulerStepperMomento<double> >(momento->clone(),true);
    feMomento->set_model(Teuchos::null);
    feMomento->set_basePoint(Teuchos::null);
    TEST_THROW(stepper->setMomento(feMomento.ptr()), std::logic_error);
  }
  {
    // Initialize a valid non-const momento:
    RCP<ForwardEulerStepperMomento<double> > momento;
    {
      RCP<SinCosModel> model = sinCosModel(false);
      Thyra::ModelEvaluatorBase::InArgs<double> model_ic = model->getNominalValues();
      RCP<ForwardEulerStepper<double> > stepper = forwardEulerStepper<double>(model);
      stepper->setInitialCondition(model_ic);
      double dt = 0.1;
      double dt_taken = stepper->takeStep(dt,STEP_TYPE_FIXED);
      TEST_ASSERT( dt == dt_taken );
      momento = Teuchos::rcp_dynamic_cast<ForwardEulerStepperMomento<double> >(stepper->getMomento()->clone(),true);
    }
    {
      // Check if isInitialized_ == true, but model_ = null or residual_vector = null
      RCP<ForwardEulerStepperMomento<double> > m = Teuchos::rcp_dynamic_cast<ForwardEulerStepperMomento<double> >(momento->clone(),true);
      RCP<ForwardEulerStepper<double> > stepper = forwardEulerStepper<double>();
      RCP<const Thyra::ModelEvaluator<double> > model = m->get_model();
      TEST_NOTHROW(stepper->setMomento(m.ptr()));
      m->set_model(Teuchos::null);
      TEST_THROW(stepper->setMomento(m.ptr()),std::logic_error);
      m->set_model(model);
      TEST_NOTHROW(stepper->setMomento(m.ptr()));
      m->set_residual_vector(Teuchos::null);
      TEST_THROW(stepper->setMomento(m.ptr()),std::logic_error);
    }
    {
      // Check if haveInitialCondition_ == true, but t_ == nan or t_old_ == nan or solution_vector == null or solution_vector_old == null
      typedef Teuchos::ScalarTraits<double> ST;
      RCP<ForwardEulerStepperMomento<double> > m = Teuchos::rcp_dynamic_cast<ForwardEulerStepperMomento<double> >(momento->clone(),true);
      RCP<ForwardEulerStepper<double> > stepper = forwardEulerStepper<double>();
      TEST_NOTHROW(stepper->setMomento(m.ptr()));
      double t = m->get_t();
      m->set_t(ST::nan());
      TEST_THROW(stepper->setMomento(m.ptr()),std::logic_error);
      m->set_t(t);
      TEST_NOTHROW(stepper->setMomento(m.ptr()));
      double t_old = m->get_t_old();
      m->set_t_old(ST::nan());
      TEST_THROW(stepper->setMomento(m.ptr()),std::logic_error);
      m->set_t_old(t_old);
      TEST_NOTHROW(stepper->setMomento(m.ptr()));
      RCP<VectorBase<double> > sol_vec = m->get_solution_vector();
      m->set_solution_vector(Teuchos::null);
      TEST_THROW(stepper->setMomento(m.ptr()),std::logic_error);
      m->set_solution_vector(sol_vec);
      TEST_NOTHROW(stepper->setMomento(m.ptr()));
      RCP<VectorBase<double> > sol_vec_old = m->get_solution_vector_old();
      m->set_solution_vector_old(Teuchos::null);
      TEST_THROW(stepper->setMomento(m.ptr()),std::logic_error);
      m->set_solution_vector_old(sol_vec_old);
      TEST_NOTHROW(stepper->setMomento(m.ptr()));
    }
    {
      // Check that if numSteps_ > 0 then isInitialized_ == true and haveInitialCondition_ == true
      RCP<ForwardEulerStepperMomento<double> > m = Teuchos::rcp_dynamic_cast<ForwardEulerStepperMomento<double> >(momento->clone(),true);
      TEST_ASSERT( m->get_numSteps() == 1 );
      RCP<ForwardEulerStepper<double> > stepper = forwardEulerStepper<double>();
      TEST_NOTHROW(stepper->setMomento(m.ptr()));
      m->set_isInitialized(false);
      TEST_THROW(stepper->setMomento(m.ptr()),std::logic_error);
      m->set_isInitialized(true);
      TEST_NOTHROW(stepper->setMomento(m.ptr()));
      m->set_haveInitialCondition(false);
      TEST_THROW(stepper->setMomento(m.ptr()),std::logic_error);
      m->set_haveInitialCondition(true);
      TEST_NOTHROW(stepper->setMomento(m.ptr()));
    }
  }
}

} // namespace Rythmos

