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
#include "Teuchos_ParameterList.hpp"
#include "Thyra_DetachedVectorView.hpp"

namespace Rythmos {

using Thyra::VectorBase;
using Thyra::VectorSpaceBase;
using Teuchos::is_null;

TEUCHOS_UNIT_TEST( Rythmos_ExplicitRKStepper, create ) {
  RCP<SinCosModel> model = sinCosModel(false);
  RCP<ExplicitRKStepper<double> > stepper = explicitRKStepper<double>(model);
  TEST_EQUALITY_CONST( is_null(stepper), false );
}

TEUCHOS_UNIT_TEST( Rythmos_ExplicitRKStepper, assertValidModel ) {
  {
    // explicit model, OK
    RCP<SinCosModel> model = sinCosModel(false);
    RCP<ExplicitRKStepper<double> > stepper = explicitRKStepper<double>();
    TEST_NOTHROW( stepper->setModel(model) );
  }
//  {
    // implicit model, throw
//    RCP<SinCosModel> model = sinCosModel(true);
//    RCP<ExplicitRKStepper<double> > stepper = explicitRKStepper<double>();
//    TEST_THROW( stepper->setModel(model), std::logic_error );
//  }
}

TEUCHOS_UNIT_TEST( Rythmos_ExplicitRKStepper, setgetRKButcherTableau ) {
  RCP<SinCosModel> model = sinCosModel(false);
  RCP<RKButcherTableauBase<double> > rkbt = createRKBT<double>("Explicit 4 Stage");
  RCP<ExplicitRKStepper<double> > stepper = explicitRKStepper<double>(model,rkbt);
  TEST_EQUALITY_CONST( is_null(stepper), false );
  RCP<const RKButcherTableauBase<double> > rkbt_out = stepper->getRKButcherTableau();
  TEST_EQUALITY_CONST( *rkbt == *rkbt_out, true );
}

TEUCHOS_UNIT_TEST( Rythmos_ExplicitRKStepper, invalidRKBT ) {
  RCP<ExplicitRKStepper<double> > stepper = explicitRKStepper<double>();
  RCP<RKButcherTableauBase<double> > rkbt;
  TEST_THROW( stepper->setRKButcherTableau(rkbt), std::logic_error ); // empty RKBT
}

TEUCHOS_UNIT_TEST( Rythmos_ExplicitRKStepper, getTimeRange ) {
  {
    RCP<SinCosModel> model = sinCosModel(false);
    RCP<ExplicitRKStepper<double> > stepper = explicitRKStepper<double>(model);
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    stepper->setInitialCondition(ic);
    TimeRange<double> tr = stepper->getTimeRange();
    TEST_EQUALITY_CONST( tr.isValid(), true );
    TEST_EQUALITY_CONST( tr.lower(), 0.0 );
    TEST_EQUALITY_CONST( tr.upper(), 0.0 );
    TEST_EQUALITY_CONST( tr.length(), 0.0 );
  }
  {
    RCP<SinCosModel> model = sinCosModel(false);
    RCP<ExplicitRKStepper<double> > stepper = explicitRKStepper<double>();
    stepper->setModel(model);
    TimeRange<double> tr;
    TEST_NOTHROW( tr = stepper->getTimeRange() );
    TEST_EQUALITY_CONST( tr.isValid(), false );
  }
} 

TEUCHOS_UNIT_TEST( Rythmos_ExplicitRKStepper, noRKBT ) {
  RCP<SinCosModel> model = sinCosModel(false);
  RCP<ExplicitRKStepper<double> > stepper = explicitRKStepper<double>();
  stepper->setModel(model);
  double dt = 1.0;
  double step_taken = 0.0;
  TEST_THROW( step_taken = stepper->takeStep(dt, STEP_TYPE_FIXED), std::logic_error ); // no RKBT defined
}

TEUCHOS_UNIT_TEST( Rythmos_ExplicitRKStepper, invalidTakeStep ) {
  RCP<SinCosModel> model = sinCosModel(false);
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<ExplicitRKStepper<double> > stepper = explicitRKStepper<double>();
  stepper->setModel(model);
  stepper->setInitialCondition(ic);
  stepper->setRKButcherTableau(createRKBT<double>("Explicit 4 Stage"));
  double dt;
#ifdef RYTHMOS_DEBUG
  TEST_THROW(dt = stepper->takeStep(0.1,STEP_TYPE_VARIABLE), std::logic_error);
#else
  dt = stepper->takeStep(0.1,STEP_TYPE_VARIABLE);
  TEST_EQUALITY_CONST( dt, -1.0 );
#endif // RYTHMOS_DEBUG
}

TEUCHOS_UNIT_TEST( Rythmos_ExplicitRKStepper, basePoint ) {
  RCP<SinCosModel> model = sinCosModel(false);
  {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set("Accept model parameters",true);
    model->setParameterList(pl);
  }
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  // t_ic
  double t_ic = 1.0; // not used
  // x_ic
  RCP<VectorBase<double> > x_ic = Thyra::createMember(*model->get_x_space());
  {
    Thyra::DetachedVectorView<double> x_ic_view( *x_ic );
    x_ic_view[0] = 5.0;
    x_ic_view[1] = 6.0;
  }
  // parameter 0 ic
  RCP<VectorBase<double> > p_ic = Thyra::createMember(*model->get_p_space(0));
  {
    Thyra::DetachedVectorView<double> p_ic_view( *p_ic );
    p_ic_view[0] = 2.0; // a
    p_ic_view[1] = 3.0; // f 
    p_ic_view[2] = 4.0; // L
  }
  ic.set_p(0,p_ic); 
  ic.set_x(x_ic);
  ic.set_t(t_ic);
  RCP<ExplicitRKStepper<double> > stepper = explicitRKStepper<double>();
  stepper->setModel(model);
  stepper->setInitialCondition(ic);
  stepper->setRKButcherTableau(createRKBT<double>("Forward Euler"));
  double dt = 0.2;
  double dt_taken;
  dt_taken = stepper->takeStep(dt,STEP_TYPE_FIXED);
  TEST_EQUALITY_CONST( dt_taken, 0.2 );
  const StepStatus<double> status = stepper->getStepStatus();
  TEST_ASSERT( !is_null(status.solution) );
  double tol = 1.0e-10;
  {
    Thyra::ConstDetachedVectorView<double> x_new_view( *(status.solution) );
    TEST_FLOATING_EQUALITY( x_new_view[0], 5.0 + 0.2*(6.0), tol );
    TEST_FLOATING_EQUALITY( x_new_view[1], 6.0 + 0.2*( (3.0/4.0)*(3.0/4.0)*(2.0-5.0) ), tol );
  }
}

} // namespace Rythmos

