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

#include "Rythmos_StepperHelpers.hpp"
#include "Rythmos_StepperBuilder.hpp"
#include "../SinCos/SinCosModel.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_ModelEvaluator.hpp"

namespace Rythmos {

// Functions to test:
// assertValidModel
TEUCHOS_UNIT_TEST( Rythmos_StepperHelpers, assertValidModel ) {
  RCP<SinCosModel> model = sinCosModel(false);
  RCP<StepperBuilder<double> > builder = stepperBuilder<double>();
  RCP<StepperBase<double> > stepper = builder->create("Implicit BDF");
  TEST_EQUALITY_CONST( is_null(stepper), false );
  // implicit Stepper and explicit model, throws
  TEST_THROW(
      assertValidModel( *stepper, *model ),
      std::logic_error
      ); 
  stepper = builder->create("Explicit RK");
  // explicit stepper and explicit model, OK
  TEST_NOTHROW(
      assertValidModel( *stepper, *model )
      );
  model = sinCosModel(true);
  // explicit stepper and implicit model, throws
//  TEST_THROW( 
//      assertValidModel( *stepper, *model ),
//      std::logic_error
//      )
  stepper = builder->create("Implicit RK");
  // implicit stepper and implicit model, OK
  TEST_NOTHROW(
      assertValidModel( *stepper, *model )
      );

}

// setDefaultInitialConditionFromNominalValues
// restart

TEUCHOS_UNIT_TEST( Rythmos_StepperHelpers, eval_model_explicit ) {
  RCP<Thyra::ModelEvaluator<double> > model;
  {
    RCP<SinCosModel> scModel = sinCosModel(false);
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set("Accept model parameters",true);
    scModel->setParameterList(pl);
    model = scModel;
  }
  Thyra::ModelEvaluatorBase::InArgs<double> basePoint = model->getNominalValues();
  {
    // Base Point values:
    // time 
    double t_ic = 1.5;
    // x
    RCP<VectorBase<double> > x_ic = Thyra::createMember(*model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_ic_view( *x_ic );
      x_ic_view[0] = 5.0;
      x_ic_view[1] = 6.0;
    }
    // Parameter 0
    RCP<VectorBase<double> > p_ic = Thyra::createMember(*model->get_p_space(0));
    {
      Thyra::DetachedVectorView<double> p_ic_view( *p_ic );
      p_ic_view[0] = 2.0; // a
      p_ic_view[1] = 3.0; // f
      p_ic_view[2] = 4.0; // L
    }
    basePoint.set_t(t_ic);
    basePoint.set_x(x_ic);
    basePoint.set_p(0,p_ic);
  }
  // Evaluation with eval_model_explicit
  double t = 2.0;
  RCP<VectorBase<double> > x = Thyra::createMember(*model->get_x_space());
  {
    Thyra::DetachedVectorView<double> x_view( *x );
    x_view[0] = 7.0;
    x_view[1] = 8.0;
  }
  RCP<VectorBase<double> > f_out = Thyra::createMember(*model->get_f_space());
  eval_model_explicit<double>(*model,basePoint,*x,t,Teuchos::outArg(*f_out));
  // Verify that our new t and x were used
  // Verify that the parameters in the basePoint were used rather than the defaults.
  // I can't test t with this model because it doesn't use it (TODO)
  double a = 2.0;
  double f = 3.0;
  double L = 4.0;
  double x0 = 7.0;
  double x1 = 8.0;
  double tol = 1.0e-10;
  {
    Thyra::ConstDetachedVectorView<double> f_out_view( *f_out );
    TEST_EQUALITY( f_out_view[0], x1 );
    TEST_FLOATING_EQUALITY( f_out_view[1], (f/L)*(f/L)*(a-x0), tol );
  }
}


} // namespace Rythmos 


