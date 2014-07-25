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

#include "Rythmos_StepperHelpers.hpp"
#include "Rythmos_StepperBuilder.hpp"
#include "../SinCos/SinCosModel.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Rythmos_UnitTestHelpers.hpp"
#include "Rythmos_LinearInterpolator.hpp"

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

TEUCHOS_UNIT_TEST( Rythmos_StepperHelpers, defaultGetPoints ) {
  using Teuchos::constOptInArg;
  using Teuchos::null;
  using Teuchos::ptr;
  double t_old = 1.0;
  RCP<VectorBase<double> > x_old = createDefaultVector(1,2.0);
  RCP<VectorBase<double> > xdot_old = createDefaultVector(1,3.0);
  double t = 2.0;
  RCP<VectorBase<double> > x = createDefaultVector(1,4.0);
  RCP<VectorBase<double> > xdot = createDefaultVector(1,5.0);

  Array<double> time_vec;
  time_vec.push_back(t_old);
  Array<RCP<const VectorBase<double> > > x_vec;
  Array<RCP<const VectorBase<double> > > xdot_vec;
  Array<double> accuracy_vec;

  Teuchos::outArg(*x_old);
  // Get the first point
  defaultGetPoints<double>(
      t_old, constOptInArg(*x_old), constOptInArg(*xdot_old),
      t, constOptInArg(*x), constOptInArg(*xdot),
      time_vec, ptr(&x_vec), ptr(&xdot_vec), ptr(&accuracy_vec),
      null
      );
  TEST_ASSERT( x_vec.length() == 1 );
  TEST_ASSERT( xdot_vec.length() == 1 );
  TEST_ASSERT( accuracy_vec.length() == 1 );
  TEST_ASSERT( x_vec[0].get() != x_old.get() );
  TEST_ASSERT( xdot_vec[0].get() != xdot_old.get() );
  {
    Thyra::ConstDetachedVectorView<double> x_view( x_vec[0] );
    Thyra::ConstDetachedVectorView<double> xdot_view( xdot_vec[0] );
    TEST_EQUALITY_CONST( x_view[0], 2.0 );
    TEST_EQUALITY_CONST( xdot_view[0], 3.0 );
    TEST_EQUALITY_CONST( accuracy_vec[0], 0.0 );
  }
  // Get the second point
  time_vec.clear();
  time_vec.push_back(t);
  defaultGetPoints<double>(
      t_old, constOptInArg(*x_old), constOptInArg(*xdot_old),
      t, constOptInArg(*x), constOptInArg(*xdot),
      time_vec, ptr(&x_vec), ptr(&xdot_vec), ptr(&accuracy_vec),
      null
      );
  TEST_ASSERT( x_vec.length() == 1 );
  TEST_ASSERT( xdot_vec.length() == 1 );
  TEST_ASSERT( accuracy_vec.length() == 1 );
  TEST_ASSERT( x_vec[0].get() != x.get() );
  TEST_ASSERT( xdot_vec[0].get() != xdot.get() );
  {
    Thyra::ConstDetachedVectorView<double> x_view( x_vec[0] );
    Thyra::ConstDetachedVectorView<double> xdot_view( xdot_vec[0] );
    TEST_EQUALITY_CONST( x_view[0], 4.0 );
    TEST_EQUALITY_CONST( xdot_view[0], 5.0 );
    TEST_EQUALITY_CONST( accuracy_vec[0], 0.0 );
  }
  // Get both points
  time_vec.clear();
  time_vec.push_back(t_old);
  time_vec.push_back(t);
  defaultGetPoints<double>(
      t_old, constOptInArg(*x_old), constOptInArg(*xdot_old),
      t, constOptInArg(*x), constOptInArg(*xdot),
      time_vec, ptr(&x_vec), ptr(&xdot_vec), ptr(&accuracy_vec),
      null
      );
  TEST_ASSERT( x_vec.length() == 2 );
  TEST_ASSERT( xdot_vec.length() == 2 );
  TEST_ASSERT( accuracy_vec.length() == 2 );
  TEST_ASSERT( x_vec[0].get() != x_old.get() );
  TEST_ASSERT( xdot_vec[0].get() != xdot_old.get() );
  TEST_ASSERT( x_vec[1].get() != x.get() );
  TEST_ASSERT( xdot_vec[1].get() != xdot.get() );
  {
    Thyra::ConstDetachedVectorView<double> x_view_0( x_vec[0] );
    Thyra::ConstDetachedVectorView<double> xdot_view_0( xdot_vec[0] );
    Thyra::ConstDetachedVectorView<double> x_view_1( x_vec[1] );
    Thyra::ConstDetachedVectorView<double> xdot_view_1( xdot_vec[1] );
    TEST_EQUALITY_CONST( x_view_0[0], 2.0 );
    TEST_EQUALITY_CONST( xdot_view_0[0], 3.0 );
    TEST_EQUALITY_CONST( x_view_1[0], 4.0 );
    TEST_EQUALITY_CONST( xdot_view_1[0], 5.0 );
    TEST_EQUALITY_CONST( accuracy_vec[0], 0.0 );
    TEST_EQUALITY_CONST( accuracy_vec[1], 0.0 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_StepperHelpers, defaultGetPoints_EdgeCases ) {
  using Teuchos::constOptInArg;
  using Teuchos::null;
  {
    // t_old > t : throw
    double t_old = 2.0;
    RCP<VectorBase<double> > x_old;
    RCP<VectorBase<double> > xdot_old;
    RCP<VectorBase<double> > x;
    RCP<VectorBase<double> > xdot;
    double t = 1.0;
    Array<double> time_vec;
    time_vec.push_back(t);
    Array<RCP<const VectorBase<double> > > x_vec;
    Array<RCP<const VectorBase<double> > > xdot_vec;
    Array<double> accuracy_vec;
    TEST_THROW(
      defaultGetPoints<double>(
          t_old,null,null,
          t,null,null,
          time_vec,ptr(&x_vec),ptr(&xdot_vec),ptr(&accuracy_vec),
          null
          ),
      std::logic_error
      );
  }
  // t_old == t : uses x_old,xdot_old TODO
  //
  // x_old = null : okay TODO
  // x = null : okay TODO
  // xdot_old = null : okay TODO
  // xdot = null : okay TODO
  // 
  // time_vec = empty : okay TODO
  {
    double t_old = 1.0;
    RCP<VectorBase<double> > x_old;
    RCP<VectorBase<double> > xdot_old;
    double t = 2.0;
    RCP<VectorBase<double> > x;
    RCP<VectorBase<double> > x_dot;
    Array<double> time_vec;
    Array<RCP<const VectorBase<double> > > x_vec;
    Array<RCP<const VectorBase<double> > > xdot_vec;
    Array<double> accuracy_vec;
    TEST_NOTHROW(
      defaultGetPoints<double>(
          t_old,null,null,
          t,null,null,
          time_vec,ptr(&x_vec),ptr(&xdot_vec),ptr(&accuracy_vec),
          null
          )
      );
  }
  {
    // time_vec = unsorted : throw
    double t_old = 1.0;
    RCP<VectorBase<double> > x_old;
    RCP<VectorBase<double> > xdot_old;
    double t = 2.0;
    RCP<VectorBase<double> > x;
    RCP<VectorBase<double> > x_dot;
    Array<double> time_vec;
    time_vec.push_back(t);
    time_vec.push_back(t_old);
    Array<RCP<const VectorBase<double> > > x_vec;
    Array<RCP<const VectorBase<double> > > xdot_vec;
    Array<double> accuracy_vec;
    TEST_THROW(
      defaultGetPoints<double>(
          t_old,null,null,
          t,null,null,
          time_vec,ptr(&x_vec),ptr(&xdot_vec),ptr(&accuracy_vec),
          null
          ),
      std::logic_error
      );
  }
  {
    // x_vec = null : okay
    // xdot_vec = null : okay
    // accuracy_vec = null : okay
    double t_old = 1.0;
    RCP<VectorBase<double> > x_old = createDefaultVector(1,2.0);
    RCP<VectorBase<double> > xdot_old = createDefaultVector(1,3.0);
    double t = 2.0;
    RCP<VectorBase<double> > x = createDefaultVector(1,4.0);
    RCP<VectorBase<double> > xdot = createDefaultVector(1,5.0);
    Array<double> time_vec;
    time_vec.push_back(t_old);
    Array<RCP<const VectorBase<double> > > x_vec;
    Array<RCP<const VectorBase<double> > > xdot_vec;
    Array<double> accuracy_vec;
    TEST_NOTHROW(
        defaultGetPoints<double>(
          t_old,constOptInArg(*x_old),constOptInArg(*xdot_old),
          t,constOptInArg(*x),constOptInArg(*xdot),
          time_vec,null,null,null,
          null
          )
        );
    TEST_NOTHROW(
        defaultGetPoints<double>(
          t_old,constOptInArg(*x_old),constOptInArg(*xdot_old),
          t,constOptInArg(*x),constOptInArg(*xdot),
          time_vec,ptr(&x_vec),null,null,
          null
          )
        );
    TEST_NOTHROW(
        defaultGetPoints<double>(
          t_old,constOptInArg(*x_old),constOptInArg(*xdot_old),
          t,constOptInArg(*x),constOptInArg(*xdot),
          time_vec,null,ptr(&xdot_vec),null,
          null
          )
        );
    TEST_NOTHROW(
        defaultGetPoints<double>(
          t_old,constOptInArg(*x_old),constOptInArg(*xdot_old),
          t,constOptInArg(*x),constOptInArg(*xdot),
          time_vec,null,null,ptr(&accuracy_vec),
          null
          )
        );
  }
  {
    // time_vec contains points before t_old  : throw
    double t_old = 1.0;
    RCP<VectorBase<double> > x_old;
    RCP<VectorBase<double> > xdot_old;
    double t = 2.0;
    RCP<VectorBase<double> > x;
    RCP<VectorBase<double> > x_dot;
    Array<double> time_vec;
    time_vec.push_back(0.0);
    time_vec.push_back(1.0);
    Array<RCP<const VectorBase<double> > > x_vec;
    Array<RCP<const VectorBase<double> > > xdot_vec;
    Array<double> accuracy_vec;
    TEST_THROW(
      defaultGetPoints<double>(
          t_old,null,null,
          t,null,null,
          time_vec,ptr(&x_vec),ptr(&xdot_vec),ptr(&accuracy_vec),
          null
          ),
      std::logic_error
      );
    // time_vec contains points after t : throw
    time_vec.clear();
    time_vec.push_back(2.0);
    time_vec.push_back(3.0);
    TEST_THROW(
      defaultGetPoints<double>(
          t_old,null,null,
          t,null,null,
          time_vec,ptr(&x_vec),ptr(&xdot_vec),ptr(&accuracy_vec),
          null
          ),
      std::logic_error
      );
    // time_vec contains points in (t_old,t) : throw
    time_vec.clear();
    time_vec.push_back(1.0);
    time_vec.push_back(1.5);
    time_vec.push_back(2.0);
    TEST_THROW(
      defaultGetPoints<double>(
          t_old,null,null,
          t,null,null,
          time_vec,ptr(&x_vec),ptr(&xdot_vec),ptr(&accuracy_vec),
          null
          ),
      std::logic_error
      );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_StepperHelpers, defaultGetPoints_interpolator ) {
  using Teuchos::outArg;
  using Teuchos::constOptInArg;
  using Teuchos::ptr;
  double t_old = 1.0;
  RCP<VectorBase<double> > x_old = createDefaultVector(1,2.0);
  RCP<VectorBase<double> > xdot_old = createDefaultVector(1,3.0);
  double t = 2.0;
  RCP<VectorBase<double> > x = createDefaultVector(1,4.0);
  RCP<VectorBase<double> > xdot = createDefaultVector(1,5.0);

  Array<double> time_vec;
  time_vec.push_back(t_old);
  time_vec.push_back(1.5);
  time_vec.push_back(t);
  Array<RCP<const VectorBase<double> > > x_vec;
  Array<RCP<const VectorBase<double> > > xdot_vec;
  Array<double> accuracy_vec;

  RCP<InterpolatorBase<double> > interp = linearInterpolator<double>();

  // Get the boundaries plus an interior point
  defaultGetPoints<double>(
      t_old, constOptInArg(*x_old), constOptInArg(*xdot_old),
      t, constOptInArg(*x), constOptInArg(*xdot),
      time_vec, ptr(&x_vec), ptr(&xdot_vec), ptr(&accuracy_vec),
      outArg(*interp)
      );
  TEST_ASSERT( x_vec.length() == 3 );
  TEST_ASSERT( xdot_vec.length() == 3 );
  TEST_ASSERT( accuracy_vec.length() == 3 );
  {
    Thyra::ConstDetachedVectorView<double> x_view( x_vec[0] );
    Thyra::ConstDetachedVectorView<double> xdot_view( xdot_vec[0] );
    TEST_EQUALITY_CONST( x_view[0], 2.0 );
    TEST_EQUALITY_CONST( xdot_view[0], 3.0 );
    TEST_EQUALITY_CONST( accuracy_vec[0], 0.0 );
  }
  {
    Thyra::ConstDetachedVectorView<double> x_view( x_vec[1] );
    Thyra::ConstDetachedVectorView<double> xdot_view( xdot_vec[1] );
    TEST_EQUALITY_CONST( x_view[0], 3.0 );
    TEST_EQUALITY_CONST( xdot_view[0], 4.0 );
    TEST_EQUALITY_CONST( accuracy_vec[1], 1.0 );
  }
  {
    Thyra::ConstDetachedVectorView<double> x_view( x_vec[2] );
    Thyra::ConstDetachedVectorView<double> xdot_view( xdot_vec[2] );
    TEST_EQUALITY_CONST( x_view[0], 4.0 );
    TEST_EQUALITY_CONST( xdot_view[0], 5.0 );
    TEST_EQUALITY_CONST( accuracy_vec[2], 0.0 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_StepperHelpers, defaultGetPoints_interpolator_EdgeCases ) {
  using Teuchos::outArg;
  using Teuchos::constOptInArg;
  using Teuchos::null;
  using Teuchos::ptr;
  {
    // xdot = null
    double t_old = 1.0;
    RCP<VectorBase<double> > x_old = createDefaultVector(1,2.0);
    RCP<VectorBase<double> > xdot_old; 
    double t = 2.0;
    RCP<VectorBase<double> > x = createDefaultVector(1,4.0);
    RCP<VectorBase<double> > xdot;

    Array<double> time_vec;
    time_vec.push_back(1.5);
    Array<RCP<const VectorBase<double> > > x_vec;
    Array<RCP<const VectorBase<double> > > xdot_vec;
    Array<double> accuracy_vec;

    RCP<InterpolatorBase<double> > interp = linearInterpolator<double>();

    // Get the boundaries plus an interior point
    defaultGetPoints<double>(
        t_old, constOptInArg(*x_old), null,
        t, constOptInArg(*x), null,
        time_vec, ptr(&x_vec), ptr(&xdot_vec), ptr(&accuracy_vec),
        outArg(*interp)
        );
    TEST_ASSERT( x_vec.length() == 1 );
    TEST_ASSERT( xdot_vec.length() == 1 );
    TEST_ASSERT( accuracy_vec.length() == 1 );
    TEST_ASSERT( is_null(xdot_vec[0]) );
    {
      Thyra::ConstDetachedVectorView<double> x_view( x_vec[0] );
      TEST_EQUALITY_CONST( x_view[0], 3.0 );
      TEST_EQUALITY_CONST( accuracy_vec[0], 1.0 );
    }
  }
  {
    // x = null
    double t_old = 1.0;
    RCP<VectorBase<double> > x_old;
    RCP<VectorBase<double> > xdot_old = createDefaultVector(1,3.0);
    double t = 2.0;
    RCP<VectorBase<double> > x;
    RCP<VectorBase<double> > xdot = createDefaultVector(1,5.0);

    Array<double> time_vec;
    time_vec.push_back(1.5);
    Array<RCP<const VectorBase<double> > > x_vec;
    Array<RCP<const VectorBase<double> > > xdot_vec;
    Array<double> accuracy_vec;

    RCP<InterpolatorBase<double> > interp = linearInterpolator<double>();

    // Get the boundaries plus an interior point
    defaultGetPoints<double>(
        t_old, null, constOptInArg(*xdot_old),
        t, null, constOptInArg(*xdot),
        time_vec, ptr(&x_vec), ptr(&xdot_vec), ptr(&accuracy_vec),
        outArg(*interp)
        );
    TEST_ASSERT( x_vec.length() == 1 );
    TEST_ASSERT( xdot_vec.length() == 1 );
    TEST_ASSERT( accuracy_vec.length() == 1 );
    TEST_ASSERT( is_null(x_vec[0]) );
    {
      Thyra::ConstDetachedVectorView<double> xdot_view( xdot_vec[0] );
      TEST_EQUALITY_CONST( xdot_view[0], 4.0 );
      TEST_EQUALITY_CONST( accuracy_vec[0], 1.0 );
    }
  }
  {
    // x = null and xdot = null
    double t_old = 1.0;
    RCP<VectorBase<double> > x_old;
    RCP<VectorBase<double> > xdot_old;
    double t = 2.0;
    RCP<VectorBase<double> > x;
    RCP<VectorBase<double> > xdot;

    Array<double> time_vec;
    time_vec.push_back(1.5);
    Array<RCP<const VectorBase<double> > > x_vec;
    Array<RCP<const VectorBase<double> > > xdot_vec;
    Array<double> accuracy_vec;

    RCP<InterpolatorBase<double> > interp = linearInterpolator<double>();

    // Get the boundaries plus an interior point
    defaultGetPoints<double>(
        t_old, null, null,
        t, null, null,
        time_vec, ptr(&x_vec), ptr(&xdot_vec), ptr(&accuracy_vec),
        outArg(*interp)
        );
    TEST_ASSERT( x_vec.length() == 1 );
    TEST_ASSERT( xdot_vec.length() == 1 );
    TEST_ASSERT( accuracy_vec.length() == 1 );
    TEST_ASSERT( is_null(x_vec[0]) );
    TEST_ASSERT( is_null(xdot_vec[0]) );
    TEST_EQUALITY_CONST( accuracy_vec[0], 1.0 );
  }
}

} // namespace Rythmos 


