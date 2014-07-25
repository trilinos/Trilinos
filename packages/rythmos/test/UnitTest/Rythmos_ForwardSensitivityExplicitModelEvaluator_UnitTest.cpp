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

#include "Rythmos_Types.hpp"
#include "Rythmos_UnitTestHelpers.hpp"
#include "Rythmos_ForwardSensitivityExplicitModelEvaluator.hpp"
#include "Rythmos_UnitTestModels.hpp"

#include "Rythmos_StepperBuilder.hpp"
#include "../SinCos/SinCosModel.hpp"

#include "Thyra_DetachedVectorView.hpp"

namespace Rythmos {


TEUCHOS_UNIT_TEST( Rythmos_ForwardSensitivityExplicitModelEvaluator, create ) { 
  RCP<ForwardSensitivityExplicitModelEvaluator<double> > sensModel =
    forwardSensitivityExplicitModelEvaluator<double>();
  TEST_ASSERT( !is_null(sensModel) );
}

TEUCHOS_UNIT_TEST( Rythmos_ForwardSensitivityExplicitModelEvaluator, args ) {
  RCP<ForwardSensitivityExplicitModelEvaluator<double> > model =
    forwardSensitivityExplicitModelEvaluator<double>();
  RCP<SinCosModel> innerModel = sinCosModel();
  {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set("Accept model parameters",true);
    pl->set("Implicit model formulation",false);
    innerModel->setParameterList(pl);
  }
  model->initializeStructure(innerModel, 0 );
  typedef Thyra::ModelEvaluatorBase MEB;
  {
    MEB::InArgs<double> inArgs = model->createInArgs();
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_x_dot), false );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_alpha), false );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_beta), true );
  }
  {
    MEB::OutArgs<double> outArgs = model->createOutArgs();
    TEST_EQUALITY_CONST( outArgs.supports(MEB::OUT_ARG_f), true );
    TEST_EQUALITY_CONST( outArgs.supports(MEB::OUT_ARG_W_op), false );
    TEST_EQUALITY_CONST( outArgs.supports(MEB::OUT_ARG_W), false );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_ForwardSensitivityExplicitModelEvaluator, spaces ) {
  RCP<ForwardSensitivityExplicitModelEvaluator<double> > model =
    forwardSensitivityExplicitModelEvaluator<double>();
  RCP<SinCosModel> innerModel = sinCosModel();
  {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set("Accept model parameters",true);
    pl->set("Implicit model formulation",false);
    innerModel->setParameterList(pl);
  }
  model->initializeStructure(innerModel, 0 );
  // x_space:
  {
    RCP<const Thyra::VectorSpaceBase<double> > x_space = model->get_x_space();
    TEST_EQUALITY_CONST( x_space->dim(), 6 );
    RCP<VectorBase<double> > x = createMember(*x_space);
    RCP<Thyra::DefaultMultiVectorProductVector<double> >
      x_bar = Teuchos::rcp_dynamic_cast<Thyra::DefaultMultiVectorProductVector<double> >(
        x, false
        );
    TEST_ASSERT( !is_null(x_bar) );
    RCP<Thyra::MultiVectorBase<double> > X = x_bar->getNonconstMultiVector();
    // Test that the vectors produced have the right number of columns
    TEST_EQUALITY_CONST( X->domain()->dim(), 3 );
    // Test that the vectors produced have the right number of rows
    TEST_EQUALITY_CONST( X->range()->dim(), 2 );
  }
  // f_space:
  {
    RCP<const Thyra::VectorSpaceBase<double> > f_space = model->get_f_space();
    TEST_EQUALITY_CONST( f_space->dim(), 6 );
    RCP<VectorBase<double> > f = createMember(*f_space);
    RCP<Thyra::DefaultMultiVectorProductVector<double> >
      f_bar = Teuchos::rcp_dynamic_cast<Thyra::DefaultMultiVectorProductVector<double> >(
        f, false
        );
    TEST_ASSERT( !is_null(f_bar) );
    RCP<Thyra::MultiVectorBase<double> > F = f_bar->getNonconstMultiVector();
    // Test that the vectors produced have the right number of columns
    TEST_EQUALITY_CONST( F->domain()->dim(), 3 );
    // Test that the vectors produced have the right number of rows
    TEST_EQUALITY_CONST( F->range()->dim(), 2 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_ForwardSensitivityExplicitModelEvaluator, evalModel ) {
  typedef Thyra::ModelEvaluatorBase MEB;
  RCP<ForwardSensitivityExplicitModelEvaluator<double> > model =
    forwardSensitivityExplicitModelEvaluator<double>();
  RCP<SinCosModel> innerModel = sinCosModel(false);
  double a = 0.4;
  double f = 1.5;
  double L = 1.6;
  {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set("Accept model parameters",true);
    pl->set("Implicit model formulation",false);
    pl->set("Coeff a", a );
    pl->set("Coeff f", f );
    pl->set("Coeff L", L );
    innerModel->setParameterList(pl);
  }
  model->initializeStructure(innerModel, 0 );
  RCP<VectorBase<double> > x;
  MEB::InArgs<double> pointInArgs;  // Used to change the solution for re-evaluation
  RCP<StepperBase<double> > stepper; // Used for initializePointState
  {
    pointInArgs = innerModel->createInArgs();
    pointInArgs.set_t(0.1);
    x = Thyra::createMember(innerModel->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view( *x );
      x_view[0] = 2.0;
      x_view[1] = 3.0;
    }
    pointInArgs.set_x(x);
    RCP<VectorBase<double> > p0 = Thyra::createMember(innerModel->get_p_space(0));
    {
      Thyra::DetachedVectorView<double> p0_view( *p0 );
      p0_view[0] = a;
      p0_view[1] = f;
      p0_view[2] = L;
    }
    pointInArgs.set_p(0,p0);
    {
      // Create a stepper with these initial conditions to use to call
      // initializePointState on this ME:
      stepper = forwardEulerStepper<double>();
      stepper->setInitialCondition(pointInArgs);
      model->initializePointState(Teuchos::inOutArg(*stepper),false);
    }
  }
  MEB::InArgs<double> inArgs = model->createInArgs();
  RCP<VectorBase<double> > x_bar = Thyra::createMember(model->get_x_space());
  RCP<Thyra::DefaultMultiVectorProductVector<double> >
    s_bar = Teuchos::rcp_dynamic_cast<Thyra::DefaultMultiVectorProductVector<double> >(
      x_bar, true
      );
  RCP<Thyra::MultiVectorBase<double> >
    S = s_bar->getNonconstMultiVector();
  // Fill S with data
  {
    TEST_EQUALITY_CONST( S->domain()->dim(), 3 );
    TEST_EQUALITY_CONST( S->range()->dim(), 2 );
    RCP<VectorBase<double> > S0 = S->col(0);
    RCP<VectorBase<double> > S1 = S->col(1);
    RCP<VectorBase<double> > S2 = S->col(2);
    TEST_EQUALITY_CONST( S0->space()->dim(), 2 );
    TEST_EQUALITY_CONST( S1->space()->dim(), 2 );
    TEST_EQUALITY_CONST( S2->space()->dim(), 2 );
    Thyra::DetachedVectorView<double> S0_view( *S0 );
    S0_view[0] = 7.0;
    S0_view[1] = 8.0;
    Thyra::DetachedVectorView<double> S1_view( *S1 );
    S1_view[0] = 9.0;
    S1_view[1] = 10.0;
    Thyra::DetachedVectorView<double> S2_view( *S2 );
    S2_view[0] = 11.0;
    S2_view[1] = 12.0;
  }
  inArgs.set_x(x_bar);
  MEB::OutArgs<double> outArgs = model->createOutArgs();
  RCP<VectorBase<double> > f_bar = Thyra::createMember(model->get_f_space());
  RCP<Thyra::DefaultMultiVectorProductVector<double> >
    f_sens = Teuchos::rcp_dynamic_cast<Thyra::DefaultMultiVectorProductVector<double> >(
      f_bar, true
      );
  RCP<Thyra::MultiVectorBase<double> >
    F_sens = f_sens->getNonconstMultiVector().assert_not_null();

  V_S(Teuchos::outArg(*f_bar),0.0);
  outArgs.set_f(f_bar);
  
  inArgs.set_t(0.1);
  model->evalModel(inArgs,outArgs);

  // Verify F_sens = df/dx*S = df/dp
  // df/dx = [ 0             1 ]
  //         [ -(f/L)*(f/L)  0 ]
  // S =   [ 7   9  11 ]    x = [ 2 ]
  //       [ 8  10  12 ]        [ 3 ]
  // df/dp = [     0             0                   0              ]
  //         [ (f/L)*(f/L) 2*f/(L*L)*(a-x_0) -2*f*f/(L*L*L)*(a-x_0) ]
  // F_sens_0 = 
  // [            8               ]
  // [ -7*(f/L)*(f/L)+(f*f)/(L*L) ]
  // F_sens_1 = 
  // [            10                    ]
  // [ -9*(f/L)*(f/L)+2*f/(L*L)*(a-x_0) ]
  // F_sens_2 = 
  // [            12                         ]
  // [ -11*(f/L)*(f/L)-2*f*f/(L*L*L)*(a-x_0) ]
  // 
  double tol = 1.0e-10;
  {
    TEST_EQUALITY_CONST( F_sens->domain()->dim(), 3 );
    TEST_EQUALITY_CONST( F_sens->range()->dim(), 2 );
    RCP<VectorBase<double> > F_sens_0 = F_sens->col(0);
    RCP<VectorBase<double> > F_sens_1 = F_sens->col(1);
    RCP<VectorBase<double> > F_sens_2 = F_sens->col(2);
    TEST_EQUALITY_CONST( F_sens_0->space()->dim(), 2 );
    TEST_EQUALITY_CONST( F_sens_1->space()->dim(), 2 );
    TEST_EQUALITY_CONST( F_sens_2->space()->dim(), 2 );

    Thyra::DetachedVectorView<double> F_sens_0_view( *F_sens_0 );
    TEST_FLOATING_EQUALITY( F_sens_0_view[0], 8.0, tol );
    TEST_FLOATING_EQUALITY( F_sens_0_view[1], -7.0*(f/L)*(f/L)+(f*f)/(L*L), tol );

    Thyra::DetachedVectorView<double> F_sens_1_view( *F_sens_1 );
    TEST_FLOATING_EQUALITY( F_sens_1_view[0], 10.0, tol );
    TEST_FLOATING_EQUALITY( F_sens_1_view[1], -9*(f/L)*(f/L)+2*f/(L*L)*(a-2.0), tol );

    Thyra::DetachedVectorView<double> F_sens_2_view( *F_sens_2 );
    TEST_FLOATING_EQUALITY( F_sens_2_view[0], 12.0, tol );
    TEST_FLOATING_EQUALITY( F_sens_2_view[1], -11*(f/L)*(f/L)-2*f*f/(L*L*L)*(a-2.0), tol );
  }

  // Now change x and evaluate again.
  {
    Thyra::DetachedVectorView<double> x_view( *x );
    x_view[0] = 20.0;
    x_view[1] = 21.0;
  }
  // We need to call initializePointState again due to the vector
  // being cloned inside.
  stepper->setInitialCondition(pointInArgs);
  model->initializePointState(Teuchos::inOutArg(*stepper),false);

  model->evalModel(inArgs,outArgs);
  {
    TEST_EQUALITY_CONST( F_sens->domain()->dim(), 3 );
    TEST_EQUALITY_CONST( F_sens->range()->dim(), 2 );
    RCP<VectorBase<double> > F_sens_0 = F_sens->col(0);
    RCP<VectorBase<double> > F_sens_1 = F_sens->col(1);
    RCP<VectorBase<double> > F_sens_2 = F_sens->col(2);
    TEST_EQUALITY_CONST( F_sens_0->space()->dim(), 2 );
    TEST_EQUALITY_CONST( F_sens_1->space()->dim(), 2 );
    TEST_EQUALITY_CONST( F_sens_2->space()->dim(), 2 );

    Thyra::DetachedVectorView<double> F_sens_0_view( *F_sens_0 );
    TEST_FLOATING_EQUALITY( F_sens_0_view[0], 8.0, tol );
    TEST_FLOATING_EQUALITY( F_sens_0_view[1], -7.0*(f/L)*(f/L)+(f*f)/(L*L), tol );

    Thyra::DetachedVectorView<double> F_sens_1_view( *F_sens_1 );
    TEST_FLOATING_EQUALITY( F_sens_1_view[0], 10.0, tol );
    TEST_FLOATING_EQUALITY( F_sens_1_view[1], -9*(f/L)*(f/L)+2*f/(L*L)*(a-20.0), tol );

    Thyra::DetachedVectorView<double> F_sens_2_view( *F_sens_2 );
    TEST_FLOATING_EQUALITY( F_sens_2_view[0], 12.0, tol );
    TEST_FLOATING_EQUALITY( F_sens_2_view[1], -11*(f/L)*(f/L)-2*f*f/(L*L*L)*(a-20.0), tol );
  }

}

} // namespace Rythmos

