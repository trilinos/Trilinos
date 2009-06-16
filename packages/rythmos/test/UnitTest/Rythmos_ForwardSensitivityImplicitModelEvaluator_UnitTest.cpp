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
#include "Rythmos_ForwardSensitivityImplicitModelEvaluator.hpp"
#include "Rythmos_UnitTestModels.hpp"

#include "Rythmos_StepperBuilder.hpp"
#include "../SinCos/SinCosModel.hpp"

#include "Thyra_DetachedVectorView.hpp"

namespace Rythmos {


TEUCHOS_UNIT_TEST( Rythmos_ForwardSensitivityImplicitModelEvaluator, create ) { 
  RCP<ForwardSensitivityImplicitModelEvaluator<double> > sensModel =
    forwardSensitivityImplicitModelEvaluator<double>();
  TEST_ASSERT( !is_null(sensModel) );
}

TEUCHOS_UNIT_TEST( Rythmos_ForwardSensitivityImplicitModelEvaluator, args ) {
  RCP<ForwardSensitivityImplicitModelEvaluator<double> > model =
    forwardSensitivityImplicitModelEvaluator<double>();
  RCP<SinCosModel> innerModel = sinCosModel();
  {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set("Accept model parameters",true);
    pl->set("Implicit model formulation",true);
    innerModel->setParameterList(pl);
  }
  model->initializeStructure(innerModel, 0 );
  typedef Thyra::ModelEvaluatorBase MEB;
  {
    MEB::InArgs<double> inArgs = model->createInArgs();
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_x_dot), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_alpha), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_beta), true );
  }
  {
    MEB::OutArgs<double> outArgs = model->createOutArgs();
    TEST_EQUALITY_CONST( outArgs.supports(MEB::OUT_ARG_f), true );
    TEST_EQUALITY_CONST( outArgs.supports(MEB::OUT_ARG_W_op), false );
    TEST_EQUALITY_CONST( outArgs.supports(MEB::OUT_ARG_W), true );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_ForwardSensitivityImplicitModelEvaluator, spaces ) {
  RCP<ForwardSensitivityImplicitModelEvaluator<double> > model =
    forwardSensitivityImplicitModelEvaluator<double>();
  RCP<SinCosModel> innerModel = sinCosModel();
  {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set("Accept model parameters",true);
    pl->set("Implicit model formulation",true);
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

/*
TEUCHOS_UNIT_TEST( Rythmos_ForwardSensitivityImplicitModelEvaluator, evalModel ) {
  typedef Thyra::ModelEvaluatorBase MEB;
  RCP<ForwardSensitivityImplicitModelEvaluator<double> > model =
    forwardSensitivityImplicitModelEvaluator<double>();
  RCP<SinCosModel> innerModel = sinCosModel();
  double a = 0.4;
  double f = 1.5;
  double L = 1.6;
  {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set("Accept model parameters",true);
    pl->set("Coeff a", a );
    pl->set("Coeff f", f );
    pl->set("Coeff L", L );
    pl->set("Implicit model formulation",true);
    innerModel->setParameterList(pl);
  }
  model->initializeStructure(innerModel, 0 );
  RCP<VectorBase<double> > x;
  RCP<VectorBase<double> > x_dot;
  {
    MEB::InArgs<double> pointInArgs = innerModel->createInArgs();
    pointInArgs.set_t(0.1);
    x = Thyra::createMember(innerModel->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view( *x );
      x_view[0] = 2.0;
      x_view[1] = 3.0;
    }
    pointInArgs.set_x(x);
    x_dot = Thyra::createMember(innerModel->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_dot_view( *x_dot );
      x_dot_view[0] = 4.0;
      x_dot_view[1] = 5.0;
    }
    pointInArgs.set_x_dot(x_dot);
    RCP<VectorBase<double> > p0 = Thyra::createMember(innerModel->get_p_space(0));
    {
      Thyra::DetachedVectorView<double> p0_view( *p0 );
      p0_view[0] = a;
      p0_view[1] = f;
      p0_view[2] = L;
    }
    pointInArgs.set_p(0,p0);
    //model->initializeState(
    //    pointInArgs,
    //    W_tilde,
    //    coeff_x_dot,
    //    coeff_x
    //    );
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
  RCP<VectorBase<double> > x_dot_bar = Thyra::createMember(model->get_x_space());
  RCP<Thyra::DefaultMultiVectorProductVector<double> >
    s_dot_bar = Teuchos::rcp_dynamic_cast<Thyra::DefaultMultiVectorProductVector<double> >(
      x_dot_bar, true
      );
  RCP<Thyra::MultiVectorBase<double> >
    S_dot = s_dot_bar->getNonconstMultiVector();
  // Fill S_dot with data
  {
    TEST_EQUALITY_CONST( S_dot->domain()->dim(), 3 );
    TEST_EQUALITY_CONST( S_dot->range()->dim(), 2 );
    RCP<VectorBase<double> > S_dot0 = S_dot->col(0);
    RCP<VectorBase<double> > S_dot1 = S_dot->col(1);
    RCP<VectorBase<double> > S_dot2 = S_dot->col(2);
    TEST_EQUALITY_CONST( S_dot0->space()->dim(), 2 );
    TEST_EQUALITY_CONST( S_dot1->space()->dim(), 2 );
    TEST_EQUALITY_CONST( S_dot2->space()->dim(), 2 );
    Thyra::DetachedVectorView<double> S_dot0_view( *S_dot0 );
    S_dot0_view[0] = 13.0;
    S_dot0_view[1] = 14.0;
    Thyra::DetachedVectorView<double> S_dot1_view( *S_dot1 );
    S_dot1_view[0] = 15.0;
    S_dot1_view[1] = 16.0;
    Thyra::DetachedVectorView<double> S_dot2_view( *S_dot2 );
    S_dot2_view[0] = 17.0;
    S_dot2_view[1] = 18.0;
  }
  inArgs.set_x_dot(x_dot_bar);
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
  
  TEST_THROW(model->evalModel(inArgs,outArgs), std::out_of_range );
  inArgs.set_t(0.1);
  model->evalModel(inArgs,outArgs);

  // Verify F_sens = df/dxdot*Sdot + df/dx*S + df/dp
  // df/dxdot = [ 1 0 ]
  //            [ 0 1 ]
  // df/dx = [ 0             -1 ]
  //         [ +(f/L)*(f/L)  0  ]
  // S =   [ 7   9  11 ]    Sdot = [ 13 15 17 ] x = [ 2 ]
  //       [ 8  10  12 ]           [ 14 16 18 ]     [ 3 ]
  // df/dp = [     0             0                   0                ]
  //         [ -(f/L)*(f/L) -2*f/(L*L)*(a-x_0) +2*f*f/(L*L*L)*(a-x_0) ]
  // F_sens_0 = 
  // [            13-8               ]
  // [ 14+7*(f/L)*(f/L)+(f*f)/(L*L)  ]
  // F_sens_1 = 
  // [            15-10                    ]
  // [ 16+9*(f/L)*(f/L)+2*f/(L*L)*(a-x_0)  ]
  // F_sens_2 = 
  // [            17-12                         ]
  // [ 18+11*(f/L)*(f/L)-2*f*f/(L*L*L)*(a-x_0)  ]
  // 
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
    TEST_EQUALITY_CONST( F_sens_0_view[0], 13.0-8.0 );
    TEST_EQUALITY      ( F_sens_0_view[1], 14.0+7.0*(f/L)*(f/L)+(f*f)/(L*L) );

    Thyra::DetachedVectorView<double> F_sens_1_view( *F_sens_1 );
    TEST_EQUALITY_CONST( F_sens_1_view[0], 15.0-10.0 );
    TEST_EQUALITY      ( F_sens_1_view[1], 16.0+9*(f/L)*(f/L)+2*f/(L*L)*(a-2.0) );

    Thyra::DetachedVectorView<double> F_sens_2_view( *F_sens_2 );
    TEST_EQUALITY_CONST( F_sens_2_view[0], 17.0-12.0 );
    TEST_EQUALITY      ( F_sens_2_view[1], 18.0+11*(f/L)*(f/L)-2*f*f/(L*L*L)*(a-2.0) );
  }

  // Now change x and evaluate again.
  {
    Thyra::DetachedVectorView<double> x_view( *x );
    x_view[0] = 20.0;
    x_view[1] = 21.0;
  }
  {
    Thyra::DetachedVectorView<double> x_dot_view( *x_dot );
    x_dot_view[0] = 30.0;
    x_dot_view[1] = 31.0;
  }
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
    TEST_EQUALITY_CONST( F_sens_0_view[0], 13.0-8.0 );
    TEST_EQUALITY      ( F_sens_0_view[1], 14.0+7.0*(f/L)*(f/L)+(f*f)/(L*L) );

    Thyra::DetachedVectorView<double> F_sens_1_view( *F_sens_1 );
    TEST_EQUALITY_CONST( F_sens_1_view[0], 15.0-10.0 );
    TEST_EQUALITY      ( F_sens_1_view[1], 16.0+9*(f/L)*(f/L)+2*f/(L*L)*(a-20.0) );

    Thyra::DetachedVectorView<double> F_sens_2_view( *F_sens_2 );
    TEST_EQUALITY_CONST( F_sens_2_view[0], 17.0-12.0 );
    TEST_EQUALITY      ( F_sens_2_view[1], 18.0+11*(f/L)*(f/L)-2*f*f/(L*L*L)*(a-20.0) );
  }

}
*/

} // namespace Rythmos

