//@HEADER
// ***********************************************************************
//
//                     Rythmos Package
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

#include "../SinCos/SinCosModel.hpp"

#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"

namespace Rythmos {

typedef ModelEvaluatorBase MEB;
using Teuchos::as;

TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, setImplicitFlag ) {
  {
    RCP<SinCosModel> explicit_model = rcp(new SinCosModel);
    MEB::InArgs<double> explicit_model_ic;
    TEST_THROW( explicit_model_ic = explicit_model->getNominalValues(), std::logic_error );
    explicit_model->setImplicitFlag(false);
    explicit_model_ic = explicit_model->getNominalValues();
    TEST_EQUALITY_CONST( explicit_model_ic.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( explicit_model_ic.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( explicit_model_ic.supports(MEB::IN_ARG_x_dot), false );
    TEST_EQUALITY_CONST( explicit_model_ic.supports(MEB::IN_ARG_alpha), false );
    TEST_EQUALITY_CONST( explicit_model_ic.supports(MEB::IN_ARG_beta), true );
  }
  {
    RCP<SinCosModel> implicit_model = rcp(new SinCosModel);
    MEB::InArgs<double> implicit_model_ic;
    TEST_THROW( implicit_model_ic = implicit_model->getNominalValues(), std::logic_error );
    implicit_model->setImplicitFlag(true);
    implicit_model_ic = implicit_model->getNominalValues();
    TEST_EQUALITY_CONST( implicit_model_ic.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( implicit_model_ic.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( implicit_model_ic.supports(MEB::IN_ARG_x_dot), true );
    TEST_EQUALITY_CONST( implicit_model_ic.supports(MEB::IN_ARG_alpha), true );
    TEST_EQUALITY_CONST( implicit_model_ic.supports(MEB::IN_ARG_beta), true );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, nominalValues ) {
  {
    RCP<SinCosModel> explicit_model = sinCosModel(false);
    MEB::InArgs<double> explicit_ic = explicit_model->getNominalValues();
    TEST_EQUALITY_CONST( explicit_ic.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( explicit_ic.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( explicit_ic.supports(MEB::IN_ARG_x_dot), false );
    TEST_EQUALITY_CONST( explicit_ic.supports(MEB::IN_ARG_alpha), false );
    TEST_EQUALITY_CONST( explicit_ic.supports(MEB::IN_ARG_beta), true );

    TEST_EQUALITY_CONST( explicit_ic.get_t(), 0.0 );
    RCP<const VectorBase<double> > explicit_ic_x = explicit_ic.get_x();
    TEST_EQUALITY_CONST( Thyra::get_ele(*explicit_ic_x,0), 0.0 );
    TEST_EQUALITY_CONST( Thyra::get_ele(*explicit_ic_x,1), 1.0 );
    TEST_EQUALITY_CONST( explicit_ic.get_beta(), 0.0 );
  }

  {
    RCP<SinCosModel> implicit_model = sinCosModel(true);
    MEB::InArgs<double> implicit_ic = implicit_model->getNominalValues();
    TEST_EQUALITY_CONST( implicit_ic.supports(MEB::IN_ARG_t), true);
    TEST_EQUALITY_CONST( implicit_ic.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( implicit_ic.supports(MEB::IN_ARG_x_dot), true );
    TEST_EQUALITY_CONST( implicit_ic.supports(MEB::IN_ARG_alpha), true );
    TEST_EQUALITY_CONST( implicit_ic.supports(MEB::IN_ARG_beta), true );

    TEST_EQUALITY_CONST( implicit_ic.get_t(), 0.0 );
    RCP<const VectorBase<double> > implicit_ic_x = implicit_ic.get_x();
    TEST_EQUALITY_CONST( Thyra::get_ele(*implicit_ic_x,0), 0.0 );
    TEST_EQUALITY_CONST( Thyra::get_ele(*implicit_ic_x,1), 1.0 );
    RCP<const VectorBase<double> > implicit_ic_x_dot = implicit_ic.get_x_dot();
    TEST_EQUALITY_CONST( Thyra::get_ele(*implicit_ic_x_dot,0), 1.0 );
    TEST_EQUALITY_CONST( Thyra::get_ele(*implicit_ic_x_dot,1), 0.0 );
    TEST_EQUALITY_CONST( implicit_ic.get_alpha(), 0.0 );
    TEST_EQUALITY_CONST( implicit_ic.get_beta(), 0.0 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, spaces ) {
  {
    RCP<SinCosModel> explicit_model = sinCosModel(false);
    RCP<const Thyra::VectorSpaceBase<double> > x_space = explicit_model->get_x_space();
    TEST_EQUALITY_CONST( x_space->dim(), 2 );
    RCP<const Thyra::VectorSpaceBase<double> > f_space = explicit_model->get_f_space();
    TEST_EQUALITY_CONST( f_space->dim(), 2 );
  }
  {
    RCP<SinCosModel> implicit_model = sinCosModel(true);
    RCP<const Thyra::VectorSpaceBase<double> > x_space = implicit_model->get_x_space();
    TEST_EQUALITY_CONST( x_space->dim(), 2 );
    RCP<const Thyra::VectorSpaceBase<double> > f_space = implicit_model->get_f_space();
    TEST_EQUALITY_CONST( f_space->dim(), 2 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, create_W_op ) {
  RCP<SinCosModel> explicit_model = sinCosModel(false);
  RCP<Thyra::LinearOpBase<double> > W_op = explicit_model->create_W_op();
  RCP<Thyra::MultiVectorBase<double> > matrix = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(W_op,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(matrix), false );
  TEST_EQUALITY_CONST( matrix->domain()->dim(), 2 );
  TEST_EQUALITY_CONST( matrix->range()->dim(), 2 );
}

TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, get_W_factory ) {
  RCP<SinCosModel> explicit_model = sinCosModel(false);
  RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > W_factory = explicit_model->get_W_factory();
  RCP<const Thyra::DefaultSerialDenseLinearOpWithSolveFactory<double> > myFactory =
    Teuchos::rcp_dynamic_cast<const Thyra::DefaultSerialDenseLinearOpWithSolveFactory<double> >(W_factory,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(myFactory), false );
}

TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, createInArgs ) {
  {
    RCP<SinCosModel> explicit_model = sinCosModel(false);
    MEB::InArgs<double> inArgs = explicit_model->createInArgs();
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_x_dot), false );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_alpha), false );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_beta), true );
  }
  {
    RCP<SinCosModel> implicit_model = sinCosModel(true);
    MEB::InArgs<double> inArgs = implicit_model->createInArgs();
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_x_dot), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_alpha), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_beta), true );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, createOutArgs ) {
  {
    RCP<SinCosModel> explicit_model = sinCosModel(false);
    MEB::OutArgs<double> outArgs = explicit_model->createOutArgs();
    TEST_EQUALITY_CONST( outArgs.supports(MEB::OUT_ARG_f), true );
    TEST_EQUALITY_CONST( outArgs.supports(MEB::OUT_ARG_W_op), true );
  }
  {
    RCP<SinCosModel> explicit_model = sinCosModel(true);
    MEB::OutArgs<double> outArgs = explicit_model->createOutArgs();
    TEST_EQUALITY_CONST( outArgs.supports(MEB::OUT_ARG_f), true );
    TEST_EQUALITY_CONST( outArgs.supports(MEB::OUT_ARG_W_op), true );
  }

}

TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, exactSolution ) {
  std::vector<double> t_values;
  int N = 10;
  double a = 25;
  for (int i=0 ; i<2*N+1 ; ++i) {
    t_values.push_back( -a + 2*a*i/(2*N) );
  }

  double tol = 1.0e-10;

  {
    RCP<SinCosModel> explicit_model = sinCosModel(false);
    MEB::InArgs<double> exact_sol = explicit_model->getExactSolution(0.0);
    TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_x_dot), false );
    TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_alpha), false );
    TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_beta), true );
    for (int i=0 ; i < as<int>(t_values.size()); ++i) {
      MEB::InArgs<double> exact_sol = explicit_model->getExactSolution(t_values[i]);
      TEST_FLOATING_EQUALITY( exact_sol.get_t(), t_values[i], tol );
      RCP<const VectorBase<double> > x = exact_sol.get_x();
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*x,0), sin(t_values[i]), tol );
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*x,1), cos(t_values[i]), tol );
      TEST_EQUALITY_CONST( exact_sol.get_beta(), 0.0 );
    }
  }
  {
    RCP<SinCosModel> implicit_model = sinCosModel(true);
    MEB::InArgs<double> exact_sol = implicit_model->getExactSolution(0.0);
    TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_x_dot), true );
    TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_alpha), true );
    TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_beta), true );
    for (int i=0 ; i < as<int>(t_values.size()); ++i) {
      MEB::InArgs<double> exact_sol = implicit_model->getExactSolution(t_values[i]);
      TEST_FLOATING_EQUALITY( exact_sol.get_t(), t_values[i], tol );
      RCP<const VectorBase<double> > x = exact_sol.get_x();
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*x,0), sin(t_values[i]), tol );
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*x,1), cos(t_values[i]), tol );
      RCP<const VectorBase<double> > x_dot = exact_sol.get_x_dot();
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*x_dot,0), cos(t_values[i]), tol );
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*x_dot,1), -sin(t_values[i]), tol );
      TEST_EQUALITY_CONST( exact_sol.get_alpha(), 0.0 );
      TEST_EQUALITY_CONST( exact_sol.get_beta(), 0.0 );
    }
  }
}

TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, evalexplicitModel ) {
  { // Explicit model, just load f.
    RCP<SinCosModel> explicit_model = sinCosModel(false);
    MEB::InArgs<double> inArgs = explicit_model->createInArgs();
    MEB::OutArgs<double> outArgs = explicit_model->createOutArgs();
    double t = 4.0; 
    RCP<VectorBase<double> > x = Thyra::createMember(explicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = 5.0;
      x_view[1] = 6.0;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);

    RCP<VectorBase<double> > f = Thyra::createMember(explicit_model->get_f_space());
    outArgs.set_f(f);

    explicit_model->evalModel(inArgs,outArgs);

    TEST_EQUALITY_CONST( Thyra::get_ele(*f,0), 6.0 );
    TEST_EQUALITY_CONST( Thyra::get_ele(*f,1), -5.0 );
  }
  { // Explicit model, load f and W_op
    RCP<SinCosModel> explicit_model = sinCosModel(false);
    MEB::InArgs<double> inArgs = explicit_model->createInArgs();
    MEB::OutArgs<double> outArgs = explicit_model->createOutArgs();
    double t = 4.0; 
    double beta = 7.0;
    RCP<VectorBase<double> > x = Thyra::createMember(explicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = 5.0;
      x_view[1] = 6.0;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);
    inArgs.set_beta(beta);

    RCP<VectorBase<double> > f = Thyra::createMember(explicit_model->get_f_space());
    RCP<Thyra::LinearOpBase<double> > W_op = explicit_model->create_W_op();
    outArgs.set_f(f);
    outArgs.set_W_op(W_op);

    explicit_model->evalModel(inArgs,outArgs);

    TEST_EQUALITY_CONST( Thyra::get_ele(*f,0), 6.0 );
    TEST_EQUALITY_CONST( Thyra::get_ele(*f,1), -5.0 );

    RCP<Thyra::MultiVectorBase<double> > matrix = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(W_op);
    Thyra::DetachedMultiVectorView<double> matrix_view(*matrix);
    TEST_EQUALITY_CONST( matrix_view(0,0), 0.0 );
    TEST_EQUALITY_CONST( matrix_view(0,1), 7.0 );
    TEST_EQUALITY_CONST( matrix_view(1,0), -7.0 );
    TEST_EQUALITY_CONST( matrix_view(1,1), 0.0 );
  }
  { // Explicit model, just load W_op
    RCP<SinCosModel> explicit_model = sinCosModel(false);
    MEB::InArgs<double> inArgs = explicit_model->createInArgs();
    MEB::OutArgs<double> outArgs = explicit_model->createOutArgs();
    double t = 4.0; 
    double beta = 7.0;
    RCP<VectorBase<double> > x = Thyra::createMember(explicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = 5.0;
      x_view[1] = 6.0;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);
    inArgs.set_beta(beta);

    RCP<Thyra::LinearOpBase<double> > W_op = explicit_model->create_W_op();
    outArgs.set_W_op(W_op);

    explicit_model->evalModel(inArgs,outArgs);

    RCP<Thyra::MultiVectorBase<double> > matrix = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(W_op);
    Thyra::DetachedMultiVectorView<double> matrix_view(*matrix);
    TEST_EQUALITY_CONST( matrix_view(0,0), 0.0 );
    TEST_EQUALITY_CONST( matrix_view(0,1), 7.0 );
    TEST_EQUALITY_CONST( matrix_view(1,0), -7.0 );
    TEST_EQUALITY_CONST( matrix_view(1,1), 0.0 );
  }

}

TEUCHOS_UNIT_TEST( Rythmos_SinCosModel, evalImplicitModel ) {
  { // Implicit model, just load f.
    RCP<SinCosModel> explicit_model = sinCosModel(true);
    MEB::InArgs<double> inArgs = explicit_model->createInArgs();
    MEB::OutArgs<double> outArgs = explicit_model->createOutArgs();
    double t = 4.0; 
    RCP<VectorBase<double> > x = Thyra::createMember(explicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = 5.0;
      x_view[1] = 6.0;
    }
    RCP<VectorBase<double> > x_dot = Thyra::createMember(explicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_dot_view(*x_dot);
      x_dot_view[0] = 8.0;
      x_dot_view[1] = 9.0;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);
    inArgs.set_x_dot(x_dot);

    RCP<VectorBase<double> > f = Thyra::createMember(explicit_model->get_f_space());
    outArgs.set_f(f);

    explicit_model->evalModel(inArgs,outArgs);

    TEST_EQUALITY_CONST( Thyra::get_ele(*f,0), 8.0 - 6.0 );
    TEST_EQUALITY_CONST( Thyra::get_ele(*f,1), 9.0 + 5.0 );
  }
  { // Implicit model, load f and W_op
    RCP<SinCosModel> explicit_model = sinCosModel(true);
    MEB::InArgs<double> inArgs = explicit_model->createInArgs();
    MEB::OutArgs<double> outArgs = explicit_model->createOutArgs();
    double t = 4.0; 
    double alpha = 11.0;
    double beta = 7.0;
    RCP<VectorBase<double> > x = Thyra::createMember(explicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = 5.0;
      x_view[1] = 6.0;
    }
    RCP<VectorBase<double> > x_dot = Thyra::createMember(explicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_dot_view(*x_dot);
      x_dot_view[0] = 8.0;
      x_dot_view[1] = 9.0;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);
    inArgs.set_x_dot(x_dot);
    inArgs.set_alpha(alpha);
    inArgs.set_beta(beta);

    RCP<VectorBase<double> > f = Thyra::createMember(explicit_model->get_f_space());
    RCP<Thyra::LinearOpBase<double> > W_op = explicit_model->create_W_op();
    outArgs.set_f(f);
    outArgs.set_W_op(W_op);

    explicit_model->evalModel(inArgs,outArgs);

    TEST_EQUALITY_CONST( Thyra::get_ele(*f,0), 8.0 - 6.0 );
    TEST_EQUALITY_CONST( Thyra::get_ele(*f,1), 9.0 + 5.0 );

    RCP<Thyra::MultiVectorBase<double> > matrix = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(W_op);
    Thyra::DetachedMultiVectorView<double> matrix_view(*matrix);
    TEST_EQUALITY_CONST( matrix_view(0,0), 11.0 );
    TEST_EQUALITY_CONST( matrix_view(0,1), -7.0 );
    TEST_EQUALITY_CONST( matrix_view(1,0), 7.0 );
    TEST_EQUALITY_CONST( matrix_view(1,1), 11.0 );
  }
  { // Implicit model, just load W_op
    RCP<SinCosModel> explicit_model = sinCosModel(true);
    MEB::InArgs<double> inArgs = explicit_model->createInArgs();
    MEB::OutArgs<double> outArgs = explicit_model->createOutArgs();
    double t = 4.0; 
    double alpha = 11.0;
    double beta = 7.0;
    RCP<VectorBase<double> > x = Thyra::createMember(explicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = 5.0;
      x_view[1] = 6.0;
    }
    RCP<VectorBase<double> > x_dot = Thyra::createMember(explicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_dot_view(*x_dot);
      x_dot_view[0] = 8.0;
      x_dot_view[1] = 9.0;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);
    inArgs.set_x_dot(x_dot);
    inArgs.set_alpha(alpha);
    inArgs.set_beta(beta);

    RCP<Thyra::LinearOpBase<double> > W_op = explicit_model->create_W_op();
    outArgs.set_W_op(W_op);

    explicit_model->evalModel(inArgs,outArgs);

    RCP<Thyra::MultiVectorBase<double> > matrix = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(W_op);
    Thyra::DetachedMultiVectorView<double> matrix_view(*matrix);
    TEST_EQUALITY_CONST( matrix_view(0,0), 11.0 );
    TEST_EQUALITY_CONST( matrix_view(0,1), -7.0 );
    TEST_EQUALITY_CONST( matrix_view(1,0), 7.0 );
    TEST_EQUALITY_CONST( matrix_view(1,1), 11.0 );
  }

}

} // namespace Rythmos



