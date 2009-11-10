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

#include "../VanderPol/VanderPolModel.hpp"
#include "../SinCos/SinCosModel.hpp"

#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"

#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "Rythmos_IntegratorBuilder.hpp"

namespace Rythmos {

typedef ModelEvaluatorBase MEB;
using Teuchos::as;


TEUCHOS_UNIT_TEST( Rythmos_VanderPolModel, setImplicitFlag ) {
  {
    RCP<VanderPolModel> explicit_model = rcp(new VanderPolModel);
    MEB::InArgs<double> explicit_model_ic;
    TEST_THROW( explicit_model_ic = explicit_model->getNominalValues(), std::logic_error );
    explicit_model->setImplicitFlag(false);
    explicit_model_ic = explicit_model->getNominalValues();
    TEST_EQUALITY_CONST( explicit_model_ic.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( explicit_model_ic.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( explicit_model_ic.supports(MEB::IN_ARG_x_dot), false );
    TEST_EQUALITY_CONST( explicit_model_ic.supports(MEB::IN_ARG_alpha), false );
    TEST_EQUALITY_CONST( explicit_model_ic.supports(MEB::IN_ARG_beta), false );
  }
  {
    RCP<VanderPolModel> implicit_model = rcp(new VanderPolModel);
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

TEUCHOS_UNIT_TEST( Rythmos_VanderPolModel, nominalValues ) {
  //double tol = 1.0e-10;
  {
    RCP<VanderPolModel> explicit_model = vanderPolModel(false);
    MEB::InArgs<double> explicit_ic = explicit_model->getNominalValues();

    TEST_EQUALITY_CONST( explicit_ic.get_t(), 0.0 );
    RCP<const VectorBase<double> > explicit_ic_x = explicit_ic.get_x();
    TEST_EQUALITY_CONST( Thyra::get_ele(*explicit_ic_x,0), 2.0 );
    TEST_EQUALITY_CONST( Thyra::get_ele(*explicit_ic_x,1), 0.0 );
    //TEST_EQUALITY_CONST( explicit_ic.get_beta(), 0.0 );
  }

  {
    RCP<VanderPolModel> explicit_model = vanderPolModel(false);
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set("Implicit model formulation",false);
    pl->set("Accept model parameters",true);
    pl->set("Coeff epsilon",10.0);
    explicit_model->setParameterList(pl);
    MEB::InArgs<double> explicit_ic = explicit_model->getNominalValues();

    TEST_EQUALITY_CONST( explicit_ic.get_t(), 0.0 );

    RCP<const VectorBase<double> > explicit_ic_x = explicit_ic.get_x();
    TEST_EQUALITY_CONST( Thyra::get_ele(*explicit_ic_x,0), 2.0 );
    TEST_EQUALITY_CONST( Thyra::get_ele(*explicit_ic_x,1), 0.0 );

    //TEST_EQUALITY_CONST( explicit_ic.get_beta(), 0.0 );

    RCP<const VectorBase<double> > explicit_ic_p = explicit_ic.get_p(0);
    TEST_EQUALITY_CONST( Thyra::get_ele(*explicit_ic_p,0), 10.0 );
  }
  
  {
    RCP<VanderPolModel> model = vanderPolModel();
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set("Implicit model formulation",false);
    pl->set("Provide nominal values",false);
    model->setParameterList(pl);
    MEB::InArgs<double> model_ic = model->getNominalValues();

    TEST_EQUALITY_CONST( model_ic.get_t(), 0.0 );
    RCP<const VectorBase<double> > x = model_ic.get_x();
    TEST_EQUALITY_CONST( is_null(x), true );
    //TEST_EQUALITY_CONST( model_ic.get_beta(), 0.0 );
  }

  {
    RCP<VanderPolModel> implicit_model = vanderPolModel(true);
    MEB::InArgs<double> implicit_ic = implicit_model->getNominalValues();

    TEST_EQUALITY_CONST( implicit_ic.get_t(), 0.0 );

    RCP<const VectorBase<double> > implicit_ic_x = implicit_ic.get_x();
    TEST_EQUALITY_CONST( Thyra::get_ele(*implicit_ic_x,0), 2.0 );
    TEST_EQUALITY_CONST( Thyra::get_ele(*implicit_ic_x,1), 0.0 );

    RCP<const VectorBase<double> > implicit_ic_x_dot = implicit_ic.get_x_dot();
    TEST_EQUALITY_CONST( Thyra::get_ele(*implicit_ic_x_dot,0), 0.0 );
    TEST_EQUALITY_CONST( Thyra::get_ele(*implicit_ic_x_dot,1), 0.0 );

    TEST_EQUALITY_CONST( implicit_ic.get_alpha(), 0.0 );
    TEST_EQUALITY_CONST( implicit_ic.get_beta(), 0.0 );
  }

  {
    RCP<VanderPolModel> model = vanderPolModel();
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set("Implicit model formulation", true);
    pl->set("Provide nominal values",false);
    model->setParameterList(pl);

    MEB::InArgs<double> model_ic = model->getNominalValues();

    TEST_EQUALITY_CONST( model_ic.get_t(), 0.0 );

    RCP<const VectorBase<double> > x = model_ic.get_x();
    TEST_EQUALITY_CONST( is_null(x), true );

    RCP<const VectorBase<double> > x_dot = model_ic.get_x_dot();
    TEST_EQUALITY_CONST( is_null(x_dot), true );

    TEST_EQUALITY_CONST( model_ic.get_alpha(), 0.0 );
    TEST_EQUALITY_CONST( model_ic.get_beta(), 0.0 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_VanderPolModel, p_names ) {
  RCP<VanderPolModel> model = vanderPolModel();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set("Accept model parameters", true);
  model->setParameterList(pl);
  RCP<const Teuchos::Array<std::string> > p_names;
#ifdef RYTHMOS_DEBUG
  TEST_THROW( p_names = model->get_p_names(1), std::logic_error );
#endif // RYTHMOS_DEBUG
  p_names = model->get_p_names(0);
  TEST_EQUALITY_CONST( Teuchos::as<int>(p_names->size()), 1 );
  TEST_EQUALITY_CONST( (*p_names)[0], "Model Coefficient:  epsilon" );
}

TEUCHOS_UNIT_TEST( Rythmos_VanderPolModel, spaces ) {
  {
    RCP<VanderPolModel> explicit_model = vanderPolModel(false);
    RCP<const Thyra::VectorSpaceBase<double> > x_space = explicit_model->get_x_space();
    TEST_EQUALITY_CONST( x_space->dim(), 2 );
    RCP<const Thyra::VectorSpaceBase<double> > f_space = explicit_model->get_f_space();
    TEST_EQUALITY_CONST( f_space->dim(), 2 );
    //RCP<const Thyra::VectorSpaceBase<double> > g_space = explicit_model->get_g_space(0);
    //TEST_EQUALITY_CONST( g_space->dim(), 0 );
#ifdef RYTHMOS_DEBUG
    //TEST_THROW( explicit_model->get_g_space(1), std::logic_error );
#endif // RYTHMOS_DEBUG
    RCP<const Thyra::VectorSpaceBase<double> > p_space = explicit_model->get_p_space(0);
    TEST_EQUALITY_CONST( is_null(p_space), true );
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set("Accept model parameters",true);
    explicit_model->setParameterList(pl);
    p_space = explicit_model->get_p_space(0);
    TEST_EQUALITY_CONST( p_space->dim(), 1 );
#ifdef RYTHMOS_DEBUG
    TEST_THROW( explicit_model->get_p_space(1), std::logic_error );
#endif // RYTHMOS_DEBUG
  }
  {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set("Implicit model formulation",true);
    RCP<VanderPolModel> implicit_model = vanderPolModel(true);
    RCP<const Thyra::VectorSpaceBase<double> > x_space = implicit_model->get_x_space();
    TEST_EQUALITY_CONST( x_space->dim(), 2 );
    RCP<const Thyra::VectorSpaceBase<double> > f_space = implicit_model->get_f_space();
    TEST_EQUALITY_CONST( f_space->dim(), 2 );
    //RCP<const Thyra::VectorSpaceBase<double> > g_space = implicit_model->get_g_space(0);
    //TEST_EQUALITY_CONST( g_space->dim(), 0 );
#ifdef RYTHMOS_DEBUG
    //TEST_THROW( implicit_model->get_g_space(1), std::logic_error );
#endif // RYTHMOS_DEBUG
    RCP<const Thyra::VectorSpaceBase<double> > p_space = implicit_model->get_p_space(0);
    TEST_EQUALITY_CONST( is_null(p_space), true );
    pl->set("Accept model parameters",true);
    implicit_model->setParameterList(pl);
    p_space = implicit_model->get_p_space(0);
    TEST_EQUALITY_CONST( p_space->dim(), 1 );
#ifdef RYTHMOS_DEBUG
    TEST_THROW( implicit_model->get_p_space(1), std::logic_error );
#endif // RYTHMOS_DEBUG
  }
}

TEUCHOS_UNIT_TEST( Rythmos_VanderPolModel, create_W_op ) {
  RCP<VanderPolModel> explicit_model = vanderPolModel(false);
  RCP<Thyra::LinearOpBase<double> > W_op = explicit_model->create_W_op();
  RCP<Thyra::MultiVectorBase<double> > matrix = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(W_op,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(matrix), false );
  TEST_EQUALITY_CONST( matrix->domain()->dim(), 2 );
  TEST_EQUALITY_CONST( matrix->range()->dim(), 2 );
}

TEUCHOS_UNIT_TEST( Rythmos_VanderPolModel, get_W_factory ) {
  RCP<VanderPolModel> explicit_model = vanderPolModel(false);
  RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > W_factory = explicit_model->get_W_factory();
  RCP<const Thyra::DefaultSerialDenseLinearOpWithSolveFactory<double> > myFactory =
    Teuchos::rcp_dynamic_cast<const Thyra::DefaultSerialDenseLinearOpWithSolveFactory<double> >(W_factory,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(myFactory), false );
}

TEUCHOS_UNIT_TEST( Rythmos_VanderPolModel, createInArgs ) {
  {
    RCP<VanderPolModel> explicit_model = vanderPolModel(false);
    MEB::InArgs<double> inArgs = explicit_model->createInArgs();
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_x_dot), false );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_alpha), false );
    //TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_beta), true );
  }
  {
    RCP<VanderPolModel> implicit_model = vanderPolModel(true);
    MEB::InArgs<double> inArgs = implicit_model->createInArgs();
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_x_dot), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_alpha), true );
    TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_beta), true );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_VanderPolModel, createOutArgs ) {
  {
    RCP<VanderPolModel> explicit_model = vanderPolModel(false);
    MEB::OutArgs<double> outArgs = explicit_model->createOutArgs();
    TEST_EQUALITY_CONST( outArgs.supports(MEB::OUT_ARG_f), true );
    TEST_EQUALITY_CONST( outArgs.supports(MEB::OUT_ARG_W_op), true );
  }
  {
    RCP<VanderPolModel> explicit_model = vanderPolModel(true);
    MEB::OutArgs<double> outArgs = explicit_model->createOutArgs();
    TEST_EQUALITY_CONST( outArgs.supports(MEB::OUT_ARG_f), true );
    TEST_EQUALITY_CONST( outArgs.supports(MEB::OUT_ARG_W_op), true );
  }

}

TEUCHOS_UNIT_TEST( Rythmos_VanderPolModel, exactSolution ) {
  std::vector<double> t_values;
  int N = 10;
  {
    double t = 25;
    for (int i=0 ; i<2*N+1 ; ++i) {
      t_values.push_back( -t + 2*t*i/(2*N) );
    }
  }

  double tol = 1.0e-10;

  double eps = 2.0; // [0.5]
  double x0 = 2.0; // [2.0]
  double x1 = 0.0; // [0.0]
  double t0 = 0.333; // [0.0]
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set("Coeff epsilon", eps);
  pl->set("IC x_0", x0);
  pl->set("IC x_1", x1);
  pl->set("IC t_0", t0);
  {
    RCP<VanderPolModel> explicit_model = vanderPolModel();
    pl->set("Implicit model formulation", false);
    explicit_model->setParameterList(pl);
    MEB::InArgs<double> exact_sol = explicit_model->getExactSolution(0.0);
    TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_x_dot), false );
    TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_alpha), false );
    //TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_beta), true );
    for (int i=0 ; i < as<int>(t_values.size()); ++i) {
      const double t = t_values[i];
      MEB::InArgs<double> exact_sol2 = explicit_model->getExactSolution(t);
      TEST_FLOATING_EQUALITY( exact_sol2.get_t(), t, tol );
      RCP<const VectorBase<double> > x = exact_sol2.get_x();
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*x,0), 2*cos(t)+eps*(0.75*sin(t)-0.25*sin(3.0*t)), tol );
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*x,1), -2*sin(t)+eps*(0.75*cos(t)-0.75*cos(3.0*t)), tol );
      //TEST_EQUALITY_CONST( exact_sol2.get_beta(), 0.0 );
    }
  }
  {
    RCP<VanderPolModel> implicit_model = vanderPolModel();
    pl->set("Implicit model formulation", true);
    implicit_model->setParameterList(pl);
    MEB::InArgs<double> exact_sol2 = implicit_model->getExactSolution(0.0);
    TEST_EQUALITY_CONST( exact_sol2.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( exact_sol2.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( exact_sol2.supports(MEB::IN_ARG_x_dot), true );
    TEST_EQUALITY_CONST( exact_sol2.supports(MEB::IN_ARG_alpha), true );
    TEST_EQUALITY_CONST( exact_sol2.supports(MEB::IN_ARG_beta), true );
    for (int i=0 ; i < as<int>(t_values.size()); ++i) {
      const double t = t_values[i];
      MEB::InArgs<double> exact_sol3 = implicit_model->getExactSolution(t);
      TEST_FLOATING_EQUALITY( exact_sol3.get_t(), t, tol );
      RCP<const VectorBase<double> > x = exact_sol3.get_x();
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*x,0),  2*cos(t)+eps*(0.75*sin(t)-0.25*sin(3.0*t)), tol );
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*x,1), -2*sin(t)+eps*(0.75*cos(t)-0.75*cos(3.0*t)), tol );
      RCP<const VectorBase<double> > x_dot = exact_sol3.get_x_dot();
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*x_dot,0), -2*sin(t)+eps*(0.75*cos(t)-0.75*cos(3.0*t)), tol );
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*x_dot,1), -2*cos(t)+eps*(-0.75*sin(t)+3.0*0.75*sin(3.0*t)), tol );
      TEST_EQUALITY_CONST( exact_sol3.get_alpha(), 0.0 );
      TEST_EQUALITY_CONST( exact_sol3.get_beta(), 0.0 );
    }
  }
}

TEUCHOS_UNIT_TEST( Rythmos_VanderPolModel, exactSensSolution ) {
  std::vector<double> t_values;
  int N = 10;
  {
    double t = 25;
    for (int i=0 ; i<2*N+1 ; ++i) {
      t_values.push_back( -t + 2*t*i/(2*N) );
    }
  }

  double tol = 1.0e-10;

  double eps = 2.0; // [0.5]
  double x0 = 2.0; // [2.0]
  double x1 = 0.0; // [0.0]
  double t0 = 0.333; // [0.0]
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set("Coeff epsilon", eps);
  pl->set("IC x_0", x0);
  pl->set("IC x_1", x1);
  pl->set("IC t_0", t0);
  {
    RCP<VanderPolModel> explicit_model = vanderPolModel();
    pl->set("Implicit model formulation", false);
    explicit_model->setParameterList(pl);
    MEB::InArgs<double> exact_sol = explicit_model->getExactSensSolution(0, 0.0);
    TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_t), true );
    TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_x), true );
    TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_x_dot), false );
    TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_alpha), false );
    //TEST_EQUALITY_CONST( exact_sol.supports(MEB::IN_ARG_beta), true );
    for (int i=0 ; i < as<int>(t_values.size()); ++i) {
      const double t = t_values[i];
      MEB::InArgs<double> exact_sol2 = explicit_model->getExactSensSolution(0,t);
      TEST_FLOATING_EQUALITY( exact_sol2.get_t(), t, tol );
      RCP<const VectorBase<double> > sens = exact_sol2.get_x();
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens,0), (0.75*sin(t)-0.25*sin(3.0*t)), tol );
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens,1), (0.75*cos(t)-0.75*cos(3.0*t)), tol );
    }
  }
  {
    RCP<VanderPolModel> implicit_model = vanderPolModel();
    pl->set("Implicit model formulation", true);
    implicit_model->setParameterList(pl);
    for (int i=0 ; i < as<int>(t_values.size()); ++i) {
      const double t = t_values[i];
      MEB::InArgs<double> exact_sol = implicit_model->getExactSensSolution(0,t);
      TEST_FLOATING_EQUALITY( exact_sol.get_t(), t, tol );
      RCP<const VectorBase<double> > sens = exact_sol.get_x();
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens,0), (0.75*sin(t)-0.25*sin(3.0*t)), tol );
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens,1), (0.75*cos(t)-0.75*cos(3.0*t)), tol );
      RCP<const VectorBase<double> > sens_dot = exact_sol.get_x_dot();
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens_dot,0), (0.75*cos(t)-3.0*0.25*cos(3.0*t)), tol );
      TEST_FLOATING_EQUALITY( Thyra::get_ele(*sens_dot,1), (-0.75*sin(t)+3.0*0.75*sin(3.0*t)), tol );
    }
  }
}

TEUCHOS_UNIT_TEST( Rythmos_VanderPolModel, evalExplicitModel ) {
  const double eps = 2.0;
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set("Implicit model formulation", false);
  pl->set("Coeff epsilon", eps);

  const double x0 = 5.0;
  const double x1 = 6.0;
  const double t = 4.0;
  //const double beta = 7.0;

  const double x0_ex = 2*cos(t)+eps*(0.75*sin(t)-0.25*sin(3.0*t));
  const double x1_ex = -2*sin(t)+eps*(0.75*cos(t)-0.75*cos(3.0*t));
  const double x1prime_ex = -2*cos(t)+eps*(-0.75*sin(t)+3.0*0.75*sin(3.0*t));
  const double forcing_term = x1prime_ex - eps*(1.0-x0_ex*x0_ex)*x1_ex+x0_ex;

  double tol = 1.0e-10;
  { // Explicit model, just load f.
    RCP<VanderPolModel> explicit_model = vanderPolModel();
    explicit_model->setParameterList(pl);
    MEB::InArgs<double> inArgs = explicit_model->createInArgs();
    MEB::OutArgs<double> outArgs = explicit_model->createOutArgs();
    RCP<VectorBase<double> > x = Thyra::createMember(explicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = x0;
      x_view[1] = x1;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);

    RCP<VectorBase<double> > f = Thyra::createMember(explicit_model->get_f_space());
    outArgs.set_f(f);

    explicit_model->evalModel(inArgs,outArgs);

    TEST_FLOATING_EQUALITY( Thyra::get_ele(*f,0), x1, tol);
    TEST_FLOATING_EQUALITY( Thyra::get_ele(*f,1), (eps*(1.0-x0*x0)*x1-x0) + forcing_term, tol);
  }
  { // Explicit model, load f and W_op
    RCP<VanderPolModel> explicit_model = vanderPolModel();
    explicit_model->setParameterList(pl);
    MEB::InArgs<double> inArgs = explicit_model->createInArgs();
    MEB::OutArgs<double> outArgs = explicit_model->createOutArgs();
    RCP<VectorBase<double> > x = Thyra::createMember(explicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = x0;
      x_view[1] = x1;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);
    //inArgs.set_beta(beta);

    RCP<VectorBase<double> > f = Thyra::createMember(explicit_model->get_f_space());
    RCP<Thyra::LinearOpBase<double> > W_op = explicit_model->create_W_op();
    outArgs.set_f(f);
    outArgs.set_W_op(W_op);

    explicit_model->evalModel(inArgs,outArgs);

    TEST_FLOATING_EQUALITY( Thyra::get_ele(*f,0), x1, tol );
    TEST_FLOATING_EQUALITY( Thyra::get_ele(*f,1), (eps*(1.0-x0*x0)*x1-x0) + forcing_term, tol );

    RCP<Thyra::MultiVectorBase<double> > matrix = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(W_op);
    Thyra::DetachedMultiVectorView<double> matrix_view(*matrix);
    TEST_EQUALITY_CONST( matrix_view(0,0), 0.0 );
    TEST_EQUALITY_CONST( matrix_view(0,1), 1.0 );
    TEST_FLOATING_EQUALITY( matrix_view(1,0), eps*(-2.0*x0)*x1-1.0, tol );
    TEST_FLOATING_EQUALITY( matrix_view(1,1), eps*(1.0-x0*x0), tol );
  }
  { // Explicit model, just load W_op
    RCP<VanderPolModel> explicit_model = vanderPolModel();
    explicit_model->setParameterList(pl);
    MEB::InArgs<double> inArgs = explicit_model->createInArgs();
    MEB::OutArgs<double> outArgs = explicit_model->createOutArgs();
    RCP<VectorBase<double> > x = Thyra::createMember(explicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = x0;
      x_view[1] = x1;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);
    //inArgs.set_beta(beta);

    RCP<Thyra::LinearOpBase<double> > W_op = explicit_model->create_W_op();
    outArgs.set_W_op(W_op);

    explicit_model->evalModel(inArgs,outArgs);

    RCP<Thyra::MultiVectorBase<double> > matrix = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(W_op);
    Thyra::DetachedMultiVectorView<double> matrix_view(*matrix);
    TEST_EQUALITY_CONST( matrix_view(0,0), 0.0 );
    TEST_EQUALITY_CONST( matrix_view(0,1), 1.0 );
    TEST_FLOATING_EQUALITY( matrix_view(1,0), eps*(-2.0*x0)*x1-1.0, tol );
    TEST_FLOATING_EQUALITY( matrix_view(1,1), eps*(1.0-x0*x0), tol );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_VanderPolModel, evalImplicitModel ) {
  const double eps = 2.0;
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set("Implicit model formulation", true);
  pl->set("Coeff epsilon", eps);

  const double x0 = 5.0;
  const double x1 = 6.0;
  const double x_dot0 = 8.0;
  const double x_dot1 = 9.0;
  const double t = 4.0;
  const double alpha = 11.0;
  const double beta = 7.0;

  const double x0_ex = 2*cos(t)+eps*(0.75*sin(t)-0.25*sin(3.0*t));
  const double x1_ex = -2*sin(t)+eps*(0.75*cos(t)-0.75*cos(3.0*t));
  const double x1prime_ex = -2*cos(t)+eps*(-0.75*sin(t)+3.0*0.75*sin(3.0*t));
  const double forcing_term = x1prime_ex - eps*(1.0-x0_ex*x0_ex)*x1_ex+x0_ex;

  double tol = 1.0e-10;
  { // Implicit model, just load f.
    RCP<VanderPolModel> implicit_model = vanderPolModel();
    implicit_model->setParameterList(pl);
    MEB::InArgs<double> inArgs = implicit_model->createInArgs();
    MEB::OutArgs<double> outArgs = implicit_model->createOutArgs();
    RCP<VectorBase<double> > x = Thyra::createMember(implicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = x0;
      x_view[1] = x1;
    }
    RCP<VectorBase<double> > x_dot = Thyra::createMember(implicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_dot_view(*x_dot);
      x_dot_view[0] = x_dot0;
      x_dot_view[1] = x_dot1;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);
    inArgs.set_x_dot(x_dot);

    RCP<VectorBase<double> > f = Thyra::createMember(implicit_model->get_f_space());
    outArgs.set_f(f);

    implicit_model->evalModel(inArgs,outArgs);

    TEST_FLOATING_EQUALITY( Thyra::get_ele(*f,0), x_dot0 - x1, tol );
    TEST_FLOATING_EQUALITY( Thyra::get_ele(*f,1), x_dot1 - (eps*(1.0-x0*x0)*x1-x0) - forcing_term, tol );
  }
  { // Implicit model, load f and W_op
    RCP<VanderPolModel> implicit_model = vanderPolModel();
    implicit_model->setParameterList(pl);
    MEB::InArgs<double> inArgs = implicit_model->createInArgs();
    MEB::OutArgs<double> outArgs = implicit_model->createOutArgs();
    RCP<VectorBase<double> > x = Thyra::createMember(implicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = x0;
      x_view[1] = x1;
    }
    RCP<VectorBase<double> > x_dot = Thyra::createMember(implicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_dot_view(*x_dot);
      x_dot_view[0] = x_dot0;
      x_dot_view[1] = x_dot1;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);
    inArgs.set_x_dot(x_dot);
    inArgs.set_alpha(alpha);
    inArgs.set_beta(beta);

    RCP<VectorBase<double> > f = Thyra::createMember(implicit_model->get_f_space());
    RCP<Thyra::LinearOpBase<double> > W_op = implicit_model->create_W_op();
    outArgs.set_f(f);
    outArgs.set_W_op(W_op);

    implicit_model->evalModel(inArgs,outArgs);

    TEST_FLOATING_EQUALITY( Thyra::get_ele(*f,0), x_dot0 - x1, tol );
    TEST_FLOATING_EQUALITY( Thyra::get_ele(*f,1), x_dot1 - (eps*(1.0-x0*x0)*x1-x0) - forcing_term, tol );

    RCP<Thyra::MultiVectorBase<double> > matrix = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(W_op);
    Thyra::DetachedMultiVectorView<double> matrix_view(*matrix);
    TEST_FLOATING_EQUALITY( matrix_view(0,0), alpha, tol );
    TEST_FLOATING_EQUALITY( matrix_view(0,1), -beta, tol );
    TEST_FLOATING_EQUALITY( matrix_view(1,0), - beta*(eps*(-2.0*x0)*x1-1.0), tol );
    TEST_FLOATING_EQUALITY( matrix_view(1,1), alpha - beta*eps*(1.0-x0*x0), tol );
  }
  { // Implicit model, just load W_op
    RCP<VanderPolModel> implicit_model = vanderPolModel();
    implicit_model->setParameterList(pl);
    MEB::InArgs<double> inArgs = implicit_model->createInArgs();
    MEB::OutArgs<double> outArgs = implicit_model->createOutArgs();
    RCP<VectorBase<double> > x = Thyra::createMember(implicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = x0;
      x_view[1] = x1;
    }
    RCP<VectorBase<double> > x_dot = Thyra::createMember(implicit_model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_dot_view(*x_dot);
      x_dot_view[0] = x_dot0;
      x_dot_view[1] = x_dot1;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);
    inArgs.set_x_dot(x_dot);
    inArgs.set_alpha(alpha);
    inArgs.set_beta(beta);

    RCP<Thyra::LinearOpBase<double> > W_op = implicit_model->create_W_op();
    outArgs.set_W_op(W_op);

    implicit_model->evalModel(inArgs,outArgs);

    RCP<Thyra::MultiVectorBase<double> > matrix = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(W_op);
    Thyra::DetachedMultiVectorView<double> matrix_view(*matrix);
    TEST_FLOATING_EQUALITY( matrix_view(0,0), alpha, tol );
    TEST_FLOATING_EQUALITY( matrix_view(0,1), -beta, tol );
    TEST_FLOATING_EQUALITY( matrix_view(1,0), - beta*(eps*(-2.0*x0)*x1-1.0), tol );
    TEST_FLOATING_EQUALITY( matrix_view(1,1), alpha - beta*eps*(1.0-x0*x0), tol );
  }

}

TEUCHOS_UNIT_TEST( Rythmos_VanderPolModel, modelParams ) {
  const double eps = 13.0;
  const double x0 = 5.0;
  const double x1 = 6.0;
  const double x_dot0 = 8.0;
  const double x_dot1 = 9.0;
  const double t = 4.0;

  const double x0_ex = 2*cos(t)+eps*(0.75*sin(t)-0.25*sin(3.0*t));
  const double x1_ex = -2*sin(t)+eps*(0.75*cos(t)-0.75*cos(3.0*t));
  const double x1prime_ex = -2*cos(t)+eps*(-0.75*sin(t)+3.0*0.75*sin(3.0*t));
  const double forcing_term = x1prime_ex - eps*(1.0-x0_ex*x0_ex)*x1_ex+x0_ex;

  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set("Accept model parameters", true);
  { // Explicit model, just load f.
    RCP<VanderPolModel> model = vanderPolModel();
    pl->set("Implicit model formulation", false);
    model->setParameterList(pl);
    MEB::InArgs<double> inArgs = model->createInArgs();
    MEB::OutArgs<double> outArgs = model->createOutArgs();
    RCP<VectorBase<double> > x = Thyra::createMember(model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = x0;
      x_view[1] = x1;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);
    RCP<VectorBase<double> > p = Thyra::createMember(model->get_p_space(0));
    {
      Thyra::DetachedVectorView<double> p_view(*p);
      p_view[0] = eps;
    }
    inArgs.set_p(0,p);

    RCP<VectorBase<double> > f = Thyra::createMember(model->get_f_space());
    outArgs.set_f(f);

    model->evalModel(inArgs,outArgs);

    double tol = 1.0e-10;
    TEST_FLOATING_EQUALITY( Thyra::get_ele(*f,0), x1, tol );
    TEST_FLOATING_EQUALITY( Thyra::get_ele(*f,1), (eps*(1.0-x0*x0)*x1-x0) + forcing_term, tol );
  }
  { // Implicit model, just load f.
    RCP<VanderPolModel> model = vanderPolModel();
    pl->set("Implicit model formulation", true);
    model->setParameterList(pl);
    MEB::InArgs<double> inArgs = model->createInArgs();
    MEB::OutArgs<double> outArgs = model->createOutArgs();
    RCP<VectorBase<double> > x = Thyra::createMember(model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = x0;
      x_view[1] = x1;
    }
    RCP<VectorBase<double> > xdot = Thyra::createMember(model->get_x_space());
    {
      Thyra::DetachedVectorView<double> xdot_view(*xdot);
      xdot_view[0] = x_dot0;
      xdot_view[1] = x_dot1;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);
    inArgs.set_x_dot(xdot);
    RCP<VectorBase<double> > p = Thyra::createMember(model->get_p_space(0));
    {
      Thyra::DetachedVectorView<double> p_view(*p);
      p_view[0] = eps;
    }
    inArgs.set_p(0,p);

    RCP<VectorBase<double> > f = Thyra::createMember(model->get_f_space());
    outArgs.set_f(f);

    model->evalModel(inArgs,outArgs);

    double tol = 1.0e-10;
    TEST_FLOATING_EQUALITY( Thyra::get_ele(*f,0), x_dot0 - x1, tol );
    TEST_FLOATING_EQUALITY( Thyra::get_ele(*f,1), x_dot1 - (eps*(1.0-x0*x0)*x1-x0) - forcing_term, tol );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_VanderPolModel, DfDp ) {
  const double eps = 13.0;
  const double x0 = 5.0;
  const double x1 = 6.0;
  const double x_dot0 = 8.0;
  const double x_dot1 = 9.0;
  const double t = 4.0;

  const double x0_ex = 2*cos(t)+eps*(0.75*sin(t)-0.25*sin(3.0*t));
  const double x1_ex = -2*sin(t)+eps*(0.75*cos(t)-0.75*cos(3.0*t));
  const double x1prime_ex = -2*cos(t)+eps*(-0.75*sin(t)+3.0*0.75*sin(3.0*t));
  const double forcing_term = x1prime_ex - eps*(1.0-x0_ex*x0_ex)*x1_ex+x0_ex;
  //const double forcing_term_d_eps = - (1.0-x0_ex*x0_ex)*x1_ex;

  double tol = 1.0e-10;

  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set("Accept model parameters", true);
  pl->set("Coeff epsilon", eps);
  pl->set("IC x_0", x0);
  pl->set("IC x_1", x1);
  { // Explicit model, load f and DfDp
    RCP<VanderPolModel> model = vanderPolModel();
    pl->set("Implicit model formulation", false);
    model->setParameterList(pl);
    MEB::InArgs<double> inArgs = model->createInArgs();
    MEB::OutArgs<double> outArgs = model->createOutArgs();
    RCP<VectorBase<double> > x = Thyra::createMember(model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = x0;
      x_view[1] = x1;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);
    RCP<VectorBase<double> > p = Thyra::createMember(model->get_p_space(0));
    {
      Thyra::DetachedVectorView<double> p_view(*p);
      p_view[0] = eps;
    }
    inArgs.set_p(0,p);

    RCP<VectorBase<double> > f = Thyra::createMember(model->get_f_space());
    outArgs.set_f(f);

    RCP<Thyra::MultiVectorBase<double> > mv = Thyra::createMembers(model->get_f_space(),1);
    MEB::Derivative<double> DfDp(mv);
    outArgs.set_DfDp(0, DfDp);

    model->evalModel(inArgs,outArgs);

    TEST_FLOATING_EQUALITY( Thyra::get_ele(*f,0), x1, tol );
    TEST_FLOATING_EQUALITY( Thyra::get_ele(*f,1), (eps*(1.0-x0*x0)*x1-x0) + forcing_term, tol );

    // ToDo: Fix this
//    // Check that DfDp is correct
//    {
//      Thyra::ConstDetachedMultiVectorView<double> mv_view( *mv );
//      TEST_EQUALITY_CONST( mv_view(0,0), 0.0 ); 
//      TEST_EQUALITY_CONST( mv_view(1,0), (1.0-x0*x0)*x1 + forcing_term_d_eps );
//    }    
  }
  { // Implicit model, load f and DfDp
    RCP<VanderPolModel> model = vanderPolModel();
    pl->set("Implicit model formulation", true);
    model->setParameterList(pl);
    MEB::InArgs<double> inArgs = model->createInArgs();
    MEB::OutArgs<double> outArgs = model->createOutArgs();
    RCP<VectorBase<double> > x = Thyra::createMember(model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_view(*x);
      x_view[0] = x0;
      x_view[1] = x1;
    }
    RCP<VectorBase<double> > x_dot = Thyra::createMember(model->get_x_space());
    {
      Thyra::DetachedVectorView<double> x_dot_view(*x_dot);
      x_dot_view[0] = x_dot0;
      x_dot_view[1] = x_dot1;
    }
    inArgs.set_t(t);
    inArgs.set_x(x);
    inArgs.set_x_dot(x_dot);
    RCP<VectorBase<double> > p = Thyra::createMember(model->get_p_space(0));
    {
      Thyra::DetachedVectorView<double> p_view(*p);
      p_view[0] = eps;
    }
    inArgs.set_p(0,p);

    RCP<VectorBase<double> > f = Thyra::createMember(model->get_f_space());
    outArgs.set_f(f);

    RCP<Thyra::MultiVectorBase<double> > mv = Thyra::createMembers(model->get_f_space(),1);
    MEB::Derivative<double> DfDp(mv);
    outArgs.set_DfDp(0, DfDp);

    model->evalModel(inArgs,outArgs);

    TEST_FLOATING_EQUALITY( Thyra::get_ele(*f,0), x_dot0 - x1, tol );
    TEST_FLOATING_EQUALITY( Thyra::get_ele(*f,1), x_dot1 - (eps*(1.0-x0*x0)*x1-x0) - forcing_term, tol );

    // ToDo: Fix this check!
//
//    // Check that DfDp is correct
//    {
//    Thyra::ConstDetachedMultiVectorView<double> mv_view( *mv );
//      TEST_EQUALITY_CONST( mv_view(0,0), 0.0 );
//      TEST_EQUALITY_CONST( mv_view(1,0), -(1.0-x0*x0)*x1 - forcing_term_d_eps );
//    }   
//


  }
}

TEUCHOS_UNIT_TEST( Rythmos_VanderPolModel, fwdIntegration ) {
  RCP<VanderPolModel> model = vanderPolModel(true);
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver = timeStepNonlinearSolver<double>();
  {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set("Default Max Iters",100);
    pl->set("Default Tol",1.0e-12);
    nlSolver->setParameterList(pl);
  }
  double dt = 0.1;
  double finalTime = 0.1;
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*ib->getValidParameters());
  pl->sublist("Stepper Settings").sublist("Stepper Selection").set("Stepper Type","Backward Euler");
  pl->sublist("Integration Control Strategy Selection").set("Integration Control Strategy Type","Simple Integration Control Strategy");
  pl->sublist("Integration Control Strategy Selection").sublist("Simple Integration Control Strategy").set("Take Variable Steps",false);
  pl->sublist("Integration Control Strategy Selection").sublist("Simple Integration Control Strategy").set("Fixed dt",dt);
  ib->setParameterList(pl);
  RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
  //integrator->setVerbLevel(Teuchos::VERB_EXTREME);
  Teuchos::Array<double> time_vec;
  Teuchos::Array<RCP<const VectorBase<double> > > x_vec;
  //time_vec.push_back(pl->sublist("Integrator Settings").get<double>("Final Time"));
  time_vec.push_back(finalTime);
  integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);
  RCP<const VectorBase<double> > x_final = x_vec[0];
  {
    double tol = 1.0e-14;
    Thyra::ConstDetachedVectorView<double> x_final_view( *x_final );
    //double eps = 0.5;
    //double t = finalTime;
    
    // These values are the outputs from running with dt=0.1
    // first order rates were observed in a mesh convergence study by varying dt
    double x0 = 1.98289662139252; //2*cos(t)+eps*(0.75*sin(t)-0.25*sin(3.0*t));
    double x1 = -0.171033786074818; //-2*sin(t)+eps*(0.75*cos(t)-0.75*cos(3.0*t));
    TEST_FLOATING_EQUALITY( x_final_view[0], x0, tol );
    TEST_FLOATING_EQUALITY( x_final_view[1], x1, tol );
  }
}

} // namespace Rythmos



