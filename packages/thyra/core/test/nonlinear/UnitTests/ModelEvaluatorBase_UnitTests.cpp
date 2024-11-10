// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Thyra_ModelEvaluatorBase.hpp"
//#include "Thyra_Simple2DModelEvaluator.hpp"
#include "Thyra_DummyTestModelEvaluator.hpp"


namespace Thyra {


using Teuchos::null;
using Teuchos::tuple;


template<typename Scalar>
RCP<const ModelEvaluator<Scalar> >
getXGTestModel(const Ordinal x_size, const Ordinal g_size)
{
  return dummyTestModelEvaluator<Scalar>(x_size, null, tuple<Ordinal>(g_size));
}



//
// MEB::Evaluation
//


TEUCHOS_UNIT_TEST( Evaluation, defaultConstruct)
{
  typedef ModelEvaluatorBase MEB;
  MEB::Evaluation<int> eval;
  TEST_ASSERT(is_null(eval));
  TEST_EQUALITY_CONST(eval.getType(), MEB::EVAL_TYPE_EXACT);
  ECHO(const RCP<int> p = eval);
  TEST_ASSERT(is_null(p));
}


TEUCHOS_UNIT_TEST( Evaluation, nullConstruct)
{
  typedef ModelEvaluatorBase MEB;
  MEB::Evaluation<int> eval(null);
  TEST_ASSERT(is_null(eval));
  TEST_EQUALITY_CONST(eval.getType(), MEB::EVAL_TYPE_EXACT);
}



TEUCHOS_UNIT_TEST( Evaluation, nullAssign)
{
  typedef ModelEvaluatorBase MEB;
  ECHO(MEB::Evaluation<int> eval = Teuchos::rcp(new int(1)));
  TEST_ASSERT(nonnull(eval));  
  ECHO(eval = null);
  TEST_EQUALITY_CONST(eval.getType(), MEB::EVAL_TYPE_EXACT);
  TEST_ASSERT(is_null(eval));
}


//
// MEB::OutArgs
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( OutArgs, set_rcp_get_rcp, Scalar )
{
  typedef ModelEvaluatorBase MEB;
  const RCP<const ModelEvaluator<Scalar> > model = getXGTestModel<Scalar>(2, 1);
  MEB::OutArgs<Scalar> outArgs = model->createOutArgs();

  out << "Test that 'f' can be set and gotten as an RCP<VB>!\n";
  const RCP<VectorBase<Scalar> > f = createMember(model->get_f_space());
  outArgs.set_f(f);
  const RCP<VectorBase<Scalar> > f_out = outArgs.get_f();
  TEST_EQUALITY(f_out, f);

  out << "Test that 'g' can be set and gotten as an RCP<VB>!\n";
  const RCP<VectorBase<Scalar> > g = createMember(model->get_g_space(0));
  outArgs.set_g(0, g);
  const RCP<VectorBase<Scalar> > g_out = outArgs.get_g(0);
  TEST_EQUALITY(g_out, g);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( OutArgs, set_rcp_get_rcp )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( OutArgs, set_eval_get_eval, Scalar )
{
  typedef ModelEvaluatorBase MEB;
  const RCP<const ModelEvaluator<Scalar> > model = getXGTestModel<Scalar>(2, 1);
  MEB::OutArgs<Scalar> outArgs = model->createOutArgs();

  out << "Test that 'f' can be set and gotten as an Evaluation<VB>!\n";
  MEB::Evaluation<VectorBase<Scalar> > f(
    createMember(model->get_f_space()), MEB::EVAL_TYPE_APPROX_DERIV);
  outArgs.set_f(f);
  MEB::Evaluation<VectorBase<Scalar> > f_out = outArgs.get_f();
  TEST_EQUALITY(f_out, f);
  TEST_EQUALITY(f_out.getType(), MEB::EVAL_TYPE_APPROX_DERIV);

  out << "Test that 'g' can be set and gotten as an Evaluation<VB>!\n";
  MEB::Evaluation<VectorBase<Scalar> > g(
    createMember(model->get_g_space(0)), MEB::EVAL_TYPE_APPROX_DERIV);
  outArgs.set_g(0, g);
  MEB::Evaluation<VectorBase<Scalar> > g_out = outArgs.get_g(0);
  TEST_EQUALITY(g_out, g);
  TEST_EQUALITY(g_out.getType(), MEB::EVAL_TYPE_APPROX_DERIV);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( OutArgs, set_eval_get_eval )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( OutArgs, set_rcp_get_eval, Scalar )
{
  typedef ModelEvaluatorBase MEB;
  const RCP<const ModelEvaluator<Scalar> > model = getXGTestModel<Scalar>(2, 1);
  MEB::OutArgs<Scalar> outArgs = model->createOutArgs();

  out << "Test that 'f' can be set as an RCP<VB> and gotten as an Evaluation<VB>!\n";
  MEB::Evaluation<VectorBase<Scalar> > f(createMember(model->get_f_space()));
  outArgs.set_f(f);
  MEB::Evaluation<VectorBase<Scalar> > f_out = outArgs.get_f();
  TEST_EQUALITY(f_out, f);
  TEST_EQUALITY(f_out.getType(), MEB::EVAL_TYPE_EXACT);

  out << "Test that 'g' can be set as an RCP<VB> and gotten as an Evaluation<VB>!\n";
  MEB::Evaluation<VectorBase<Scalar> > g(createMember(model->get_g_space(0)));
  outArgs.set_g(0,g);
  MEB::Evaluation<VectorBase<Scalar> > g_out = outArgs.get_g(0);
  TEST_EQUALITY(g_out, g);
  TEST_EQUALITY(g_out.getType(), MEB::EVAL_TYPE_EXACT);

  out << "NOTE: When set as an RCP<> object, we get the right default eval type of EXACT!\n";
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( OutArgs, set_rcp_get_eval )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( OutArgs, set_eval_get_rcp, Scalar )
{
  typedef ModelEvaluatorBase MEB;
  const RCP<const ModelEvaluator<Scalar> > model = getXGTestModel<Scalar>(2, 1);
  MEB::OutArgs<Scalar> outArgs = model->createOutArgs();

  out << "Test that 'f' can be set as an Evaluation <VB> and gotten as an RCP<VB> !\n";
  MEB::Evaluation<VectorBase<Scalar> > f(
    createMember(model->get_f_space()), MEB::EVAL_TYPE_APPROX_DERIV);
  outArgs.set_f(f);
  RCP<VectorBase<Scalar> > f_out = outArgs.get_f();
  TEST_EQUALITY(f_out, f);

  out << "Test that 'g' can be set as an Evaluation <VB> and gotten as an RCP<VB> !\n";
  MEB::Evaluation<VectorBase<Scalar> > g(
    createMember(model->get_g_space(0)), MEB::EVAL_TYPE_APPROX_DERIV);
  outArgs.set_g(0, g);
  RCP<VectorBase<Scalar> > g_out = outArgs.get_g(0);
  TEST_EQUALITY(g_out, g);

  out << "NOTE: We loose the Evaluation type when we get this back as an RCP<> object!\n";
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( OutArgs, set_eval_get_rcp )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( OutArgs, setArgs, Scalar )
{
  typedef ModelEvaluatorBase MEB;
  const RCP<const ModelEvaluator<Scalar> > model = getXGTestModel<Scalar>(2, 1);
  MEB::OutArgs<Scalar> outArgs = model->createOutArgs();

  MEB::Evaluation<VectorBase<Scalar> > f(
    createMember(model->get_f_space()), MEB::EVAL_TYPE_APPROX_DERIV);
  outArgs.set_f(f);

  MEB::Evaluation<VectorBase<Scalar> > g(
    createMember(model->get_g_space(0)), MEB::EVAL_TYPE_VERY_APPROX_DERIV);
  outArgs.set_g(0, g);

  MEB::OutArgs<Scalar> outArgs2 = model->createOutArgs();
  outArgs2.setArgs(outArgs);

  MEB::Evaluation<VectorBase<Scalar> > f_out = outArgs2.get_f();
  TEST_EQUALITY(f_out, f);
  TEST_EQUALITY(f_out.getType(), MEB::EVAL_TYPE_APPROX_DERIV);

  MEB::Evaluation<VectorBase<Scalar> > g_out = outArgs2.get_g(0);
  TEST_EQUALITY(g_out, g);
  TEST_EQUALITY(g_out.getType(), MEB::EVAL_TYPE_VERY_APPROX_DERIV);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( OutArgs, setArgs )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( OutArgs, get_g_names, Scalar )
{
  const RCP<const ModelEvaluator<Scalar> > model = getXGTestModel<Scalar>(2, 1);
  TEST_ASSERT( model->get_g_names(0).size() == 0);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( OutArgs, get_g_names )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( OutArgs, set_DgDx_get_DgDx, Scalar )
{
  typedef ModelEvaluatorBase MEB;
  const Ordinal x_size = 2;
  const Ordinal g_size = 2;
  const RCP<const ModelEvaluator<Scalar> > model =
    dummyTestModelEvaluator<Scalar>(x_size, null, Teuchos::tuple<Ordinal>(1,1), false, false, true, true, true);
  MEB::OutArgs<Scalar> outArgs = model->createOutArgs();

  out << "Test that 'DgDx' can be set and gotten as a DerivativeMultiVector<Scalar>!\n";
  for (int j=0; j<g_size; ++j)
  {
    auto g_space = model->get_g_space(j);
    auto dgdx = Thyra::createMembers(g_space, g_size);
    MEB::DerivativeMultiVector<Scalar> DgDx(dgdx);

    outArgs.set_DgDx(j, DgDx);
    const MEB::Derivative<Scalar> DgDx_out = outArgs.get_DgDx(j);
    TEST_EQUALITY(DgDx_out.getDerivativeMultiVector().getMultiVector(), DgDx.getMultiVector());
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( OutArgs, set_DgDx_get_DgDx )

#ifdef Thyra_BUILD_HESSIAN_SUPPORT

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( OutArgs, set_hess_vec_prod_f_get_hess_vec_prod_f, Scalar )
{
  typedef ModelEvaluatorBase MEB;
  const Ordinal x_size = 2;
  const Ordinal g_size = 2;
  const Ordinal p_size = 3;
  const RCP<const ModelEvaluator<Scalar> > model =
    dummyTestModelEvaluator<Scalar>(x_size, Teuchos::tuple<Ordinal>(1,1,1), Teuchos::tuple<Ordinal>(1,1), false, false, true, true, true);
  MEB::OutArgs<Scalar> outArgs = model->createOutArgs();

  out << "Test that 'hess_vec_prod_f' can be set and gotten as a DerivativeMultiVector<Scalar>!\n";

  auto g_space = model->get_g_space(0);

  auto hess_vec_prod_f_xx_tmp = Thyra::createMembers(g_space, g_size);

  RCP<MultiVectorBase<Scalar> > hess_vec_prod_f_xx(hess_vec_prod_f_xx_tmp);

  outArgs.set_hess_vec_prod_f_xx(hess_vec_prod_f_xx);

  RCP<MultiVectorBase<Scalar> > hess_vec_prod_f_xx_out = outArgs.get_hess_vec_prod_f_xx();

  TEST_EQUALITY(hess_vec_prod_f_xx_out, hess_vec_prod_f_xx);

  for (int l1=0; l1<p_size; ++l1)
  {
    auto p_space_l1 = model->get_p_space(l1);

    auto hess_vec_prod_f_xp_tmp = Thyra::createMembers(g_space, g_size);
    auto hess_vec_prod_f_px_tmp = Thyra::createMembers(p_space_l1, p_size);

    RCP<MultiVectorBase<Scalar> > hess_vec_prod_f_xp(hess_vec_prod_f_xp_tmp);
    RCP<MultiVectorBase<Scalar> > hess_vec_prod_f_px(hess_vec_prod_f_px_tmp);

    outArgs.set_hess_vec_prod_f_xp(l1, hess_vec_prod_f_xp);
    outArgs.set_hess_vec_prod_f_px(l1, hess_vec_prod_f_px);

    RCP<MultiVectorBase<Scalar> > hess_vec_prod_f_xp_out = outArgs.get_hess_vec_prod_f_xp(l1);
    RCP<MultiVectorBase<Scalar> > hess_vec_prod_f_px_out = outArgs.get_hess_vec_prod_f_px(l1);

    TEST_EQUALITY(hess_vec_prod_f_xp_out, hess_vec_prod_f_xp);
    TEST_EQUALITY(hess_vec_prod_f_px_out, hess_vec_prod_f_px);

    for (int l2=0; l2<p_size; ++l2)
    {
      auto hess_vec_prod_f_pp_tmp = Thyra::createMembers(p_space_l1, p_size);
      RCP<MultiVectorBase<Scalar> > hess_vec_prod_f_pp(hess_vec_prod_f_pp_tmp);

      outArgs.set_hess_vec_prod_f_pp(l1, l2, hess_vec_prod_f_pp);

      RCP<MultiVectorBase<Scalar> > hess_vec_prod_f_pp_out = outArgs.get_hess_vec_prod_f_pp(l1, l2);

      TEST_EQUALITY(hess_vec_prod_f_pp_out, hess_vec_prod_f_pp);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( OutArgs, set_hess_vec_prod_f_get_hess_vec_prod_f )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( OutArgs, set_hess_vec_prod_g_get_hess_vec_prod_g, Scalar )
{
  typedef ModelEvaluatorBase MEB;
  const Ordinal x_size = 2;
  const Ordinal g_size = 2;
  const Ordinal p_size = 3;
  const RCP<const ModelEvaluator<Scalar> > model =
    dummyTestModelEvaluator<Scalar>(x_size, Teuchos::tuple<Ordinal>(1,1,1), Teuchos::tuple<Ordinal>(1,1), false, false, true, true, true);
  MEB::OutArgs<Scalar> outArgs = model->createOutArgs();

  out << "Test that 'hess_vec_prod_g' can be set and gotten as a DerivativeMultiVector<Scalar>!\n";
  for (int j=0; j<g_size; ++j)
  {
    auto g_space = model->get_g_space(j);

    auto hess_vec_prod_g_xx_tmp = Thyra::createMembers(g_space, g_size);

    RCP<MultiVectorBase<Scalar> > hess_vec_prod_g_xx(hess_vec_prod_g_xx_tmp);

    outArgs.set_hess_vec_prod_g_xx(j, hess_vec_prod_g_xx);

    RCP<MultiVectorBase<Scalar> > hess_vec_prod_g_xx_out = outArgs.get_hess_vec_prod_g_xx(j);

    TEST_EQUALITY(hess_vec_prod_g_xx_out, hess_vec_prod_g_xx);

    for (int l1=0; l1<p_size; ++l1)
    {
      auto p_space_l1 = model->get_p_space(l1);

      auto hess_vec_prod_g_xp_tmp = Thyra::createMembers(g_space, g_size);
      auto hess_vec_prod_g_px_tmp = Thyra::createMembers(p_space_l1, p_size);

      RCP<MultiVectorBase<Scalar> > hess_vec_prod_g_xp(hess_vec_prod_g_xp_tmp);
      RCP<MultiVectorBase<Scalar> > hess_vec_prod_g_px(hess_vec_prod_g_px_tmp);

      outArgs.set_hess_vec_prod_g_xp(j, l1, hess_vec_prod_g_xp);
      outArgs.set_hess_vec_prod_g_px(j, l1, hess_vec_prod_g_px);

      RCP<MultiVectorBase<Scalar> > hess_vec_prod_g_xp_out = outArgs.get_hess_vec_prod_g_xp(j, l1);
      RCP<MultiVectorBase<Scalar> > hess_vec_prod_g_px_out = outArgs.get_hess_vec_prod_g_px(j, l1);

      TEST_EQUALITY(hess_vec_prod_g_xp_out, hess_vec_prod_g_xp);
      TEST_EQUALITY(hess_vec_prod_g_px_out, hess_vec_prod_g_px);

      for (int l2=0; l2<p_size; ++l2)
      {
        auto hess_vec_prod_g_pp_tmp = Thyra::createMembers(p_space_l1, p_size);
        RCP<MultiVectorBase<Scalar> > hess_vec_prod_g_pp(hess_vec_prod_g_pp_tmp);

        outArgs.set_hess_vec_prod_g_pp(j, l1, l2, hess_vec_prod_g_pp);

        RCP<MultiVectorBase<Scalar> > hess_vec_prod_g_pp_out = outArgs.get_hess_vec_prod_g_pp(j, l1, l2);

        TEST_EQUALITY(hess_vec_prod_g_pp_out, hess_vec_prod_g_pp);
      }
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( OutArgs, set_hess_vec_prod_g_get_hess_vec_prod_g )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( OutArgs, set_hess_f_get_hess_f, Scalar )
{
  typedef ModelEvaluatorBase MEB;
  const Ordinal x_size = 2;
  const Ordinal g_size = 2;
  const Ordinal p_size = 3;
  const RCP<const ModelEvaluator<Scalar> > model =
    dummyTestModelEvaluator<Scalar>(x_size, Teuchos::tuple<Ordinal>(1,1,1), Teuchos::tuple<Ordinal>(1,1), false, false, true, true, true);
  MEB::OutArgs<Scalar> outArgs = model->createOutArgs();

  out << "Test that 'hess_f' can be set and gotten as a RCP<LinearOpBase<Scalar> >!\n";

  auto g_space = model->get_g_space(0);

  auto hess_f_xx_tmp = Thyra::createMembers(g_space, g_size);

  RCP<LinearOpBase<Scalar> > hess_f_xx(hess_f_xx_tmp);

  outArgs.set_hess_f_xx(hess_f_xx);

  RCP<LinearOpBase<Scalar> > hess_f_xx_out = outArgs.get_hess_f_xx();

  TEST_EQUALITY(hess_f_xx_out, hess_f_xx);

  for (int l1=0; l1<p_size; ++l1)
  {
    auto p_space_l1 = model->get_p_space(l1);

    auto hess_f_xp_tmp = Thyra::createMembers(g_space, g_size);

    RCP<LinearOpBase<Scalar> > hess_f_xp(hess_f_xp_tmp);

    outArgs.set_hess_f_xp(l1, hess_f_xp);

    RCP<LinearOpBase<Scalar> > hess_f_xp_out = outArgs.get_hess_f_xp(l1);

    TEST_EQUALITY(hess_f_xp_out, hess_f_xp);

    for (int l2=0; l2<p_size; ++l2)
    {
      auto hess_f_pp_tmp = Thyra::createMembers(p_space_l1, p_size);
      RCP<LinearOpBase<Scalar> > hess_f_pp(hess_f_pp_tmp);

      outArgs.set_hess_f_pp(l1, l2, hess_f_pp);

      RCP<LinearOpBase<Scalar> > hess_f_pp_out = outArgs.get_hess_f_pp(l1, l2);

      TEST_EQUALITY(hess_f_pp_out, hess_f_pp);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( OutArgs, set_hess_f_get_hess_f )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( OutArgs, set_hess_g_get_hess_g, Scalar )
{
  typedef ModelEvaluatorBase MEB;
  const Ordinal x_size = 2;
  const Ordinal g_size = 2;
  const Ordinal p_size = 3;
  const RCP<const ModelEvaluator<Scalar> > model =
    dummyTestModelEvaluator<Scalar>(x_size, Teuchos::tuple<Ordinal>(1,1,1), Teuchos::tuple<Ordinal>(1,1), false, false, true, true, true);
  MEB::OutArgs<Scalar> outArgs = model->createOutArgs();

  out << "Test that 'hess_g' can be set and gotten as a RCP<LinearOpBase<Scalar> >!\n";
  for (int j=0; j<g_size; ++j)
  {
    auto g_space = model->get_g_space(j);

    auto hess_g_xx_tmp = Thyra::createMembers(g_space, g_size);

    RCP<LinearOpBase<Scalar> > hess_g_xx(hess_g_xx_tmp);

    outArgs.set_hess_g_xx(j, hess_g_xx);

    RCP<LinearOpBase<Scalar> > hess_g_xx_out = outArgs.get_hess_g_xx(j);

    TEST_EQUALITY(hess_g_xx_out, hess_g_xx);

    for (int l1=0; l1<p_size; ++l1)
    {
      auto p_space_l1 = model->get_p_space(l1);

      auto hess_g_xp_tmp = Thyra::createMembers(g_space, g_size);

      RCP<LinearOpBase<Scalar> > hess_g_xp(hess_g_xp_tmp);

      outArgs.set_hess_g_xp(j, l1, hess_g_xp);

      RCP<LinearOpBase<Scalar> > hess_g_xp_out = outArgs.get_hess_g_xp(j, l1);

      TEST_EQUALITY(hess_g_xp_out, hess_g_xp);

      for (int l2=0; l2<p_size; ++l2)
      {
        auto hess_g_pp_tmp = Thyra::createMembers(p_space_l1, p_size);
        RCP<LinearOpBase<Scalar> > hess_g_pp(hess_g_pp_tmp);

        outArgs.set_hess_g_pp(j, l1, l2, hess_g_pp);

        RCP<LinearOpBase<Scalar> > hess_g_pp_out = outArgs.get_hess_g_pp(j, l1, l2);

        TEST_EQUALITY(hess_g_pp_out, hess_g_pp);
      }
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( OutArgs, set_hess_g_get_hess_g )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( OutArgs, set_H_get_H, Scalar )
{
  typedef ModelEvaluatorBase MEB;
  const Ordinal x_size = 2;
  const Ordinal g_size = 2;
  const Ordinal p_size = 3;
  const RCP<const ModelEvaluator<Scalar> > model =
    dummyTestModelEvaluator<Scalar>(x_size, Teuchos::tuple<Ordinal>(1,1,1), Teuchos::tuple<Ordinal>(1,1), false, false, true, true, true);
  MEB::OutArgs<Scalar> outArgs = model->createOutArgs();

  out << "Test that 'H' can be set and gotten as a RCP<LinearOpBase<Scalar> >!\n";

  auto g_space = model->get_g_space(0);

  auto H_xx_tmp = Thyra::createMembers(g_space, g_size);

  RCP<LinearOpBase<Scalar> > H_xx(H_xx_tmp);

  outArgs.set_H_xx(H_xx);

  RCP<LinearOpBase<Scalar> > H_xx_out = outArgs.get_H_xx();

  TEST_EQUALITY(H_xx_out, H_xx);

  for (int l1=0; l1<p_size; ++l1)
  {
    auto p_space_l1 = model->get_p_space(l1);

    auto H_xp_tmp = Thyra::createMembers(g_space, g_size);

    RCP<LinearOpBase<Scalar> > H_xp(H_xp_tmp);

    outArgs.set_H_xp(l1, H_xp);

    RCP<LinearOpBase<Scalar> > H_xp_out = outArgs.get_H_xp(l1);

    TEST_EQUALITY(H_xp_out, H_xp);

    for (int l2=0; l2<p_size; ++l2)
    {
      auto H_pp_tmp = Thyra::createMembers(p_space_l1, p_size);
      RCP<LinearOpBase<Scalar> > H_pp(H_pp_tmp);

      outArgs.set_H_pp(l1, l2, H_pp);

      RCP<LinearOpBase<Scalar> > H_pp_out = outArgs.get_H_pp(l1, l2);

      TEST_EQUALITY(H_pp_out, H_pp);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( OutArgs, set_H_get_H )

#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT

//
// MEB::InArgs
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( InArgs, setSolutionArgs, Scalar )
{
  //const RCP<const ModelEvaluator<Scalar> > model = getXGTestModel<Scalar>(2, 1);
  const RCP<const ModelEvaluator<Scalar> > model = 
    dummyTestModelEvaluator<Scalar>(2, Teuchos::tuple<Ordinal>(1), Teuchos::tuple<Ordinal>(1), true, true);

  auto inArgs = model->createInArgs();

  RCP<VectorBase<Scalar>> x = createMember(model->get_x_space());
  inArgs.set_x(x);

  RCP<VectorBase<Scalar>> x_dot = createMember(model->get_x_space());
  inArgs.set_x_dot(x_dot);

  RCP<VectorBase<Scalar>> x_dot_dot = createMember(model->get_x_space());
  inArgs.set_x_dot_dot(x_dot_dot);

  RCP<VectorBase<Scalar>> x_direction = createMember(model->get_x_space());
  inArgs.set_x_direction(x_direction);

  RCP<VectorBase<Scalar>> p_direction = createMember(model->get_p_space(0));
  inArgs.set_p_direction(0, p_direction);

  RCP<VectorBase<Scalar>> f_multiplier = createMember(model->get_x_space());
  inArgs.set_f_multiplier(f_multiplier);

  RCP<VectorBase<Scalar>> g_multiplier = createMember(model->get_x_space());
  inArgs.set_g_multiplier(0, g_multiplier);

  auto inArgs2 = model->createInArgs();
  inArgs2.setArgs(inArgs);

  auto x_out = inArgs2.get_x();
  TEST_EQUALITY(x_out, x);

  auto x_dot_out = inArgs2.get_x_dot();
  TEST_EQUALITY(x_dot_out, x_dot);

  auto x_dot_dot_out = inArgs2.get_x_dot_dot();
  TEST_EQUALITY(x_dot_dot_out, x_dot_dot);

  auto x_direction_out = inArgs2.get_x_direction();
  TEST_EQUALITY(x_direction_out, x_direction);

  auto p_direction_out = inArgs2.get_p_direction(0);
  TEST_EQUALITY(p_direction_out, p_direction);

  auto f_multiplier_out = inArgs2.get_f_multiplier();
  TEST_EQUALITY(f_multiplier_out, f_multiplier);

  auto g_multiplier_out = inArgs2.get_g_multiplier(0);
  TEST_EQUALITY(g_multiplier_out, g_multiplier);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( InArgs, setSolutionArgs )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( InArgs, extendedInArgs, Scalar )
{
  const RCP<const ModelEvaluator<Scalar> > model = 
    dummyTestModelEvaluator<Scalar>(2, null, Teuchos::tuple<Ordinal>(1), false, false);

  // Create the extended inarg data
  RCP<Thyra::MockExtendedInArgs<Scalar>> my_data =
    Teuchos::rcp(new Thyra::MockExtendedInArgs<Scalar>);
  my_data->a = createMember(model->get_x_space());

  auto inArgs = model->createInArgs();

  // Check that the ME supports the extended type
  TEST_ASSERT(inArgs.template supports<Thyra::MockExtendedInArgs<Scalar>>());
  TEST_ASSERT(!inArgs.template supports<double>()); // unsupported type

  // Check set
  RCP<const Thyra::MockExtendedInArgs<Scalar>> const_my_data = my_data;
  inArgs.set(const_my_data);
  
  // Check copy operation
  auto inArgs2 = model->createInArgs();
  inArgs2.setArgs(inArgs);

  // Check get
  auto my_data_2 = inArgs2.template get<const Thyra::MockExtendedInArgs<Scalar>>();
  TEST_EQUALITY(my_data->a, my_data_2->a);

  // Make sure get() throws for unsupported type
  TEST_THROW(inArgs.template get<const double>(),std::runtime_error);
  TEST_THROW(inArgs2.template get<const double>(),std::runtime_error);

  // Disable extended support (tests setting supports to false)
  const RCP<const ModelEvaluator<Scalar> > unsupported_model = 
    dummyTestModelEvaluator<Scalar>(2, null, Teuchos::tuple<Ordinal>(1), false, false,false,false);
  auto inArgs3 = unsupported_model->createInArgs();
  TEST_ASSERT(!inArgs3.template supports<Thyra::MockExtendedInArgs<Scalar>>());
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( InArgs, extendedInArgs )

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( OutArgs, extendedOutArgs, Scalar )
{
  const RCP<const ModelEvaluator<Scalar> > model = 
    dummyTestModelEvaluator<Scalar>(2, null, Teuchos::tuple<Ordinal>(1), false, false);

  // Create the extended inarg data
  RCP<Thyra::MockExtendedOutArgs<Scalar>> my_data = 
    Teuchos::rcp(new Thyra::MockExtendedOutArgs<Scalar>);
  my_data->b = createMember(model->get_x_space());

  auto outArgs = model->createOutArgs();

  // Check that the ME supports the extended type
  TEST_ASSERT(outArgs.template supports<Thyra::MockExtendedOutArgs<Scalar>>());
  TEST_ASSERT(!outArgs.template supports<double>()); // unsupported type

  // Check set
  RCP<const Thyra::MockExtendedOutArgs<Scalar>> const_my_data = my_data;
  outArgs.set(const_my_data);
  
  // Check copy operation
  auto outArgs2 = model->createOutArgs();
  outArgs2.setArgs(outArgs);

  // Check get
  auto my_data_2 = outArgs2.template get<const Thyra::MockExtendedOutArgs<Scalar>>();
  TEST_EQUALITY(my_data->b, my_data_2->b);

  // Make sure get() throws for unsupported type
  TEST_THROW(outArgs.template get<const double>(),std::runtime_error);
  TEST_THROW(outArgs2.template get<const double>(),std::runtime_error);

  // Disable extended support (tests setting supports to false)
  const RCP<const ModelEvaluator<Scalar> > unsupported_model = 
    dummyTestModelEvaluator<Scalar>(2, null, Teuchos::tuple<Ordinal>(1), false, false,false,false);
  auto outArgs3 = unsupported_model->createOutArgs();
  TEST_ASSERT(!outArgs3.template supports<Thyra::MockExtendedOutArgs<Scalar>>());
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( OutArgs, extendedOutArgs )

} // namespace Thyra
