/*
// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER
*/


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

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Responses, get_g_names, Scalar )
{
  typedef ModelEvaluatorBase MEB;
  const RCP<const ModelEvaluator<Scalar> > model = getXGTestModel<Scalar>(2, 1);
  TEST_ASSERT( model->get_g_names(0) == Teuchos::null);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( OutArgs, setArgs )


} // namespace Thyra
