// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_Simple2DModelEvaluator.hpp"
#include "Thyra_SimpleDenseLinearOp.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace Thyra {


using Teuchos::null;
typedef ModelEvaluatorBase MEB;


//
// Unit tests
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Simple2DModelEvaluator, construct, Scalar )
{
  RCP<Simple2DModelEvaluator<Scalar> > model = simple2DModelEvaluator<Scalar>();
  TEST_ASSERT(model != null);
  TEST_EQUALITY(model->Np(), 0);
  TEST_EQUALITY(model->Ng(), 0);
  TEST_ASSERT(model->get_x_space() != null);
  TEST_EQUALITY(model->get_x_space()->dim(), 2);
  TEST_ASSERT(model->get_f_space() != null);
  TEST_EQUALITY(model->get_f_space()->dim(), 2);
  // ToDo: Test getNominalValues()
  TEST_ASSERT(model->create_W_op() != null);
  TEST_ASSERT(model->get_W_factory() != null);
  TEST_ASSERT(model->create_W_prec() != null);
  MEB::InArgs<Scalar> inArgs = model->createInArgs();
  TEST_ASSERT(inArgs.supports(MEB::IN_ARG_x));
  TEST_EQUALITY(inArgs.Np(), 0);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES(
  Simple2DModelEvaluator, construct )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Simple2DModelEvaluator, eval, Scalar )
{
  using Teuchos::as; using Teuchos::rcp_dynamic_cast;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;

  out << "A) Create the Simple2DModelEvaluator and objects  ...\n";
  
  RCP<Simple2DModelEvaluator<Scalar> > model = simple2DModelEvaluator<Scalar>();

  ModelEvaluatorBase::InArgs<Scalar> in_args = model->getNominalValues();
  ModelEvaluatorBase::OutArgs<Scalar> out_args = model->createOutArgs();

  const RCP<const VectorSpaceBase<Scalar> > f_space = model->get_f_space();
  const RCP<VectorBase<Scalar> > f = createMember(f_space);
  const RCP<LinearOpBase<Scalar> > W_op = model->create_W_op();
  const RCP<PreconditionerBase<Scalar> > W_prec = model->create_W_prec();

  out << "B) Get the concrete output objects ...\n";

  const RCP<SimpleDenseLinearOp<Scalar> > W_sdlo = 
    rcp_dynamic_cast<SimpleDenseLinearOp<Scalar> >(W_op, true);
  const RCP<MultiVectorBase<Scalar> > W_mv =  W_sdlo->getNonconstMultiVector();
  
  const RCP<SimpleDenseLinearOp<Scalar> > W_prec_sdlo =
    rcp_dynamic_cast<SimpleDenseLinearOp<Scalar> >(W_prec->getNonconstUnspecifiedPrecOp(), true);
  const RCP<MultiVectorBase<Scalar> > W_prec_mv =  W_prec_sdlo->getNonconstMultiVector();

  out << "C) Fill the output objects with the wrong values ...\n";

  const Scalar threethree = 3.3;
  V_S(f.ptr(), threethree);
  assign(W_mv.ptr(), threethree);
  assign(W_prec_mv.ptr(), threethree);

  out << "D) Make sure all output objects are initialized to the wrong values ...\n";

  const ScalarMag tol = as<ScalarMag>(10.0) * ST::eps();

  {
    const ConstDetachedVectorView<Scalar> f_dv(f);
    TEST_FLOATING_EQUALITY(f_dv[0], threethree, tol);
    TEST_FLOATING_EQUALITY(f_dv[1], threethree, tol);  

    const ConstDetachedMultiVectorView<Scalar> W_dv(*W_mv);
    TEST_FLOATING_EQUALITY(W_dv(0,0), threethree, tol);
    TEST_FLOATING_EQUALITY(W_dv(0,1), threethree, tol);
    TEST_FLOATING_EQUALITY(W_dv(1,0), threethree, tol);
    TEST_FLOATING_EQUALITY(W_dv(1,1), threethree, tol);

    const ConstDetachedMultiVectorView<Scalar> W_prec_dv(*W_mv);
    TEST_FLOATING_EQUALITY(W_prec_dv(0,0), threethree, tol);
    TEST_FLOATING_EQUALITY(W_prec_dv(0,1), threethree, tol);
    TEST_FLOATING_EQUALITY(W_prec_dv(1,0), threethree, tol);
    TEST_FLOATING_EQUALITY(W_prec_dv(1,1), threethree, tol);
  }

  out << "E) Evaluate the output objects ...\n";

  out_args.set_f(f);
  out_args.set_W_op(W_op);
  out_args.set_W_prec(W_prec);

  model->evalModel(in_args, out_args);

  out << "F) Check the values of the output objects ...\n";

  // Based on nominalValue settings x0=1, x1=1, p0=2, p1=0, d=10

  const ScalarMag zero = ST::zero();

  {
    const ConstDetachedVectorView<Scalar> f_dv(f);
    TEST_FLOATING_EQUALITY(f_dv[0], as<Scalar>(1.0+1.0*1.0-2.0), tol);
    TEST_FLOATING_EQUALITY(f_dv[0], as<Scalar>(10.0*(1.0*1.0-1.0-0.0)), tol);

    const ConstDetachedMultiVectorView<Scalar> W_dv(*W_mv);
    TEST_FLOATING_EQUALITY(W_dv(0,0), as<Scalar>(1.0), tol);
    TEST_FLOATING_EQUALITY(W_dv(0,1), as<Scalar>(2.0), tol);
    TEST_FLOATING_EQUALITY(W_dv(1,0), as<Scalar>(10.0 * 2.0 * 1.0), tol);
    TEST_FLOATING_EQUALITY(W_dv(1,1), as<Scalar>(-10.0), tol);

    const ConstDetachedMultiVectorView<Scalar> W_prec_dv(*W_prec_mv);
    TEST_FLOATING_EQUALITY(W_prec_dv(0,0), as<Scalar>(1.0), tol);
    TEST_FLOATING_EQUALITY(W_prec_dv(0,1), zero, tol);
    TEST_FLOATING_EQUALITY(W_prec_dv(1,0), zero, tol);
    TEST_FLOATING_EQUALITY(W_prec_dv(1,1), as<Scalar>(-0.1), tol);
  }
    
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES(
  Simple2DModelEvaluator, eval )


} // namespace Thyra
