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
