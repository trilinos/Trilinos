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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/


#include "Thyra_Simple2DModelEvaluator.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace {


using Teuchos::null;
using Teuchos::RCP;
typedef Thyra::ModelEvaluatorBase MEB;
using Thyra::Simple2DModelEvaluator;
using Thyra::simple2DModelEvaluator;


//
// Unit tests
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SimpleModelEvaluator, construct, Scalar )
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
  MEB::InArgs<Scalar> inArgs = model->createInArgs();
  TEST_ASSERT(inArgs.supports(MEB::IN_ARG_x));
  TEST_EQUALITY(inArgs.Np(), 0);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES(
  SimpleModelEvaluator, construct )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SimpleModelEvaluator, eval, Scalar )
{
  RCP<Simple2DModelEvaluator<Scalar> > model = simple2DModelEvaluator<Scalar>();

  Thyra::ModelEvaluatorBase::InArgs<Scalar> in_args = model->getNominalValues();
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> out_args = model->createOutArgs();

  RCP<const Thyra::VectorSpaceBase<Scalar> > f_space = model->get_f_space();
  RCP<Thyra::VectorBase<Scalar> > f = createMember(f_space);
  RCP<Thyra::LinearOpBase<Scalar> > W_op = model->create_W_op();

  Thyra::V_S(f.ptr(),Teuchos::ScalarTraits<Scalar>::zero());

  RCP<Thyra::MultiVectorBase<Scalar> > M = 
    Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(W_op);
  TEUCHOS_ASSERT(Teuchos::nonnull(M));
  Thyra::DetachedMultiVectorView<Scalar> M_dv(*M);

  Thyra::assign(M.ptr(),Teuchos::ScalarTraits<Scalar>::zero());

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
  ScalarMag tol = 
    Teuchos::as<ScalarMag>(10.0) * Teuchos::ScalarTraits<Scalar>::eps();

  // Make sure all entries zeroed out
  TEST_FLOATING_EQUALITY(Thyra::get_ele(*f,0), 0.0, tol);
  TEST_FLOATING_EQUALITY(Thyra::get_ele(*f,1), 0.0, tol);  

  TEST_FLOATING_EQUALITY(M_dv(0,0), 0.0, tol);
  TEST_FLOATING_EQUALITY(M_dv(0,1), 0.0, tol);
  TEST_FLOATING_EQUALITY(M_dv(1,0), 0.0, tol);
  TEST_FLOATING_EQUALITY(M_dv(1,1), 0.0, tol);

  out_args.set_f(f);
  out_args.set_W_op(W_op);

  model->evalModel(in_args, out_args);

  // Based on nominalValue settings x0=1, x1=1, p0=2, p1=0, d=10
  TEST_FLOATING_EQUALITY(Thyra::get_ele(*f,0), 1.0+1.0*1.0-2.0, tol);
  TEST_FLOATING_EQUALITY(Thyra::get_ele(*f,1), 10.0*(1.0*1.0-1.0-0.0), tol);

  TEST_FLOATING_EQUALITY(M_dv(0,0), 1.0, tol);
  TEST_FLOATING_EQUALITY(M_dv(0,1), 2.0, tol);
  TEST_FLOATING_EQUALITY(M_dv(1,0), 10.0 * 2.0 * 1.0, tol);
  TEST_FLOATING_EQUALITY(M_dv(1,1), -10.0, tol);
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES(
  SimpleModelEvaluator, eval )


} // namespace
