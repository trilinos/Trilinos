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


#include "Thyra_SimpleDenseLinearOp.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultSpmdMultiVector.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

#include "Teuchos_UnitTestHarness.hpp"
#include "Thyra_UnitTestHelpers.hpp"



//
// Helper code and declarations
//


bool g_dumpAll = false;
const int g_dim = 4;


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "dump-all", "no-dump-all", &g_dumpAll,
    "Dump lots of data" );
}



//
// Unit Tests
//


namespace Thyra {


template<class Scalar>
RCP<MultiVectorBase<Scalar> >
createSerialMultiVector(const Ordinal rowDim, const Ordinal colDim,
  const Scalar constVal)
{ 
  const RCP<MultiVectorBase<Scalar> > mv =
    createMembers<Scalar>(defaultSpmdVectorSpace<Scalar>(rowDim), colDim);
  assign<Scalar>(mv.ptr(), constVal);
  return mv;
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SimpleDenseLinearOp, basic,
  Scalar )
{

  using Teuchos::rcp_dynamic_cast;
  typedef ScalarTraits<Scalar> ST;

  const RCP<MultiVectorBase<Scalar> > mv =
    createSerialMultiVector<Scalar>(g_dim, g_dim/2, ST::one());
  const RCP<LinearOpBase<Scalar> > op =
    createNonconstSimpleDenseLinearOp<Scalar>(mv);

  TEST_EQUALITY(
    mv,
    rcp_dynamic_cast<SimpleDenseLinearOp<Scalar> >(op)->getNonconstMultiVector()
    );

  Thyra::LinearOpTester<Scalar> linearOpTester;
  linearOpTester.dump_all(g_dumpAll);
  TEST_ASSERT(linearOpTester.check(*op, ptrFromRef(out)));

}

THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( SimpleDenseLinearOp,
  basic )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SimpleDenseLinearOp, scaleLeft,
  Scalar )
{

  using Teuchos::as;
  using Teuchos::rcp_dynamic_cast;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef ScalarTraits<ScalarMag> SMT;

  // Set up the op as all 1's

  const Ordinal rowDim = g_dim, colDim = g_dim/2;

  const RCP<MultiVectorBase<Scalar> > mv =
    createSerialMultiVector<Scalar>(rowDim, colDim, ST::one());
  const RCP<LinearOpBase<Scalar> > op =
    createNonconstSimpleDenseLinearOp<Scalar>(mv);

  // Scale the op from the left (row scaling)

  const RCP<VectorBase<Scalar> > row_scaling = createMember<Scalar>(mv->range());
  const DetachedVectorView<Scalar> row_scaling_dvv(*row_scaling);
  for (Ordinal i = 0; i < rowDim; ++i) {
    row_scaling_dvv(i) = as<Scalar>(i+1);
  }

  rcp_dynamic_cast<ScaledLinearOpBase<Scalar> >(op)->scaleLeft(*row_scaling);

  // Test that resulting left scaling

  const Scalar two = 2.0;

  const RCP<VectorBase<Scalar> > rhs_vec = createMember<Scalar>(mv->domain());
  assign<Scalar>(rhs_vec.ptr(), two);

  const RCP<VectorBase<Scalar> > lhs_vec = createMember<Scalar>(mv->range());

  apply<Scalar>(*op, NOTRANS, *rhs_vec, lhs_vec.ptr());

  if (g_dumpAll) {
    out << "op = " << *op;
    out << "row_scaling = " << *row_scaling;
    out << "rhs_vec = " << *rhs_vec;
    out << "lhs_vec = " << *lhs_vec;
  }


  TEST_FLOATING_EQUALITY(
    sum<Scalar>(*lhs_vec),
    as<Scalar>(sum<Scalar>(*row_scaling) * as<Scalar>(colDim) * two),
    as<ScalarMag>(10.0 * SMT::eps())
    );

}

THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( SimpleDenseLinearOp,
  scaleLeft )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SimpleDenseLinearOp, scaleRight,
  Scalar )
{
  using Teuchos::as;
  using Teuchos::rcp_dynamic_cast;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef ScalarTraits<ScalarMag> SMT;

  // Set up the op as all 1's

  const Ordinal rowDim = g_dim, colDim = g_dim/2;

  const RCP<MultiVectorBase<Scalar> > mv =
    createSerialMultiVector<Scalar>(rowDim, colDim, ST::one());
  const RCP<LinearOpBase<Scalar> > op =
    createNonconstSimpleDenseLinearOp<Scalar>(mv);

  // Scale the op from the right (col scaling)

  const RCP<VectorBase<Scalar> > col_scaling = createMember<Scalar>(mv->domain());
  const DetachedVectorView<Scalar> col_scaling_dvv(*col_scaling);
  for (Ordinal i = 0; i < colDim; ++i) {
    col_scaling_dvv(i) = as<Scalar>(i+1);
  }

  rcp_dynamic_cast<ScaledLinearOpBase<Scalar> >(op, true)->scaleRight(*col_scaling);

  // Test that resulting right scaling

  const Scalar two = 2.0;

  const RCP<VectorBase<Scalar> > rhs_vec = createMember<Scalar>(mv->domain());
  assign<Scalar>(rhs_vec.ptr(), two);

  const RCP<VectorBase<Scalar> > lhs_vec = createMember<Scalar>(mv->range());

  apply<Scalar>(*op, NOTRANS, *rhs_vec, lhs_vec.ptr());

  if (g_dumpAll) {
    out << "op = " << *op;
    out << "col_scaling = " << *col_scaling;
    out << "rhs_vec = " << *rhs_vec;
    out << "lhs_vec = " << *lhs_vec;
  }

  TEST_FLOATING_EQUALITY(
    sum<Scalar>(*lhs_vec),
    as<Scalar>((sum<Scalar>(*col_scaling) * two) * as<Scalar>(rowDim)),
    as<ScalarMag>(10.0 * SMT::eps())
    );

}

THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( SimpleDenseLinearOp,
  scaleRight )


} // namespace Thyra
