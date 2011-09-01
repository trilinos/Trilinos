/*
// @HEADER
// ***********************************************************************
// 
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//                Copyright (2006) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/


#include "RTOpPack_ROpMaxIndexLessThanBound.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"
#include "Teuchos_implicit_cast.hpp"

// Must come last!
#include "supportUnitTestsHelpers.hpp"


namespace {


template<class Scalar>
void basicTest(const int stride, FancyOStream &out, bool &success)
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;

  SubVectorView<Scalar> sv = newStridedSubVectorView<Scalar>(4, stride, as<Scalar>(0.0));

  const Scalar three = as<Scalar>(3.0);
  const Scalar four = as<Scalar>(4.0);
  const Scalar five = as<Scalar>(5.0);
  const Scalar bound = as<Scalar>(4.2);

  sv(0) = ST::zero();
  sv(1) = five;
  sv(2) = four;
  sv(3) = three;

  RTOpPack::ROpMaxIndexLessThanBound<Scalar> maxIndexLessThanBoundOp(bound);
  RCP<RTOpPack::ReductTarget> maxIndex = maxIndexLessThanBoundOp.reduct_obj_create();
  maxIndexLessThanBoundOp.apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv)(),
    Teuchos::null,
    maxIndex.ptr()
    );

  const ScalarIndex<Scalar> maxIndex_vals = maxIndexLessThanBoundOp(*maxIndex);
  TEST_EQUALITY( maxIndex_vals.scalar, four );
  TEST_EQUALITY_CONST( maxIndex_vals.index, 2 );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpMaxIndexLessThanBound, unitStride, Scalar )
{
  basicTest<Scalar>(1, out, success);
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpMaxIndexLessThanBound, nonunitStride, Scalar )
{
  basicTest<Scalar>(3, out, success);
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpMaxIndexLessThanBound, reduct, Scalar )
{
  using Teuchos::dyn_cast;
  typedef ScalarTraits<Scalar> ST;

  const Scalar two = as<Scalar>(2.0);
  const Scalar four = as<Scalar>(4.0);

  RTOpPack::ROpMaxIndexLessThanBound<Scalar> maxIndexLessThanBoundOp;

  RCP<ReductTarget> reduct1 = maxIndexLessThanBoundOp.reduct_obj_create();
  RCP<ReductTarget> reduct2 = maxIndexLessThanBoundOp.reduct_obj_create();

  DefaultReductTarget<ScalarIndex<Scalar> > &scalarReduct1 =
    dyn_cast<DefaultReductTarget<ScalarIndex<Scalar> > >(*reduct1); 

  scalarReduct1.set(ScalarIndex<Scalar>(two, 10));
  maxIndexLessThanBoundOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(ScalarIndex<Scalar>(four, 2));
  maxIndexLessThanBoundOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  const ScalarIndex<Scalar> maxIndex_vals = maxIndexLessThanBoundOp(*reduct2);
  TEST_EQUALITY( maxIndex_vals.scalar, four );
  TEST_EQUALITY_CONST( maxIndex_vals.index, 2 );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpMaxIndexLessThanBound, reductTie_1, Scalar )
{
  using Teuchos::dyn_cast;
  typedef ScalarTraits<Scalar> ST;

  const Scalar two = as<Scalar>(2.0);

  RTOpPack::ROpMaxIndexLessThanBound<Scalar> maxIndexLessThanBoundOp;

  RCP<ReductTarget> reduct1 = maxIndexLessThanBoundOp.reduct_obj_create();
  RCP<ReductTarget> reduct2 = maxIndexLessThanBoundOp.reduct_obj_create();

  DefaultReductTarget<ScalarIndex<Scalar> > &scalarReduct1 =
    dyn_cast<DefaultReductTarget<ScalarIndex<Scalar> > >(*reduct1); 

  scalarReduct1.set(ScalarIndex<Scalar>(two, 4));
  maxIndexLessThanBoundOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(ScalarIndex<Scalar>(two, 10));
  maxIndexLessThanBoundOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  const ScalarIndex<Scalar> maxIndex_vals = maxIndexLessThanBoundOp(*reduct2);
  TEST_EQUALITY( maxIndex_vals.scalar, two );
  TEST_EQUALITY_CONST( maxIndex_vals.index, 4 );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpMaxIndexLessThanBound, reductTie_2, Scalar )
{
  using Teuchos::dyn_cast;
  typedef ScalarTraits<Scalar> ST;

  const Scalar two = as<Scalar>(2.0);

  RTOpPack::ROpMaxIndexLessThanBound<Scalar> maxIndexLessThanBoundOp;

  RCP<ReductTarget> reduct1 = maxIndexLessThanBoundOp.reduct_obj_create();
  RCP<ReductTarget> reduct2 = maxIndexLessThanBoundOp.reduct_obj_create();

  DefaultReductTarget<ScalarIndex<Scalar> > &scalarReduct1 =
    dyn_cast<DefaultReductTarget<ScalarIndex<Scalar> > >(*reduct1); 

  scalarReduct1.set(ScalarIndex<Scalar>(two, 10));
  maxIndexLessThanBoundOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(ScalarIndex<Scalar>(two, 4));
  maxIndexLessThanBoundOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  const ScalarIndex<Scalar> maxIndex_vals = maxIndexLessThanBoundOp(*reduct2);
  TEST_EQUALITY( maxIndex_vals.scalar, two );
  TEST_EQUALITY_CONST( maxIndex_vals.index, 4 );

}


#define INSTANT_UNIT_TESTS(SCALAR) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ROpMaxIndexLessThanBound, unitStride, SCALAR) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ROpMaxIndexLessThanBound, nonunitStride, SCALAR) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ROpMaxIndexLessThanBound, reduct, SCALAR) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ROpMaxIndexLessThanBound, reductTie_1, SCALAR) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ROpMaxIndexLessThanBound, reductTie_2, SCALAR)


TEUCHOS_MACRO_TEMPLATE_INSTANT_REAL_SCALAR_TYPES(INSTANT_UNIT_TESTS)


} // namespace
