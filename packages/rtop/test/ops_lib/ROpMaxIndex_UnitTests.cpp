// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RTOpPack_ROpMaxIndex.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"
#include "Teuchos_implicit_cast.hpp"

// Must come last!
#include "supportUnitTestsHelpers.hpp"


namespace {


template<class Scalar>
void basicTest(const int stride, FancyOStream &out, bool &success)
{
  using Teuchos::as;
  // typedef ScalarTraits<Scalar> ST; // unused

  SubVectorView<Scalar> sv = newStridedSubVectorView<Scalar>(3, stride, as<Scalar>(0.0));

  const Scalar three = as<Scalar>(3.0);
  const Scalar two = as<Scalar>(2.0);
  const Scalar four = as<Scalar>(4.0);

  sv(0) = three;
  sv(1) = four;
  sv(2) = two;

  RTOpPack::ROpMaxIndex<Scalar> maxIndexOp;
  RCP<RTOpPack::ReductTarget> maxIndex = maxIndexOp.reduct_obj_create();
  maxIndexOp.apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv)(),
    Teuchos::null,
    maxIndex.ptr()
    );

  const ScalarIndex<Scalar> maxIndex_vals = maxIndexOp(*maxIndex);
  TEST_EQUALITY( maxIndex_vals.scalar, four );
  TEST_EQUALITY_CONST( maxIndex_vals.index, 1 );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpMaxIndex, unitStride, Scalar )
{
  basicTest<Scalar>(1, out, success);
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpMaxIndex, nonunitStride, Scalar )
{
  basicTest<Scalar>(3, out, success);
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpMaxIndex, reduct, Scalar )
{
  using Teuchos::dyn_cast;
  // typedef ScalarTraits<Scalar> ST; // unused

  const Scalar two = as<Scalar>(2.0);
  const Scalar four = as<Scalar>(4.0);
  const Scalar three = as<Scalar>(3.0);

  RTOpPack::ROpMaxIndex<Scalar> maxIndexOp;

  RCP<ReductTarget> reduct1 = maxIndexOp.reduct_obj_create();
  RCP<ReductTarget> reduct2 = maxIndexOp.reduct_obj_create();

  DefaultReductTarget<ScalarIndex<Scalar> > &scalarReduct1 =
    dyn_cast<DefaultReductTarget<ScalarIndex<Scalar> > >(*reduct1);

  scalarReduct1.set(ScalarIndex<Scalar>(three, 3));
  maxIndexOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(ScalarIndex<Scalar>(four, 10));
  maxIndexOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(ScalarIndex<Scalar>(two, 2));
  maxIndexOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  const ScalarIndex<Scalar> maxIndex_vals = maxIndexOp(*reduct2);
  TEST_EQUALITY( maxIndex_vals.scalar, four );
  TEST_EQUALITY_CONST( maxIndex_vals.index, 10 );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpMaxIndex, reductTie_1, Scalar )
{
  using Teuchos::dyn_cast;
  // typedef ScalarTraits<Scalar> ST; // unused

  const Scalar two = as<Scalar>(2.0);
  const Scalar three = as<Scalar>(3.0);

  RTOpPack::ROpMaxIndex<Scalar> maxIndexOp;

  RCP<ReductTarget> reduct1 = maxIndexOp.reduct_obj_create();
  RCP<ReductTarget> reduct2 = maxIndexOp.reduct_obj_create();

  DefaultReductTarget<ScalarIndex<Scalar> > &scalarReduct1 =
    dyn_cast<DefaultReductTarget<ScalarIndex<Scalar> > >(*reduct1);

  scalarReduct1.set(ScalarIndex<Scalar>(two, 2));
  maxIndexOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(ScalarIndex<Scalar>(three, 4));
  maxIndexOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(ScalarIndex<Scalar>(three, 10));
  maxIndexOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  const ScalarIndex<Scalar> maxIndex_vals = maxIndexOp(*reduct2);
  TEST_EQUALITY( maxIndex_vals.scalar, three );
  TEST_EQUALITY_CONST( maxIndex_vals.index, 4 );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpMaxIndex, reductTie_2, Scalar )
{
  using Teuchos::dyn_cast;
  // typedef ScalarTraits<Scalar> ST; // unused

  const Scalar two = as<Scalar>(2.0);
  const Scalar three = as<Scalar>(3.0);

  RTOpPack::ROpMaxIndex<Scalar> maxIndexOp;

  RCP<ReductTarget> reduct1 = maxIndexOp.reduct_obj_create();
  RCP<ReductTarget> reduct2 = maxIndexOp.reduct_obj_create();

  DefaultReductTarget<ScalarIndex<Scalar> > &scalarReduct1 =
    dyn_cast<DefaultReductTarget<ScalarIndex<Scalar> > >(*reduct1);

  scalarReduct1.set(ScalarIndex<Scalar>(three, 10));
  maxIndexOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(ScalarIndex<Scalar>(two, 5));
  maxIndexOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(ScalarIndex<Scalar>(three, 4));
  maxIndexOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  const ScalarIndex<Scalar> maxIndex_vals = maxIndexOp(*reduct2);
  TEST_EQUALITY( maxIndex_vals.scalar, three );
  TEST_EQUALITY_CONST( maxIndex_vals.index, 4 );

}


#define INSTANT_UNIT_TESTS(SCALAR) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ROpMaxIndex, unitStride, SCALAR) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ROpMaxIndex, nonunitStride, SCALAR) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ROpMaxIndex, reduct, SCALAR) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ROpMaxIndex, reductTie_1, SCALAR) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ROpMaxIndex, reductTie_2, SCALAR)


TEUCHOS_MACRO_TEMPLATE_INSTANT_REAL_SCALAR_TYPES(INSTANT_UNIT_TESTS)


} // namespace
