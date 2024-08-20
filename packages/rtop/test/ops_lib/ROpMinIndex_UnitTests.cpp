// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RTOpPack_ROpMinIndex.hpp"
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

  SubVectorView<Scalar> sv = newStridedSubVectorView<Scalar>(3, stride, as<Scalar>(0.0));

  const Scalar three = as<Scalar>(3.0);
  const Scalar two = as<Scalar>(2.0);
  const Scalar four = as<Scalar>(4.0);

  sv(0) = three;
  sv(1) = two;
  sv(2) = four;

  RTOpPack::ROpMinIndex<Scalar> minIndexOp;
  RCP<RTOpPack::ReductTarget> minIndex = minIndexOp.reduct_obj_create();
  minIndexOp.apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv)(),
    Teuchos::null,
    minIndex.ptr()
    );

  const ScalarIndex<Scalar> minIndex_vals = minIndexOp(*minIndex);
  TEST_EQUALITY( minIndex_vals.scalar, two );
  TEST_EQUALITY_CONST( minIndex_vals.index, 1 );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpMinIndex, unitStride, Scalar )
{
  basicTest<Scalar>(1, out, success);
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpMinIndex, nonunitStride, Scalar )
{
  basicTest<Scalar>(3, out, success);
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpMinIndex, reduct, Scalar )
{
  using Teuchos::dyn_cast;
  typedef ScalarTraits<Scalar> ST;

  const Scalar two = as<Scalar>(2.0);
  const Scalar four = as<Scalar>(4.0);

  RTOpPack::ROpMinIndex<Scalar> minIndexOp;

  RCP<ReductTarget> reduct1 = minIndexOp.reduct_obj_create();
  RCP<ReductTarget> reduct2 = minIndexOp.reduct_obj_create();

  DefaultReductTarget<ScalarIndex<Scalar> > &scalarReduct1 =
    dyn_cast<DefaultReductTarget<ScalarIndex<Scalar> > >(*reduct1); 

  scalarReduct1.set(ScalarIndex<Scalar>(two, 10));
  minIndexOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(ScalarIndex<Scalar>(four, 2));
  minIndexOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  const ScalarIndex<Scalar> minIndex_vals = minIndexOp(*reduct2);
  TEST_EQUALITY( minIndex_vals.scalar, two );
  TEST_EQUALITY_CONST( minIndex_vals.index, 10 );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpMinIndex, reductTie_1, Scalar )
{
  using Teuchos::dyn_cast;
  typedef ScalarTraits<Scalar> ST;

  const Scalar two = as<Scalar>(2.0);

  RTOpPack::ROpMinIndex<Scalar> minIndexOp;

  RCP<ReductTarget> reduct1 = minIndexOp.reduct_obj_create();
  RCP<ReductTarget> reduct2 = minIndexOp.reduct_obj_create();

  DefaultReductTarget<ScalarIndex<Scalar> > &scalarReduct1 =
    dyn_cast<DefaultReductTarget<ScalarIndex<Scalar> > >(*reduct1); 

  scalarReduct1.set(ScalarIndex<Scalar>(two, 4));
  minIndexOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(ScalarIndex<Scalar>(two, 10));
  minIndexOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  const ScalarIndex<Scalar> minIndex_vals = minIndexOp(*reduct2);
  TEST_EQUALITY( minIndex_vals.scalar, two );
  TEST_EQUALITY_CONST( minIndex_vals.index, 4 );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpMinIndex, reductTie_2, Scalar )
{
  using Teuchos::dyn_cast;
  typedef ScalarTraits<Scalar> ST;

  const Scalar two = as<Scalar>(2.0);

  RTOpPack::ROpMinIndex<Scalar> minIndexOp;

  RCP<ReductTarget> reduct1 = minIndexOp.reduct_obj_create();
  RCP<ReductTarget> reduct2 = minIndexOp.reduct_obj_create();

  DefaultReductTarget<ScalarIndex<Scalar> > &scalarReduct1 =
    dyn_cast<DefaultReductTarget<ScalarIndex<Scalar> > >(*reduct1); 

  scalarReduct1.set(ScalarIndex<Scalar>(two, 10));
  minIndexOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(ScalarIndex<Scalar>(two, 4));
  minIndexOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  const ScalarIndex<Scalar> minIndex_vals = minIndexOp(*reduct2);
  TEST_EQUALITY( minIndex_vals.scalar, two );
  TEST_EQUALITY_CONST( minIndex_vals.index, 4 );

}


#define INSTANT_UNIT_TESTS(SCALAR) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ROpMinIndex, unitStride, SCALAR) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ROpMinIndex, nonunitStride, SCALAR) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ROpMinIndex, reduct, SCALAR) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ROpMinIndex, reductTie_1, SCALAR) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ROpMinIndex, reductTie_2, SCALAR)


TEUCHOS_MACRO_TEMPLATE_INSTANT_REAL_SCALAR_TYPES(INSTANT_UNIT_TESTS)


} // namespace
