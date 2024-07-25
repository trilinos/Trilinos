// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RTOpPack_ROpMin.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"
#include "supportUnitTestsHelpers.hpp"


namespace {


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpMin, unitStride, Scalar )
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;

  SubVectorView<Scalar> sv = newSubVectorView<Scalar>(3, as<Scalar>(0.0));

  const ScalarMag three = as<ScalarMag>(3.0);
  const ScalarMag two = as<ScalarMag>(2.0);
  const ScalarMag four = as<ScalarMag>(4.0);

  sv(0) = three;
  sv(1) = two;
  sv(2) = four;

  RTOpPack::ROpMin<Scalar> minOp;
  RCP<RTOpPack::ReductTarget> min = minOp.reduct_obj_create();
  minOp.apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv)(),
    Teuchos::null,
    min.ptr()
    );

  TEST_EQUALITY( minOp(*min), two );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpMin, nonunitStride, Scalar )
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;

  SubVectorView<Scalar> sv = newStridedSubVectorView<Scalar>(3, 3, as<Scalar>(0.0));

  const ScalarMag three = as<ScalarMag>(3.0);
  const ScalarMag two = as<ScalarMag>(2.0);
  const ScalarMag four = as<ScalarMag>(4.0);

  sv(0) = three;
  sv(1) = two;
  sv(2) = four;

  RTOpPack::ROpMin<Scalar> minOp;
  RCP<RTOpPack::ReductTarget> min = minOp.reduct_obj_create();
  minOp.apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv)(),
    Teuchos::null,
    min.ptr()
    );

  TEST_EQUALITY( minOp(*min), two );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpMin, reduct, Scalar )
{
  using Teuchos::as;
  using Teuchos::dyn_cast;
  using RTOpPack::DefaultReductTarget;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef ScalarTraits<ScalarMag> SMT;

  const ScalarMag two = as<ScalarMag>(2.0);
  const ScalarMag four = as<ScalarMag>(4.0);

  RTOpPack::ROpMin<Scalar> minOp;

  RCP<RTOpPack::ReductTarget> reduct1 = minOp.reduct_obj_create();
  RCP<RTOpPack::ReductTarget> reduct2 = minOp.reduct_obj_create();

  DefaultReductTarget<ScalarMag> &scalarReduct1 =
    dyn_cast<DefaultReductTarget<ScalarMag> >(*reduct1); 

  scalarReduct1.set(two);
  minOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(four);
  minOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  TEST_EQUALITY( minOp(*reduct2), two );

}


#define INSTANT_UNIT_TESTS(SCALAR) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ROpMin, unitStride, SCALAR) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ROpMin, nonunitStride, SCALAR) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ROpMin, reduct, SCALAR)


TEUCHOS_MACRO_TEMPLATE_INSTANT_REAL_SCALAR_TYPES(INSTANT_UNIT_TESTS)


} // namespace
