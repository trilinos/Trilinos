// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RTOpPack_ROpMax.hpp"
#include "supportUnitTestsHelpers.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"


namespace {


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpMax, unitStride, Scalar )
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;

  SubVectorView<Scalar> sv = newSubVectorView<Scalar>(3, as<Scalar>(0.0));

  const ScalarMag m_three = as<ScalarMag>(-3.0);
  const ScalarMag m_two = as<ScalarMag>(-2.0);
  const ScalarMag m_four = as<ScalarMag>(-4.0);

  sv(0) = m_three;
  sv(1) = m_two;
  sv(2) = m_four;

  RTOpPack::ROpMax<Scalar> maxOp;
  RCP<RTOpPack::ReductTarget> max = maxOp.reduct_obj_create();
  maxOp.apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv)(),
    Teuchos::null,
    max.ptr()
    );

  TEST_EQUALITY( maxOp(*max), m_two );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpMax, nonunitStride, Scalar )
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;

  SubVectorView<Scalar> sv = newStridedSubVectorView<Scalar>(3, 3, as<Scalar>(0.0));

  const ScalarMag m_three = as<ScalarMag>(-3.0);
  const ScalarMag m_two = as<ScalarMag>(-2.0);
  const ScalarMag m_four = as<ScalarMag>(-4.0);

  sv(0) = m_three;
  sv(1) = m_two;
  sv(2) = m_four;

  RTOpPack::ROpMax<Scalar> maxOp;
  RCP<RTOpPack::ReductTarget> max = maxOp.reduct_obj_create();
  maxOp.apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv)(),
    Teuchos::null,
    max.ptr()
    );

  TEST_EQUALITY( maxOp(*max), m_two );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpMax, reduct, Scalar )
{
  using Teuchos::as;
  using Teuchos::dyn_cast;
  using RTOpPack::DefaultReductTarget;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  //typedef ScalarTraits<ScalarMag> SMT; // unused


  const ScalarMag m_two = as<ScalarMag>(-2.0);
  const ScalarMag m_four = as<ScalarMag>(-4.0);

  RTOpPack::ROpMax<Scalar> maxOp;

  RCP<RTOpPack::ReductTarget> reduct1 = maxOp.reduct_obj_create();
  RCP<RTOpPack::ReductTarget> reduct2 = maxOp.reduct_obj_create();

  DefaultReductTarget<ScalarMag> &scalarReduct1 =
    dyn_cast<DefaultReductTarget<ScalarMag> >(*reduct1);

  scalarReduct1.set(m_two);
  maxOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(m_four);
  maxOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  TEST_EQUALITY( maxOp(*reduct2), m_two );

}


#define INSTANT_UNIT_TESTS(SCALAR) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ROpMax, unitStride, SCALAR) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ROpMax, nonunitStride, SCALAR) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(ROpMax, reduct, SCALAR)


TEUCHOS_MACRO_TEMPLATE_INSTANT_REAL_SCALAR_TYPES(INSTANT_UNIT_TESTS)


} // namespace
