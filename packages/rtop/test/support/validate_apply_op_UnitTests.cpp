// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RTOpPack_ROpSum.hpp"
#include "RTOpPack_ROpDotProd.hpp"
#include "RTOpPack_TOpAXPY.hpp"

#include "supportUnitTestsHelpers.hpp"


namespace {


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( validate_apply_op, noVectors, Scalar )
{
  using Teuchos::null;
  RTOpPack::ROpSum<Scalar> sumOp;
  RCP<RTOpPack::ReductTarget> sum = sumOp.reduct_obj_create();
  TEST_THROW(
    RTOpPack::validate_apply_op<Scalar>( sumOp, 1, 0, true, null, null, sum.ptr() ),
    RTOpPack::InvalidNumVecs
    );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( validate_apply_op, noVectors )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( validate_apply_op, wrongNumberOfSubVecs, Scalar )
{

  using Teuchos::null;
  using Teuchos::tuple;
  using Teuchos::as;
  
  typedef ScalarTraits<Scalar> ST;

  RTOpPack::ROpSum<Scalar> sumOp;

  ConstSubVectorView<Scalar> sv = newSubVectorView<Scalar>(n, ST::random());
  RCP<RTOpPack::ReductTarget> sum = sumOp.reduct_obj_create();
  TEST_THROW(
    RTOpPack::validate_apply_op<Scalar>( sumOp, 1, 0, true,
      tuple(sv, sv)(), null, sum.ptr() ),
    RTOpPack::InvalidNumVecs
    );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( validate_apply_op,
  wrongNumberOfSubVecs )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( validate_apply_op, wrongNumberOfTargSubVecs, Scalar )
{

  using Teuchos::null;
  using Teuchos::tuple;
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;

  RTOpPack::ROpSum<Scalar> sumOp;

  SubVectorView<Scalar> sv = newSubVectorView<Scalar>(n, ST::random());
  ConstSubVectorView<Scalar> csv = sv;
  RCP<RTOpPack::ReductTarget> sum = sumOp.reduct_obj_create();
  TEST_THROW(
    RTOpPack::validate_apply_op<Scalar>( sumOp, 1, 0, true,
      tuple(csv)(), tuple(sv)(), sum.ptr() ),
    RTOpPack::InvalidNumTargVecs
    );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( validate_apply_op,
  wrongNumberOfTargSubVecs)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( validate_apply_op, inCompatibleSubVecs_1_1_a, Scalar )
{

  using Teuchos::null;
  using Teuchos::tuple;
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;

  RTOpPack::TOpAXPY<Scalar> axpyOp(ST::one());

  ConstSubVectorView<Scalar> sv = newSubVectorView<Scalar>(n, ST::random());
  SubVectorView<Scalar> tsv = newSubVectorView<Scalar>(n+1, ST::random());
  TEST_THROW(
    RTOpPack::validate_apply_op<Scalar>( axpyOp, 1, 1, false,
      tuple(sv)(), tuple(tsv)(), null ),
    RTOpPack::IncompatibleVecs
    );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( validate_apply_op,
  inCompatibleSubVecs_1_1_a)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( validate_apply_op, inCompatibleSubVecs_1_1_b, Scalar )
{

  using Teuchos::null;
  using Teuchos::tuple;
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;

  RTOpPack::TOpAXPY<Scalar> axpyOp(ST::one());

  ConstSubVectorView<Scalar> sv = newSubVectorView<Scalar>(n, ST::random());
  SubVectorView<Scalar> tsv = newSubVectorView<Scalar>(n, ST::random());
  tsv.setGlobalOffset(1);
  TEST_THROW(
    RTOpPack::validate_apply_op<Scalar>( axpyOp, 1, 1, false,
      tuple(sv)(), tuple(tsv)(), null ),
    RTOpPack::IncompatibleVecs
    );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( validate_apply_op,
  inCompatibleSubVecs_1_1_b)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( validate_apply_op, inCompatibleSubVecs_2_0_a, Scalar )
{

  using Teuchos::null;
  using Teuchos::tuple;
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;

  RTOpPack::TOpAXPY<Scalar> axpyOp(ST::one());

  ConstSubVectorView<Scalar> sv1 = newSubVectorView<Scalar>(n, ST::random());
  ConstSubVectorView<Scalar> sv2 = newSubVectorView<Scalar>(n+1, ST::random());
  TEST_THROW(
    RTOpPack::validate_apply_op<Scalar>( axpyOp, 2, 0, false,
      tuple(sv1, sv2)(), null, null ),
    RTOpPack::IncompatibleVecs
    );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( validate_apply_op,
  inCompatibleSubVecs_2_0_a)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( validate_apply_op, inCompatibleSubVecs_2_0_b, Scalar )
{

  using Teuchos::null;
  using Teuchos::tuple;
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;

  RTOpPack::TOpAXPY<Scalar> axpyOp(ST::one());

  ConstSubVectorView<Scalar> sv1 = newSubVectorView<Scalar>(n, ST::random());
  ConstSubVectorView<Scalar> sv2 = newSubVectorView<Scalar>(n, ST::random());
  sv2.setGlobalOffset(1);
  TEST_THROW(
    RTOpPack::validate_apply_op<Scalar>( axpyOp, 2, 0, false,
      tuple(sv1, sv2)(), null, null ),
    RTOpPack::IncompatibleVecs
    );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( validate_apply_op,
  inCompatibleSubVecs_2_0_b)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( validate_apply_op, missingReductObj, Scalar )
{

  using Teuchos::null;
  using Teuchos::tuple;
  typedef ScalarTraits<Scalar> ST;

  RTOpPack::ROpSum<Scalar> sumOp;

  ConstSubVectorView<Scalar> sv = newSubVectorView<Scalar>(n, ST::random());
  RCP<RTOpPack::ReductTarget> sum = sumOp.reduct_obj_create();
  TEST_THROW(
    RTOpPack::validate_apply_op<Scalar>( sumOp, 1, 0, true,
      tuple(sv)(), null, null ),
    RTOpPack::IncompatibleReductObj
    );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( validate_apply_op,
  missingReductObj )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( validate_apply_op, incompatibleReductObj, Scalar )
{

  using Teuchos::null;
  using Teuchos::tuple;
  typedef ScalarTraits<Scalar> ST;

  RTOpPack::ROpSum<Scalar> sumOp;

  ConstSubVectorView<Scalar> sv = newSubVectorView<Scalar>(n, ST::random());
  RCP<RTOpPack::ReductTarget> sum =
    Teuchos::rcp(new RTOpPack::DefaultReductTarget<RTOpPack::index_type>(0));
  TEST_THROW(
    RTOpPack::validate_apply_op<Scalar>( sumOp, 1, 0, true,
      tuple(sv)(), null, sum.ptr() ),
    RTOpPack::IncompatibleReductObj
    );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( validate_apply_op,
  incompatibleReductObj )


} // namespace
