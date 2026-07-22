// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RTOpPack_ROpCountNanInf.hpp"
#include "supportUnitTestsHelpers.hpp"


namespace {


using TestingSupportHelpers::print;


template<class Scalar>
void basicTest(const int stride, FancyOStream &out, bool &success)
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;

  SubVectorView<Scalar> sv = newStridedSubVectorView<Scalar>(4, stride, as<Scalar>(0.0));

  sv(0) = ST::nan();
  sv(1) = ST::random();
  sv(2) = ST::one()/ST::zero();
  sv(3) = ST::zero();

  TEST_ASSERT(ST::isnaninf(sv(0)));
  TEST_ASSERT(!ST::isnaninf(sv(1)));
  TEST_ASSERT(ST::isnaninf(sv(2)));
  TEST_ASSERT(!ST::isnaninf(sv(3)));

  print(sv, "sv", out);

  RTOpPack::ROpCountNanInf<Scalar> countNanInfOp;
  RCP<RTOpPack::ReductTarget> countNanInf = countNanInfOp.reduct_obj_create();
  Teuchos::implicit_ref_cast<RTOpPack::RTOpT<Scalar> >(countNanInfOp).apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv)(),
    Teuchos::null,
    countNanInf.ptr()
    );

  const index_type countNanInf_val = countNanInfOp(*countNanInf);
  TEST_EQUALITY( countNanInf_val, 2 );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpCountNanInf, unitStride, Scalar )
{
  basicTest<Scalar>(1, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpCountNanInf, unitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpCountNanInf, nonunitStride, Scalar )
{
  basicTest<Scalar>(4, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpCountNanInf, nonunitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpCountNanInf, reduct, Scalar )
{
  using Teuchos::as;
  using Teuchos::dyn_cast;
  using RTOpPack::DefaultReductTarget;
  //typedef ScalarTraits<Scalar> ST; // unused

  RTOpPack::ROpCountNanInf<Scalar> countNanInfOp;

  RCP<ReductTarget> reduct1 = countNanInfOp.reduct_obj_create();
  RCP<ReductTarget> reduct2 = countNanInfOp.reduct_obj_create();

  DefaultReductTarget<index_type> &scalarReduct1 =
    dyn_cast<DefaultReductTarget<index_type> >(*reduct1);

  scalarReduct1.set(1);
  countNanInfOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(0);
  countNanInfOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(2);
  countNanInfOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  TEST_EQUALITY( countNanInfOp(*reduct2), as<index_type>(3) )

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpCountNanInf, reduct )


} // namespace
