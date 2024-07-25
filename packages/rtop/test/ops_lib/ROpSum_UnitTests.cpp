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
#include "supportUnitTestsHelpers.hpp"


namespace {


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpSum, unitStride, Scalar )
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  const Scalar v = - ST::random();
  out << "n="<<n<<"\n";
  out << "v="<<v<<"\n";
  SubVectorView<Scalar> sv = newSubVectorView<Scalar>(n, v);
  RTOpPack::ROpSum<Scalar> sumOp;
  RCP<RTOpPack::ReductTarget> sum = sumOp.reduct_obj_create();
  sumOp.apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv)(),
    Teuchos::null,
    sum.ptr()
    );
  TEST_FLOATING_EQUALITY( sumOp(*sum), as<Scalar>(as<Scalar>(n)*v),
    as<ScalarMag>(ST::eps() * errorTolSlack * n) );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpSum, unitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpSum, reduct, Scalar )
{
  using Teuchos::as;
  using Teuchos::dyn_cast;
  using RTOpPack::DefaultReductTarget;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef ScalarTraits<ScalarMag> SMT;

  const Scalar v1 = - ST::random();
  const Scalar v2 = + ST::random();
  const Scalar v3 = + ST::random();

  RTOpPack::ROpSum<Scalar> sumOp;

  RCP<RTOpPack::ReductTarget> reduct1 = sumOp.reduct_obj_create();
  RCP<RTOpPack::ReductTarget> reduct2 = sumOp.reduct_obj_create();

  DefaultReductTarget<Scalar> &scalarReduct1 =
    dyn_cast<DefaultReductTarget<Scalar> >(*reduct1); 
  DefaultReductTarget<Scalar> &scalarReduct2 =
    dyn_cast<DefaultReductTarget<Scalar> >(*reduct2); 

  scalarReduct1.set(v1);
  scalarReduct2.set(v2);
  sumOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(v3);
  sumOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  TEST_FLOATING_EQUALITY( sumOp(*reduct2), v1+v2+v3,
    as<ScalarMag>(ST::eps() * errorTolSlack) );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpSum, reduct )


} // namespace
