// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RTOpPack_ROpNorm1.hpp"
#include "supportUnitTestsHelpers.hpp"


namespace {


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpNorm1, unitStride, Scalar )
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  const Scalar v = ST::random();
  out << "v="<<v<<"\n";
  SubVectorView<Scalar> sv = newSubVectorView<Scalar>(n, v);
  RTOpPack::ROpNorm1<Scalar> norm1Op;
  RCP<RTOpPack::ReductTarget> norm1 = norm1Op.reduct_obj_create();
  norm1Op.apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv)(),
    Teuchos::null,
    norm1.ptr()
    );
  TEST_FLOATING_EQUALITY( norm1Op(*norm1), as<ScalarMag>(n*ST::magnitude(v)),
    as<ScalarMag>(ST::eps() * errorTolSlack) );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpNorm1, unitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpNorm1, nonunitStride, Scalar )
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  const Scalar v = ST::random();
  out << "v="<<v<<"\n";
  SubVectorView<Scalar> sv = newStridedSubVectorView<Scalar>(n, 3, v);
  RTOpPack::ROpNorm1<Scalar> norm1Op;
  RCP<RTOpPack::ReductTarget> norm1 = norm1Op.reduct_obj_create();
  norm1Op.apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv)(),
    Teuchos::null,
    norm1.ptr()
    );
  TEST_FLOATING_EQUALITY( norm1Op(*norm1), as<ScalarMag>(n*ST::magnitude(v)),
    as<ScalarMag>(ST::eps() * errorTolSlack) );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpNorm1, nonunitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpNorm1, reduct, Scalar )
{
  using Teuchos::as;
  using Teuchos::dyn_cast;
  using RTOpPack::DefaultReductTarget;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef ScalarTraits<ScalarMag> SMT;

  const ScalarMag three = as<ScalarMag>(3.0);
  const ScalarMag four = as<ScalarMag>(4.0);
  const ScalarMag two = as<ScalarMag>(2.0);

  RTOpPack::ROpNorm1<Scalar> norm1Op;

  RCP<RTOpPack::ReductTarget> reduct1 = norm1Op.reduct_obj_create();
  RCP<RTOpPack::ReductTarget> reduct2 = norm1Op.reduct_obj_create();

  DefaultReductTarget<ScalarMag> &scalarReduct1 =
    dyn_cast<DefaultReductTarget<ScalarMag> >(*reduct1); 
  DefaultReductTarget<ScalarMag> &scalarReduct2 =
    dyn_cast<DefaultReductTarget<ScalarMag> >(*reduct2); 

  scalarReduct1.set(three);
  scalarReduct2.set(four);
  norm1Op.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(two);
  norm1Op.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  TEST_FLOATING_EQUALITY( norm1Op(*reduct2), three+four+two,
    as<ScalarMag>(ST::eps() * errorTolSlack) );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpNorm1, reduct )


} // namespace
