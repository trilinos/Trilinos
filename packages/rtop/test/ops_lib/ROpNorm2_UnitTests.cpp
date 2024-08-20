// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RTOpPack_ROpNorm2.hpp"
#include "supportUnitTestsHelpers.hpp"


namespace {


template<class Scalar>
void basicTest(const int stride, FancyOStream &out, bool &success)
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef ScalarTraits<ScalarMag> SMT;
  const Scalar v = ST::random();
  out << "v="<<v<<"\n";
  ConstSubVectorView<Scalar> sv = newStridedSubVectorView<Scalar>(n, stride, v);
  RTOpPack::ROpNorm2<Scalar> norm2Op;
  RCP<RTOpPack::ReductTarget> norm2 = norm2Op.reduct_obj_create();
  norm2Op.apply_op(
    tuple(sv)(),
    Teuchos::null,
    norm2.ptr()
    );
  TEST_FLOATING_EQUALITY( norm2Op(*norm2),
    SMT::squareroot(n)*ST::magnitude(v),
    as<ScalarMag>(ST::eps() * errorTolSlack) );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpNorm2, unitStride, Scalar )
{
  basicTest<Scalar>(1, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpNorm2, unitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpNorm2, nonunitStride, Scalar )
{
  basicTest<Scalar>(4, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpNorm2, nonunitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpNorm2, reduct, Scalar )
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

  RTOpPack::ROpNorm2<Scalar> norm2Op;

  RCP<RTOpPack::ReductTarget> reduct1 = norm2Op.reduct_obj_create();
  RCP<RTOpPack::ReductTarget> reduct2 = norm2Op.reduct_obj_create();

  DefaultReductTarget<Scalar> &scalarReduct1 =
    dyn_cast<DefaultReductTarget<Scalar> >(*reduct1); 
  DefaultReductTarget<Scalar> &scalarReduct2 =
    dyn_cast<DefaultReductTarget<Scalar> >(*reduct2); 

  scalarReduct1.set(three);
  scalarReduct2.set(four);
  norm2Op.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(two);
  norm2Op.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  TEST_FLOATING_EQUALITY( norm2Op(*reduct2), SMT::squareroot(three+four+two),
    as<ScalarMag>(ST::eps() * errorTolSlack) );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpNorm2, reduct )


} // namespace
