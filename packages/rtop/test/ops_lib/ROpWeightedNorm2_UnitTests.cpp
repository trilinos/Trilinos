// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RTOpPack_ROpWeightedNorm2.hpp"
#include "Teuchos_implicit_cast.hpp"

#include "supportUnitTestsHelpers.hpp"


namespace {


template<class Scalar>
void basicTest(const int stride, FancyOStream &out, bool &success)
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef ScalarTraits<ScalarMag> SMT;
  const ScalarMag w = as<ScalarMag>(2.0);
  const Scalar v = ST::random();
  out << "w="<<w<<"\n";
  out << "v="<<v<<"\n";
  ConstSubVectorView<Scalar> svw = newStridedSubVectorView<Scalar>(n, stride, w);
  ConstSubVectorView<Scalar> svv = newStridedSubVectorView<Scalar>(n, stride, v);
  RTOpPack::ROpWeightedNorm2<Scalar> weightedNorm2Op;
  RCP<RTOpPack::ReductTarget> weightedNorm2 = weightedNorm2Op.reduct_obj_create();
  Teuchos::implicit_ref_cast<RTOpPack::RTOpT<Scalar> >(weightedNorm2Op).apply_op(
    tuple(svw, svv)(),
    Teuchos::null,
    weightedNorm2.ptr()
    );
  TEST_FLOATING_EQUALITY( weightedNorm2Op(*weightedNorm2),
    SMT::squareroot(n)*SMT::squareroot(w)*ST::magnitude(v),
    as<ScalarMag>(ST::eps() * errorTolSlack / n) );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpWeightedNorm2, unitStride, Scalar )
{
  basicTest<Scalar>(1, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpWeightedNorm2, unitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpWeightedNorm2, nonunitStride, Scalar )
{
  basicTest<Scalar>(4, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpWeightedNorm2, nonunitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpWeightedNorm2, reduct, Scalar )
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

  RTOpPack::ROpWeightedNorm2<Scalar> weightedNorm2Op;

  RCP<RTOpPack::ReductTarget> reduct1 = weightedNorm2Op.reduct_obj_create();
  RCP<RTOpPack::ReductTarget> reduct2 = weightedNorm2Op.reduct_obj_create();

  DefaultReductTarget<Scalar> &scalarReduct1 =
    dyn_cast<DefaultReductTarget<Scalar> >(*reduct1); 
  DefaultReductTarget<Scalar> &scalarReduct2 =
    dyn_cast<DefaultReductTarget<Scalar> >(*reduct2); 

  scalarReduct1.set(three);
  scalarReduct2.set(four);
  weightedNorm2Op.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(two);
  weightedNorm2Op.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  TEST_FLOATING_EQUALITY( weightedNorm2Op(*reduct2), SMT::squareroot(three+four+two),
    as<ScalarMag>(ST::eps() * errorTolSlack) );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpWeightedNorm2, reduct )


} // namespace
