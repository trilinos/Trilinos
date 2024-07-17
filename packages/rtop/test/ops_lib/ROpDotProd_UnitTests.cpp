// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RTOpPack_ROpDotProd.hpp"
#include "supportUnitTestsHelpers.hpp"


namespace {


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpDotProd, unitStride, Scalar )
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  const Scalar v1 = ST::random();
  const Scalar v2 = ST::random();
  out << "v1="<<v1<<", v2="<<v2<<"\n";
  SubVectorView<Scalar> sv1 = newSubVectorView<Scalar>(n, v1);
  SubVectorView<Scalar> sv2 = newSubVectorView<Scalar>(n, v2);
  RTOpPack::ROpDotProd<Scalar> dotProdOp;
  RCP<RTOpPack::ReductTarget> dotProd = dotProdOp.reduct_obj_create();
  dotProdOp.apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv1, sv2)(),
    Teuchos::null,
    dotProd.ptr()
    );
  TEST_FLOATING_EQUALITY( dotProdOp(*dotProd), ST::conjugate(v1)*v2*as<Scalar>(n),
    as<ScalarMag>(ST::eps() * errorTolSlack) );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpDotProd, unitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpDotProd, nonunitStride, Scalar )
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  const Scalar v1 = ST::random();
  const Scalar v2 = ST::random();
  out << "v1="<<v1<<", v2="<<v2<<"\n";
  SubVectorView<Scalar> sv1 = newStridedSubVectorView<Scalar>(n, 2, v1);
  SubVectorView<Scalar> sv2 = newStridedSubVectorView<Scalar>(n, 3, v2);
  RTOpPack::ROpDotProd<Scalar> dotProdOp;
  RCP<RTOpPack::ReductTarget> dotProd = dotProdOp.reduct_obj_create();
  dotProdOp.apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv1, sv2)(),
    Teuchos::null,
    dotProd.ptr()
    );
  TEST_FLOATING_EQUALITY( dotProdOp(*dotProd), ST::conjugate(v1)*v2*as<Scalar>(n),
    as<ScalarMag>(ST::eps() * errorTolSlack) );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpDotProd, nonunitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpDotProd, reduct, Scalar )
{
  using Teuchos::as;
  using Teuchos::dyn_cast;
  using RTOpPack::DefaultReductTarget;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;

  const Scalar v1 = ST::random();
  const Scalar v2 = ST::random();
  out << "v1="<<v1<<", v2="<<v2<<"\n";

  RTOpPack::ROpDotProd<Scalar> dotProdOp;

  RCP<RTOpPack::ReductTarget> reduct1 = dotProdOp.reduct_obj_create();
  RCP<RTOpPack::ReductTarget> reduct2 = dotProdOp.reduct_obj_create();

  DefaultReductTarget<Scalar> &scalarReduct1 =
    dyn_cast<DefaultReductTarget<Scalar> >(*reduct1); 
  DefaultReductTarget<Scalar> &scalarReduct2 =
    dyn_cast<DefaultReductTarget<Scalar> >(*reduct2); 

  scalarReduct1.set(v1);
  scalarReduct2.set(v2);

  dotProdOp.reduce_reduct_objs( *reduct1, reduct2.ptr() );

  TEST_FLOATING_EQUALITY( scalarReduct2.get(), v1+v2, as<ScalarMag>(ST::eps()*errorTolSlack) );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpDotProd, reduct )


} // namespace
