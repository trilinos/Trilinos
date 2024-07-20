// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RTOpPack_TOpEleWiseProd.hpp"

#include "supportUnitTestsHelpers.hpp"


namespace {


template<class Scalar>
void basicTest(const int stride, FancyOStream &out, bool &success)
{

  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;

  ConstSubVectorView<Scalar> v0 =
    newStridedRandomSubVectorView<Scalar>(n, stride);

  ConstSubVectorView<Scalar> v1 =
    newStridedRandomSubVectorView<Scalar>(n, stride);

  SubVectorView<Scalar> orig_z0 =
    newStridedRandomSubVectorView<Scalar>(n, stride);
  
  SubVectorView<Scalar> z0 =
    newStridedSubVectorView<Scalar>(n, stride, ST::nan());
  RTOpPack::assign_entries<Scalar>( Teuchos::outArg(z0), orig_z0 );

  const Scalar alpha = ST::random();

  RTOpPack::TOpEleWiseProd<Scalar> op(alpha);
  op.apply_op( tuple(v0,v1), tuple(z0)(), null );

  SubVectorView<Scalar> expected_z0
    = newStridedSubVectorView<Scalar>(n, stride, ST::nan());
  for (int k = 0; k < v0.subDim(); ++k)
    expected_z0[k] = orig_z0[k] + alpha * v0[k] * v1[k];

  if (verbose) {
    out << "alpha = " << alpha << "\n";
    dumpSubVectorView(v0, "v0", out);
    dumpSubVectorView(v1, "v1", out);
    dumpSubVectorView(orig_z0, "orig_z0", out);
    dumpSubVectorView(z0, "z0", out);
    dumpSubVectorView(expected_z0, "expected_z0", out);
  }

  TEST_COMPARE_FLOATING_ARRAYS( constSubVectorViewAsArray(z0),
    constSubVectorViewAsArray(expected_z0),
    as<ScalarMag>(ST::eps() * errorTolSlack)
    );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpEleWiseProd, unitStride, Scalar )
{
  basicTest<Scalar>(1, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpEleWiseProd, unitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpEleWiseProd, nonunitStride, Scalar )
{
  basicTest<Scalar>(4, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpEleWiseProd, nonunitStride )


} // namespace
