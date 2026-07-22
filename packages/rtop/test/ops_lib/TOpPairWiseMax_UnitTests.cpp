// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RTOpPack_TOpPairWiseMax.hpp"

#include "supportUnitTestsHelpers.hpp"


namespace {


template<class Scalar>
void basicTest(const int stride, FancyOStream &out, bool &success)
{

  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;

  const int local_n = 4;

  SubVectorView<Scalar> v0 =
    newStridedSubVectorView<Scalar>(local_n, stride, ST::nan());

  SubVectorView<Scalar> v1 =
    newStridedSubVectorView<Scalar>(local_n, stride, ST::nan());

  const Scalar alpha = 1.5*ST::one();

  SubVectorView<Scalar> expected_z0 =
    newStridedSubVectorView<Scalar>(local_n, stride, ST::nan());

  v0[0] = 1.0;  v1[0] = 2.0;  expected_z0[0] = 3.0;
  v0[1] = -1.0; v1[1] = -2.0; expected_z0[1] = -1.5;
  v0[2] = 1.0;  v1[2] = -2.0; expected_z0[2] = 1.5;
  v0[3] = 2.0;  v1[3] = 1.0;  expected_z0[3] = 3.0;
  
  SubVectorView<Scalar> z0 =
    newStridedSubVectorView<Scalar>(n, stride, ST::nan());

  RTOpPack::TOpPairWiseMax<Scalar> op(alpha);
  op.apply_op( tuple<ConstSubVectorView<Scalar>>(v0,v1), tuple(z0)(), null );

  if (verbose) {
    out << "alpha = " << alpha << "\n";
    dumpSubVectorView(v0, "v0", out);
    dumpSubVectorView(v1, "v1", out);
    dumpSubVectorView(z0, "z0", out);
    dumpSubVectorView(expected_z0, "expected_z0", out);
  }

  TEST_COMPARE_FLOATING_ARRAYS( constSubVectorViewAsArray(z0),
    constSubVectorViewAsArray(expected_z0),
    as<ScalarMag>(ST::eps() * errorTolSlack)
    );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpPairWiseMax, unitStride, Scalar )
{
  basicTest<Scalar>(1, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( TOpPairWiseMax, unitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpPairWiseMax, nonunitStride, Scalar )
{
  basicTest<Scalar>(4, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( TOpPairWiseMax, nonunitStride )


} // namespace
