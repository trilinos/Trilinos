// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RTOpPack_TOpAssignVectors.hpp"

#include "supportUnitTestsHelpers.hpp"


namespace {


template<class Scalar>
void basicTest(const int stride, FancyOStream &out, bool &success)
{

  typedef ScalarTraits<Scalar> ST;

  ConstSubVectorView<Scalar> v0 =
    newStridedRandomSubVectorView<Scalar>(n, stride);

  SubVectorView<Scalar> z0 =
    newStridedSubVectorView<Scalar>(n, stride, ST::nan());
  
  RTOpPack::TOpAssignVectors<Scalar> op;
  op.apply_op( tuple(v0), tuple(z0)(), null );

  SubVectorView<Scalar> expected_z0
    = newStridedSubVectorView<Scalar>(n, stride, ST::nan());
  for (int k = 0; k < v0.subDim(); ++k)
    expected_z0[k] = v0[k];

  if (verbose) {
    dumpSubVectorView(v0, "v0", out);
    dumpSubVectorView(z0, "z0", out);
    dumpSubVectorView(expected_z0, "expected_z0", out);
  }

  TEST_COMPARE_ARRAYS( constSubVectorViewAsArray(z0),
    constSubVectorViewAsArray(expected_z0) );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpAssignVectors, unitStride, Scalar )
{
  basicTest<Scalar>(1, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpAssignVectors, unitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpAssignVectors, nonunitStride, Scalar )
{
  basicTest<Scalar>(4, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpAssignVectors, nonunitStride )


} // namespace
