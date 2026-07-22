// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RTOpPack_TOpAssignScalar.hpp"

#include "supportUnitTestsHelpers.hpp"


namespace {


template<class Scalar>
void basicTest(const int stride, FancyOStream &out, bool &success)
{

  typedef ScalarTraits<Scalar> ST;

  SubVectorView<Scalar> sv =
    newStridedSubVectorView<Scalar>(n, stride, ST::nan());

  const Scalar alpha = ST::random();
  
  RTOpPack::TOpAssignScalar<Scalar> addScalarOp(alpha);
  addScalarOp.apply_op( null, tuple(sv)(), null );

  SubVectorView<Scalar> expectedSv
    = newStridedSubVectorView<Scalar>(n, stride, alpha);

  if (verbose) {
    out << "alpha = " << alpha << "\n";
    dumpSubVectorView(sv, "sv", out);
    dumpSubVectorView(expectedSv, "expectedSv", out);
  }

  TEST_COMPARE_ARRAYS( constSubVectorViewAsArray(sv),
    constSubVectorViewAsArray(expectedSv) );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpAssignScalar, unitStride, Scalar )
{
  basicTest<Scalar>(1, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpAssignScalar, unitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpAssignScalar, nonunitStride, Scalar )
{
  basicTest<Scalar>(4, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpAssignScalar, nonunitStride )


} // namespace
