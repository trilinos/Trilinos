// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RTOpPack_TOpScaleVector.hpp"

#include "supportUnitTestsHelpers.hpp"


namespace {


template<class Scalar>
void basicTest(const int stride, FancyOStream &out, bool &success)
{

  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;

  ConstSubVectorView<Scalar> origSv =
    newStridedRandomSubVectorView<Scalar>(n, stride);
  SubVectorView<Scalar> sv =
    newStridedSubVectorView<Scalar>(n, stride, ST::nan());
  RTOpPack::assign_entries<Scalar>( outArg(sv), origSv );

  const Scalar alpha = ST::random();
  
  RTOpPack::TOpScaleVector<Scalar> addScalarOp(alpha);
  addScalarOp.apply_op( null, tuple(sv)(), null );

  SubVectorView<Scalar> expectedSv
    = newStridedSubVectorView<Scalar>(n, stride, ST::nan());
  for (int k = 0; k < sv.subDim(); ++k)
    expectedSv[k] = alpha * origSv[k];

  if (verbose) {
    dumpSubVectorView(origSv, "origSv", out);
    out << "alpha = " << alpha << "\n";
    dumpSubVectorView(sv, "sv", out);
    dumpSubVectorView(expectedSv, "expectedSv", out);
  }

  TEST_COMPARE_FLOATING_ARRAYS( constSubVectorViewAsArray(sv),
    constSubVectorViewAsArray(expectedSv),
    as<ScalarMag>(ST::eps() * errorTolSlack)
    );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpScaleVector, unitStride, Scalar )
{
  basicTest<Scalar>(1, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpScaleVector, unitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpScaleVector, nonunitStride, Scalar )
{
  basicTest<Scalar>(4, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpScaleVector, nonunitStride )


} // namespace
