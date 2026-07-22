// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RTOpPack_TOpRandomize.hpp"

#include "supportUnitTestsHelpers.hpp"


namespace {


template<class Scalar>
void basicTest(const int stride, FancyOStream &out, bool &success)
{

  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  //typedef typename ST::magnitudeType ScalarMag; // unused

  SubVectorView<Scalar> z0 =
    newStridedSubVectorView<Scalar>(n, stride, ST::nan());

  const Scalar l(1.5), u(3.0);

  RTOpPack::TOpRandomize<Scalar> op(l, u);
  op.apply_op( null, tuple(z0)(), null );

  if (verbose) {
    dumpSubVectorView(z0, "z0", out);
  }

  for (int k = 0; k < z0.subDim(); ++k)
  {
    TEST_COMPARE( ST::magnitude(l), <=, ST::magnitude(z0[k]) );
    TEST_COMPARE( ST::magnitude(z0[k]), <=, ST::magnitude(u) );
    // 2008/07/25: The standard rand() function generates some duplicate
    // numbers so we can't check this!
    //if (k > 0)
    //  TEST_COMPARE( z0[k], !=, z0[k-1] );
  }

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpRandomize, unitStride, Scalar )
{
  basicTest<Scalar>(1, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpRandomize, unitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpRandomize, nonunitStride, Scalar )
{
  basicTest<Scalar>(4, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpRandomize, nonunitStride )


} // namespace
