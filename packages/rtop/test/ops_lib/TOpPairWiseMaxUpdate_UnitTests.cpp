// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RTOpPack_TOpPairWiseMaxUpdate.hpp"

#include "supportUnitTestsHelpers.hpp"


namespace {


template<class Scalar>
void basicTest(const int stride, FancyOStream &out, bool &success)
{

  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;

  const int local_n = 4;

  SubVectorView<Scalar> z0 =
    newStridedSubVectorView<Scalar>(local_n, stride, ST::nan());

  SubVectorView<Scalar> v0 =
    newStridedSubVectorView<Scalar>(local_n, stride, ST::nan());

  const Scalar alpha = 1.5*ST::one();

  SubVectorView<Scalar> expected_z0 =
    newStridedSubVectorView<Scalar>(local_n, stride, ST::nan());

  z0[0] = 1.0;  v0[0] = 2.0;  expected_z0[0] = 3.0;
  z0[1] = -1.0; v0[1] = -2.0; expected_z0[1] = -1.5;
  z0[2] = 1.0;  v0[2] = -2.0; expected_z0[2] = 1.5;
  z0[3] = 2.0;  v0[3] = 1.0;  expected_z0[3] = 3.0;

  SubVectorView<Scalar> orig_z0 =
    newStridedSubVectorView<Scalar>(local_n, stride, ST::nan());
  RTOpPack::assign_entries<Scalar>(outArg(orig_z0), z0);

  RTOpPack::TOpPairWiseMaxUpdate<Scalar> op(alpha);
  op.apply_op( tuple<ConstSubVectorView<Scalar>>(v0), tuple(z0)(), null );

  if (verbose) {
    out << "alpha = " << alpha << "\n";
    dumpSubVectorView(v0, "v0", out);
    dumpSubVectorView(z0, "z0", out);
    dumpSubVectorView(orig_z0, "orig)z0", out);
    dumpSubVectorView(expected_z0, "expected_z0", out);
  }

  TEST_COMPARE_FLOATING_ARRAYS( constSubVectorViewAsArray(z0),
    constSubVectorViewAsArray(expected_z0),
    as<ScalarMag>(ST::eps() * errorTolSlack)
    );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpPairWiseMaxUpdate, unitStride, Scalar )
{
  basicTest<Scalar>(1, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( TOpPairWiseMaxUpdate, unitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpPairWiseMaxUpdate, nonunitStride, Scalar )
{
  basicTest<Scalar>(4, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( TOpPairWiseMaxUpdate, nonunitStride )


} // namespace
