// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RTOpPack_TOpEleWiseScale.hpp"

#include "supportUnitTestsHelpers.hpp"


namespace RTOp {


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpEleWiseScale, basic, Scalar )
{

  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;

  const int stride=1;

  ConstSubVectorView<Scalar> v0 =
    newStridedRandomSubVectorView<Scalar>(n, stride);

  SubVectorView<Scalar> orig_z0 =
    newStridedRandomSubVectorView<Scalar>(n, stride);
  
  SubVectorView<Scalar> z0 =
    newStridedSubVectorView<Scalar>(n, stride, ST::nan());
  RTOpPack::assign_entries<Scalar>( Teuchos::outArg(z0), orig_z0 );

  RTOpPack::TOpEleWiseScale<Scalar> op;
  op.apply_op( tuple(v0), tuple(z0)(), null );

  SubVectorView<Scalar> expected_z0 = 
    newStridedSubVectorView<Scalar>(n, stride, ST::nan());
  for (int k = 0; k < v0.subDim(); ++k)
    expected_z0[k] = orig_z0[k] * v0[k];

  if (verbose) {
    dumpSubVectorView(v0, "v0", out);
    dumpSubVectorView(orig_z0, "orig_z0", out);
    dumpSubVectorView(z0, "z0", out);
    dumpSubVectorView(expected_z0, "expected_z0", out);
  }

  TEST_COMPARE_FLOATING_ARRAYS( constSubVectorViewAsArray(z0),
    constSubVectorViewAsArray(expected_z0),
    as<ScalarMag>(ST::eps() * errorTolSlack)
    );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpEleWiseScale, basic )


} // namespace RTOp
