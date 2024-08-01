// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RTOpPack_TOpSetSubVector.hpp"

#include "supportUnitTestsHelpers.hpp"


namespace {


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpSetSubVector, dense, Scalar )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;

  SubVectorView<Scalar> input_sv = newSubVectorView<Scalar>(n, ST::zero());
  for (index_type k = 0; k < n; ++k)
    input_sv(k) = ST::random();
  
  RTOpPack::TOpSetSubVector<Scalar> setSubVectorOp;
  setSubVectorOp.set_sub_vec(input_sv);
  // 2008/07/23: rabartl: Above: For some reason, I have to use the default
  // constructor for this class object and then set the subvector or gcc 3.4.6
  // issues an error and thinks that setSubVectorOp is an operator function.
  // This must be a compiler bug!

  SubVectorView<Scalar> sv = newSubVectorView<Scalar>(n, ST::nan());
  setSubVectorOp.apply_op( null, tuple(sv)(), null );

  if (verbose) {
    dumpSubVectorView(input_sv, "input_sv", out);
    dumpSubVectorView(sv, "sv", out);
  }

  TEST_COMPARE_ARRAYS( constSubVectorViewAsArray(input_sv),
    constSubVectorViewAsArray(sv) );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpSetSubVector, dense )


} // namespace
