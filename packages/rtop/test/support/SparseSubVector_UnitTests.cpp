// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include "RTOpPack_SparseSubVectorT.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_as.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace RTOpPack {

//
// Test default constructor
//
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SparseSubVector,
  stridedConstruct, Scalar )
{
  using Teuchos::tuple;
  using Teuchos::arcpFromArray;
  using Teuchos::as;
  using Teuchos::Array;

  Array<Scalar> values = tuple<Scalar>(1.1, -1.0, 2.2, -2.0)().getConst();
  Array<Ordinal> indices = tuple<Ordinal>(1, -1, -2, 3, -3, -4)().getConst();

  SparseSubVectorT<Scalar> ssv(
    1,  // glboalOfset_in
    6,  // subDim_in
    2, // subNz_in
    arcpFromArray(values), // values_in
    2, // valuesStride_in
    arcpFromArray(indices), // indices_in
    3, // indiciesStride_in
    7, // localOffset_in
    true // isSorted_in
    );
  TEST_EQUALITY_CONST(ssv.globalOffset(), as<Ordinal>(1));
  TEST_EQUALITY_CONST(ssv.subDim(), as<Ordinal>(6));
  TEST_EQUALITY_CONST(ssv.subNz(), as<Ordinal>(2));
  TEST_COMPARE_ARRAYS(ssv.values(), values);
  TEST_EQUALITY_CONST(ssv.valuesStride(), as<Ordinal>(2));
  TEST_COMPARE_ARRAYS(ssv.indices(), indices);
  TEST_EQUALITY_CONST(ssv.indicesStride(), as<Ordinal>(3));
  TEST_EQUALITY_CONST(ssv.localOffset(), as<Ordinal>(7));
  TEST_EQUALITY_CONST(ssv.isSorted(), true);
  
  // ToDo: Add more checks ...
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( SparseSubVector,
  stridedConstruct )

} // namespace RTOpPack
