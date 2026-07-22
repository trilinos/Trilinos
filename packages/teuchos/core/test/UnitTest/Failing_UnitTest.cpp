// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"


namespace {


TEUCHOS_UNIT_TEST( Int,  BadAssignment )
{
  int i1 = 4;
  int i2 = i1 + 1;
  TEST_EQUALITY( i2, i1 );
}


TEUCHOS_UNIT_TEST( VectorInt, OutOfRangeAt )
{
  const size_t n = 1;
  std::vector<int> v(n);
  const int i = v.at(n); // Should throw std::out_of_range!
  TEST_EQUALITY_CONST( i, 10 ); // Will never be executed!
}


} // namespace
