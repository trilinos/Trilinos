// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_as.hpp"


template<class T>
T someFunc( const T &t )
{
  return t*t;
}


namespace {


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( someFunc, test1, T )
{
  using Teuchos::as;
  T t1 = 5;
  TEST_EQUALITY_CONST( t1, as<T>(5) );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( someFunc, test1)


} // namespace
