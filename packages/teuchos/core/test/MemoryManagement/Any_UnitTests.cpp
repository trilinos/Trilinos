// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_any.hpp"

#include <sstream>

namespace {


TEUCHOS_UNIT_TEST( any, noThrowComparePrintDouble )
{
  double value = 25.0;
  auto a = Teuchos::any(value);
  auto b = Teuchos::any(value);
  TEST_NOTHROW(a == b);
  TEST_EQUALITY_CONST(true, b.same(a));
  std::stringstream ss;
  TEST_NOTHROW(ss << a);
}

TEUCHOS_UNIT_TEST( any, throwPrintVector )
{
  std::vector<double> value;
  auto a = Teuchos::any(value);
  auto b = Teuchos::any(value);
  TEST_NOTHROW(a == b);
  TEST_EQUALITY_CONST(true, b.same(a));
  std::stringstream ss;
  TEST_THROW(ss << a, std::runtime_error);
  TEST_THROW(ss << b, std::runtime_error);
}

struct NotComparableOrPrintable {
  int x;
};

TEUCHOS_UNIT_TEST( any, throwComparePrintStruct )
{
  NotComparableOrPrintable value;
  value.x = 15;
  auto a = Teuchos::any(value);
  auto b = Teuchos::any(value);
  TEST_THROW(a == b, std::runtime_error);
  std::stringstream ss;
  TEST_THROW(ss << a, std::runtime_error);
  TEST_THROW(ss << b, std::runtime_error);
}

} // namespace
