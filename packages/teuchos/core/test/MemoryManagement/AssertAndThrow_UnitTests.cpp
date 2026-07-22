// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_Assert.hpp"


namespace {


TEUCHOS_UNIT_TEST( TEUCHOS_TEST_FOR_EXCEPTION, noThrowExcept )
{
  const int i = 5;
  TEST_NOTHROW(
    TEUCHOS_TEST_FOR_EXCEPTION(i != 5, std::logic_error, "Blah blah blah")
    );
}


TEUCHOS_UNIT_TEST( TEUCHOS_TEST_FOR_EXCEPTION, throwExcept )
{
  const int i = 5;
  TEST_THROW(
    TEUCHOS_TEST_FOR_EXCEPTION(i == 5, std::logic_error, "Blah blah blah"),
    std::logic_error
    );
}


TEUCHOS_UNIT_TEST( TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC, noThrowExcept )
{
  const int i = 5;
  std::string tfecfFuncName("someMethod");
  TEST_NOTHROW(
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(i != 5, std::logic_error, "Blah blah blah")
    );
}


TEUCHOS_UNIT_TEST( TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC, throwExcept )
{
  const int i = 5;
  std::string tfecfFuncName("someMethod");
  TEST_THROW(
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(i == 5, std::logic_error, "Blah blah blah"),
    std::logic_error);
}


TEUCHOS_UNIT_TEST( TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG, noThrowExcept )
{
  const int i = 5;
  TEST_NOTHROW(
    TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(i != 5, std::logic_error, "Blah blah blah")
    );
}


TEUCHOS_UNIT_TEST( TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG, throwExcept )
{
  const int i = 5;
  TEST_THROW(
    TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(i == 5, std::logic_error, "Blah blah blah"),
    std::logic_error);
}


TEUCHOS_UNIT_TEST( TEUCHOS_TEST_FOR_EXCEPT, noThrowExcept )
{
  const int i = 5;
  TEST_NOTHROW(
    TEUCHOS_TEST_FOR_EXCEPT(i != 5)
    );
}


TEUCHOS_UNIT_TEST( TEUCHOS_TEST_FOR_EXCEPT, throwExcept )
{
  const int i = 5;
  TEST_THROW(
    TEUCHOS_TEST_FOR_EXCEPT(i == 5),
    std::logic_error);
}


TEUCHOS_UNIT_TEST( TEUCHOS_TEST_FOR_EXCEPT_MSG, noThrowExcept )
{
  const int i = 5;
  TEST_NOTHROW(
    TEUCHOS_TEST_FOR_EXCEPT_MSG(i != 5, "Blah blah blah")
    );
}


TEUCHOS_UNIT_TEST( TEUCHOS_TEST_FOR_EXCEPT_MSG, throwExcept )
{
  const int i = 5;
  TEST_THROW(
    TEUCHOS_TEST_FOR_EXCEPT_MSG(i == 5, "Blah blah blah"),
    std::logic_error);
}


TEUCHOS_UNIT_TEST( TEUCHOS_TEST_FOR_EXCEPTION_PRINT, noThrowExcept )
{
  const int i = 5;
  TEST_NOTHROW(
    TEUCHOS_TEST_FOR_EXCEPTION_PRINT(i != 5, std::logic_error, "Blah blah blah", &out)
    );
}


TEUCHOS_UNIT_TEST( TEUCHOS_TEST_FOR_EXCEPTION_PRINT, throwExcept )
{
  const int i = 5;
  std::ostringstream oss;
  TEST_THROW(
    TEUCHOS_TEST_FOR_EXCEPTION_PRINT(i == 5, std::logic_error, "Blah blah blah", &oss),
    std::logic_error);
  TEST_INEQUALITY(oss.str().find("Throwing an std::exception of type"), std::string::npos);
  TEST_INEQUALITY(oss.str().find("logic_error"), std::string::npos);
  TEST_INEQUALITY(oss.str().find("AssertAndThrow_UnitTests.cpp"), std::string::npos);
  TEST_INEQUALITY(oss.str().find("Throw number ="), std::string::npos);
  TEST_INEQUALITY(oss.str().find("Throw test that evaluated to true"), std::string::npos);
  TEST_INEQUALITY(oss.str().find("i == 5"), std::string::npos);
  TEST_INEQUALITY(oss.str().find("Blah blah blah"), std::string::npos);
  // NOTE: The above test asserts that the exception message is being built
  // somewhat correctly!
}


TEUCHOS_UNIT_TEST( TEUCHOS_TEST_FOR_EXCEPT_PRINT, noThrowExcept )
{
  const int i = 5;
  TEST_NOTHROW(
    TEUCHOS_TEST_FOR_EXCEPT_PRINT(i != 5, &out)
    );
}


TEUCHOS_UNIT_TEST( TEUCHOS_TEST_FOR_EXCEPT_PRINT, throwExcept )
{
  const int i = 5;
  std::ostringstream oss;
  TEST_THROW(
    TEUCHOS_TEST_FOR_EXCEPT_PRINT(i == 5, &oss),
    std::logic_error);
  TEST_INEQUALITY(oss.str().find("AssertAndThrow_UnitTests.cpp"), std::string::npos);
  TEST_INEQUALITY(oss.str().find("i == 5"), std::string::npos);
}


TEUCHOS_UNIT_TEST( TEUCHOS_TRACE, basic )
{
  std::logic_error localExcept("Blah blah blah");
  TEST_THROW(
    TEUCHOS_TRACE(localExcept),
    std::runtime_error);
}


} // namespace
