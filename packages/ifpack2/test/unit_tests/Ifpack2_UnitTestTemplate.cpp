// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/*! \file Ifpack2_UnitTestTemplate.cpp

\brief Ifpack2 Unit testing template.

This file demonstrates how you create a unit test for template code.

*/


#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

template<class T>
T my_trivial_function(T in)
{
  T out = in*in;
  return out;
}

//this macro declares the test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(Ifpack2Group0, Ifpack2Test0, T)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  std::string version = Ifpack2::Version();
  bool empty_version = version.empty();
  TEST_EQUALITY(empty_version, false);

  T input = 5;
  T result = my_trivial_function(input);
  T expected_result = input*input;

  TEST_EQUALITY(result, expected_result);
}

//this macro instantiates and registers the test:
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES(Ifpack2Group0, Ifpack2Test0)

