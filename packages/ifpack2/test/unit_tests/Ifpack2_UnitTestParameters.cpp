// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/*! \file Ifpack2_UnitTestParameters.cpp

\brief Ifpack2 Parameters testing.

*/

#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <iostream>
#include <Ifpack2_Parameters.hpp>

TEUCHOS_UNIT_TEST(Ifpack2Parameters, Test0)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  Teuchos::ParameterList params;

  params.set("fact: iluk level-of-fill", (int) 2);

  Teuchos::ParameterList validparams;

  TEST_NOTHROW(Ifpack2::getValidParameters(validparams));

  params.validateParameters(validparams);

  int level_of_fill = 0;

  //call getParameter with a wrong name:
  Ifpack2::getParameter(params, "level-of-fill", level_of_fill);
  TEST_EQUALITY(level_of_fill, 0)

  //call getParameter with a valid name:
  Ifpack2::getParameter(params, "fact: iluk level-of-fill", level_of_fill);
  TEST_EQUALITY(level_of_fill, 2)
}

