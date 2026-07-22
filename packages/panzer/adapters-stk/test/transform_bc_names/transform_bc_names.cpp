// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_STK_TransformBCNameForIOSS.hpp"
#include "Teuchos_UnitTestHarness.hpp"

TEUCHOS_UNIT_TEST(TransformBCNameForIOSS, basic)
{
  // Cubit converts whitespace in sideset names to underscore when
  // writing to exodus.  Prior to 05/2026, IOSS/STK converted capital
  // letters to lower case when reading exodus. It now supports upper
  // case. We allow for supporting both as IOSS behavior can be
  // changed through properties.

  // Explicit false
  {
    std::string name = "My Awesome BC";
    const std::string return_name = panzer_stk::transformBCNameForIOSS(name,false);
    const std::string gold_preserving_case = "My_Awesome_BC";
    TEST_EQUALITY(name, gold_preserving_case);
    TEST_EQUALITY(return_name, gold_preserving_case);
  }

  // Explicit true
  {
    std::string name = "My Awesome BC";
    const std::string return_name = panzer_stk::transformBCNameForIOSS(name,true);
    const std::string gold_lower_case = "my_awesome_bc";
    TEST_EQUALITY(name, gold_lower_case);
    TEST_EQUALITY(return_name, gold_lower_case);
  }
  
  // Default (false)
  {
    std::string name = "My Awesome BC";
    const std::string return_name = panzer_stk::transformBCNameForIOSS(name);
    const std::string gold_preserving_case = "My_Awesome_BC";
    TEST_EQUALITY(name, gold_preserving_case);
    TEST_EQUALITY(return_name, gold_preserving_case);
  }
}
