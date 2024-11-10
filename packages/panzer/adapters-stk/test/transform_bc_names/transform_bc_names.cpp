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
  // Cubit converts whitespace to underscore when writing to exodus.
  // STK converts capital letters to lower case when reading exodus.
  std::string name = "My Awesome BC";
  const std::string return_name = panzer_stk::transformBCNameForIOSS(name);
  const std::string gold = "my_awesome_bc";
  TEST_EQUALITY(name, gold);
  TEST_EQUALITY(return_name, gold);
}
