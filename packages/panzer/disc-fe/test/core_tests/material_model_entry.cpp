// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Panzer_MaterialModelEntry.hpp"

#include <sstream>

namespace panzer {

  TEUCHOS_UNIT_TEST(material_model_entry, no_params)
  {
    
    std::string factory_name_m1 = "one";
    std::string factory_name_m2 = "two";
    
    panzer::MaterialModelEntry m1(factory_name_m1);
    panzer::MaterialModelEntry m2(factory_name_m2);

    panzer::MaterialModelEntry copy_m1(m1);

    panzer::MaterialModelEntry copy_m2;
    copy_m2 = m2;

    TEST_EQUALITY(m1.factoryName(), factory_name_m1);
    TEST_EQUALITY(m2.factoryName(), factory_name_m2);

    TEST_EQUALITY(copy_m1.factoryName(), factory_name_m1);
    TEST_EQUALITY(copy_m2.factoryName(), factory_name_m2);

    TEST_ASSERT(m1 == copy_m1);
    TEST_ASSERT(m2 == copy_m2);
    TEST_ASSERT(m1 != m2);
    TEST_ASSERT(copy_m1 != copy_m2);

    std::stringstream s;
    s << m1;
    TEST_EQUALITY(s.str(),"Material Model Entry: one"); 
  }

  TEUCHOS_UNIT_TEST(material_model_entry, with_params)
  {
    
    std::string factory_name = "one";
    Teuchos::ParameterList p1;
    p1.set<double>("value", 1.0);
    Teuchos::ParameterList p2;
    p2.set<int>("value", 1);
    

    panzer::MaterialModelEntry m1(factory_name,p1);
    panzer::MaterialModelEntry m2(factory_name,p2);

    panzer::MaterialModelEntry copy_m1(m1);

    panzer::MaterialModelEntry copy_m2;
    copy_m2 = m2;

    TEST_EQUALITY(m1.factoryName(), factory_name);
    TEST_EQUALITY(m2.factoryName(), factory_name);

    TEST_EQUALITY(copy_m1.factoryName(), factory_name);
    TEST_EQUALITY(copy_m2.factoryName(), factory_name);

    TEST_EQUALITY(m1.params(), p1);
    TEST_EQUALITY(m2.params(), p2);
    TEST_INEQUALITY(m1.params(), p2);

    TEST_EQUALITY(copy_m1.params(), p1);
    TEST_EQUALITY(copy_m2.params(), p2);

    TEST_ASSERT(m1 == copy_m1);
    TEST_ASSERT(m2 == copy_m2);
    TEST_ASSERT(m1 != m2);
    TEST_ASSERT(copy_m1 != copy_m2);

  }


}
