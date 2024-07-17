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

#include "Panzer_CellData.hpp"
#include "user_app_EquationSetFactory_Physics1.hpp"
#include "user_app_EquationSetFactory_Physics2.hpp"
#include "Panzer_EquationSet_Factory_Composite.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(equation_set, composite_factory)
  {    

    const int num_cells = 20;
    const int default_int_order = 2;

    panzer::CellData cell_data(num_cells,Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >())));

    Teuchos::RCP<panzer::GlobalData> global_data = panzer::createGlobalData();

    Teuchos::RCP<Teuchos::ParameterList> p1 = Teuchos::parameterList();
    p1->set("Type","Energy 1");
    p1->set("Prefix","ION_");
    p1->set("Model ID","solid");
    p1->set("Basis Type","HGrad");
    p1->set("Basis Order",2);

    Teuchos::RCP<Teuchos::ParameterList> p2 = Teuchos::parameterList();
    p2->set("Type","Energy 2");
    p2->set("Prefix","ION_");
    p2->set("Model ID","solid");
    p2->set("Basis Type","HGrad");
    p2->set("Basis Order",2);

    Teuchos::RCP<Teuchos::ParameterList> p_bad = Teuchos::parameterList();
    p_bad->set("Type","UNKNOWN");
    p_bad->set("Prefix","ION_");
    p_bad->set("Model ID","solid");
    p_bad->set("Basis Type","HGrad");
    p_bad->set("Basis Order",2);

    Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > eq_set;
  
    // Check that eqf1 can build Energy 1 and fails to build Energy 2 
    {
      user_app::MyFactory1 my_factory1(true);
      TEST_NOTHROW(eq_set = my_factory1.buildEquationSet(p1, default_int_order, cell_data, global_data, false));
      TEST_ASSERT(nonnull(eq_set));
      TEST_THROW(eq_set = my_factory1.buildEquationSet(p2, default_int_order, cell_data, global_data, false),std::logic_error);
    }

    // Check that eqf2 can build Energy 2 and fails to build Energy 1 
    {
      eq_set = Teuchos::null;
      user_app::MyFactory2 my_factory2(true);
      TEST_NOTHROW(eq_set = my_factory2.buildEquationSet(p2, default_int_order, cell_data, global_data, false));
      TEST_ASSERT(nonnull(eq_set));
      TEST_THROW(eq_set = my_factory2.buildEquationSet(p1, default_int_order, cell_data, global_data, false),std::logic_error);
    }

    Teuchos::RCP<panzer::EquationSetFactory> eqsf1 = Teuchos::rcp(new user_app::MyFactory1(false));
    Teuchos::RCP<panzer::EquationSetFactory> eqsf2 = Teuchos::rcp(new user_app::MyFactory2(false));

    std::vector<Teuchos::RCP<panzer::EquationSetFactory> > factories;
    factories.push_back(eqsf1);
    factories.push_back(eqsf2);

    Teuchos::RCP<panzer::EquationSet_FactoryComposite> composite_factory = 
      Teuchos::rcp(new panzer::EquationSet_FactoryComposite(factories));

    // Build each physics
    eq_set = Teuchos::null;
    TEST_NOTHROW(eq_set = composite_factory->buildEquationSet(p1, default_int_order, cell_data, global_data, false));
    TEST_ASSERT(nonnull(eq_set));

    eq_set = Teuchos::null;
    TEST_NOTHROW(eq_set = composite_factory->buildEquationSet(p2, default_int_order, cell_data, global_data, false));
    TEST_ASSERT(nonnull(eq_set));
    
    // try bad eq set
    eq_set = Teuchos::null;
    TEST_THROW(eq_set = composite_factory->buildEquationSet(p_bad, default_int_order, cell_data, global_data, false),std::logic_error);
    TEST_ASSERT(is_null(eq_set));
    
  }

}
