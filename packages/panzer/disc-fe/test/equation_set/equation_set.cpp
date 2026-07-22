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
#include "user_app_EquationSetFactory.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(equation_set, steady_state)
  {

    Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
    p->set("Type","Energy");
    p->set("Prefix","ION_");
    p->set("Model ID","solid");
    p->set("Basis Type","HGrad");
    p->set("Basis Order",1);

    int default_integration_order = 1;    
    int num_cells = 20;
    panzer::CellData cell_data(num_cells,Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >())));

    Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > eq_set;
  
    user_app::MyFactory my_factory;
    Teuchos::RCP<panzer::GlobalData> global_data = panzer::createGlobalData();
    TEST_NOTHROW(eq_set = my_factory.buildEquationSet(p, default_integration_order, cell_data, global_data, false));
  }

  TEUCHOS_UNIT_TEST(equation_set, transient)
  {

    Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
    p->set("Type","Energy");
    p->set("Prefix","ION_");
    p->set("Model ID","solid");
    p->set("Basis Type","HGrad");
    p->set("Basis Order",1);

    int default_integration_order = 1; 
    int num_cells = 20;
    panzer::CellData cell_data(num_cells,
                        Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >())));

    Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > eq_set;
  
    user_app::MyFactory my_factory;
    Teuchos::RCP<panzer::GlobalData> global_data = panzer::createGlobalData();
    TEST_NOTHROW(eq_set = my_factory.buildEquationSet(p, default_integration_order, cell_data, global_data, true));
  }

}
