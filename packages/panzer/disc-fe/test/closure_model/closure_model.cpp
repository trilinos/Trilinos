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
#include <Teuchos_ParameterList.hpp>

#include "Panzer_CellData.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_FieldLibrary.hpp"
#include "Phalanx_FieldManager.hpp"
#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_ClosureModel_Factory.hpp"
#include <iostream>
#include <vector>

namespace panzer {

  TEUCHOS_UNIT_TEST(evaluator_factory, basic_construction)
  {

    panzer::FieldLayoutLibrary fl;
    Teuchos::RCP<panzer::IntegrationRule> ir;
    {
      Teuchos::RCP<shards::CellTopology> topo = 
         Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<8> >()));
    
      const int num_cells = 20;
      const panzer::CellData cell_data(num_cells,topo);
      const int cubature_degree = 2;      
      ir = Teuchos::rcp(new panzer::IntegrationRule(cubature_degree, cell_data));
      Teuchos::RCP<panzer::BasisIRLayout> basis = Teuchos::rcp(new panzer::BasisIRLayout("Q1",0,*ir));
      
      fl.addFieldAndLayout("Ux",basis);
    }

    std::string model_id = "fluid model";

    Teuchos::ParameterList eqset_params; 

    Teuchos::ParameterList p("Closure Models");
    {
      p.sublist("fluid model").sublist("Density").set<double>("Value",1.0);
      p.sublist("fluid model").sublist("Viscosity").set<double>("Value",1.0);
      p.sublist("fluid model").sublist("Heat Capacity").set<double>("Value",1.0);
      p.sublist("fluid model").sublist("Thermal Conductivity").set<double>("Value",1.0);
    }

    user_app::MyModelFactory<panzer::Traits::Residual> mf;

    Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > evaluators;

    Teuchos::ParameterList user_data("User Data");

    Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

    PHX::FieldManager<panzer::Traits> fm;

    evaluators = mf.buildClosureModels(model_id, p, fl, ir, eqset_params, user_data, gd, fm);

    TEST_EQUALITY(evaluators->size(), 8);

    user_app::MyModelFactory_TemplateBuilder builder;
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> model_factory;
    model_factory.buildObjects(builder);
    evaluators = model_factory.getAsObject<panzer::Traits::Residual>()->buildClosureModels(model_id, p, fl, ir, eqset_params, user_data, gd, fm);

    TEST_EQUALITY(evaluators->size(), 8);

    // Add an unsupported type
    p.sublist("fluid model").sublist("garbage").set<std::string>("Value","help!");
    
    TEST_THROW(model_factory.getAsObject<panzer::Traits::Residual>()->buildClosureModels(model_id, p, fl, ir, eqset_params, user_data, gd, fm), std::logic_error);

  }

}
