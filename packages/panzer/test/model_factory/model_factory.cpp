#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_ParameterList.hpp>

#include "Panzer_CellData.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Traits.hpp"
#include "user_app_ModelFactory.hpp"
#include "Panzer_ModelFactory_TemplateManager.hpp"
#include "user_app_ModelFactory_TemplateBuilder.hpp"
#include <iostream>
#include <vector>

namespace panzer {

  TEUCHOS_UNIT_TEST(evaluator_factory, basic_construction)
  {
    panzer::InputEquationSet ies;
    {
      ies.name = "Momentum";
      ies.basis = "Q2";
      ies.integration_order = 1;
      ies.model_id = 6;
      ies.model_factory = "rf";
      ies.prefix = "";
      ies.params.set<int>("junk", 1);
    }

    Teuchos::ParameterList default_params; 
    {
      const int num_cells = 20;
      const int base_cell_dimension = 3;
      const panzer::CellData cell_data(num_cells, base_cell_dimension);
      const int cubature_degree = 2;      
      Teuchos::RCP<panzer::IntegrationRule> ir = 
	Teuchos::rcp(new panzer::IntegrationRule(cubature_degree, cell_data));
      default_params.set("IR",ir);
    }

    Teuchos::ParameterList p;
    {
      p.set("Object Type","Constant");
      p.set("Name","UX");
      p.set("Value",4.5);
    }

    std::vector<Teuchos::ParameterList> evaluators_to_build;
    evaluators_to_build.push_back(p);
    evaluators_to_build.push_back(p);
    evaluators_to_build.push_back(p);

    TEST_EQUALITY(evaluators_to_build.size(), 3);
    

    user_app::MyModelFactory<panzer::Traits::Residual> mf;

    Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > evaluators;

    evaluators = mf.buildModels(ies, evaluators_to_build, default_params);

    TEST_EQUALITY(evaluators->size(), 3);

    // Now using template manger
    evaluators_to_build.push_back(p);  // add one to diferentiate from above

    user_app::MyModelFactory_TemplateBuilder builder;
    panzer::ModelFactory_TemplateManager<panzer::Traits> model_factory;
    model_factory.buildObjects(builder);
    evaluators = model_factory.getAsObject<panzer::Traits::Residual>()->buildModels(ies, evaluators_to_build, default_params);

    TEST_EQUALITY(evaluators->size(), 4);
  }

}
