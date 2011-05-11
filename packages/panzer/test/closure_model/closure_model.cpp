#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_ParameterList.hpp>

#include "Panzer_CellData.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"
#include "Phalanx_FieldManager.hpp"
#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_ClosureModel_Factory.hpp"
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
      ies.model_id = "fluid model";
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
      Teuchos::RCP<panzer::Basis> basis = 
	Teuchos::rcp(new panzer::Basis("Q1", *ir));
      default_params.set("Basis",basis);
      
    }

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

    PHX::FieldManager<panzer::Traits> fm;

    evaluators = mf.buildClosureModels(ies.model_id, ies, p, default_params, user_data, fm);

    TEST_EQUALITY(evaluators->size(), 8);

    user_app::MyModelFactory_TemplateBuilder builder;
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> model_factory;
    model_factory.buildObjects(builder);
    evaluators = model_factory.getAsObject<panzer::Traits::Residual>()->buildClosureModels(ies.model_id, ies, p, default_params, user_data, fm);

    TEST_EQUALITY(evaluators->size(), 8);

    // Add an unsupported type
    p.sublist("fluid model").sublist("garbage").set<std::string>("Value","help!");
    
    TEST_THROW(model_factory.getAsObject<panzer::Traits::Residual>()->buildClosureModels(ies.model_id, ies, p, default_params, user_data, fm), std::logic_error);

  }

}
