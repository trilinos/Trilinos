#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_ParameterList.hpp>

#include "Panzer_Traits.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_BCStrategy.hpp"
#include "user_app_BCStrategy_Factory.hpp"
#include <iostream>

namespace panzer {

  TEUCHOS_UNIT_TEST(bcstrategy, basic_construction)
  {

    std::size_t bc_id = 0;
    panzer::BCType neumann = BCT_Dirichlet;
    std::string sideset_id = "4";
    std::string element_block_id = "fluid";
    std::string dof_name = "UX";
    std::string strategy = "Constant";
    double value = 5.0;
    Teuchos::ParameterList p;
    p.set("Value",value);
    panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		  strategy, p);

    Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> > bcs;
  
    user_app::MyFactory my_factory;
    bcs = my_factory.buildBCStrategy(bc);

  }

}
