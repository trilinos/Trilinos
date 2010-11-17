#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_ParameterList.hpp>

#include "Panzer_BC.hpp"
#include <iostream>
#include <sstream>
#include <map>

namespace panzer {

  TEUCHOS_UNIT_TEST(bc, neumann_no_param_list)
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

    TEST_EQUALITY(bc.bcID(), bc_id);
    TEST_EQUALITY(bc.bcType(), neumann);
    TEST_EQUALITY(bc.sidesetID(), sideset_id);
    TEST_EQUALITY(bc.elementBlockID(), element_block_id);
    TEST_EQUALITY(bc.equationSetName(), dof_name);

    std::stringstream s;
    s << bc << std::endl;
  }

  TEUCHOS_UNIT_TEST(bc, dirichlet_with_param_list)
  {
    std::size_t bc_id = 0;
    panzer::BCType dirichlet = BCT_Dirichlet;
    std::string sideset_id = "4";
    std::string element_block_id = "fluid";
    std::string dof_name = "UX";
    std::string strategy = "Constant";
    double value = 5.0;
    Teuchos::ParameterList p;
    p.set("Value",value);
    panzer::BC bc(bc_id, dirichlet, sideset_id, element_block_id, dof_name, 
		  strategy, p);

    TEST_EQUALITY(bc.bcID(), bc_id);
    TEST_EQUALITY(bc.bcType(), dirichlet);
    TEST_EQUALITY(bc.sidesetID(), sideset_id);
    TEST_EQUALITY(bc.elementBlockID(), element_block_id);
    TEST_EQUALITY(bc.equationSetName(), dof_name);

    std::stringstream s;
    s << bc << std::endl;
  }

  TEUCHOS_UNIT_TEST(bc, map_comparitor)
  {
    using panzer::BC;

    BC bc1(0,BCT_Dirichlet,"3","fluid","VELOCITY","Constant");
    BC bc2(1,BCT_Dirichlet,"3","fluid","VELOCITY","Constant");
    BC bc3(2,BCT_Dirichlet,"3","fluid","VELOCITY","Constant");

    std::map<BC,int,panzer::LessBC> my_bcs;

    my_bcs[bc1] = 8;
    my_bcs[bc2] = 2;
    my_bcs[bc3] = 11;

    TEST_EQUALITY(my_bcs[bc1], 8);
    TEST_EQUALITY(my_bcs[bc2], 2);
    TEST_EQUALITY(my_bcs[bc3], 11);

    TEST_INEQUALITY(my_bcs[bc3], 4);    
  }
}
