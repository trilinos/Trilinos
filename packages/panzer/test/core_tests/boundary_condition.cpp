#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Panzer_BoundaryCondition.hpp"
#include <iostream>
#include <sstream>

namespace panzer {


  TEUCHOS_UNIT_TEST(bc, dirichlet_constant)
  {
    panzer::BCType dirichlet = BCT_Dirichlet;
    panzer::BCType neumann = BCT_Dirichlet;

    int bc_id = 0;
    int sideset_id = 4;
    int element_block_id = 7;
    std::string dof_name = "UX";
    double value = 5.0;
    panzer::BoundaryCondition bc(bc_id, dirichlet, sideset_id, 
				 element_block_id, dof_name, value);

    TEST_ASSERT(bc.bcID() == bc_id);
    TEST_ASSERT(bc.bcType() == dirichlet);
    TEST_ASSERT(bc.sidesetID() == sideset_id);
    TEST_ASSERT(bc.elementBlockID() == element_block_id);
    TEST_ASSERT(bc.equationSetName() == dof_name);
    TEST_ASSERT(bc.constantValue() == value);

    TEST_ASSERT(bc.isConstant());
    TEST_ASSERT(!bc.isStrategy());
    TEST_ASSERT(!bc.isMethodFunction());

    
    std::stringstream s;
    bc.print(s);

    std::string compare = "  BoundaryCondition ID =0\n  Type = Dirichlet\n  Side Set ID = 4\n  Element Block ID =7\n  Function Type = Constant\n  Variable Name = UX\n  Value = 5\n";

    std::cout << s.str() << std::endl;
    std::cout << compare << std::endl;

    //TEST_ASSERT(s.str() == compare);

  }

}
