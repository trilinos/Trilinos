#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Panzer_CellData.hpp"
#include "Panzer_IntegrationRule.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(integration_rule, volume)
  {
    const int num_cells = 20;
    const int base_cell_dimension = 3;
    const panzer::CellData cell_data(num_cells, base_cell_dimension);
    const int cubature_degree = 2;
    
    panzer::IntegrationRule int_rule(cubature_degree, cell_data);
    
    TEST_ASSERT(cubature_degree == int_rule.cubature_degree);
    TEST_ASSERT(num_cells == int_rule.workset_size);
    TEST_ASSERT(int_rule.num_points == 8);
    TEST_ASSERT(base_cell_dimension == int_rule.spatial_dimension);
    TEST_ASSERT(!int_rule.isSide());
    TEST_ASSERT(int_rule.side == -1);
  }

  TEUCHOS_UNIT_TEST(integration_rule, side)
  {
    const int num_cells = 20;
    const int base_cell_dimension = 3;
    const int cell_local_side_id = 1;
    const panzer::CellData cell_data(num_cells, base_cell_dimension,
				     cell_local_side_id);
    const int cubature_degree = 2;
    
    panzer::IntegrationRule int_rule(cubature_degree, cell_data);
    
    TEST_ASSERT(cubature_degree == int_rule.cubature_degree);
    TEST_ASSERT(num_cells == int_rule.workset_size);
    TEST_ASSERT(int_rule.num_points == 4);
    TEST_ASSERT(base_cell_dimension == int_rule.spatial_dimension);
    TEST_ASSERT(int_rule.isSide());
    TEST_ASSERT(int_rule.side == 1);
  }

}
