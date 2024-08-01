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
#include "Panzer_IntegrationRule.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(integration_rule, volume)
  {
    
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<8> >()));
    
    const int num_cells = 20;
    const int base_cell_dimension = 3;
    const panzer::CellData cell_data(num_cells,topo);
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
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<8> >()));

    const int num_cells = 20;
    const int base_cell_dimension = 3;
    const int cell_local_side_id = 1;
    const panzer::CellData cell_data(num_cells,cell_local_side_id,topo);
    const int cubature_degree = 2;
    
    panzer::IntegrationRule int_rule(cubature_degree, cell_data);
    
    TEST_ASSERT(cubature_degree == int_rule.cubature_degree);
    TEST_ASSERT(num_cells == int_rule.workset_size);
    TEST_ASSERT(int_rule.num_points == 4);
    TEST_ASSERT(base_cell_dimension == int_rule.spatial_dimension);
    TEST_ASSERT(int_rule.isSide());
    TEST_ASSERT(int_rule.side == 1);
  }

  TEUCHOS_UNIT_TEST(integration_rule, cv_vol)
  {
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    const int num_cells = 20;
    const int base_cell_dimension = 2;
    const panzer::CellData cell_data(num_cells,topo);
    std::string cv_type = "volume";
    
    panzer::IntegrationRule int_rule(cell_data,cv_type);
    
    TEST_ASSERT(num_cells == int_rule.workset_size);
    TEST_ASSERT(int_rule.num_points == 4);
    TEST_ASSERT(base_cell_dimension == int_rule.spatial_dimension);
    TEST_ASSERT(cv_type == int_rule.cv_type);
  }

  TEUCHOS_UNIT_TEST(integration_rule, cv_side)
  {
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    const int num_cells = 20;
    const int base_cell_dimension = 2;
    const int cell_local_side_id = 1;
    const panzer::CellData cell_data(num_cells,cell_local_side_id,topo);
    std::string cv_type = "side";
    
    panzer::IntegrationRule int_rule(cell_data,cv_type);
    
    TEST_ASSERT(num_cells == int_rule.workset_size);
    TEST_ASSERT(int_rule.num_points == 4);
    TEST_ASSERT(base_cell_dimension == int_rule.spatial_dimension);
    TEST_ASSERT(cv_type == int_rule.cv_type);
  }

  TEUCHOS_UNIT_TEST(integration_rule, cv_bc)
  {
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    const int num_cells = 20;
    const int base_cell_dimension = 2;
    const int cell_local_side_id = 1;
    const panzer::CellData cell_data(num_cells,cell_local_side_id,topo);
    std::string cv_type = "boundary";
    
    panzer::IntegrationRule int_rule(cell_data,cv_type);
    
    TEST_ASSERT(num_cells == int_rule.workset_size);
    TEST_ASSERT(int_rule.num_points == 2);
    TEST_ASSERT(base_cell_dimension == int_rule.spatial_dimension);
    TEST_ASSERT(cv_type == int_rule.cv_type);
  }

}
