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

#include "Panzer_CellTopologyInfo.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(cell_topology_info, quad_test)
  {
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));
    
    const int num_cells = 20;
    const int num_dims = 2;
    
    panzer::CellTopologyInfo cellTopoInfo(num_cells, topo);

    TEST_ASSERT(cellTopoInfo.getNumCells() == 20);
    TEST_ASSERT(cellTopoInfo.getDimension() == num_dims);
    TEST_ASSERT(cellTopoInfo.getNumEdges() == 4);
    TEST_ASSERT(cellTopoInfo.getCellName() == "Quadrilateral_4");
    
    Teuchos::RCP<const shards::CellTopology> topology = cellTopoInfo.getCellTopology();
    TEST_ASSERT(!Teuchos::is_null(topology));
    
    const int num_edges = cellTopoInfo.getNumEdges();
    
    TEST_ASSERT(static_cast<int>(cellTopoInfo.edge_scalar->size()) == num_cells*num_edges);
    TEST_ASSERT(static_cast<int>(cellTopoInfo.edge_vector->size()) == num_cells*num_edges*num_dims);
  }
  
  TEUCHOS_UNIT_TEST(cell_topology_info, tri_test)
  {
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Triangle<3> >()));
    
    const int num_cells = 20;
    const int num_dims = 2;
    
    panzer::CellTopologyInfo cellTopoInfo(num_cells, topo);

    TEST_ASSERT(cellTopoInfo.getNumCells() == 20);
    TEST_ASSERT(cellTopoInfo.getDimension() == num_dims);
    TEST_ASSERT(cellTopoInfo.getNumEdges() == 3);
    TEST_ASSERT(cellTopoInfo.getCellName() == "Triangle_3");
    
    Teuchos::RCP<const shards::CellTopology> topology = cellTopoInfo.getCellTopology();
    TEST_ASSERT(!Teuchos::is_null(topology));
    
    const int num_edges = cellTopoInfo.getNumEdges();

    TEST_ASSERT(static_cast<int>(cellTopoInfo.edge_scalar->size()) == num_cells*num_edges);
    TEST_ASSERT(static_cast<int>(cellTopoInfo.edge_vector->size()) == num_cells*num_edges*num_dims);
  }
  
  TEUCHOS_UNIT_TEST(cell_topology_info, tet_test)
  {
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<4> >()));
    
    const int num_cells = 20;
    const int num_dims = 3;
    
    panzer::CellTopologyInfo cellTopoInfo(num_cells, topo);

    TEST_ASSERT(cellTopoInfo.getNumCells() == 20);
    TEST_ASSERT(cellTopoInfo.getDimension() == num_dims);
    TEST_ASSERT(cellTopoInfo.getNumEdges() == 6);
    TEST_ASSERT(cellTopoInfo.getCellName() == "Tetrahedron_4");
    
    Teuchos::RCP<const shards::CellTopology> topology = cellTopoInfo.getCellTopology();
    TEST_ASSERT(!Teuchos::is_null(topology));
    
    const int num_edges = cellTopoInfo.getNumEdges();

    TEST_ASSERT(static_cast<int>(cellTopoInfo.edge_scalar->size()) == num_cells*num_edges);
    TEST_ASSERT(static_cast<int>(cellTopoInfo.edge_vector->size()) == num_cells*num_edges*num_dims);
  }

  TEUCHOS_UNIT_TEST(cell_topology_info, hex_test)
  {
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<8> >()));
    
    const int num_cells = 20;
    const int num_dims = 3;
    
    panzer::CellTopologyInfo cellTopoInfo(num_cells, topo);

    TEST_ASSERT(cellTopoInfo.getNumCells() == 20);
    TEST_ASSERT(cellTopoInfo.getDimension() == num_dims);
    TEST_ASSERT(cellTopoInfo.getNumEdges() == 12);
    TEST_ASSERT(cellTopoInfo.getCellName() == "Hexahedron_8");
    
    Teuchos::RCP<const shards::CellTopology> topology = cellTopoInfo.getCellTopology();
    TEST_ASSERT(!Teuchos::is_null(topology));
    
    const int num_edges = cellTopoInfo.getNumEdges();

    TEST_ASSERT(static_cast<int>(cellTopoInfo.edge_scalar->size()) == num_cells*num_edges);
    TEST_ASSERT(static_cast<int>(cellTopoInfo.edge_vector->size()) == num_cells*num_edges*num_dims);
  }

}
