// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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
    
    TEST_ASSERT(cellTopoInfo.edge_scalar->size() == num_cells*num_edges);
    TEST_ASSERT(cellTopoInfo.edge_vector->size() == num_cells*num_edges*num_dims);
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

    TEST_ASSERT(cellTopoInfo.edge_scalar->size() == num_cells*num_edges);
    TEST_ASSERT(cellTopoInfo.edge_vector->size() == num_cells*num_edges*num_dims);
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

    TEST_ASSERT(cellTopoInfo.edge_scalar->size() == num_cells*num_edges);
    TEST_ASSERT(cellTopoInfo.edge_vector->size() == num_cells*num_edges*num_dims);
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

    TEST_ASSERT(cellTopoInfo.edge_scalar->size() == num_cells*num_edges);
    TEST_ASSERT(cellTopoInfo.edge_vector->size() == num_cells*num_edges*num_dims);
  }

}
