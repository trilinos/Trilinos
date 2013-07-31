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

}
