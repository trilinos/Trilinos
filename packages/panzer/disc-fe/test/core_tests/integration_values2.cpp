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
#include "Panzer_IntegrationValues2.hpp"
#include "Panzer_ArrayTraits.hpp"
#include "Panzer_CommonArrayFactories.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using panzer::IntegrationRule;

namespace panzer {

  TEUCHOS_UNIT_TEST(integration_values, volume)
  {    
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    const int num_cells = 20;
    const int base_cell_dimension = 2;
    const panzer::CellData cell_data(num_cells,topo);

    const int cubature_degree = 2;    
    RCP<IntegrationRule> int_rule = 
      rcp(new IntegrationRule(cubature_degree, cell_data));
    
    panzer::IntegrationValues2<double> int_values("prefix_",true);
    panzer::MDFieldArrayFactory af("prefix_",true);

    int_values.setupArrays(int_rule);

    const int num_vertices = int_rule->topology->getNodeCount();
    PHX::MDField<double,Cell,NODE,Dim> node_coordinates 
        = af.buildStaticArray<double,Cell,NODE,Dim>("nc",num_cells, num_vertices, base_cell_dimension);

    // Set up node coordinates.  Here we assume the following
    // ordering.  This needs to be consistent with shards topology,
    // otherwise we will get negative determinates

    // 3(0,1)---2(1,1)
    //   |    0  |
    //   |       |
    // 0(0,0)---1(1,0)

    typedef panzer::ArrayTraits<double,PHX::MDField<double> >::size_type size_type;
    const size_type x = 0;
    const size_type y = 1;
    for (size_type cell = 0; cell < node_coordinates.extent(0); ++cell) {
      node_coordinates(cell,0,x) = 0.0;
      node_coordinates(cell,0,y) = 0.0;
      node_coordinates(cell,1,x) = 1.0;
      node_coordinates(cell,1,y) = 0.0;
      node_coordinates(cell,2,x) = 1.0;
      node_coordinates(cell,2,y) = 1.0;
      node_coordinates(cell,3,x) = 0.0;
      node_coordinates(cell,3,y) = 1.0;
    }

    int_values.evaluateValues(node_coordinates);
    
    TEST_EQUALITY(int_values.ip_coordinates.extent(1), 4);
    double realspace_x_coord = (1.0/std::sqrt(3.0) + 1.0) / 2.0;
    double realspace_y_coord = (1.0/std::sqrt(3.0) + 1.0) / 2.0;
    TEST_FLOATING_EQUALITY(int_values.ip_coordinates(0,0,0), 
                           realspace_x_coord, 1.0e-8);
    TEST_FLOATING_EQUALITY(int_values.ip_coordinates(0,0,1), 
                           realspace_y_coord, 1.0e-8);

  }

  TEUCHOS_UNIT_TEST(integration_values, control_volume)
  {    
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    const int num_cells = 20;
    const int base_cell_dimension = 2;
    const panzer::CellData cell_data(num_cells,topo);

    std::string cv_type = "volume";
    RCP<IntegrationRule> int_rule_vol = 
      rcp(new IntegrationRule(cell_data, cv_type));
    
    panzer::IntegrationValues2<double> int_values_vol("prefix_",true);
    panzer::MDFieldArrayFactory af("prefix_",true);

    int_values_vol.setupArrays(int_rule_vol);

    cv_type = "side";
    RCP<IntegrationRule> int_rule_side = 
      rcp(new IntegrationRule(cell_data, cv_type));
    
    panzer::IntegrationValues2<double> int_values_side("prefix_",true);

    int_values_side.setupArrays(int_rule_side);

    const int num_vertices = int_rule_vol->topology->getNodeCount();
    PHX::MDField<double,Cell,NODE,Dim> node_coordinates 
        = af.buildStaticArray<double,Cell,NODE,Dim>("nc",num_cells, num_vertices, base_cell_dimension);

    // Set up node coordinates.  Here we assume the following
    // ordering.  This needs to be consistent with shards topology,
    // otherwise we will get negative determinates

    // 3(0,1)---2(1,1)
    //   |    0  |
    //   |       |
    // 0(0,0)---1(1,0)

    typedef panzer::ArrayTraits<double,PHX::MDField<double> >::size_type size_type;
    const size_type x = 0;
    const size_type y = 1;
    for (size_type cell = 0; cell < node_coordinates.extent(0); ++cell) {
      node_coordinates(cell,0,x) = 0.0;
      node_coordinates(cell,0,y) = 0.0;
      node_coordinates(cell,1,x) = 1.0;
      node_coordinates(cell,1,y) = 0.0;
      node_coordinates(cell,2,x) = 1.0;
      node_coordinates(cell,2,y) = 1.0;
      node_coordinates(cell,3,x) = 0.0;
      node_coordinates(cell,3,y) = 1.0;
    }

    int_values_vol.evaluateValues(node_coordinates);
    int_values_side.evaluateValues(node_coordinates);
    
    TEST_EQUALITY(int_values_vol.ip_coordinates.extent(1), 4);
    TEST_EQUALITY(int_values_side.ip_coordinates.extent(1), 4);
    TEST_EQUALITY(int_values_side.weighted_normals.extent(1), 4);
    double realspace_x_coord_1 = 0.25;
    double realspace_y_coord_1 = 0.25;
    TEST_FLOATING_EQUALITY(int_values_vol.ip_coordinates(0,0,0), 
                           realspace_x_coord_1, 1.0e-8);
    TEST_FLOATING_EQUALITY(int_values_vol.ip_coordinates(0,0,1), 
                           realspace_y_coord_1, 1.0e-8);
    double realspace_x_coord_2 = 0.5;
    double realspace_y_coord_2 = 0.25;
    TEST_FLOATING_EQUALITY(int_values_side.ip_coordinates(0,0,0), 
                           realspace_x_coord_2, 1.0e-8);
    TEST_FLOATING_EQUALITY(int_values_side.ip_coordinates(0,0,1), 
                           realspace_y_coord_2, 1.0e-8);

  }

  TEUCHOS_UNIT_TEST(integration_values, control_volume_boundary)
  {
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    const int num_cells = 2;
    const int base_cell_dimension = 2;
    const int cell_side = 1;
    const panzer::CellData cell_data(num_cells,cell_side,topo);

    std::string cv_type = "boundary";
    RCP<IntegrationRule> int_rule_bc = 
      rcp(new IntegrationRule(cell_data, cv_type));

    panzer::IntegrationValues2<double> int_values_bc("prefix_",true);
    panzer::MDFieldArrayFactory af("prefix_",true);

    int_values_bc.setupArrays(int_rule_bc);

    const int num_vertices = int_rule_bc->topology->getNodeCount();
    PHX::MDField<double,Cell,NODE,Dim> node_coordinates 
        = af.buildStaticArray<double,Cell,NODE,Dim>("nc",num_cells, num_vertices, base_cell_dimension);

    // Set up node coordinates.  Here we assume the following
    // ordering.  This needs to be consistent with shards topology,
    // otherwise we will get negative determinates

    // 3(0,1)---2(1,1)
    //   |    0  |
    //   |       |
    // 0(0,0)---1(1,0)

    typedef panzer::ArrayTraits<double,PHX::MDField<double> >::size_type size_type;
    const size_type x = 0;
    const size_type y = 1;
    for (size_type cell = 0; cell < node_coordinates.extent(0); ++cell) {
      node_coordinates(cell,0,x) = 0.0;
      node_coordinates(cell,0,y) = 0.0;
      node_coordinates(cell,1,x) = 1.0;
      node_coordinates(cell,1,y) = 0.0;
      node_coordinates(cell,2,x) = 1.0;
      node_coordinates(cell,2,y) = 1.0;
      node_coordinates(cell,3,x) = 0.0;
      node_coordinates(cell,3,y) = 1.0;
    }

    int_values_bc.evaluateValues(node_coordinates);
    
    TEST_EQUALITY(int_values_bc.ip_coordinates.extent(1), 2);
    double realspace_x_coord_1 = 1.0;
    double realspace_y_coord_1 = 0.25;
    TEST_FLOATING_EQUALITY(int_values_bc.ip_coordinates(0,0,0), 
                           realspace_x_coord_1, 1.0e-8);
    TEST_FLOATING_EQUALITY(int_values_bc.ip_coordinates(0,0,1), 
                           realspace_y_coord_1, 1.0e-8);

  }
}
