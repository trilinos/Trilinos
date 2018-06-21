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

  TEUCHOS_UNIT_TEST(integration_values, coord_ordering_1d)  
  {
    typedef IntegrationValues2<double> IV;
    MDFieldArrayFactory af("",true);

    int num_faces = 2;
    int num_points_per_face = 4;
    std::vector<int> order(num_points_per_face,-1);
  
    IV::Array_CellIPDim coords = af.template buildStaticArray<double,Cell,IP,Dim>("coord",2,num_faces*num_points_per_face,1); 
    coords(0,0,0) = -1.0; coords(0,1,0) =  0.0; coords(0,2,0) =  1.0; coords(0,3,0) = 2.0;
    coords(0,4,0) =  1.0; coords(0,5,0) =  2.0; coords(0,6,0) = -1.0; coords(0,7,0) = 0.0;
    coords(1,0,0) = -1.0; coords(1,1,0) =  0.0; coords(1,2,0) =  1.0; coords(1,3,0) = 2.0;
    coords(1,4,0) =  1.0; coords(1,5,0) =  2.0; coords(1,6,0) = -1.0; coords(1,7,0) = 0.0;

    int cell = 0, offset = 0;

    cell = 0; offset = 0;
    IV::uniqueCoordOrdering(coords,cell,offset,order);
    TEST_ASSERT(order==std::vector<int>({0,1,2,3}));

    cell = 0; offset = 4;
    IV::uniqueCoordOrdering(coords,cell,offset,order);
    TEST_ASSERT(order==std::vector<int>({2,3,0,1}));

    cell = 1; offset = 0;
    IV::uniqueCoordOrdering(coords,cell,offset,order);
    TEST_ASSERT(order==std::vector<int>({0,1,2,3}));

    cell = 1; offset = 4;
    IV::uniqueCoordOrdering(coords,cell,offset,order);
    TEST_ASSERT(order==std::vector<int>({2,3,0,1}));
  }

  TEUCHOS_UNIT_TEST(integration_values, coord_ordering_2d)  
  {
    typedef IntegrationValues2<double> IV;
    MDFieldArrayFactory af("",true);

    {
      int num_faces = 2;
      int num_points_per_face = 4;
      std::vector<int> order(num_points_per_face,-1);
  
      IV::Array_CellIPDim coords = af.template buildStaticArray<double,Cell,IP,Dim>("coord",2,num_faces*num_points_per_face,2); 

      // cell 0
      coords(0,0,0) = -1.0; coords(0,1,0) =  0.0; coords(0,2,0) =  1.0; coords(0,3,0) = 2.0;
      coords(0,0,1) =  0.0; coords(0,1,1) =  0.0; coords(0,2,1) =  0.0; coords(0,3,1) = 0.0;

      coords(0,4,0) =  1.0; coords(0,5,0) =  2.0; coords(0,6,0) = -1.0; coords(0,7,0) = 0.0;
      coords(0,4,1) =  0.0; coords(0,5,1) =  0.0; coords(0,6,1) =  0.0; coords(0,7,1) = 0.0;

      // cell 1
      coords(1,0,0) =  2.0; coords(1,1,0) =  2.0; coords(1,2,0) =  2.0; coords(1,3,0) = 2.0;
      coords(1,0,0) = -1.1; coords(1,1,1) =  0.0; coords(1,2,1) =  1.0; coords(1,3,1) = 2.0;

      coords(1,4,0) =  2.0; coords(1,5,0) =  2.0; coords(1,6,0) =  2.0; coords(1,7,0) = 2.0;
      coords(1,4,1) =  1.0; coords(1,5,1) =  2.0; coords(1,6,1) = -1.0; coords(1,7,1) = 0.0;

      int cell = 0, offset = 0;

      cell = 0; offset = 0;
      IV::uniqueCoordOrdering(coords,cell,offset,order);
      TEST_ASSERT(order==std::vector<int>({0,1,2,3}));

      cell = 0; offset = 4;
      IV::uniqueCoordOrdering(coords,cell,offset,order);
      TEST_ASSERT(order==std::vector<int>({2,3,0,1}));

      cell = 1; offset = 0;
      IV::uniqueCoordOrdering(coords,cell,offset,order);
      TEST_ASSERT(order==std::vector<int>({0,1,2,3}));

      cell = 1; offset = 4;
      IV::uniqueCoordOrdering(coords,cell,offset,order);
      TEST_ASSERT(order==std::vector<int>({2,3,0,1}));
    }

    // this was the original failing case
    {
      std::vector<int> order(2,-1);
  
      IV::Array_CellIPDim coords = af.template buildStaticArray<double,Cell,IP,Dim>("coord",4,2,2); 

      coords(0,0,0) = 0.0; coords(0,0,1) = 4.4088517374119224e-01;
      coords(0,1,0) = 0.0; coords(0,1,1) = 1.1813482625880771e-01; 

      coords(1,0,0) = 2.5; coords(1,0,1) = 4.4088517374119224e-01;
      coords(1,1,0) = 2.5; coords(1,1,1) = 1.1813482625880771e-01; 

      coords(2,0,0) = 0.0; coords(2,0,1) = 1.1813482625880771e-01; 
      coords(2,1,0) = 0.0; coords(2,1,1) = 4.4088517374119224e-01;

      coords(3,0,0) = 2.5; coords(3,0,1) = 1.1813482625880771e-01; 
      coords(3,1,0) = 2.5; coords(3,1,1) = 4.4088517374119224e-01;

      int cell = 0, offset = 0;

      cell = 0; 
      IV::uniqueCoordOrdering(coords,cell,offset,order);
      TEST_ASSERT(order==std::vector<int>({1,0}));

      cell = 1; 
      IV::uniqueCoordOrdering(coords,cell,offset,order);
      TEST_ASSERT(order==std::vector<int>({1,0}));

      cell = 2;
      IV::uniqueCoordOrdering(coords,cell,offset,order);
      TEST_ASSERT(order==std::vector<int>({0,1}));

      cell = 3;
      IV::uniqueCoordOrdering(coords,cell,offset,order);
      TEST_ASSERT(order==std::vector<int>({0,1}));
    }
  }

  TEUCHOS_UNIT_TEST(integration_values, quadpt_swap)
  {    
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    const int in_num_cells = 20;
    const int base_cell_dimension = 2;
    const panzer::CellData cell_data(in_num_cells,topo);

    const int cubature_degree = 2;    
    RCP<IntegrationRule> int_rule = 
      rcp(new IntegrationRule(cubature_degree, cell_data));
    
    panzer::IntegrationValues2<double> int_values("prefix_",true);
    panzer::MDFieldArrayFactory af("prefix_",true);

    int_values.setupArrays(int_rule);

    const int num_vertices = int_rule->topology->getNodeCount();
    PHX::MDField<double,Cell,NODE,Dim> node_coordinates 
        = af.buildStaticArray<double,Cell,NODE,Dim>("nc",in_num_cells, num_vertices, base_cell_dimension);

    // Set up node coordinates.  Here we assume the following
    // ordering.  This needs to be consistent with shards topology,
    // otherwise we will get negative determinates

    // 3(0,1)---2(1.5,1.5)
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
      node_coordinates(cell,2,x) = 1.5;
      node_coordinates(cell,2,y) = 1.5;
      node_coordinates(cell,3,x) = 0.0;
      node_coordinates(cell,3,y) = 1.0;
    }

    int_values.evaluateValues(node_coordinates);

    // pull out arrays to look at
    
    auto weighted_measure   = int_values.weighted_measure;
    auto jac_det            = int_values.jac_det;
    auto ref_ip_coordinates = int_values.ref_ip_coordinates;
    auto ip_coordinates     = int_values.ip_coordinates;
    auto jac                = int_values.jac;
    auto jac_inv            = int_values.jac_inv;

    int num_cells  = ip_coordinates.extent(0); 
    int num_points = ip_coordinates.extent(1); 
    int num_dim    = ip_coordinates.extent(2); 

    TEST_EQUALITY(num_cells,in_num_cells);
    TEST_EQUALITY(num_points,4);
    TEST_EQUALITY(num_dim,2);

    int ref_cell = 0;

    {
      int tst_cell = 1;
      int org_pt = 0;
      int new_pt = 1;

      int_values.swapQuadraturePoints(tst_cell,org_pt,new_pt);
   
      TEST_EQUALITY(    weighted_measure(ref_cell,org_pt),     weighted_measure(tst_cell,new_pt));
      TEST_EQUALITY(             jac_det(ref_cell,org_pt),              jac_det(tst_cell,new_pt));
      TEST_EQUALITY(ref_ip_coordinates(ref_cell,org_pt,0), ref_ip_coordinates(tst_cell,new_pt,0));
      TEST_EQUALITY(ref_ip_coordinates(ref_cell,org_pt,1), ref_ip_coordinates(tst_cell,new_pt,1));
      TEST_EQUALITY(    ip_coordinates(ref_cell,org_pt,0),     ip_coordinates(tst_cell,new_pt,0));
      TEST_EQUALITY(    ip_coordinates(ref_cell,org_pt,1),     ip_coordinates(tst_cell,new_pt,1));
      TEST_EQUALITY(    jac(ref_cell,org_pt,0,0),     jac(tst_cell,new_pt,0,0));
      TEST_EQUALITY(    jac(ref_cell,org_pt,0,1),     jac(tst_cell,new_pt,0,1));
      TEST_EQUALITY(    jac(ref_cell,org_pt,1,0),     jac(tst_cell,new_pt,1,0));
      TEST_EQUALITY(    jac(ref_cell,org_pt,1,1),     jac(tst_cell,new_pt,1,1));
      TEST_EQUALITY(    jac_inv(ref_cell,org_pt,0,0),     jac_inv(tst_cell,new_pt,0,0));
      TEST_EQUALITY(    jac_inv(ref_cell,org_pt,0,1),     jac_inv(tst_cell,new_pt,0,1));
      TEST_EQUALITY(    jac_inv(ref_cell,org_pt,1,0),     jac_inv(tst_cell,new_pt,1,0));
      TEST_EQUALITY(    jac_inv(ref_cell,org_pt,1,1),     jac_inv(tst_cell,new_pt,1,1));
    }

    {
      int tst_cell = 7;
      int org_pt = 1;
      int new_pt = 3;

      int_values.swapQuadraturePoints(tst_cell,org_pt,new_pt);
   
      TEST_EQUALITY(    weighted_measure(ref_cell,org_pt),     weighted_measure(tst_cell,new_pt));
      TEST_EQUALITY(             jac_det(ref_cell,org_pt),              jac_det(tst_cell,new_pt));
      TEST_EQUALITY(ref_ip_coordinates(ref_cell,org_pt,0), ref_ip_coordinates(tst_cell,new_pt,0));
      TEST_EQUALITY(ref_ip_coordinates(ref_cell,org_pt,1), ref_ip_coordinates(tst_cell,new_pt,1));
      TEST_EQUALITY(    ip_coordinates(ref_cell,org_pt,0),     ip_coordinates(tst_cell,new_pt,0));
      TEST_EQUALITY(    ip_coordinates(ref_cell,org_pt,1),     ip_coordinates(tst_cell,new_pt,1));
      TEST_EQUALITY(    jac(ref_cell,org_pt,0,0),     jac(tst_cell,new_pt,0,0));
      TEST_EQUALITY(    jac(ref_cell,org_pt,0,1),     jac(tst_cell,new_pt,0,1));
      TEST_EQUALITY(    jac(ref_cell,org_pt,1,0),     jac(tst_cell,new_pt,1,0));
      TEST_EQUALITY(    jac(ref_cell,org_pt,1,1),     jac(tst_cell,new_pt,1,1));
      TEST_EQUALITY(    jac_inv(ref_cell,org_pt,0,0),     jac_inv(tst_cell,new_pt,0,0));
      TEST_EQUALITY(    jac_inv(ref_cell,org_pt,0,1),     jac_inv(tst_cell,new_pt,0,1));
      TEST_EQUALITY(    jac_inv(ref_cell,org_pt,1,0),     jac_inv(tst_cell,new_pt,1,0));
      TEST_EQUALITY(    jac_inv(ref_cell,org_pt,1,1),     jac_inv(tst_cell,new_pt,1,1));
    }
  }
}
