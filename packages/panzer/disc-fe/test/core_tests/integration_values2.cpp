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
#include "Panzer_IntegrationValues2_impl.hpp" // for unit testing wihtout CUDA RDC
#include "Panzer_ArrayTraits.hpp"
#include "Panzer_CommonArrayFactories.hpp"
#include "Panzer_SubcellConnectivity.hpp"
#include "Panzer_LocalMeshInfo.hpp"

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
    auto node_coordinates_k = node_coordinates.get_view();
    Kokkos::parallel_for("initialize node coords",node_coordinates.extent(0),
                         KOKKOS_LAMBDA (const int cell) {
      node_coordinates_k(cell,0,x) = 0.0;
      node_coordinates_k(cell,0,y) = 0.0;
      node_coordinates_k(cell,1,x) = 1.0;
      node_coordinates_k(cell,1,y) = 0.0;
      node_coordinates_k(cell,2,x) = 1.0;
      node_coordinates_k(cell,2,y) = 1.0;
      node_coordinates_k(cell,3,x) = 0.0;
      node_coordinates_k(cell,3,y) = 1.0;
    });

    int_values.evaluateValues(node_coordinates);
    
    TEST_EQUALITY(int_values.ip_coordinates.extent(1), 4);
    double realspace_x_coord = (1.0/std::sqrt(3.0) + 1.0) / 2.0;
    double realspace_y_coord = (1.0/std::sqrt(3.0) + 1.0) / 2.0;
    auto tmp_coords = int_values.ip_coordinates;
    auto tmp_coords_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),tmp_coords.get_view());
    TEST_FLOATING_EQUALITY(tmp_coords_host(0,0,0), 
                           realspace_x_coord, 1.0e-8);
    TEST_FLOATING_EQUALITY(tmp_coords_host(0,0,1), 
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
    auto node_coordinates_k = node_coordinates.get_view();
    Kokkos::parallel_for("initialize node coords",node_coordinates.extent(0),
                         KOKKOS_LAMBDA (const int cell) {
      node_coordinates_k(cell,0,x) = 0.0;
      node_coordinates_k(cell,0,y) = 0.0;
      node_coordinates_k(cell,1,x) = 1.0;
      node_coordinates_k(cell,1,y) = 0.0;
      node_coordinates_k(cell,2,x) = 1.0;
      node_coordinates_k(cell,2,y) = 1.0;
      node_coordinates_k(cell,3,x) = 0.0;
      node_coordinates_k(cell,3,y) = 1.0;
    });

    int_values_vol.evaluateValues(node_coordinates);
    int_values_side.evaluateValues(node_coordinates);
    
    TEST_EQUALITY(int_values_vol.ip_coordinates.extent(1), 4);
    TEST_EQUALITY(int_values_side.ip_coordinates.extent(1), 4);
    TEST_EQUALITY(int_values_side.weighted_normals.extent(1), 4);
    double realspace_x_coord_1 = 0.25;
    double realspace_y_coord_1 = 0.25;
    auto tmp_coords = int_values_vol.ip_coordinates;
    auto tmp_coords_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),tmp_coords.get_view());
    TEST_FLOATING_EQUALITY(tmp_coords_host(0,0,0), 
                           realspace_x_coord_1, 1.0e-8);
    TEST_FLOATING_EQUALITY(tmp_coords_host(0,0,1), 
                           realspace_y_coord_1, 1.0e-8);
    double realspace_x_coord_2 = 0.5;
    double realspace_y_coord_2 = 0.25;
    auto tmp_side_coords = int_values_side.ip_coordinates;
    auto tmp_side_coords_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),tmp_side_coords.get_view());
    TEST_FLOATING_EQUALITY(tmp_side_coords_host(0,0,0), 
                           realspace_x_coord_2, 1.0e-8);
    TEST_FLOATING_EQUALITY(tmp_side_coords_host(0,0,1), 
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
    auto node_coordinates_k = node_coordinates.get_view();
    Kokkos::parallel_for("initialize node coords",node_coordinates.extent(0),
                         KOKKOS_LAMBDA (const int cell) {
      node_coordinates_k(cell,0,x) = 0.0;
      node_coordinates_k(cell,0,y) = 0.0;
      node_coordinates_k(cell,1,x) = 1.0;
      node_coordinates_k(cell,1,y) = 0.0;
      node_coordinates_k(cell,2,x) = 1.0;
      node_coordinates_k(cell,2,y) = 1.0;
      node_coordinates_k(cell,3,x) = 0.0;
      node_coordinates_k(cell,3,y) = 1.0;
    });

    int_values_bc.evaluateValues(node_coordinates);

    auto ip_coords_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),int_values_bc.ip_coordinates.get_view());
    
    TEST_EQUALITY(int_values_bc.ip_coordinates.extent(1), 2);
    double realspace_x_coord_1 = 1.0;
    double realspace_y_coord_1 = 0.25;
    TEST_FLOATING_EQUALITY(ip_coords_host(0,0,0), 
                           realspace_x_coord_1, 1.0e-8);
    TEST_FLOATING_EQUALITY(ip_coords_host(0,0,1), 
                           realspace_y_coord_1, 1.0e-8);

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
    auto node_coordinates_k = node_coordinates.get_view();
    Kokkos::parallel_for("initialize node coords",node_coordinates.extent(0),
                         KOKKOS_LAMBDA (const int cell) {
      node_coordinates_k(cell,0,x) = 0.0;
      node_coordinates_k(cell,0,y) = 0.0;
      node_coordinates_k(cell,1,x) = 1.0;
      node_coordinates_k(cell,1,y) = 0.0;
      node_coordinates_k(cell,2,x) = 1.5;
      node_coordinates_k(cell,2,y) = 1.5;
      node_coordinates_k(cell,3,x) = 0.0;
      node_coordinates_k(cell,3,y) = 1.0;
    });

    int_values.evaluateValues(node_coordinates);

    // pull out arrays to look at
    
    auto weighted_measure   = Kokkos::create_mirror_view(int_values.weighted_measure.get_view());
    auto jac_det            = Kokkos::create_mirror_view(int_values.jac_det.get_view());
    auto ref_ip_coordinates = Kokkos::create_mirror_view(int_values.ref_ip_coordinates.get_view());
    auto ip_coordinates     = Kokkos::create_mirror_view(int_values.ip_coordinates.get_view());
    auto jac                = Kokkos::create_mirror_view(int_values.jac.get_view());
    auto jac_inv            = Kokkos::create_mirror_view(int_values.jac_inv.get_view());

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

      auto ref_ip_coordinates_k = int_values.ref_ip_coordinates.get_view();
      auto ip_coordinates_k = int_values.ip_coordinates.get_view();
      auto weighted_measure_k = int_values.weighted_measure.get_view();
      auto jac_k = int_values.jac.get_view();
      auto jac_det_k = int_values.jac_det.get_view();
      auto jac_inv_k = int_values.jac_inv.get_view();
      auto surface_normals_k = int_values.surface_normals.get_view();
      auto surface_rotation_matrices_k = int_values.surface_rotation_matrices.get_view();
      Kokkos::parallel_for("swap qp",1,KOKKOS_LAMBDA (const int ) {
        panzer::swapQuadraturePoints<double>(tst_cell,org_pt,new_pt,
                                             ref_ip_coordinates_k,
                                             ip_coordinates_k,
                                             weighted_measure_k,
                                             jac_k,
                                             jac_det_k,
                                             jac_inv_k,
                                             surface_normals_k,
                                             surface_rotation_matrices_k);
      });

      Kokkos::deep_copy(weighted_measure,int_values.weighted_measure.get_view());
      Kokkos::deep_copy(jac_det,int_values.jac_det.get_view());
      Kokkos::deep_copy(ref_ip_coordinates,int_values.ref_ip_coordinates.get_view());
      Kokkos::deep_copy(ip_coordinates,int_values.ip_coordinates.get_view());
      Kokkos::deep_copy(jac,int_values.jac.get_view());
      Kokkos::deep_copy(jac_inv,int_values.jac_inv.get_view());
      
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

      auto ref_ip_coordinates_k = int_values.ref_ip_coordinates.get_view();
      auto ip_coordinates_k = int_values.ip_coordinates.get_view();
      auto weighted_measure_k = int_values.weighted_measure.get_view();
      auto jac_k = int_values.jac.get_view();
      auto jac_det_k = int_values.jac_det.get_view();
      auto jac_inv_k = int_values.jac_inv.get_view();
      auto surface_normals_k = int_values.surface_normals.get_view();
      auto surface_rotation_matrices_k = int_values.surface_rotation_matrices.get_view();
      // Kokkos::parallel_for("swap qp",1,DoSwapOnDevice<double>(tst_cell,org_pt,new_pt,int_values));
      Kokkos::parallel_for("swap qp",1,KOKKOS_LAMBDA (const int ) {
        panzer::swapQuadraturePoints<double>(tst_cell,org_pt,new_pt,
                                             ref_ip_coordinates_k,
                                             ip_coordinates_k,
                                             weighted_measure_k,
                                             jac_k,
                                             jac_det_k,
                                             jac_inv_k,
                                             surface_normals_k,
                                             surface_rotation_matrices_k);
      });

      Kokkos::deep_copy(weighted_measure,int_values.weighted_measure.get_view());
      Kokkos::deep_copy(jac_det,int_values.jac_det.get_view());
      Kokkos::deep_copy(ref_ip_coordinates,int_values.ref_ip_coordinates.get_view());
      Kokkos::deep_copy(ip_coordinates,int_values.ip_coordinates.get_view());
      Kokkos::deep_copy(jac,int_values.jac.get_view());
      Kokkos::deep_copy(jac_inv,int_values.jac_inv.get_view());

   
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

  TEUCHOS_UNIT_TEST(integration_values, surface_quadrature){

    /*

      We are going to apply a periodic boundary condition set at an angle to make sure we can align points across the boundary.
      The mesh takes the form:

      0-----1-----2
       \    |    /
        \   |   /
         3--4--5
 
      Cell 0: 0,3,4,1
      Cell 1: 1,4,5,2

      Face connectivity will be:
      face ID, Cell 0, Cell 1, (Node 0, Node 1)
      0:  0,  1 (0,3) = (5,2)
      1:  0, -1 (3,4)
      2:  0,  1 (4,1) = (1,4)
      3:  0, -1 (1,0)
      4:  1, -1 (4,5)
      5:  1, -1 (2,1)
       
     */

    // First we build the mesh
    panzer::LocalMeshPartition mesh;
    {
      mesh.num_owned_cells = 2;
      mesh.num_ghstd_cells = 0;
      mesh.num_virtual_cells = 0;
      mesh.cell_topology = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));
      mesh.local_cells = PHX::View<panzer::LocalOrdinal*>("local_cells",2);
      auto local_cells_host = Kokkos::create_mirror_view(mesh.local_cells);
      local_cells_host(0) = 0;
      local_cells_host(1) = 1;
      Kokkos::deep_copy(mesh.local_cells,local_cells_host);
      mesh.cell_vertices = PHX::View<double***>("cell_vertices",2,4,2);

      PHX::View<double**> coordinates("coordinates",6,2);
      auto coordinates_host = Kokkos::create_mirror_view(coordinates);
      coordinates_host(0,0) = 1.; coordinates_host(0,1) = 3.;
      coordinates_host(1,0) = 3.; coordinates_host(1,1) = 3.;
      coordinates_host(2,0) = 5.; coordinates_host(2,1) = 3.;
      coordinates_host(3,0) = 2.; coordinates_host(3,1) = 1.;
      coordinates_host(4,0) = 3.; coordinates_host(4,1) = 1.;
      coordinates_host(5,0) = 4.; coordinates_host(5,1) = 1.;
      Kokkos::deep_copy(coordinates,coordinates_host);

      auto cell_vertices_host = Kokkos::create_mirror_view(mesh.cell_vertices);
#define SET_NC(cell,node,vertex) cell_vertices_host(cell,node,0) = coordinates_host(vertex,0); cell_vertices_host(cell,node,1) = coordinates_host(vertex,1);
      SET_NC(0,0, 0);
      SET_NC(0,1, 3);
      SET_NC(0,2, 4);
      SET_NC(0,3, 1);
      SET_NC(1,0, 1);
      SET_NC(1,1, 4);
      SET_NC(1,2, 5);
      SET_NC(1,3, 2);
#undef SET_NC
      Kokkos::deep_copy(mesh.cell_vertices,cell_vertices_host);

      mesh.face_to_cells = PHX::View<panzer::LocalOrdinal*[2]>("face_to_cells",6);
      auto face_to_cells_host = Kokkos::create_mirror_view(mesh.face_to_cells);
      face_to_cells_host(0,0) = 0; face_to_cells_host(0,1) =  1;
      face_to_cells_host(1,0) = 0; face_to_cells_host(1,1) = -1;
      face_to_cells_host(2,0) = 0; face_to_cells_host(2,1) =  1;
      face_to_cells_host(3,0) = 0; face_to_cells_host(3,1) = -1;
      face_to_cells_host(4,0) = 1; face_to_cells_host(4,1) = -1;
      face_to_cells_host(5,0) = 1; face_to_cells_host(5,1) = -1;
      Kokkos::deep_copy(mesh.face_to_cells,face_to_cells_host);

      mesh.face_to_lidx = PHX::View<panzer::LocalOrdinal*[2]>("face_to_lidx",6);
      auto face_to_lidx_host = Kokkos::create_mirror_view(mesh.face_to_lidx);
      face_to_lidx_host(0,0) = 0; face_to_lidx_host(0,1) =  2;
      face_to_lidx_host(1,0) = 1; face_to_lidx_host(1,1) = -1;
      face_to_lidx_host(2,0) = 2; face_to_lidx_host(2,1) =  0;
      face_to_lidx_host(3,0) = 3; face_to_lidx_host(3,1) = -1;
      face_to_lidx_host(4,0) = 1; face_to_lidx_host(4,1) = -1;
      face_to_lidx_host(5,0) = 3; face_to_lidx_host(5,1) = -1;
      Kokkos::deep_copy(mesh.face_to_lidx,face_to_lidx_host);

      mesh.cell_to_faces = PHX::View<panzer::LocalOrdinal**>("cell_to_faces",2,4);
      auto cell_to_faces_host = Kokkos::create_mirror_view(mesh.cell_to_faces);
      cell_to_faces_host(0,0) = 0; cell_to_faces_host(0,1) = 1; cell_to_faces_host(0,2) = 2; cell_to_faces_host(0,3) = 3;
      cell_to_faces_host(1,0) = 2; cell_to_faces_host(1,1) = 4; cell_to_faces_host(1,2) = 0; cell_to_faces_host(1,3) = 5;
      Kokkos::deep_copy(mesh.cell_to_faces,cell_to_faces_host);
    }
   
    const auto id = panzer::IntegrationDescriptor(2, panzer::IntegrationDescriptor::SURFACE);
    auto int_rule = Teuchos::rcp(new IntegrationRule(id, mesh.cell_topology, 2, 6));
    panzer::IntegrationValues2<double> int_values("prefix_",true);
    int_values.setupArrays(int_rule);

    // Fill node coordinates
    panzer::MDFieldArrayFactory af("prefix_",true);
    auto node_coordinates = af.buildStaticArray<double,Cell,NODE,Dim>("node_coordinates",2,4,2);
    {
      auto node_coordinates_host = Kokkos::create_mirror_view(node_coordinates.get_static_view());
      auto cell_vertices_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),mesh.cell_vertices);
      for(int i=0; i<mesh.cell_vertices.extent_int(0); ++i)
        for(int j=0; j<mesh.cell_vertices.extent_int(1); ++j)
          for(int k=0; k<mesh.cell_vertices.extent_int(2); ++k)
            node_coordinates_host(i,j,k) = cell_vertices_host(i,j,k);
      Kokkos::deep_copy(node_coordinates.get_static_view(),node_coordinates_host);
    }

    // Make sure the integration values fail to construct without the subcell connectivity
    TEST_THROW(int_values.evaluateValues(node_coordinates), std::logic_error);

    // Build the subcell connectivity
    auto connectivity = Teuchos::rcp(new panzer::FaceConnectivity);
    connectivity->setup(mesh);

    // Now we try again
    int_values.evaluateValues(node_coordinates,-1,connectivity);
    
    // This test will focus on normals and points
    const auto normals = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),int_values.surface_normals.get_static_view());
    const auto points = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),int_values.ip_coordinates.get_static_view());
    const int num_points_per_cell = points.extent_int(1);
    const int num_points_per_face = num_points_per_cell / 4;

    for(int cell=0; cell<2; ++cell){
      out << "CELL " << cell << ":\n";
      for(int face=0; face<4; ++face){
        out << "\tFACE " << face << ":\n";
        for(int face_point=0; face_point<num_points_per_face; ++face_point){
          const int point = face*num_points_per_face + face_point;
          out << "\t\tPOINT  "<<face_point<<": ("<<points(cell,point,0)<<", "<<points(cell,point,1)<<")\n";
          out << "\t\tNORMAL "<<face_point<<": ("<<normals(cell,point,0)<<", "<<normals(cell,point,1)<<")\n";
        }
      }
    }

    // We need to make sure that face 0 and 2 have aligned quadrature points

    const double tolerance = 1.e-14;
    
    // Face 0
    {
      const int cell_0 = connectivity->cellForSubcellHost(0,0);
      const int cell_1 = connectivity->cellForSubcellHost(0,1);

      const int lidx_0 = connectivity->localSubcellForSubcellHost(0,0);
      const int lidx_1 = connectivity->localSubcellForSubcellHost(0,1);

      TEST_EQUALITY(cell_0,0);
      TEST_EQUALITY(cell_1,1);

      TEST_EQUALITY(lidx_0,0);
      TEST_EQUALITY(lidx_1,2);

      const double sqrt5 = std::sqrt(5.);
      const double normal_0[2] = {-2./sqrt5,-1./sqrt5};
      const double normal_1[2] = { 2./sqrt5,-1./sqrt5};

      // Note that y will be equal, but x will be different
      
      for(int face_point=0; face_point<num_points_per_face; ++face_point){
        const int point_0 = lidx_0*num_points_per_face+face_point;
        const int point_1 = lidx_1*num_points_per_face+face_point;

        TEST_ASSERT(std::fabs(points(cell_0,point_0,0) - points(cell_1,point_1,0)) > 0.5);
        TEST_FLOATING_EQUALITY(points(cell_0,point_0,1), points(cell_1,point_1,1), tolerance);

	TEST_FLOATING_EQUALITY(normals(cell_0,point_0,0), normal_0[0], tolerance);
        TEST_FLOATING_EQUALITY(normals(cell_0,point_0,1), normal_0[1], tolerance);

	TEST_FLOATING_EQUALITY(normals(cell_1,point_1,0), normal_1[0], tolerance);
        TEST_FLOATING_EQUALITY(normals(cell_1,point_1,1), normal_1[1], tolerance);

      }
    }
    
    // Face 2
    {
      const int cell_0 = connectivity->cellForSubcellHost(2,0);
      const int cell_1 = connectivity->cellForSubcellHost(2,1);

      const int lidx_0 = connectivity->localSubcellForSubcellHost(2,0);
      const int lidx_1 = connectivity->localSubcellForSubcellHost(2,1);

      TEST_EQUALITY(cell_0,0);
      TEST_EQUALITY(cell_1,1);

      TEST_EQUALITY(lidx_0,2);
      TEST_EQUALITY(lidx_1,0);

      const double normal_0[2] = { 1.,0.};
      const double normal_1[2] = {-1.,0.};
      
      for(int face_point=0; face_point<num_points_per_face; ++face_point){
        const int point_0 = lidx_0*num_points_per_face+face_point;
        const int point_1 = lidx_1*num_points_per_face+face_point;

        TEST_FLOATING_EQUALITY(points(cell_0,point_0,0), points(cell_1,point_1,0), tolerance);
        TEST_FLOATING_EQUALITY(points(cell_0,point_0,1), points(cell_1,point_1,1), tolerance);

        TEST_FLOATING_EQUALITY(normals(cell_0,point_0,0), normal_0[0], tolerance);
        TEST_FLOATING_EQUALITY(normals(cell_0,point_0,1), normal_0[1], tolerance);
        
        TEST_FLOATING_EQUALITY(normals(cell_1,point_1,0), normal_1[0], tolerance);
        TEST_FLOATING_EQUALITY(normals(cell_1,point_1,1), normal_1[1], tolerance);

      }
    }
  }
}
