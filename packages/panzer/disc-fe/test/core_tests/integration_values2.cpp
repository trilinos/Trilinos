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
#include "Panzer_IntegrationValues2.hpp"
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
    
    panzer::MDFieldArrayFactory af("prefix_",true);

    const int num_nodes = int_rule->topology->getNodeCount();
    PHX::MDField<double,Cell,NODE,Dim> node_coordinates 
        = af.buildStaticArray<double,Cell,NODE,Dim>("nc",num_cells, num_nodes, base_cell_dimension);

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

    panzer::IntegrationValues2<double> int_values("prefix_");
    int_values.setup(int_rule, node_coordinates);
    
    TEST_EQUALITY(int_values.getCubaturePoints().extent(1), 4);
    double realspace_x_coord = (1.0/std::sqrt(3.0) + 1.0) / 2.0;
    double realspace_y_coord = (1.0/std::sqrt(3.0) + 1.0) / 2.0;
    auto tmp_coords = int_values.getCubaturePoints();
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

    panzer::MDFieldArrayFactory af("prefix_",true);

    cv_type = "side";
    RCP<IntegrationRule> int_rule_side =
      rcp(new IntegrationRule(cell_data, cv_type));

    const int num_nodes = int_rule_vol->topology->getNodeCount();
    PHX::MDField<double,Cell,NODE,Dim> node_coordinates
        = af.buildStaticArray<double,Cell,NODE,Dim>("nc",num_cells, num_nodes, base_cell_dimension);

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


    panzer::IntegrationValues2<double> int_values_vol("prefix_");
    int_values_vol.setup(int_rule_vol, node_coordinates);

    panzer::IntegrationValues2<double> int_values_side("prefix_");
    int_values_side.setup(int_rule_side, node_coordinates);

    TEST_EQUALITY(int_values_vol.getCubaturePoints().extent(1), 4);
    TEST_EQUALITY(int_values_side.getCubaturePoints().extent(1), 4);
    TEST_EQUALITY(int_values_side.getWeightedNormals().extent(1), 4);
    double realspace_x_coord_1 = 0.25;
    double realspace_y_coord_1 = 0.25;
    auto tmp_coords = int_values_vol.getCubaturePoints();
    auto tmp_coords_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),tmp_coords.get_view());
    TEST_FLOATING_EQUALITY(tmp_coords_host(0,0,0),
                           realspace_x_coord_1, 1.0e-8);
    TEST_FLOATING_EQUALITY(tmp_coords_host(0,0,1),
                           realspace_y_coord_1, 1.0e-8);
    double realspace_x_coord_2 = 0.5;
    double realspace_y_coord_2 = 0.25;
    auto tmp_side_coords = int_values_side.getCubaturePoints();
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

    panzer::MDFieldArrayFactory af("prefix_",true);

    const int num_nodes = int_rule_bc->topology->getNodeCount();
    PHX::MDField<double,Cell,NODE,Dim> node_coordinates
        = af.buildStaticArray<double,Cell,NODE,Dim>("nc",num_cells, num_nodes, base_cell_dimension);

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


    panzer::IntegrationValues2<double> int_values_bc("prefix_");
    int_values_bc.setup(int_rule_bc, node_coordinates);

    auto ip_coords_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),int_values_bc.getCubaturePoints().get_view());

    out << "Points:\n";
    for(int i=0; i<ip_coords_host.extent_int(0); ++i){
      for(int j=0; j<ip_coords_host.extent_int(1); ++j){
        out << "("<<i<<", "<<j<<"): ";
        for(int k=0; k<ip_coords_host.extent_int(2); ++k)
          out << ip_coords_host(i,j,k) << " ";
        out << "\n";
      }
    }

    TEST_EQUALITY(int_values_bc.getCubaturePoints().extent(1), 2);
    double realspace_x_coord_1 = 1.0;
    double realspace_y_coord_1 = 0.25;
    TEST_FLOATING_EQUALITY(ip_coords_host(0,0,0),
                           realspace_x_coord_1, 1.0e-8);
    TEST_FLOATING_EQUALITY(ip_coords_host(0,0,1),
                           realspace_y_coord_1, 1.0e-8);

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
      mesh.cell_nodes = PHX::View<double***>("cell_nodes",2,4,2);

      PHX::View<double**> coordinates("coordinates",6,2);
      auto coordinates_host = Kokkos::create_mirror_view(coordinates);
      coordinates_host(0,0) = 1.; coordinates_host(0,1) = 3.;
      coordinates_host(1,0) = 3.; coordinates_host(1,1) = 3.;
      coordinates_host(2,0) = 5.; coordinates_host(2,1) = 3.;
      coordinates_host(3,0) = 2.; coordinates_host(3,1) = 1.;
      coordinates_host(4,0) = 3.; coordinates_host(4,1) = 1.;
      coordinates_host(5,0) = 4.; coordinates_host(5,1) = 1.;
      Kokkos::deep_copy(coordinates,coordinates_host);

      auto cell_nodes_host = Kokkos::create_mirror_view(mesh.cell_nodes);
#define SET_NC(cell,node,vertex) cell_nodes_host(cell,node,0) = coordinates_host(vertex,0); cell_nodes_host(cell,node,1) = coordinates_host(vertex,1);
      SET_NC(0,0, 0);
      SET_NC(0,1, 3);
      SET_NC(0,2, 4);
      SET_NC(0,3, 1);
      SET_NC(1,0, 1);
      SET_NC(1,1, 4);
      SET_NC(1,2, 5);
      SET_NC(1,3, 2);
#undef SET_NC
      Kokkos::deep_copy(mesh.cell_nodes,cell_nodes_host);

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

    // Fill node coordinates
    panzer::MDFieldArrayFactory af("prefix_",true);
    auto node_coordinates = af.buildStaticArray<double,Cell,NODE,Dim>("node_coordinates",2,4,2);
    {
      auto node_coordinates_host = Kokkos::create_mirror_view(node_coordinates.get_static_view());
      auto cell_nodes_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),mesh.cell_nodes);
      for(int i=0; i<mesh.cell_nodes.extent_int(0); ++i)
        for(int j=0; j<mesh.cell_nodes.extent_int(1); ++j)
          for(int k=0; k<mesh.cell_nodes.extent_int(2); ++k)
            node_coordinates_host(i,j,k) = cell_nodes_host(i,j,k);
      Kokkos::deep_copy(node_coordinates.get_static_view(),node_coordinates_host);
    }

    panzer::IntegrationValues2<double> int_values("prefix_");
    int_values.setup(int_rule, node_coordinates);

    // Make sure the integration values fail to construct without the subcell connectivity
    TEST_THROW(int_values.evaluateValues(node_coordinates), std::logic_error);

    // Build the subcell connectivity
    auto connectivity = Teuchos::rcp(new panzer::FaceConnectivity);
    connectivity->setup(mesh);

    // Now we try again
    int_values.evaluateValues(node_coordinates,-1,connectivity,0);

    // This test will focus on normals and points
    const auto normals = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),int_values.getSurfaceNormals().get_static_view());
    const auto points = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),int_values.getCubaturePoints().get_static_view());
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
      out << "Comparing cell_0 "<<cell_0<<" to cell_1 "<<cell_1<<"\n";

      for(int face_point=0; face_point<num_points_per_face; ++face_point){
        const int point_0 = lidx_0*num_points_per_face+face_point;
        const int point_1 = lidx_1*num_points_per_face+face_point;

        out << "Comparing point_0 "<<point_0<<" to point_1 "<<point_1<<"\n";

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

      out << "Comparing cell_0 "<<cell_0<<" to cell_1 "<<cell_1<<"\n";

      for(int face_point=0; face_point<num_points_per_face; ++face_point){
        const int point_0 = lidx_0*num_points_per_face+face_point;
        const int point_1 = lidx_1*num_points_per_face+face_point;

        out << "Comparing point_0 "<<point_0<<" to point_1 "<<point_1<<"\n";

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
