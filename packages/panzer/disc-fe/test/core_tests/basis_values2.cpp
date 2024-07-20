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
#include "Panzer_BasisValues2.hpp"
#include "Panzer_PointValues2.hpp"
#include "Panzer_CommonArrayFactories.hpp"
#include "Panzer_Traits.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using panzer::IntegrationRule;

namespace panzer {

  TEUCHOS_UNIT_TEST(basis_values, const_basis)
  {
    typedef panzer::ArrayTraits<double,PHX::MDField<double> >::size_type size_type;

    Teuchos::RCP<shards::CellTopology> topo =
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    const int num_cells = 4;
    const int base_cell_dimension = 2;
    const panzer::CellData cell_data(num_cells,topo);

    const int cubature_degree = 4;
    RCP<IntegrationRule> int_rule =
      rcp(new IntegrationRule(cubature_degree, cell_data));
    const int num_qp = int_rule->num_points;

    panzer::IntegrationValues2<double> int_values("prefix_",true);
    int_values.setupArrays(int_rule);

    const int num_vertices = int_rule->topology->getNodeCount();
    panzer::MDFieldArrayFactory af("prefix_",true);
    PHX::MDField<double,Cell,NODE,Dim> node_coordinates = af.buildStaticArray<double,Cell,NODE,Dim>("nc",num_cells, num_vertices,
                                            base_cell_dimension);


    // Set up node coordinates.  Here we assume the following
    // ordering.  This needs to be consistent with shards topology,
    // otherwise we will get negative determinates

    // 3(0,1)---2(1,1)
    //   |    0  |
    //   |       |
    // 0(0,0)---1(1,0)

    const size_type x = 0;
    const size_type y = 1;
    auto node_coordinates_host = Kokkos::create_mirror_view(node_coordinates.get_static_view());
    for (size_type cell = 0; cell < node_coordinates.extent(0); ++cell) {
      int xleft = cell % 2;
      int yleft = int(cell/2);

      node_coordinates_host(cell,0,x) = xleft*0.5;
      node_coordinates_host(cell,0,y) = yleft*0.5;

      node_coordinates_host(cell,1,x) = (xleft+1)*0.5;
      node_coordinates_host(cell,1,y) = yleft*0.5;

      node_coordinates_host(cell,2,x) = (xleft+1)*0.5;
      node_coordinates_host(cell,2,y) = (yleft+1)*0.5;

      node_coordinates_host(cell,3,x) = xleft*0.5;
      node_coordinates_host(cell,3,y) = (yleft+1)*0.5;

      out << "Cell " << cell << " = ";
      for(int i=0;i<4;i++)
         out << "(" << node_coordinates_host(cell,i,x) << ", "
                    << node_coordinates_host(cell,i,y) << ") ";
      out << std::endl;
    }
    Kokkos::deep_copy(node_coordinates.get_static_view(),node_coordinates_host);

    int_values.evaluateValues(node_coordinates);

    const std::string basis_type = "Const";

    RCP<panzer::BasisIRLayout> basis = rcp(new panzer::BasisIRLayout(basis_type, 0, *int_rule));
    TEST_EQUALITY(basis->getIntrepid2Basis()->requireOrientation(), false);

    panzer::BasisValues2<double> basis_values("",true,true);

    out << "setupArrays" << std::endl;
    basis_values.setupArrays(basis);

    out << "evaluateValues" << std::endl;
    basis_values.evaluateValues(int_values.cub_points,
                                int_values.jac,
                                int_values.jac_det,
                                int_values.jac_inv,
                                int_values.weighted_measure,
                                node_coordinates);

    out << "check output" << std::endl;
    double relCellVol = 0.25*0.25; // this is the relative (to the reference cell) volume
    auto cub_weights_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),int_values.cub_weights.get_static_view());
    auto cub_weighted_measure_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),int_values.weighted_measure.get_view()); 
    auto basis_ref_scalar_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),basis_values.basis_ref_scalar.get_static_view());
    auto jac_det_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),int_values.jac_det.get_static_view());
    auto weighted_basis_scalar_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),basis_values.weighted_basis_scalar.get_static_view());
    auto basis_scalar_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),basis_values.basis_scalar.get_static_view());
    for(int i=0;i<num_qp;i++) {
       // double x = int_values.cub_points(i,0);
       // double y = int_values.cub_points(i,1);
       double weight = cub_weights_host(i);

       // check reference values
       TEST_EQUALITY(basis_ref_scalar_host(0,i),1.0);

       // check basis values
       for(int cell=0;cell<num_cells;cell++) {

          out << "wm = " << cub_weighted_measure_host(cell,i) << " " << std::endl;;
          TEST_EQUALITY(jac_det_host(cell,i),relCellVol);

          // check out basis on transformed elemented
          TEST_EQUALITY(basis_ref_scalar_host(0,i),basis_scalar_host(cell,0,i));
          TEST_EQUALITY(weighted_basis_scalar_host(cell,0,i),relCellVol*weight*basis_scalar_host(cell,0,i));
       }
       out << std::endl;
    }
  }

  TEUCHOS_UNIT_TEST(basis_values, grad_quad)
  {
    typedef panzer::ArrayTraits<double,PHX::MDField<double> >::size_type size_type;

    Teuchos::RCP<shards::CellTopology> topo =
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    const int num_cells = 4;
    const int base_cell_dimension = 2;
    const panzer::CellData cell_data(num_cells,topo);

    const int cubature_degree = 4;
    RCP<IntegrationRule> int_rule =
      rcp(new IntegrationRule(cubature_degree, cell_data));
    const int num_qp = int_rule->num_points;

    panzer::IntegrationValues2<double> int_values("prefix_",true);
    int_values.setupArrays(int_rule);

    const int num_vertices = int_rule->topology->getNodeCount();
    panzer::MDFieldArrayFactory af("prefix_",true);
    PHX::MDField<double,Cell,NODE,Dim> node_coordinates = af.buildStaticArray<double,Cell,NODE,Dim>("nc",num_cells, num_vertices,
                                            base_cell_dimension);



    // Set up node coordinates.  Here we assume the following
    // ordering.  This needs to be consistent with shards topology,
    // otherwise we will get negative determinates

    // 3(0,1)---2(1,1)
    //   |    0  |
    //   |       |
    // 0(0,0)---1(1,0)
    {
      const size_type x = 0;
      const size_type y = 1;
      auto node_coordinates_host = Kokkos::create_mirror_view(node_coordinates.get_view());
      for (size_type cell = 0; cell < node_coordinates.extent(0); ++cell) {
        int xleft = cell % 2;
        int yleft = int(cell/2);

        node_coordinates_host(cell,0,x) = xleft*0.5;
        node_coordinates_host(cell,0,y) = yleft*0.5;

        node_coordinates_host(cell,1,x) = (xleft+1)*0.5;
        node_coordinates_host(cell,1,y) = yleft*0.5;

        node_coordinates_host(cell,2,x) = (xleft+1)*0.5;
        node_coordinates_host(cell,2,y) = (yleft+1)*0.5;

        node_coordinates_host(cell,3,x) = xleft*0.5;
        node_coordinates_host(cell,3,y) = (yleft+1)*0.5;

        out << "Cell " << cell << " = ";
        for(int i=0;i<4;i++)
          out << "(" << node_coordinates_host(cell,i,x) << ", "
              << node_coordinates_host(cell,i,y) << ") ";
        out << std::endl;
      }
      Kokkos::deep_copy(node_coordinates.get_view(),node_coordinates_host);
    }

    int_values.evaluateValues(node_coordinates);

    const std::string basis_type = "Q1";

    RCP<panzer::BasisIRLayout> basis = rcp(new panzer::BasisIRLayout(basis_type, 0, *int_rule));
    TEST_EQUALITY(basis->getIntrepid2Basis()->requireOrientation(), false);

    panzer::BasisValues2<double> basis_values("",true,true);
    basis_values.setupArrays(basis);
    basis_values.evaluateValues(int_values.cub_points,
                                int_values.jac,
                                int_values.jac_det,
                                int_values.jac_inv,
				int_values.weighted_measure,
                                node_coordinates);

    double relCellVol = 0.25*0.25; // this is the relative (to the reference cell) volume
    auto cub_points_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),int_values.cub_points.get_view());
    auto cub_weights_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),int_values.cub_weights.get_view());
    auto basis_ref_scalar_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),basis_values.basis_ref_scalar.get_view());
    auto basis_scalar_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),basis_values.basis_scalar.get_view());
    auto grad_basis_ref_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),basis_values.grad_basis_ref.get_view());
    auto weighted_basis_scalar_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),basis_values.weighted_basis_scalar.get_view());
    auto grad_basis_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),basis_values.grad_basis.get_view());
    auto weighted_grad_basis_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),basis_values.weighted_grad_basis.get_view());
    auto jac_det_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),int_values.jac_det.get_view());
    for(int i=0;i<num_qp;i++) {
       double x = cub_points_host(i,0);
       double y = cub_points_host(i,1);
       double weight = cub_weights_host(i);

       // check reference values
       TEST_EQUALITY(basis_ref_scalar_host(0,i),0.25*(x-1.0)*(y-1.0));
       TEST_EQUALITY(grad_basis_ref_host(0,i,0),0.25*(y-1.0));
       TEST_EQUALITY(grad_basis_ref_host(0,i,1),0.25*(x-1.0));

       // check basis values
       for(int cell=0;cell<num_cells;cell++) {

          TEST_EQUALITY(jac_det_host(cell,i),relCellVol);

          // check out basis on transformed elemented
          TEST_EQUALITY(basis_ref_scalar_host(0,i),basis_scalar_host(cell,0,i));
          TEST_EQUALITY(weighted_basis_scalar_host(cell,0,i),relCellVol*weight*basis_scalar_host(cell,0,i));

          TEST_EQUALITY(grad_basis_host(cell,0,i,0),4.0*grad_basis_ref_host(0,i,0));
          TEST_EQUALITY(grad_basis_host(cell,0,i,1),4.0*grad_basis_ref_host(0,i,1));

          TEST_EQUALITY(weighted_grad_basis_host(cell,0,i,0),relCellVol*weight*grad_basis_host(cell,0,i,0));
          TEST_EQUALITY(weighted_grad_basis_host(cell,0,i,1),relCellVol*weight*grad_basis_host(cell,0,i,1));
       }
    }
  }

  TEUCHOS_UNIT_TEST(basis_values, hcurl_basis)
  {
    typedef panzer::ArrayTraits<double,PHX::MDField<double> >::size_type size_type;

    Teuchos::RCP<shards::CellTopology> topo =
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    const int num_cells = 4;
    const int base_cell_dimension = 2;
    const panzer::CellData cell_data(num_cells,topo);

    const int cubature_degree = 4;
    RCP<IntegrationRule> int_rule =
      rcp(new IntegrationRule(cubature_degree, cell_data));
    const unsigned int num_qp = Teuchos::as<unsigned int>(int_rule->num_points);

    panzer::IntegrationValues2<double> int_values("prefix_",true);
    int_values.setupArrays(int_rule);

    const int num_vertices = int_rule->topology->getNodeCount();
    panzer::MDFieldArrayFactory af("prefix_",true);
    PHX::MDField<double,Cell,NODE,Dim> node_coordinates
        = af.buildStaticArray<double,Cell,NODE,Dim>("nc",num_cells, num_vertices, base_cell_dimension);
    // const int num_edges = int_rule->topology->getEdgeCount();
    // FieldContainer<double> edge_orientation(num_cells, num_edges);

    // Set up node coordinates.  Here we assume the following
    // ordering.  This needs to be consistent with shards topology,
    // otherwise we will get negative determinates

    // 3(0,1)---2(1,1)
    //   |    0  |
    //   |       |
    // 0(0,0)---1(1,0)

    // and now the edges

    //   +---2---+
    //   |       |
    //   3       1
    //   |       |
    //   +---0---+
    {
      const size_type x = 0;
      const size_type y = 1;
      auto node_coordinates_host = Kokkos::create_mirror_view(node_coordinates.get_view());
      for (size_type cell = 0; cell < node_coordinates.extent(0); ++cell) {
        int xleft = cell % 2;
        int yleft = int(cell/2);

        node_coordinates_host(cell,0,x) = xleft*0.5;
        node_coordinates_host(cell,0,y) = yleft*0.5;

        node_coordinates_host(cell,1,x) = (xleft+1)*0.5;
        node_coordinates_host(cell,1,y) = yleft*0.5;

        node_coordinates_host(cell,2,x) = (xleft+1)*0.5;
        node_coordinates_host(cell,2,y) = (yleft+1)*0.5;

        node_coordinates_host(cell,3,x) = xleft*0.5;
        node_coordinates_host(cell,3,y) = (yleft+1)*0.5;

        out << "Cell " << cell << " = ";
        for(int i=0;i<4;i++)
          out << "(" << node_coordinates_host(cell,i,x) << ", "
              << node_coordinates_host(cell,i,y) << ") ";
        out << std::endl;
      }
      Kokkos::deep_copy(node_coordinates.get_view(),node_coordinates_host);
    }

    int_values.evaluateValues(node_coordinates);

    const std::string basis_type = "HCurl";

    Teuchos::RCP<PureBasis> basis = Teuchos::rcp(new PureBasis(basis_type,1,cell_data));
    RCP<panzer::BasisIRLayout> basisIRLayout = rcp(new panzer::BasisIRLayout(basis, *int_rule));
    TEST_EQUALITY(basis->getIntrepid2Basis()->requireOrientation(), true);

    panzer::BasisValues2<double> basis_values("",true,true);

    basis_values.setupArrays(basisIRLayout);
    basis_values.evaluateValues(int_values.cub_points,
                                int_values.jac,
                                int_values.jac_det,
                                int_values.jac_inv,
				int_values.weighted_measure,
                                node_coordinates);

    TEST_EQUALITY(basis_values.basis_ref_vector.extent(0),4);
    TEST_EQUALITY(basis_values.basis_ref_vector.extent(1),num_qp);
    TEST_EQUALITY(basis_values.basis_ref_vector.extent(2),2);

    TEST_EQUALITY(basis_values.grad_basis_ref.size(),0);

    TEST_EQUALITY(basis_values.curl_basis_ref_scalar.extent(0),4);
    TEST_EQUALITY(basis_values.curl_basis_ref_scalar.extent(1),num_qp);

    double relCellVol = 0.25*0.25; // this is the relative (to the reference cell) volume
    auto cub_points_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),int_values.cub_points.get_view());
    auto cub_weights_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),int_values.cub_weights.get_view());
    auto jac_det_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),int_values.jac_det.get_view());
    auto basis_ref_vector_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),basis_values.basis_ref_vector.get_view());
    auto basis_vector_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),basis_values.basis_vector.get_view());
    auto curl_basis_ref_scalar_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),basis_values.curl_basis_ref_scalar.get_view());
    auto curl_basis_scalar_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),basis_values.curl_basis_scalar.get_view());
    auto weighted_basis_vector_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),basis_values.weighted_basis_vector.get_view());
    auto weighted_curl_basis_scalar_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),basis_values.weighted_curl_basis_scalar.get_view());
    for(unsigned int i=0;i<num_qp;i++) {
       // double x = int_values.cub_points(i,0);
       double y = cub_points_host(i,1);
       double weight = cub_weights_host(i);

       // check reference values
       TEST_EQUALITY(basis_ref_vector_host(0,i,0),-0.5*(y-1.0));
       TEST_EQUALITY(basis_ref_vector_host(0,i,1),0.0);
       TEST_EQUALITY(curl_basis_ref_scalar_host(0,i),0.5);

       // check basis values
       for(int cell=0;cell<num_cells;cell++) {

          TEST_EQUALITY(jac_det_host(cell,i),relCellVol);

          // check out basis on transformed elemented
          TEST_EQUALITY(basis_ref_vector_host(0,i,0),0.25*basis_vector_host(cell,0,i,0));
          TEST_EQUALITY(basis_ref_vector_host(0,i,1),0.25*basis_vector_host(cell,0,i,1));
          TEST_EQUALITY(curl_basis_ref_scalar_host(0,i),relCellVol*curl_basis_scalar_host(cell,0,i));

          TEST_EQUALITY(basis_ref_vector_host(1,i,0),0.25*basis_vector_host(cell,1,i,0));
          TEST_EQUALITY(basis_ref_vector_host(1,i,1),0.25*basis_vector_host(cell,1,i,1));
          TEST_EQUALITY(curl_basis_ref_scalar_host(1,i),relCellVol*curl_basis_scalar_host(cell,1,i));

          TEST_EQUALITY(weighted_basis_vector_host(cell,0,i,0),relCellVol*weight*basis_vector_host(cell,0,i,0));
          TEST_EQUALITY(weighted_basis_vector_host(cell,0,i,1),relCellVol*weight*basis_vector_host(cell,0,i,1));

          TEST_EQUALITY(weighted_basis_vector_host(cell,1,i,0),relCellVol*weight*basis_vector_host(cell,1,i,0));
          TEST_EQUALITY(weighted_basis_vector_host(cell,1,i,1),relCellVol*weight*basis_vector_host(cell,1,i,1));

          TEST_EQUALITY(weighted_curl_basis_scalar_host(cell,0,i),relCellVol*weight*curl_basis_scalar_host(cell,0,i));
          TEST_EQUALITY(weighted_curl_basis_scalar_host(cell,1,i),relCellVol*weight*curl_basis_scalar_host(cell,1,i));
       }
    }
  }

  TEUCHOS_UNIT_TEST(basis_values, hdiv_basis)
  {
    typedef panzer::ArrayTraits<double,PHX::MDField<double> >::size_type size_type;

    Teuchos::RCP<shards::CellTopology> topo =
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<4> >()));

    const int num_cells = 4;
    const int base_cell_dimension = 3;
    const panzer::CellData cell_data(num_cells,topo);

    const int cubature_degree = 4;
    RCP<IntegrationRule> int_rule =
      rcp(new IntegrationRule(cubature_degree, cell_data));
    const unsigned int num_qp = Teuchos::as<unsigned int>(int_rule->num_points);

    panzer::IntegrationValues2<double> int_values("prefix_",true);
    int_values.setupArrays(int_rule);

    const int num_vertices = int_rule->topology->getNodeCount();
    panzer::MDFieldArrayFactory af("prefix_",true);
    PHX::MDField<double,Cell,NODE,Dim> node_coordinates
        = af.buildStaticArray<double,Cell,NODE,Dim>("nc",num_cells, num_vertices, base_cell_dimension);
    // const int num_edges = int_rule->topology->getEdgeCount();
    // FieldContainer<double> edge_orientation(num_cells, num_edges);

    {
      const size_type x = 0;
      const size_type y = 1;
      const size_type z = 2;
      auto node_coordinates_host = Kokkos::create_mirror_view(node_coordinates.get_view());
      for (size_type cell = 0; cell < node_coordinates.extent(0); ++cell) {
        int xleft = cell % 2;
        int yleft = int(cell/2);

        node_coordinates_host(cell,0,x) = xleft*0.5;
        node_coordinates_host(cell,0,y) = yleft*0.5;
        node_coordinates_host(cell,0,z) = 0.0;

        node_coordinates_host(cell,1,x) = (xleft+1)*0.5;
        node_coordinates_host(cell,1,y) = yleft*0.5;
        node_coordinates_host(cell,1,z) = 0.0;

        node_coordinates_host(cell,2,x) = xleft*0.5;
        node_coordinates_host(cell,2,y) = (yleft+1)*0.5;
        node_coordinates_host(cell,2,z) = 0.0;

        node_coordinates_host(cell,3,x) = xleft*0.5;
        node_coordinates_host(cell,3,y) = yleft*0.5;
        node_coordinates_host(cell,3,z) = 1.0/3.0;

        out << "Cell " << cell << " = ";
        for(int i=0;i<4;i++)
          out << "(" << node_coordinates_host(cell,i,x) << ", "
              << node_coordinates_host(cell,i,y) << ", "
              << node_coordinates_host(cell,i,z) << ") ";
        out << std::endl;
      }
      Kokkos::deep_copy(node_coordinates.get_view(),node_coordinates_host);
    }

    int_values.evaluateValues(node_coordinates);

    const std::string basis_type = "HDiv";

    Teuchos::RCP<PureBasis> basis = Teuchos::rcp(new PureBasis(basis_type,1,cell_data));
    RCP<panzer::BasisIRLayout> basisIRLayout = rcp(new panzer::BasisIRLayout(basis, *int_rule));
    TEST_EQUALITY(basis->getIntrepid2Basis()->requireOrientation(), true);

    panzer::BasisValues2<double> basis_values("",true,true);

    basis_values.setupArrays(basisIRLayout);
    basis_values.evaluateValues(int_values.cub_points,
                                int_values.jac,
                                int_values.jac_det,
                                int_values.jac_inv,
				int_values.weighted_measure,
                                node_coordinates);

    TEST_EQUALITY(basis_values.basis_ref_vector.extent(0),4);
    TEST_EQUALITY(basis_values.basis_ref_vector.extent(1),num_qp);
    TEST_EQUALITY(basis_values.basis_ref_vector.extent_int(2),base_cell_dimension);

    TEST_EQUALITY(basis_values.grad_basis_ref.size(),0);

    TEST_EQUALITY(basis_values.div_basis_ref.extent(0),4);
    TEST_EQUALITY(basis_values.div_basis_ref.extent(1),num_qp);

    double relCellVol = (1.0/72) / (1.0/6.0);
    auto cub_points_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),int_values.cub_points.get_view());
    auto cub_weights_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),int_values.cub_weights.get_view());
    auto jac_det_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),int_values.jac_det.get_view());
    auto basis_ref_vector_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),basis_values.basis_ref_vector.get_view());
    auto basis_vector_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),basis_values.basis_vector.get_view());
    auto div_basis_ref_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),basis_values.div_basis_ref.get_view());
    auto div_basis_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),basis_values.div_basis.get_view());
    auto weighted_basis_vector_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),basis_values.weighted_basis_vector.get_view());
    auto weighted_div_basis_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),basis_values.weighted_div_basis.get_view());
    for(unsigned int i=0;i<num_qp;i++) {
       double x = cub_points_host(i,0);
       double y = cub_points_host(i,1);
       double z = cub_points_host(i,2);
       double weight = cub_weights_host(i);

       // check reference values
       TEST_EQUALITY(basis_ref_vector_host(0,i,0),x);
       TEST_EQUALITY(basis_ref_vector_host(0,i,1),(y-1.0));
       TEST_EQUALITY(basis_ref_vector_host(0,i,2),z);

       TEST_EQUALITY(div_basis_ref_host(0,i),3.0);

       // check basis values
       for(int cell=0;cell<num_cells;cell++) {

          TEST_EQUALITY(jac_det_host(cell,i),relCellVol);

          // check out basis on transformed elemented
          TEST_EQUALITY(0.5*basis_ref_vector_host(0,i,0)/relCellVol,basis_vector_host(cell,0,i,0));
          TEST_EQUALITY(0.5*basis_ref_vector_host(0,i,1)/relCellVol,basis_vector_host(cell,0,i,1));
          TEST_EQUALITY((1.0/3.0)*basis_ref_vector_host(0,i,2)/relCellVol,basis_vector_host(cell,0,i,2));

          TEST_EQUALITY(div_basis_ref_host(0,i)/relCellVol,div_basis_host(cell,0,i));
          TEST_EQUALITY(weighted_basis_vector_host(cell,0,i,0),relCellVol*weight*basis_vector_host(cell,0,i,0));
          TEST_EQUALITY(weighted_basis_vector_host(cell,0,i,1),relCellVol*weight*basis_vector_host(cell,0,i,1));
          TEST_EQUALITY(weighted_basis_vector_host(cell,0,i,2),relCellVol*weight*basis_vector_host(cell,0,i,2));

          TEST_EQUALITY(weighted_basis_vector_host(cell,1,i,0),relCellVol*weight*basis_vector_host(cell,1,i,0));
          TEST_EQUALITY(weighted_basis_vector_host(cell,1,i,1),relCellVol*weight*basis_vector_host(cell,1,i,1));
          TEST_EQUALITY(weighted_basis_vector_host(cell,1,i,2),relCellVol*weight*basis_vector_host(cell,1,i,2));

          TEST_EQUALITY(weighted_div_basis_host(cell,0,i),relCellVol*weight*div_basis_host(cell,0,i));
          TEST_EQUALITY(weighted_div_basis_host(cell,1,i),relCellVol*weight*div_basis_host(cell,1,i));
       }
    }

    std::vector<Intrepid2::Orientation> orientations;

    int faceOrts[4][4], cnt = 0;
    for(int cell=0;cell<num_cells;cell++) {
      faceOrts[cell][0] = cnt++ % 6;
      faceOrts[cell][1] = cnt++ % 6;
      faceOrts[cell][2] = cnt++ % 6;
      faceOrts[cell][3] = cnt++ % 6;

      Intrepid2::Orientation ort;
      ort.setFaceOrientation(4, faceOrts[cell]);
      orientations.push_back(ort);
    }

    for(int i=0;i<4;i++) {
      out << "Orientations: Cell " << i << " = ";
      out << faceOrts[i][0] << ", "
          << faceOrts[i][1] << ", "
          << faceOrts[i][2] << ", "
          << faceOrts[i][3] << std::endl;
    }
    out << std::endl;

    basis_values.applyOrientations(orientations);

    Kokkos::deep_copy(basis_vector_host,basis_values.basis_vector.get_view());
    Kokkos::deep_copy(div_basis_host,basis_values.div_basis.get_view());
    Kokkos::deep_copy(weighted_basis_vector_host,basis_values.weighted_basis_vector.get_view());
    Kokkos::deep_copy(weighted_div_basis_host,basis_values.weighted_div_basis.get_view());

    // check with orientations
    const double ortVal[6] = {  1.0,  1.0,  1.0,
                               -1.0, -1.0, -1.0 };
    for(unsigned int i=0;i<num_qp;i++) {
       double weight = cub_weights_host(i);

       // check basis values
       for(int cell=0;cell<num_cells;cell++) {
         TEST_EQUALITY(jac_det_host(cell,i),relCellVol);

         for(int b=0;b<4;b++) {
           // check out basis on transformed elemented
           TEST_EQUALITY(ortVal[faceOrts[cell][b]]*0.5      *basis_ref_vector_host(b,i,0)/relCellVol,  basis_vector_host(cell,b,i,0));
           TEST_EQUALITY(ortVal[faceOrts[cell][b]]*0.5      *basis_ref_vector_host(b,i,1)/relCellVol,  basis_vector_host(cell,b,i,1));
           TEST_EQUALITY(ortVal[faceOrts[cell][b]]*(1.0/3.0)*basis_ref_vector_host(b,i,2)/relCellVol,  basis_vector_host(cell,b,i,2));

           TEST_EQUALITY(ortVal[faceOrts[cell][b]]*div_basis_ref_host(b,i)/relCellVol,div_basis_host(cell,b,i));

           TEST_EQUALITY(weighted_basis_vector_host(cell,b,i,0),relCellVol*weight*basis_vector_host(cell,b,i,0));
           TEST_EQUALITY(weighted_basis_vector_host(cell,b,i,1),relCellVol*weight*basis_vector_host(cell,b,i,1));
           TEST_EQUALITY(weighted_basis_vector_host(cell,b,i,2),relCellVol*weight*basis_vector_host(cell,b,i,2));

           TEST_EQUALITY(weighted_div_basis_host(cell,b,i),relCellVol*weight*div_basis_host(cell,b,i));
         }
       }
    }
  }


  TEUCHOS_UNIT_TEST(basis_values, hcurl_basis_3d_pv)
  {
    typedef panzer::ArrayTraits<double,PHX::MDField<double> >::size_type size_type;

    Teuchos::RCP<shards::CellTopology> topo =
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<8> >()));

    const int num_cells = 4;
    const int base_cell_dimension = 3;
    const panzer::CellData cell_data(num_cells,topo);

    const int cubature_degree = 4;
    RCP<IntegrationRule> int_rule =
      rcp(new IntegrationRule(cubature_degree, cell_data));
    const unsigned int num_qp = Teuchos::as<unsigned int>(int_rule->num_points);

    panzer::IntegrationValues2<double> int_values("prefix_",true);
    int_values.setupArrays(int_rule);

    const int num_vertices = int_rule->topology->getNodeCount();
    panzer::MDFieldArrayFactory af("prefix_",true);
    PHX::MDField<double,Cell,NODE,Dim> node_coordinates
        = af.buildStaticArray<double,Cell,NODE,Dim>("nc",num_cells, num_vertices, base_cell_dimension);
    // const int num_edges = int_rule->topology->getEdgeCount();
    // FieldContainer<double> edge_orientation(num_cells, num_edges);

    // Set up node coordinates.  Here we assume the following
    // ordering.  This needs to be consistent with shards topology,
    // otherwise we will get negative determinates

    // 3(0,1)---2(1,1)
    //   |    0  |
    //   |       |
    // 0(0,0)---1(1,0)

    // and now the edges

    //   +---2---+
    //   |       |
    //   3       1
    //   |       |
    //   +---0---+

    const size_type x = 0;
    const size_type y = 1;
    const size_type z = 2;
    auto node_coordinates_host = Kokkos::create_mirror_view(node_coordinates.get_view());
    for (size_type cell = 0; cell < node_coordinates.extent(0); ++cell) {
      int znum = cell % 4;
      int xleft = znum % 2;
      int yleft = int(znum/2);
      int zleft = int(cell/4);

      node_coordinates_host(cell,0,x) = xleft*0.5;
      node_coordinates_host(cell,0,y) = yleft*0.5;
      node_coordinates_host(cell,0,z) = zleft*0.5;

      node_coordinates_host(cell,1,x) = (xleft+1)*0.5;
      node_coordinates_host(cell,1,y) = yleft*0.5;
      node_coordinates_host(cell,1,z) = zleft*0.5;

      node_coordinates_host(cell,2,x) = (xleft+1)*0.5;
      node_coordinates_host(cell,2,y) = (yleft+1)*0.5;
      node_coordinates_host(cell,2,z) = zleft*0.5;

      node_coordinates_host(cell,3,x) = xleft*0.5;
      node_coordinates_host(cell,3,y) = (yleft+1)*0.5;
      node_coordinates_host(cell,3,z) = zleft*0.5;

      node_coordinates_host(cell,4,x) = xleft*0.5;
      node_coordinates_host(cell,4,y) = yleft*0.5;
      node_coordinates_host(cell,4,z) = (zleft+1)*0.5;

      node_coordinates_host(cell,5,x) = (xleft+1)*0.5;
      node_coordinates_host(cell,5,y) = yleft*0.5;
      node_coordinates_host(cell,5,z) = (zleft+1)*0.5;

      node_coordinates_host(cell,6,x) = (xleft+1)*0.5;
      node_coordinates_host(cell,6,y) = (yleft+1)*0.5;
      node_coordinates_host(cell,6,z) = (zleft+1)*0.5;

      node_coordinates_host(cell,7,x) = xleft*0.5;
      node_coordinates_host(cell,7,y) = (yleft+1)*0.5;
      node_coordinates_host(cell,7,z) = (zleft+1)*0.5;

      out << "Cell " << cell << " = ";
      for(int i=0;i<8;i++)
         out << "(" << node_coordinates_host(cell,i,x) << ", "
                    << node_coordinates_host(cell,i,y) << ", "
                    << node_coordinates_host(cell,i,z) << ") ";
      out << std::endl;
    }
    Kokkos::deep_copy(node_coordinates.get_view(),node_coordinates_host);

    int_values.evaluateValues(node_coordinates);

    const std::string basis_type = "QEdge1";

    Teuchos::RCP<PureBasis> basis = Teuchos::rcp(new PureBasis(basis_type,1,cell_data));
    RCP<panzer::BasisIRLayout> basisIRLayout = rcp(new panzer::BasisIRLayout(basis, *int_rule));
    TEST_EQUALITY(basis->getIntrepid2Basis()->requireOrientation(), true);

    panzer::BasisValues2<double> basis_values("prefix_",true,true);

    basis_values.setupArrays(basisIRLayout);
    basis_values.evaluateValues(int_values.cub_points,
                                int_values.jac,
                                int_values.jac_det,
                                int_values.jac_inv);

    TEST_EQUALITY(basis_values.basis_ref_vector.extent(0),12);
    TEST_EQUALITY(basis_values.basis_ref_vector.extent(1),num_qp);
    TEST_EQUALITY(basis_values.basis_ref_vector.extent(2),3);
    TEST_EQUALITY(basis_values.weighted_basis_vector.extent_int(0),num_cells);
    TEST_EQUALITY(basis_values.weighted_basis_vector.extent(1),12);
    TEST_EQUALITY(basis_values.weighted_basis_vector.extent(2),num_qp);
    TEST_EQUALITY(basis_values.weighted_basis_vector.extent(3),3);

    TEST_EQUALITY(basis_values.grad_basis_ref.size(),0);
    TEST_EQUALITY(basis_values.weighted_grad_basis.size(),0);
    TEST_EQUALITY(basis_values.curl_basis_ref_scalar.size(),0);
    TEST_EQUALITY(basis_values.curl_basis_scalar.size(),0);
    TEST_EQUALITY(basis_values.weighted_curl_basis_scalar.size(),0);

    TEST_EQUALITY(basis_values.curl_basis_ref_vector.extent(0),12);
    TEST_EQUALITY(basis_values.curl_basis_ref_vector.extent(1),num_qp);
    TEST_EQUALITY(basis_values.curl_basis_ref_vector.extent(2),3);
    TEST_EQUALITY(basis_values.curl_basis_vector.extent_int(0),num_cells);
    TEST_EQUALITY(basis_values.curl_basis_vector.extent(1),12);
    TEST_EQUALITY(basis_values.curl_basis_vector.extent(2),num_qp);
    TEST_EQUALITY(basis_values.curl_basis_vector.extent(3),3);
    TEST_EQUALITY(basis_values.weighted_curl_basis_vector.extent_int(0),num_cells);
    TEST_EQUALITY(basis_values.weighted_curl_basis_vector.extent(1),12);
    TEST_EQUALITY(basis_values.weighted_curl_basis_vector.extent(2),num_qp);
    TEST_EQUALITY(basis_values.weighted_curl_basis_vector.extent(3),3);

  }

  TEUCHOS_UNIT_TEST(basis_values, md_field_setup)
  {

    Teuchos::RCP<shards::CellTopology> topo =
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    const int num_cells = 20;
    const panzer::CellData cell_data(num_cells,topo);
    const std::string basis_type = "Q2";
    const int cubature_degree = 2;

    RCP<IntegrationRule> int_rule =
      rcp(new IntegrationRule(cubature_degree, cell_data));
    const unsigned int num_qp = Teuchos::as<unsigned int>(int_rule->num_points);

    RCP<PureBasis> basis = Teuchos::rcp(new PureBasis(basis_type,2,cell_data));
    RCP<panzer::BasisIRLayout> basisPtLayout = rcp(new panzer::BasisIRLayout(basis, *int_rule));

    panzer::BasisValues2<double> basis_values("prefix_",true,true);

    basis_values.setupArrays(basisPtLayout);

    // check to make sure all data layouts and field names are as
    // expected. In a simulation environment the field manager will
    // build these values.

    // check basis
    TEST_EQUALITY(basis_values.basis_ref_scalar.fieldTag().dataLayout().rank(),2);
    TEST_EQUALITY(basis_values.basis_ref_scalar.fieldTag().dataLayout().extent(0),9);
    TEST_EQUALITY(basis_values.basis_ref_scalar.fieldTag().dataLayout().extent(1),num_qp);
    TEST_EQUALITY(basis_values.basis_ref_scalar.fieldTag().name(),"prefix_basis_ref");

    TEST_EQUALITY(basis_values.basis_scalar.fieldTag().dataLayout().rank(),3);
    TEST_EQUALITY(basis_values.basis_scalar.fieldTag().dataLayout().extent_int(0),num_cells);
    TEST_EQUALITY(basis_values.basis_scalar.fieldTag().dataLayout().extent(1),9);
    TEST_EQUALITY(basis_values.basis_scalar.fieldTag().dataLayout().extent(2),num_qp);
    TEST_EQUALITY(basis_values.basis_scalar.fieldTag().name(),"prefix_basis");

    TEST_EQUALITY(basis_values.weighted_basis_scalar.fieldTag().dataLayout().rank(),3);
    TEST_EQUALITY(basis_values.weighted_basis_scalar.fieldTag().dataLayout().extent_int(0),num_cells);
    TEST_EQUALITY(basis_values.weighted_basis_scalar.fieldTag().dataLayout().extent(1),9);
    TEST_EQUALITY(basis_values.weighted_basis_scalar.fieldTag().dataLayout().extent(2),num_qp);
    TEST_EQUALITY(basis_values.weighted_basis_scalar.fieldTag().name(),"prefix_weighted_basis");

    // check gradients
    TEST_EQUALITY(basis_values.grad_basis_ref.fieldTag().dataLayout().rank(),3);
    TEST_EQUALITY(basis_values.grad_basis_ref.fieldTag().dataLayout().extent(0),9);
    TEST_EQUALITY(basis_values.grad_basis_ref.fieldTag().dataLayout().extent(1),num_qp);
    TEST_EQUALITY(basis_values.grad_basis_ref.fieldTag().dataLayout().extent(2),2);
    TEST_EQUALITY(basis_values.grad_basis_ref.fieldTag().name(),"prefix_grad_basis_ref");

    TEST_EQUALITY(basis_values.grad_basis.fieldTag().dataLayout().rank(),4);
    TEST_EQUALITY(basis_values.grad_basis.fieldTag().dataLayout().extent_int(0),num_cells);
    TEST_EQUALITY(basis_values.grad_basis.fieldTag().dataLayout().extent(1),9);
    TEST_EQUALITY(basis_values.grad_basis.fieldTag().dataLayout().extent(2),num_qp);
    TEST_EQUALITY(basis_values.grad_basis.fieldTag().dataLayout().extent(3),2);
    TEST_EQUALITY(basis_values.grad_basis.fieldTag().name(),"prefix_grad_basis");

    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().dataLayout().rank(),4);
    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().dataLayout().extent_int(0),num_cells);
    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().dataLayout().extent(1),9);
    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().dataLayout().extent(2),num_qp);
    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().dataLayout().extent(3),2);
    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().name(),"prefix_weighted_grad_basis");

    // check coordinates
    TEST_EQUALITY(basis_values.basis_coordinates_ref.fieldTag().dataLayout().rank(),2);
    TEST_EQUALITY(basis_values.basis_coordinates_ref.fieldTag().dataLayout().extent(0),9);
    TEST_EQUALITY(basis_values.basis_coordinates_ref.fieldTag().dataLayout().extent(1),2);
    TEST_EQUALITY(basis_values.basis_coordinates_ref.fieldTag().name(),"prefix_basis_coordinates_ref");

    TEST_EQUALITY(basis_values.basis_coordinates.fieldTag().dataLayout().rank(),3);
    TEST_EQUALITY(basis_values.basis_coordinates.fieldTag().dataLayout().extent_int(0),num_cells);
    TEST_EQUALITY(basis_values.basis_coordinates.fieldTag().dataLayout().extent(1),9);
    TEST_EQUALITY(basis_values.basis_coordinates.fieldTag().dataLayout().extent(2),2);
    TEST_EQUALITY(basis_values.basis_coordinates.fieldTag().name(),"prefix_basis_coordinates");
  }

  TEUCHOS_UNIT_TEST(basis_values, md_field_setup_fad_disabled)
  {
    Teuchos::RCP<shards::CellTopology> topo =
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    const int num_cells = 20;
    const panzer::CellData cell_data(num_cells,topo);
    const std::string basis_type = "Q2";
    const int cubature_degree = 2;

    RCP<IntegrationRule> int_rule =
      rcp(new IntegrationRule(cubature_degree, cell_data));
    const unsigned int num_qp = Teuchos::as<unsigned int>(int_rule->num_points);

    RCP<PureBasis> basis = Teuchos::rcp(new PureBasis(basis_type,2,cell_data));
    RCP<panzer::BasisIRLayout> basisPtLayout = rcp(new panzer::BasisIRLayout(basis, *int_rule));

    // NOTE: BasisValues2 with Fad scalar types failing at high order
    // in Intrepid2. See ticket #3340 on github for details. Changing
    // this test to double for now.
    //panzer::BasisValues2<ScalarType> basis_values("prefix_",true,true);
    panzer::BasisValues2<double> basis_values("prefix_",true,true);

    basis_values.setupArrays(basisPtLayout);

    // check to make sure all data layouts and field names are as
    // expected. In a simulation environment the field manager will
    // build these values.

    // check basis
    TEST_EQUALITY(basis_values.basis_ref_scalar.fieldTag().dataLayout().rank(),2);
    TEST_EQUALITY(basis_values.basis_ref_scalar.fieldTag().dataLayout().extent(0),9);
    TEST_EQUALITY(basis_values.basis_ref_scalar.fieldTag().dataLayout().extent(1),num_qp);
    TEST_EQUALITY(basis_values.basis_ref_scalar.fieldTag().name(),"prefix_basis_ref");

    TEST_EQUALITY(basis_values.basis_scalar.fieldTag().dataLayout().rank(),3);
    TEST_EQUALITY(basis_values.basis_scalar.fieldTag().dataLayout().extent_int(0),num_cells);
    TEST_EQUALITY(basis_values.basis_scalar.fieldTag().dataLayout().extent(1),9);
    TEST_EQUALITY(basis_values.basis_scalar.fieldTag().dataLayout().extent(2),num_qp);
    TEST_EQUALITY(basis_values.basis_scalar.fieldTag().name(),"prefix_basis");

    TEST_EQUALITY(basis_values.weighted_basis_scalar.fieldTag().dataLayout().rank(),3);
    TEST_EQUALITY(basis_values.weighted_basis_scalar.fieldTag().dataLayout().extent_int(0),num_cells);
    TEST_EQUALITY(basis_values.weighted_basis_scalar.fieldTag().dataLayout().extent(1),9);
    TEST_EQUALITY(basis_values.weighted_basis_scalar.fieldTag().dataLayout().extent(2),num_qp);
    TEST_EQUALITY(basis_values.weighted_basis_scalar.fieldTag().name(),"prefix_weighted_basis");

    // check gradients
    TEST_EQUALITY(basis_values.grad_basis_ref.fieldTag().dataLayout().rank(),3);
    TEST_EQUALITY(basis_values.grad_basis_ref.fieldTag().dataLayout().extent(0),9);
    TEST_EQUALITY(basis_values.grad_basis_ref.fieldTag().dataLayout().extent(1),num_qp);
    TEST_EQUALITY(basis_values.grad_basis_ref.fieldTag().dataLayout().extent(2),2);
    TEST_EQUALITY(basis_values.grad_basis_ref.fieldTag().name(),"prefix_grad_basis_ref");

    TEST_EQUALITY(basis_values.grad_basis.fieldTag().dataLayout().rank(),4);
    TEST_EQUALITY(basis_values.grad_basis.fieldTag().dataLayout().extent_int(0),num_cells);
    TEST_EQUALITY(basis_values.grad_basis.fieldTag().dataLayout().extent(1),9);
    TEST_EQUALITY(basis_values.grad_basis.fieldTag().dataLayout().extent(2),num_qp);
    TEST_EQUALITY(basis_values.grad_basis.fieldTag().dataLayout().extent(3),2);
    TEST_EQUALITY(basis_values.grad_basis.fieldTag().name(),"prefix_grad_basis");

    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().dataLayout().rank(),4);
    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().dataLayout().extent_int(0),num_cells);
    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().dataLayout().extent(1),9);
    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().dataLayout().extent(2),num_qp);
    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().dataLayout().extent(3),2);
    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().name(),"prefix_weighted_grad_basis");

    // check coordinates
    TEST_EQUALITY(basis_values.basis_coordinates_ref.fieldTag().dataLayout().rank(),2);
    TEST_EQUALITY(basis_values.basis_coordinates_ref.fieldTag().dataLayout().extent(0),9);
    TEST_EQUALITY(basis_values.basis_coordinates_ref.fieldTag().dataLayout().extent(1),2);
    TEST_EQUALITY(basis_values.basis_coordinates_ref.fieldTag().name(),"prefix_basis_coordinates_ref");

    TEST_EQUALITY(basis_values.basis_coordinates.fieldTag().dataLayout().rank(),3);
    TEST_EQUALITY(basis_values.basis_coordinates.fieldTag().dataLayout().extent_int(0),num_cells);
    TEST_EQUALITY(basis_values.basis_coordinates.fieldTag().dataLayout().extent(1),9);
    TEST_EQUALITY(basis_values.basis_coordinates.fieldTag().dataLayout().extent(2),2);
    TEST_EQUALITY(basis_values.basis_coordinates.fieldTag().name(),"prefix_basis_coordinates");
  }

  TEUCHOS_UNIT_TEST(basis_values, control_vol_hgrad)
  {
    typedef panzer::ArrayTraits<double,PHX::MDField<double> >::size_type size_type;

    Teuchos::RCP<shards::CellTopology> topo =
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    const int num_cells = 4;
    const int base_cell_dimension = 2;
    const panzer::CellData cell_data(num_cells,topo);

    std::string cv_type = "volume";
    RCP<IntegrationRule> int_rule_vol =
      rcp(new IntegrationRule(cell_data, cv_type));
    const unsigned int num_qp = Teuchos::as<unsigned int>(int_rule_vol->num_points);

    panzer::IntegrationValues2<double> int_values_vol("prefix_",true);
    int_values_vol.setupArrays(int_rule_vol);

    const int num_vertices = int_rule_vol->topology->getNodeCount();
    panzer::MDFieldArrayFactory af("prefix_",true);
    PHX::MDField<double,Cell,NODE,Dim> node_coordinates
        = af.buildStaticArray<double,Cell,NODE,Dim>("nc",num_cells, num_vertices, base_cell_dimension);

    // Set up node coordinates.  Here we assume the following
    // ordering.  This needs to be consistent with shards topology,
    // otherwise we will get negative determinates

    // 3(0,1)---2(1,1)
    //   |    0  |
    //   |       |
    // 0(0,0)---1(1,0)
    {
      const size_type x = 0;
      const size_type y = 1;
      auto node_coordinates_host = Kokkos::create_mirror_view(node_coordinates.get_view());
      for (size_type cell = 0; cell < node_coordinates.extent(0); ++cell) {
        int xleft = cell % 2;
        int yleft = int(cell/2);

        node_coordinates_host(cell,0,x) = xleft*0.5;
        node_coordinates_host(cell,0,y) = yleft*0.5;

        node_coordinates_host(cell,1,x) = (xleft+1)*0.5;
        node_coordinates_host(cell,1,y) = yleft*0.5;

        node_coordinates_host(cell,2,x) = (xleft+1)*0.5;
        node_coordinates_host(cell,2,y) = (yleft+1)*0.5;

        node_coordinates_host(cell,3,x) = xleft*0.5;
        node_coordinates_host(cell,3,y) = (yleft+1)*0.5;

        out << "Cell " << cell << " = ";
        for(int i=0;i<4;i++)
          out << "(" << node_coordinates_host(cell,i,x) << ", "
              << node_coordinates_host(cell,i,y) << ") ";
        out << std::endl;
      }
      Kokkos::deep_copy(node_coordinates.get_view(),node_coordinates_host);
    }

    int_values_vol.evaluateValues(node_coordinates);

    const std::string basis_type = "HGrad";

    Teuchos::RCP<PureBasis> basis = Teuchos::rcp(new PureBasis(basis_type,1,cell_data));
    RCP<panzer::BasisIRLayout> basisIRLayout = rcp(new panzer::BasisIRLayout(basis, *int_rule_vol));

    panzer::BasisValues2<double> basis_values("",true,true);

    basis_values.setupArrays(basisIRLayout);
    basis_values.evaluateValuesCV(int_values_vol.ref_ip_coordinates,
                                  int_values_vol.jac,
                                  int_values_vol.jac_det,
                                  int_values_vol.jac_inv);

    TEST_EQUALITY(basis_values.basis_scalar.extent_int(0),num_cells);
    TEST_EQUALITY(basis_values.basis_scalar.extent(1),4);
    TEST_EQUALITY(basis_values.basis_scalar.extent(2),num_qp);

    // check basis values, control volume integration points are defined on physical cells
    auto ref_ip_coordinates_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),int_values_vol.ref_ip_coordinates.get_view());
    auto basis_scalar_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),basis_values.basis_scalar.get_view());
    for(int cell=0;cell<num_cells;cell++) {

       for(unsigned int i=0;i<num_qp;i++) {
          double x = ref_ip_coordinates_host(cell,i,0);
          double y = ref_ip_coordinates_host(cell,i,1);

          // check values
          TEST_EQUALITY(basis_scalar_host(cell,0,i),0.25*(1.0-x)*(1.0-y));
          TEST_EQUALITY(basis_scalar_host(cell,1,i),0.25*(1.0+x)*(1.0-y));
          TEST_EQUALITY(basis_scalar_host(cell,2,i),0.25*(1.0+x)*(1.0+y));
          TEST_EQUALITY(basis_scalar_host(cell,3,i),0.25*(1.0-x)*(1.0+y));
       }
    }
  }
}
