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

#include "Panzer_Traits.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_PointRule.hpp"
#include "Panzer_PointValues2.hpp"
#include "Panzer_CommonArrayFactories.hpp"

#include "Kokkos_DynRankView.hpp"

#include "Phalanx_KokkosViewFactory.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer {

  TEUCHOS_UNIT_TEST(point_values2, md_field_setup)
  {
    Teuchos::RCP<shards::CellTopology> topo =
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    const int num_cells = 20;
    const int base_cell_dimension = 2;
    const panzer::CellData cell_data(num_cells,topo);
    unsigned int num_points = 3;

    RCP<PointRule> point_rule = rcp(new PointRule("RandomPoints",num_points, cell_data));

    TEST_EQUALITY(point_rule->num_points,Teuchos::as<int>(num_points));

    panzer::PointValues2<double > point_values2("prefix_");
    point_values2.setupArrays(point_rule);

    // check to make sure all data layouts and field names are as
    // expected. In a simulation environment the field manager will
    // build these values.

    // check basis
    TEST_EQUALITY(point_values2.getRefCoordinates().fieldTag().dataLayout().rank(),2);
    TEST_EQUALITY(point_values2.getRefCoordinates().fieldTag().dataLayout().extent(0),num_points);
    TEST_EQUALITY(point_values2.getRefCoordinates().fieldTag().dataLayout().extent_int(1),base_cell_dimension);
    TEST_EQUALITY(point_values2.getRefCoordinates().fieldTag().name(),"prefix_coords_ref");

    TEST_EQUALITY(point_values2.getVertexCoordinates().fieldTag().dataLayout().rank(),3);
    TEST_EQUALITY(point_values2.getVertexCoordinates().fieldTag().dataLayout().extent_int(0),num_cells);
    TEST_EQUALITY(point_values2.getVertexCoordinates().fieldTag().dataLayout().extent(1),4);
    TEST_EQUALITY(point_values2.getVertexCoordinates().fieldTag().dataLayout().extent_int(2),base_cell_dimension);
    TEST_EQUALITY(point_values2.getVertexCoordinates().fieldTag().name(),"prefix_node_coordinates");

    TEST_EQUALITY(point_values2.point_coords.fieldTag().dataLayout().rank(),3);
    TEST_EQUALITY(point_values2.point_coords.fieldTag().dataLayout().extent_int(0),num_cells);
    TEST_EQUALITY(point_values2.point_coords.fieldTag().dataLayout().extent(1),num_points);
    TEST_EQUALITY(point_values2.point_coords.fieldTag().dataLayout().extent_int(2),base_cell_dimension);
    TEST_EQUALITY(point_values2.point_coords.fieldTag().name(),"prefix_point_coords");

    TEST_EQUALITY(point_values2.jac.fieldTag().dataLayout().rank(),4);
    TEST_EQUALITY(point_values2.jac.fieldTag().dataLayout().extent_int(0),num_cells);
    TEST_EQUALITY(point_values2.jac.fieldTag().dataLayout().extent(1),num_points);
    TEST_EQUALITY(point_values2.jac.fieldTag().dataLayout().extent_int(2),base_cell_dimension);
    TEST_EQUALITY(point_values2.jac.fieldTag().dataLayout().extent_int(3),base_cell_dimension);
    TEST_EQUALITY(point_values2.jac.fieldTag().name(),"prefix_jac");

    TEST_EQUALITY(point_values2.jac_inv.fieldTag().dataLayout().rank(),4);
    TEST_EQUALITY(point_values2.jac_inv.fieldTag().dataLayout().extent_int(0),num_cells);
    TEST_EQUALITY(point_values2.jac_inv.fieldTag().dataLayout().extent(1),num_points);
    TEST_EQUALITY(point_values2.jac_inv.fieldTag().dataLayout().extent_int(2),base_cell_dimension);
    TEST_EQUALITY(point_values2.jac_inv.fieldTag().dataLayout().extent_int(3),base_cell_dimension);
    TEST_EQUALITY(point_values2.jac_inv.fieldTag().name(),"prefix_jac_inv");

    TEST_EQUALITY(point_values2.jac_det.fieldTag().dataLayout().rank(),2);
    TEST_EQUALITY(point_values2.jac_det.fieldTag().dataLayout().extent_int(0),num_cells);
    TEST_EQUALITY(point_values2.jac_det.fieldTag().dataLayout().extent(1),num_points);
    TEST_EQUALITY(point_values2.jac_det.fieldTag().name(),"prefix_jac_det");
  }

  TEUCHOS_UNIT_TEST(point_values2, md_field_evaluate_double)
  {
    using ScalarType = double;
    using ArrayType = PHX::MDField<ScalarType>;
    using ViewFactory = PHX::KokkosViewFactory<ScalarType,typename PHX::DevLayout<ScalarType>::type,PHX::Device>;
    using size_type = PHX::MDField<double>::size_type;

    Teuchos::RCP<shards::CellTopology> topo =
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    const int num_cells = 4;
    const int base_cell_dimension = 2;
    const panzer::CellData cell_data(num_cells,topo);
    int num_points = 3;

    RCP<PointRule> point_rule = rcp(new PointRule("RandomPoints",num_points, cell_data));

    TEST_EQUALITY(point_rule->num_points,num_points);

    const size_type derivative_dim = 4;
    const std::vector<PHX::index_size_type> ddims(1,derivative_dim);
    panzer::MDFieldArrayFactory af("prefix_",ddims);

    panzer::PointValues2<ScalarType> point_values2("prefix_",ddims);
    point_values2.setupArrays(point_rule);

    // Set up node coordinates.  Here we assume the following
    // ordering.  This needs to be consistent with shards topology,
    // otherwise we will get negative determinates

    // 3(0,1)---2(1,1)
    //   |    0  |
    //   |       |
    // 0(0,0)---1(1,0)

    const int num_vertices = point_rule->topology->getNodeCount();
    ArrayType node_coordinates = af.buildArray<ScalarType,Cell,NODE,Dim>("node_coordinates",num_cells, num_vertices, base_cell_dimension);
    node_coordinates.setFieldData(ViewFactory::buildView(node_coordinates.fieldTag(),ddims));
    auto node_coordinates_host = Kokkos::create_mirror_view(node_coordinates.get_static_view());
    {
      const size_type x = 0;
      const size_type y = 1;
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
    }

    // Build the evaluation points

    ArrayType point_coordinates =  af.buildArray<ScalarType,IP,Dim>("points",num_points, base_cell_dimension);
    point_coordinates.setFieldData(ViewFactory::buildView(point_coordinates.fieldTag(),ddims));
    auto point_coordinates_host = Kokkos::create_mirror_view(point_coordinates.get_static_view());
    point_coordinates_host(0,0) =  0.0; point_coordinates_host(0,1) = 0.0; // mid point
    point_coordinates_host(1,0) =  0.5; point_coordinates_host(1,1) = 0.5; // mid point of upper left quadrant
    point_coordinates_host(2,0) = -0.5; point_coordinates_host(2,1) = 0.0; // mid point of line from center to left side
    Kokkos::deep_copy(point_coordinates.get_static_view(),point_coordinates_host);

    const_cast<PHX::MDField<ScalarType,IP,Dim>&>(point_values2.getRefCoordinates()).setFieldData(ViewFactory::buildView(point_values2.getRefCoordinates().fieldTag(),ddims));
    const_cast<PHX::MDField<ScalarType,Cell,NODE,Dim>&>(point_values2.getVertexCoordinates()).setFieldData(ViewFactory::buildView(point_values2.getVertexCoordinates().fieldTag(),ddims));

    point_values2.point_coords.setFieldData(ViewFactory::buildView(point_values2.point_coords.fieldTag(),ddims));
    point_values2.jac.setFieldData(ViewFactory::buildView(point_values2.jac.fieldTag(),ddims));
    point_values2.jac_inv.setFieldData(ViewFactory::buildView(point_values2.jac_inv.fieldTag(),ddims));
    point_values2.jac_det.setFieldData(ViewFactory::buildView(point_values2.jac_det.fieldTag(),ddims));

    point_values2.evaluateValues(node_coordinates,point_coordinates);

    // check the reference values (ensure copying)
    auto ref_coords_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),point_values2.getRefCoordinates().get_static_view());
    for(int p=0;p<num_points;p++)
       for(size_type d=0;d<static_cast<size_type>(base_cell_dimension);d++)
          TEST_EQUALITY(ref_coords_host(p,d),point_coordinates_host(p,d));

    // check the shifted values (ensure physical mapping)
    auto pv2_point_coords_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),point_values2.point_coords.get_static_view());
    for(int c=0;c<num_cells;c++) {
       double dx = 0.5;
       double dy = 0.5;
       for(int p=0;p<num_points;p++) {
          double x = dx*(point_coordinates_host(p,0)+1.0)/2.0 + node_coordinates_host(c,0,0);
          double y = dy*(point_coordinates_host(p,1)+1.0)/2.0 + node_coordinates_host(c,0,1);
          TEST_FLOATING_EQUALITY(pv2_point_coords_host(c,p,0),x,1e-10);
          TEST_FLOATING_EQUALITY(pv2_point_coords_host(c,p,1),y,1e-10);
       }
    }

    // check the jacobian
    auto pv2_jac_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),point_values2.jac.get_static_view());
    for(int c=0;c<num_cells;c++) {
       double dx = 0.5;
       double dy = 0.5;
       for(int p=0;p<num_points;p++) {
          TEST_FLOATING_EQUALITY(pv2_jac_host(c,p,0,0),dx/2.0,1e-10);
          TEST_FLOATING_EQUALITY(pv2_jac_host(c,p,0,1),0.0,1e-10);
          TEST_FLOATING_EQUALITY(pv2_jac_host(c,p,1,0),0.0,1e-10);
          TEST_FLOATING_EQUALITY(pv2_jac_host(c,p,1,1),dy/2.0,1e-10);
       }
    }
    auto pv2_jac_det_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),point_values2.jac_det.get_static_view());
    for(int c=0;c<num_cells;c++) {
       double dx = 0.5;
       double dy = 0.5;
       for(int p=0;p<num_points;p++) {
          TEST_FLOATING_EQUALITY(pv2_jac_det_host(c,p),dy*dx/4.0,1e-10);
       }
    }

    // check the inverse jacobian
    auto pv2_jac_inv_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),point_values2.jac_inv.get_static_view());
    for(int c=0;c<num_cells;c++) {
       double dx = 0.5;
       double dy = 0.5;
       for(int p=0;p<num_points;p++) {
          TEST_FLOATING_EQUALITY(pv2_jac_inv_host(c,p,0,0),2.0/dx,1e-10);
          TEST_FLOATING_EQUALITY(pv2_jac_inv_host(c,p,0,1),0.0,1e-10);
          TEST_FLOATING_EQUALITY(pv2_jac_inv_host(c,p,1,0),0.0,1e-10);
          TEST_FLOATING_EQUALITY(pv2_jac_inv_host(c,p,1,1),2.0/dy,1e-10);
       }
    }
  }

  // Disabled FAD support due to long build times on cuda (in debug
  // mode it takes multiple hours on some platforms). If we need
  // sensitivities wrt coordinates, we can reenable.
  /*
  TEUCHOS_UNIT_TEST(point_values2, md_field_evaluate)
  {
    typedef panzer::Traits::FadType ScalarType;
    typedef PHX::MDField<ScalarType> ArrayType;
    typedef PHX::KokkosViewFactory<ScalarType,typename PHX::DevLayout<ScalarType>::type,PHX::Device> ViewFactory;
    typedef PHX::MDField<double>::size_type size_type;


    Teuchos::RCP<shards::CellTopology> topo =
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    const int num_cells = 4;
    const int base_cell_dimension = 2;
    const panzer::CellData cell_data(num_cells,topo);
    int num_points = 3;

    RCP<PointRule> point_rule = rcp(new PointRule("RandomPoints",num_points, cell_data));

    TEST_EQUALITY(point_rule->num_points,num_points);

    const size_type derivative_dim = 4;
    const std::vector<PHX::index_size_type> ddims(1,derivative_dim);
    panzer::MDFieldArrayFactory af("prefix_",ddims);

    panzer::PointValues2<ScalarType> point_values2("prefix_",ddims);
    point_values2.setupArrays(point_rule);

    // Set up node coordinates.  Here we assume the following
    // ordering.  This needs to be consistent with shards topology,
    // otherwise we will get negative determinates

    // 3(0,1)---2(1,1)
    //   |    0  |
    //   |       |
    // 0(0,0)---1(1,0)

    const int num_vertices = point_rule->topology->getNodeCount();
    ArrayType node_coordinates = af.buildArray<ScalarType,Cell,NODE,Dim>("node_coordinates",num_cells, num_vertices, base_cell_dimension);
    node_coordinates.setFieldData(ViewFactory::buildView(node_coordinates.fieldTag(),ddims));
    {
      const size_type x = 0;
      const size_type y = 1;
      for (size_type cell = 0; cell < node_coordinates.extent(0); ++cell) {
        int xleft = cell % 2;
        int yleft = int(cell/2);

        node_coordinates(cell,0,x) = xleft*0.5;
        node_coordinates(cell,0,y) = yleft*0.5;

        node_coordinates(cell,1,x) = (xleft+1)*0.5;
        node_coordinates(cell,1,y) = yleft*0.5;

        node_coordinates(cell,2,x) = (xleft+1)*0.5;
        node_coordinates(cell,2,y) = (yleft+1)*0.5;

        node_coordinates(cell,3,x) = xleft*0.5;
        node_coordinates(cell,3,y) = (yleft+1)*0.5;

        out << "Cell " << cell << " = ";
        for(int i=0;i<4;i++)
          out << "(" << node_coordinates(cell,i,x) << ", "
              << node_coordinates(cell,i,y) << ") ";
        out << std::endl;
      }
    }

    // Build the evaluation points

    ArrayType point_coordinates =  af.buildArray<ScalarType,IP,Dim>("points",num_points, base_cell_dimension);
    point_coordinates.setFieldData(ViewFactory::buildView(point_coordinates.fieldTag(),ddims));
    point_coordinates(0,0) =  0.0; point_coordinates(0,1) = 0.0; // mid point
    point_coordinates(1,0) =  0.5; point_coordinates(1,1) = 0.5; // mid point of upper left quadrant
    point_coordinates(2,0) = -0.5; point_coordinates(2,1) = 0.0; // mid point of line from center to left side

    point_values2.getRefCoordinates().setFieldData(ViewFactory::buildView(point_values2.getRefCoordinates().fieldTag(),ddims));
    point_values2.getVertexCoordinates().setFieldData(ViewFactory::buildView(point_values2.getVertexCoordinates().fieldTag(),ddims));

    point_values2.point_coords.setFieldData(ViewFactory::buildView(point_values2.point_coords.fieldTag(),ddims));
    point_values2.jac.setFieldData(ViewFactory::buildView(point_values2.jac.fieldTag(),ddims));
    point_values2.jac_inv.setFieldData(ViewFactory::buildView(point_values2.jac_inv.fieldTag(),ddims));
    point_values2.jac_det.setFieldData(ViewFactory::buildView(point_values2.jac_det.fieldTag(),ddims));

    point_values2.evaluateValues(node_coordinates,point_coordinates);

    // check the reference values (ensure copying)
    for(int p=0;p<num_points;p++)
       for(size_type d=0;d<static_cast<size_type>(base_cell_dimension);d++)
          TEST_EQUALITY(point_values2.getRefCoordinates()(p,d).val(),point_coordinates(p,d).val());

    // check the shifted values (ensure physical mapping)
    for(int c=0;c<num_cells;c++) {
       double dx = 0.5;
       double dy = 0.5;
       for(int p=0;p<num_points;p++) {
          double x = dx*(point_coordinates(p,0).val()+1.0)/2.0 + node_coordinates(c,0,0).val();
          double y = dy*(point_coordinates(p,1).val()+1.0)/2.0 + node_coordinates(c,0,1).val();
          TEST_FLOATING_EQUALITY(point_values2.point_coords(c,p,0).val(),x,1e-10);
          TEST_FLOATING_EQUALITY(point_values2.point_coords(c,p,1).val(),y,1e-10);
       }
    }

    // check the jacobian
    for(int c=0;c<num_cells;c++) {
       double dx = 0.5;
       double dy = 0.5;
       for(int p=0;p<num_points;p++) {
          TEST_FLOATING_EQUALITY(point_values2.jac(c,p,0,0).val(),dx/2.0,1e-10);
          TEST_FLOATING_EQUALITY(point_values2.jac(c,p,0,1).val(),0.0,1e-10);
          TEST_FLOATING_EQUALITY(point_values2.jac(c,p,1,0).val(),0.0,1e-10);
          TEST_FLOATING_EQUALITY(point_values2.jac(c,p,1,1).val(),dy/2.0,1e-10);
       }
    }
    for(int c=0;c<num_cells;c++) {
       double dx = 0.5;
       double dy = 0.5;
       for(int p=0;p<num_points;p++) {
          TEST_FLOATING_EQUALITY(point_values2.jac_det(c,p).val(),dy*dx/4.0,1e-10);
       }
    }

    // check the inverse jacobian
    for(int c=0;c<num_cells;c++) {
       double dx = 0.5;
       double dy = 0.5;
       for(int p=0;p<num_points;p++) {
          TEST_FLOATING_EQUALITY(point_values2.jac_inv(c,p,0,0).val(),2.0/dx,1e-10);
          TEST_FLOATING_EQUALITY(point_values2.jac_inv(c,p,0,1).val(),0.0,1e-10);
          TEST_FLOATING_EQUALITY(point_values2.jac_inv(c,p,1,0).val(),0.0,1e-10);
          TEST_FLOATING_EQUALITY(point_values2.jac_inv(c,p,1,1).val(),2.0/dy,1e-10);
       }
    }
  }
  */

}
