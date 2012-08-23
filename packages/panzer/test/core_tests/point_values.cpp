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

#include "Panzer_Traits.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_ArrayTraits.hpp"
#include "Panzer_PointRule.hpp"
#include "Panzer_PointValues.hpp"
#include "Panzer_CommonArrayFactories.hpp"

#include "Intrepid_FieldContainer.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Intrepid::FieldContainer;

namespace panzer {
  TEUCHOS_UNIT_TEST(point_values, intrepid_container)
  {
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    const int num_cells = 4;
    const int base_cell_dimension = 2;
    const panzer::CellData cell_data(num_cells, base_cell_dimension,topo);
    int num_points = 3;

    RCP<PointRule> point_rule = rcp(new PointRule("RandomPoints",num_points, cell_data));

    TEST_EQUALITY(point_rule->num_points,num_points);
  
    panzer::PointValues<double,Intrepid::FieldContainer<double> > point_values;
    panzer::IntrepidFieldContainerFactory af;

    point_values.setupArrays(point_rule,af);

    // Set up node coordinates.  Here we assume the following
    // ordering.  This needs to be consistent with shards topology,
    // otherwise we will get negative determinates

    // 3(0,1)---2(1,1)
    //   |    0  |
    //   |       |
    // 0(0,0)---1(1,0)

    const int num_vertices = point_rule->topology->getNodeCount();
    Intrepid::FieldContainer<double> node_coordinates(num_cells, num_vertices,
	 				              base_cell_dimension);


    typedef panzer::ArrayTraits<double,FieldContainer<double> >::size_type size_type;
    const size_type x = 0;
    const size_type y = 1;
    for (size_type cell = 0; cell < node_coordinates.dimension(0); ++cell) {
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

    // Build the evaluation points

    Intrepid::FieldContainer<double> point_coordinates(num_points, base_cell_dimension);
    point_coordinates(0,0) =  0.0; point_coordinates(0,1) = 0.0; // mid point
    point_coordinates(1,0) =  0.5; point_coordinates(1,1) = 0.5; // mid point of upper left quadrant
    point_coordinates(2,0) = -0.5; point_coordinates(2,1) = 0.0; // mid point of line from center to left side
    
    point_values.evaluateValues(node_coordinates,point_coordinates);

    for(size_type p=0;p<num_points;p++)
       for(size_type d=0;d<base_cell_dimension;d++)
          TEST_EQUALITY(point_values.coords_ref(p,d),point_coordinates(p,d));

    for(size_type c=0;c<num_cells;c++) {
       double dx = 0.5;
       double dy = 0.5;
       for(size_type p=0;p<num_points;p++) {
          double x = dx*(point_coordinates(p,0)+1.0)/2.0 + node_coordinates(c,0,0); 
          double y = dy*(point_coordinates(p,1)+1.0)/2.0 + node_coordinates(c,0,1);
          TEST_FLOATING_EQUALITY(point_values.point_coords(c,p,0),x,1e-10);
          TEST_FLOATING_EQUALITY(point_values.point_coords(c,p,1),y,1e-10);
       }
    }

  }

  TEUCHOS_UNIT_TEST(point_values, intrepid_container_dfad)
  {
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    const int num_cells = 4;
    const int base_cell_dimension = 2;
    const panzer::CellData cell_data(num_cells, base_cell_dimension,topo);
    int num_points = 3;

    RCP<PointRule> point_rule = rcp(new PointRule("RandomPoints",num_points, cell_data));

    TEST_EQUALITY(point_rule->num_points,num_points);
  
    typedef panzer::Traits::FadType ScalarType;
    panzer::PointValues<ScalarType,Intrepid::FieldContainer<ScalarType> > point_values;
    panzer::IntrepidFieldContainerFactory af;

    point_values.setupArrays(point_rule,af);

    // Set up node coordinates.  Here we assume the following
    // ordering.  This needs to be consistent with shards topology,
    // otherwise we will get negative determinates

    // 3(0,1)---2(1,1)
    //   |    0  |
    //   |       |
    // 0(0,0)---1(1,0)

    const int num_vertices = point_rule->topology->getNodeCount();
    Intrepid::FieldContainer<ScalarType> node_coordinates(num_cells, num_vertices,
	 				              base_cell_dimension);

    typedef panzer::ArrayTraits<ScalarType,FieldContainer<ScalarType> >::size_type size_type;
    const size_type x = 0;
    const size_type y = 1;
    for (size_type cell = 0; cell < node_coordinates.dimension(0); ++cell) {
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

    // Build the evaluation points

    Intrepid::FieldContainer<ScalarType> point_coordinates(num_points, base_cell_dimension);
    point_coordinates(0,0) =  0.0; point_coordinates(0,1) = 0.0; // mid point
    point_coordinates(1,0) =  0.5; point_coordinates(1,1) = 0.5; // mid point of upper left quadrant
    point_coordinates(2,0) = -0.5; point_coordinates(2,1) = 0.0; // mid point of line from center to left side
    
    point_values.evaluateValues(node_coordinates,point_coordinates);
  }

  TEUCHOS_UNIT_TEST(point_values, md_field_setup)
  {
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    const int num_cells = 20;
    const int base_cell_dimension = 2;
    const panzer::CellData cell_data(num_cells, base_cell_dimension,topo);
    int num_points = 3;

    RCP<PointRule> point_rule = rcp(new PointRule("RandomPoints",num_points, cell_data));

    TEST_EQUALITY(point_rule->num_points,num_points);
  
    panzer::PointValues<double,PHX::MDField<double> > point_values;
    panzer::MDFieldArrayFactory af("prefix_");

    point_values.setupArrays(point_rule,af);

    // check to make sure all data layouts and field names are as 
    // expected. In a simulation environment the field manager will
    // build these values.

    // check basis
    TEST_EQUALITY(point_values.coords_ref.fieldTag().dataLayout().rank(),2);
    TEST_EQUALITY(point_values.coords_ref.fieldTag().dataLayout().dimension(0),num_points);
    TEST_EQUALITY(point_values.coords_ref.fieldTag().dataLayout().dimension(1),base_cell_dimension);
    TEST_EQUALITY(point_values.coords_ref.fieldTag().name(),"prefix_coords_ref");

    TEST_EQUALITY(point_values.node_coordinates.fieldTag().dataLayout().rank(),3);
    TEST_EQUALITY(point_values.node_coordinates.fieldTag().dataLayout().dimension(0),num_cells);
    TEST_EQUALITY(point_values.node_coordinates.fieldTag().dataLayout().dimension(1),4);
    TEST_EQUALITY(point_values.node_coordinates.fieldTag().dataLayout().dimension(2),base_cell_dimension);
    TEST_EQUALITY(point_values.node_coordinates.fieldTag().name(),"prefix_node_coordinates");

    TEST_EQUALITY(point_values.point_coords.fieldTag().dataLayout().rank(),3);
    TEST_EQUALITY(point_values.point_coords.fieldTag().dataLayout().dimension(0),num_cells);
    TEST_EQUALITY(point_values.point_coords.fieldTag().dataLayout().dimension(1),num_points);
    TEST_EQUALITY(point_values.point_coords.fieldTag().dataLayout().dimension(2),base_cell_dimension);
    TEST_EQUALITY(point_values.point_coords.fieldTag().name(),"prefix_point_coords");

    TEST_EQUALITY(point_values.jac.fieldTag().dataLayout().rank(),4);
    TEST_EQUALITY(point_values.jac.fieldTag().dataLayout().dimension(0),num_cells);
    TEST_EQUALITY(point_values.jac.fieldTag().dataLayout().dimension(1),num_points);
    TEST_EQUALITY(point_values.jac.fieldTag().dataLayout().dimension(2),base_cell_dimension);
    TEST_EQUALITY(point_values.jac.fieldTag().dataLayout().dimension(3),base_cell_dimension);
    TEST_EQUALITY(point_values.jac.fieldTag().name(),"prefix_jac");

    TEST_EQUALITY(point_values.jac_inv.fieldTag().dataLayout().rank(),4);
    TEST_EQUALITY(point_values.jac_inv.fieldTag().dataLayout().dimension(0),num_cells);
    TEST_EQUALITY(point_values.jac_inv.fieldTag().dataLayout().dimension(1),num_points);
    TEST_EQUALITY(point_values.jac_inv.fieldTag().dataLayout().dimension(2),base_cell_dimension);
    TEST_EQUALITY(point_values.jac_inv.fieldTag().dataLayout().dimension(3),base_cell_dimension);
    TEST_EQUALITY(point_values.jac_inv.fieldTag().name(),"prefix_jac_inv");

    TEST_EQUALITY(point_values.jac_det.fieldTag().dataLayout().rank(),2);
    TEST_EQUALITY(point_values.jac_det.fieldTag().dataLayout().dimension(0),num_cells);
    TEST_EQUALITY(point_values.jac_det.fieldTag().dataLayout().dimension(1),num_points);
    TEST_EQUALITY(point_values.jac_det.fieldTag().name(),"prefix_jac_det");
  }
}
