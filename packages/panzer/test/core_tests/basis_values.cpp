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
#include "Panzer_IntegrationValues.hpp"
#include "Panzer_ArrayTraits.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Panzer_BasisValues.hpp"
#include "Panzer_CommonArrayFactories.hpp"
#include "Panzer_Traits.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using panzer::IntegrationRule;
using Intrepid::FieldContainer;

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
    
    panzer::IntegrationValues<double,Intrepid::FieldContainer<double> > 
      int_values;

    int_values.setupArrays(int_rule);

    const int num_vertices = int_rule->topology->getNodeCount();
    FieldContainer<double> node_coordinates(num_cells, num_vertices,
					    base_cell_dimension);



    // Set up node coordinates.  Here we assume the following
    // ordering.  This needs to be consistent with shards topology,
    // otherwise we will get negative determinates

    // 3(0,1)---2(1,1)
    //   |    0  |
    //   |       |
    // 0(0,0)---1(1,0)

    typedef panzer::ArrayTraits<double,FieldContainer<double> >::size_type size_type;
    const size_type x = 0;
    const size_type y = 1;
    for (size_type cell = 0; cell < node_coordinates.dimension(0); ++cell) {
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
    
    const std::string basis_type = "Q2";
  
    RCP<panzer::BasisIRLayout> basis = rcp(new panzer::BasisIRLayout(basis_type, 0, *int_rule));

    panzer::BasisValues<double,Intrepid::FieldContainer<double> > basis_values;
    panzer::IntrepidFieldContainerFactory af;

    basis_values.setupArrays(basis,af);
    
    basis_values.evaluateValues(int_values.cub_points,
				int_values.jac,
				int_values.jac_det,
				int_values.jac_inv,
				int_values.weighted_measure,
				node_coordinates);

    TEST_FLOATING_EQUALITY(basis_values.basis_coordinates(0,0,x),
			   0.0, 1.0e-8);
    TEST_FLOATING_EQUALITY(basis_values.basis_coordinates(0,0,y),
			   0.0, 1.0e-8);

    TEST_FLOATING_EQUALITY(basis_values.basis_coordinates(0,1,x),
			   1.0, 1.0e-8);
    TEST_FLOATING_EQUALITY(basis_values.basis_coordinates(0,1,y),
			   0.0, 1.0e-8);

    TEST_FLOATING_EQUALITY(basis_values.basis_coordinates(0,2,x),
			   1.0, 1.0e-8);
    TEST_FLOATING_EQUALITY(basis_values.basis_coordinates(0,2,y),
			   1.0, 1.0e-8);

    TEST_FLOATING_EQUALITY(basis_values.basis_coordinates(0,3,x),
			   0.0, 1.0e-8);
    TEST_FLOATING_EQUALITY(basis_values.basis_coordinates(0,3,y),
			   1.0, 1.0e-8);

  }

  TEUCHOS_UNIT_TEST(basis_values, grad_quad)
  {
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    const int num_cells = 4;
    const int base_cell_dimension = 2;
    const panzer::CellData cell_data(num_cells,topo);

    const int cubature_degree = 4;    
    RCP<IntegrationRule> int_rule = 
      rcp(new IntegrationRule(cubature_degree, cell_data));
    const int num_qp = int_rule->num_points;
    
    panzer::IntegrationValues<double,Intrepid::FieldContainer<double> > 
      int_values;

    int_values.setupArrays(int_rule);

    const int num_vertices = int_rule->topology->getNodeCount();
    FieldContainer<double> node_coordinates(num_cells, num_vertices,
					    base_cell_dimension);



    // Set up node coordinates.  Here we assume the following
    // ordering.  This needs to be consistent with shards topology,
    // otherwise we will get negative determinates

    // 3(0,1)---2(1,1)
    //   |    0  |
    //   |       |
    // 0(0,0)---1(1,0)

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

    int_values.evaluateValues(node_coordinates);
    
    const std::string basis_type = "Q1";
  
    RCP<panzer::BasisIRLayout> basis = rcp(new panzer::BasisIRLayout(basis_type, 0, *int_rule));

    panzer::BasisValues<double,Intrepid::FieldContainer<double> > basis_values;
    panzer::IntrepidFieldContainerFactory af;

    basis_values.setupArrays(basis,af);
    
    basis_values.evaluateValues(int_values.cub_points,
				int_values.jac,
				int_values.jac_det,
				int_values.jac_inv,
				int_values.weighted_measure,
				node_coordinates);

    double relCellVol = 0.25*0.25; // this is the relative (to the reference cell) volume
    for(int i=0;i<num_qp;i++) {
       double x = int_values.cub_points(i,0);
       double y = int_values.cub_points(i,1);
       double weight = int_values.cub_weights(i);

       // check reference values
       TEST_EQUALITY(basis_values.basis_ref(0,i),0.25*(x-1.0)*(y-1.0));
       TEST_EQUALITY(basis_values.grad_basis_ref(0,i,0),0.25*(y-1.0));
       TEST_EQUALITY(basis_values.grad_basis_ref(0,i,1),0.25*(x-1.0));

       // check basis values
       for(int cell=0;cell<num_cells;cell++) {

          TEST_EQUALITY(int_values.jac_det(cell,i),relCellVol);

          // check out basis on transformed elemented
          TEST_EQUALITY(basis_values.basis_ref(0,i),basis_values.basis(cell,0,i));
          TEST_EQUALITY(basis_values.weighted_basis(cell,0,i),relCellVol*weight*basis_values.basis(cell,0,i));

          TEST_EQUALITY(basis_values.grad_basis(cell,0,i,0),4.0*basis_values.grad_basis_ref(0,i,0));
          TEST_EQUALITY(basis_values.grad_basis(cell,0,i,1),4.0*basis_values.grad_basis_ref(0,i,1));

          TEST_EQUALITY(basis_values.weighted_grad_basis(cell,0,i,0),relCellVol*weight*basis_values.grad_basis(cell,0,i,0));
          TEST_EQUALITY(basis_values.weighted_grad_basis(cell,0,i,1),relCellVol*weight*basis_values.grad_basis(cell,0,i,1));
       }
    }
  }

  TEUCHOS_UNIT_TEST(basis_values, hcurl_basis)
  {
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    const int num_cells = 4;
    const int base_cell_dimension = 2;
    const panzer::CellData cell_data(num_cells,topo);

    const int cubature_degree = 4;    
    RCP<IntegrationRule> int_rule = 
      rcp(new IntegrationRule(cubature_degree, cell_data));
    const int num_qp = int_rule->num_points;
    
    panzer::IntegrationValues<double,Intrepid::FieldContainer<double> > int_values;
    int_values.setupArrays(int_rule);

    const int num_vertices = int_rule->topology->getNodeCount();
    FieldContainer<double> node_coordinates(num_cells, num_vertices,
					    base_cell_dimension);
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

      // edge_orientation(cell,0) =  1.0;
      // edge_orientation(cell,1) = -1.0;
      // edge_orientation(cell,2) =  1.0;
      // edge_orientation(cell,3) = -1.0;

      out << "Cell " << cell << " = ";
      for(int i=0;i<4;i++)
         out << "(" << node_coordinates(cell,i,x) << ", "
                    << node_coordinates(cell,i,y) << ") ";
      out << std::endl;
    }

    int_values.evaluateValues(node_coordinates);
    
    const std::string basis_type = "QEdge1";
  
    Teuchos::RCP<PureBasis> basis = Teuchos::rcp(new PureBasis(basis_type,1,cell_data));
    RCP<panzer::BasisIRLayout> basisIRLayout = rcp(new panzer::BasisIRLayout(basis, *int_rule));

    panzer::BasisValues<double,Intrepid::FieldContainer<double> > basis_values;
    panzer::IntrepidFieldContainerFactory af;

    basis_values.setupArrays(basisIRLayout,af);
    basis_values.evaluateValues(int_values.cub_points,
				int_values.jac,
				int_values.jac_det,
				int_values.jac_inv,
				int_values.weighted_measure,
				node_coordinates);

    TEST_EQUALITY(basis_values.basis_ref.dimension(0),4);
    TEST_EQUALITY(basis_values.basis_ref.dimension(1),num_qp);
    TEST_EQUALITY(basis_values.basis_ref.dimension(2),2);
    TEST_EQUALITY(basis_values.weighted_basis.dimension(0),num_cells);
    TEST_EQUALITY(basis_values.weighted_basis.dimension(1),4);
    TEST_EQUALITY(basis_values.weighted_basis.dimension(2),num_qp);
    TEST_EQUALITY(basis_values.weighted_basis.dimension(3),2);

    TEST_EQUALITY(basis_values.grad_basis_ref.size(),0);
    TEST_EQUALITY(basis_values.weighted_grad_basis.size(),0);

    TEST_EQUALITY(basis_values.curl_basis_ref.dimension(0),4);
    TEST_EQUALITY(basis_values.curl_basis_ref.dimension(1),num_qp);
    TEST_EQUALITY(basis_values.weighted_curl_basis.dimension(0),num_cells);
    TEST_EQUALITY(basis_values.weighted_curl_basis.dimension(1),4);
    TEST_EQUALITY(basis_values.weighted_curl_basis.dimension(2),num_qp);

    double relCellVol = 0.25*0.25; // this is the relative (to the reference cell) volume
    for(int i=0;i<num_qp;i++) {
       // double x = int_values.cub_points(i,0);
       double y = int_values.cub_points(i,1);
       double weight = int_values.cub_weights(i);

       // check reference values
       TEST_EQUALITY(basis_values.basis_ref(0,i,0),-0.25*(y-1.0));
       TEST_EQUALITY(basis_values.basis_ref(0,i,1),0.0);
       TEST_EQUALITY(basis_values.curl_basis_ref(0,i),0.25);

       // check basis values
       for(int cell=0;cell<num_cells;cell++) {

          TEST_EQUALITY(int_values.jac_det(cell,i),relCellVol);

          // check out basis on transformed elemented
          TEST_EQUALITY(basis_values.basis_ref(0,i,0),0.25*basis_values.basis(cell,0,i,0));
          TEST_EQUALITY(basis_values.basis_ref(0,i,1),0.25*basis_values.basis(cell,0,i,1));
          TEST_EQUALITY(basis_values.curl_basis_ref(0,i),relCellVol*basis_values.curl_basis(cell,0,i));

          TEST_EQUALITY(basis_values.basis_ref(1,i,0),0.25*basis_values.basis(cell,1,i,0));
          TEST_EQUALITY(basis_values.basis_ref(1,i,1),0.25*basis_values.basis(cell,1,i,1));
          TEST_EQUALITY(basis_values.curl_basis_ref(1,i),relCellVol*basis_values.curl_basis(cell,1,i));

          TEST_EQUALITY(basis_values.weighted_basis(cell,0,i,0),relCellVol*weight*basis_values.basis(cell,0,i,0));
          TEST_EQUALITY(basis_values.weighted_basis(cell,0,i,1),relCellVol*weight*basis_values.basis(cell,0,i,1));

          TEST_EQUALITY(basis_values.weighted_basis(cell,1,i,0),relCellVol*weight*basis_values.basis(cell,1,i,0));
          TEST_EQUALITY(basis_values.weighted_basis(cell,1,i,1),relCellVol*weight*basis_values.basis(cell,1,i,1));

          TEST_EQUALITY(basis_values.weighted_curl_basis(cell,0,i),relCellVol*weight*basis_values.curl_basis(cell,0,i));
          TEST_EQUALITY(basis_values.weighted_curl_basis(cell,1,i),relCellVol*weight*basis_values.curl_basis(cell,1,i));
       }
    }
  }

  TEUCHOS_UNIT_TEST(basis_values, hcurl_basis_3d)
  {
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<8> >()));

    const int num_cells = 4;
    const int base_cell_dimension = 3;
    const panzer::CellData cell_data(num_cells,topo);

    const int cubature_degree = 4;    
    RCP<IntegrationRule> int_rule = 
      rcp(new IntegrationRule(cubature_degree, cell_data));
    const int num_qp = int_rule->num_points;
    
    panzer::IntegrationValues<double,Intrepid::FieldContainer<double> > int_values;
    int_values.setupArrays(int_rule);

    const int num_vertices = int_rule->topology->getNodeCount();
    FieldContainer<double> node_coordinates(num_cells, num_vertices,
					    base_cell_dimension);
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

    typedef panzer::ArrayTraits<double,FieldContainer<double> >::size_type size_type;
    const size_type x = 0;
    const size_type y = 1;
    const size_type z = 2;
    for (size_type cell = 0; cell < node_coordinates.dimension(0); ++cell) {
      int znum = cell % 4;
      int xleft = znum % 2;
      int yleft = int(znum/2);
      int zleft = int(cell/4);

      node_coordinates(cell,0,x) = xleft*0.5;
      node_coordinates(cell,0,y) = yleft*0.5;
      node_coordinates(cell,0,z) = zleft*0.5;

      node_coordinates(cell,1,x) = (xleft+1)*0.5;
      node_coordinates(cell,1,y) = yleft*0.5; 
      node_coordinates(cell,1,z) = zleft*0.5;

      node_coordinates(cell,2,x) = (xleft+1)*0.5;
      node_coordinates(cell,2,y) = (yleft+1)*0.5;
      node_coordinates(cell,2,z) = zleft*0.5;

      node_coordinates(cell,3,x) = xleft*0.5;
      node_coordinates(cell,3,y) = (yleft+1)*0.5;
      node_coordinates(cell,3,z) = zleft*0.5;

      node_coordinates(cell,4,x) = xleft*0.5;
      node_coordinates(cell,4,y) = yleft*0.5;
      node_coordinates(cell,4,z) = (zleft+1)*0.5;

      node_coordinates(cell,5,x) = (xleft+1)*0.5;
      node_coordinates(cell,5,y) = yleft*0.5; 
      node_coordinates(cell,5,z) = (zleft+1)*0.5;

      node_coordinates(cell,6,x) = (xleft+1)*0.5;
      node_coordinates(cell,6,y) = (yleft+1)*0.5;
      node_coordinates(cell,6,z) = (zleft+1)*0.5;

      node_coordinates(cell,7,x) = xleft*0.5;
      node_coordinates(cell,7,y) = (yleft+1)*0.5;
      node_coordinates(cell,7,z) = (zleft+1)*0.5;

      out << "Cell " << cell << " = ";
      for(int i=0;i<8;i++)
         out << "(" << node_coordinates(cell,i,x) << ", "
                    << node_coordinates(cell,i,y) << ", "
                    << node_coordinates(cell,i,z) << ") ";
      out << std::endl;
    }

    int_values.evaluateValues(node_coordinates);
    
    const std::string basis_type = "QEdge1";
  
    Teuchos::RCP<PureBasis> basis = Teuchos::rcp(new PureBasis(basis_type,1,cell_data));
    RCP<panzer::BasisIRLayout> basisIRLayout = rcp(new panzer::BasisIRLayout(basis, *int_rule));

    panzer::BasisValues<double,Intrepid::FieldContainer<double> > basis_values;
    panzer::IntrepidFieldContainerFactory af;

    basis_values.setupArrays(basisIRLayout,af);
    basis_values.evaluateValues(int_values.cub_points,
				int_values.jac,
				int_values.jac_det,
				int_values.jac_inv,
				int_values.weighted_measure,
				node_coordinates);

    TEST_EQUALITY(basis_values.basis_ref.dimension(0),12);
    TEST_EQUALITY(basis_values.basis_ref.dimension(1),num_qp);
    TEST_EQUALITY(basis_values.basis_ref.dimension(2),3);
    TEST_EQUALITY(basis_values.weighted_basis.dimension(0),num_cells);
    TEST_EQUALITY(basis_values.weighted_basis.dimension(1),12);
    TEST_EQUALITY(basis_values.weighted_basis.dimension(2),num_qp);
    TEST_EQUALITY(basis_values.weighted_basis.dimension(3),3);

    TEST_EQUALITY(basis_values.grad_basis_ref.size(),0);
    TEST_EQUALITY(basis_values.weighted_grad_basis.size(),0);

    TEST_EQUALITY(basis_values.curl_basis_ref.dimension(0),12);
    TEST_EQUALITY(basis_values.curl_basis_ref.dimension(1),num_qp);
    TEST_EQUALITY(basis_values.curl_basis_ref.dimension(2),3);
    TEST_EQUALITY(basis_values.weighted_curl_basis.dimension(0),num_cells);
    TEST_EQUALITY(basis_values.weighted_curl_basis.dimension(1),12);
    TEST_EQUALITY(basis_values.weighted_curl_basis.dimension(2),num_qp);
    TEST_EQUALITY(basis_values.weighted_curl_basis.dimension(3),3);
  }

  TEUCHOS_UNIT_TEST(basis_values, md_field_setup)
  {
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    const int num_cells = 20;
    const int base_cell_dimension = 2;
    const panzer::CellData cell_data(num_cells,topo);
    const std::string basis_type = "Q2";
    const int cubature_degree = 2;    

    RCP<IntegrationRule> int_rule = 
      rcp(new IntegrationRule(cubature_degree, cell_data));
    const int num_qp = int_rule->num_points;
  
    RCP<PureBasis> basis = Teuchos::rcp(new PureBasis(basis_type,2,cell_data));
    RCP<panzer::BasisIRLayout> basisPtLayout = rcp(new panzer::BasisIRLayout(basis, *int_rule));

    panzer::BasisValues<double,PHX::MDField<double> > basis_values;
    panzer::MDFieldArrayFactory af("prefix_");

    basis_values.setupArrays(basisPtLayout,af);

    // check to make sure all data layouts and field names are as 
    // expected. In a simulation environment the field manager will
    // build these values.

    // check basis
    TEST_EQUALITY(basis_values.basis_ref.fieldTag().dataLayout().rank(),2);
    TEST_EQUALITY(basis_values.basis_ref.fieldTag().dataLayout().dimension(0),9);
    TEST_EQUALITY(basis_values.basis_ref.fieldTag().dataLayout().dimension(1),num_qp);
    TEST_EQUALITY(basis_values.basis_ref.fieldTag().name(),"prefix_basis_ref");

    TEST_EQUALITY(basis_values.basis.fieldTag().dataLayout().rank(),3);
    TEST_EQUALITY(basis_values.basis.fieldTag().dataLayout().dimension(0),num_cells);
    TEST_EQUALITY(basis_values.basis.fieldTag().dataLayout().dimension(1),9);
    TEST_EQUALITY(basis_values.basis.fieldTag().dataLayout().dimension(2),num_qp);
    TEST_EQUALITY(basis_values.basis.fieldTag().name(),"prefix_basis");

    TEST_EQUALITY(basis_values.weighted_basis.fieldTag().dataLayout().rank(),3);
    TEST_EQUALITY(basis_values.weighted_basis.fieldTag().dataLayout().dimension(0),num_cells);
    TEST_EQUALITY(basis_values.weighted_basis.fieldTag().dataLayout().dimension(1),9);
    TEST_EQUALITY(basis_values.weighted_basis.fieldTag().dataLayout().dimension(2),num_qp);
    TEST_EQUALITY(basis_values.weighted_basis.fieldTag().name(),"prefix_weighted_basis");

    // check gradients
    TEST_EQUALITY(basis_values.grad_basis_ref.fieldTag().dataLayout().rank(),3);
    TEST_EQUALITY(basis_values.grad_basis_ref.fieldTag().dataLayout().dimension(0),9);
    TEST_EQUALITY(basis_values.grad_basis_ref.fieldTag().dataLayout().dimension(1),num_qp);
    TEST_EQUALITY(basis_values.grad_basis_ref.fieldTag().dataLayout().dimension(2),2);
    TEST_EQUALITY(basis_values.grad_basis_ref.fieldTag().name(),"prefix_grad_basis_ref");

    TEST_EQUALITY(basis_values.grad_basis.fieldTag().dataLayout().rank(),4);
    TEST_EQUALITY(basis_values.grad_basis.fieldTag().dataLayout().dimension(0),num_cells);
    TEST_EQUALITY(basis_values.grad_basis.fieldTag().dataLayout().dimension(1),9);
    TEST_EQUALITY(basis_values.grad_basis.fieldTag().dataLayout().dimension(2),num_qp);
    TEST_EQUALITY(basis_values.grad_basis.fieldTag().dataLayout().dimension(3),2);
    TEST_EQUALITY(basis_values.grad_basis.fieldTag().name(),"prefix_grad_basis");

    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().dataLayout().rank(),4);
    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().dataLayout().dimension(0),num_cells);
    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().dataLayout().dimension(1),9);
    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().dataLayout().dimension(2),num_qp);
    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().dataLayout().dimension(3),2);
    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().name(),"prefix_weighted_grad_basis");

    // check coordinates
    TEST_EQUALITY(basis_values.basis_coordinates_ref.fieldTag().dataLayout().rank(),2);
    TEST_EQUALITY(basis_values.basis_coordinates_ref.fieldTag().dataLayout().dimension(0),9);
    TEST_EQUALITY(basis_values.basis_coordinates_ref.fieldTag().dataLayout().dimension(1),2);
    TEST_EQUALITY(basis_values.basis_coordinates_ref.fieldTag().name(),"prefix_basis_coordinates_ref");

    TEST_EQUALITY(basis_values.basis_coordinates.fieldTag().dataLayout().rank(),3);
    TEST_EQUALITY(basis_values.basis_coordinates.fieldTag().dataLayout().dimension(0),num_cells);
    TEST_EQUALITY(basis_values.basis_coordinates.fieldTag().dataLayout().dimension(1),9);
    TEST_EQUALITY(basis_values.basis_coordinates.fieldTag().dataLayout().dimension(2),2);
    TEST_EQUALITY(basis_values.basis_coordinates.fieldTag().name(),"prefix_basis_coordinates");
  }

  TEUCHOS_UNIT_TEST(basis_values, md_field_setup_fad)
  {
    typedef panzer::Traits::FadType ScalarType;

    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

    const int num_cells = 20;
    const int base_cell_dimension = 2;
    const panzer::CellData cell_data(num_cells,topo);
    const std::string basis_type = "Q2";
    const int cubature_degree = 2;    

    RCP<IntegrationRule> int_rule = 
      rcp(new IntegrationRule(cubature_degree, cell_data));
    const int num_qp = int_rule->num_points;
  
    RCP<PureBasis> basis = Teuchos::rcp(new PureBasis(basis_type,2,cell_data));
    RCP<panzer::BasisIRLayout> basisPtLayout = rcp(new panzer::BasisIRLayout(basis, *int_rule));

    panzer::BasisValues<ScalarType,PHX::MDField<ScalarType> > basis_values;
    panzer::MDFieldArrayFactory af("prefix_");

    basis_values.setupArrays(basisPtLayout,af);

    // check to make sure all data layouts and field names are as 
    // expected. In a simulation environment the field manager will
    // build these values.

    // check basis
    TEST_EQUALITY(basis_values.basis_ref.fieldTag().dataLayout().rank(),2);
    TEST_EQUALITY(basis_values.basis_ref.fieldTag().dataLayout().dimension(0),9);
    TEST_EQUALITY(basis_values.basis_ref.fieldTag().dataLayout().dimension(1),num_qp);
    TEST_EQUALITY(basis_values.basis_ref.fieldTag().name(),"prefix_basis_ref");

    TEST_EQUALITY(basis_values.basis.fieldTag().dataLayout().rank(),3);
    TEST_EQUALITY(basis_values.basis.fieldTag().dataLayout().dimension(0),num_cells);
    TEST_EQUALITY(basis_values.basis.fieldTag().dataLayout().dimension(1),9);
    TEST_EQUALITY(basis_values.basis.fieldTag().dataLayout().dimension(2),num_qp);
    TEST_EQUALITY(basis_values.basis.fieldTag().name(),"prefix_basis");

    TEST_EQUALITY(basis_values.weighted_basis.fieldTag().dataLayout().rank(),3);
    TEST_EQUALITY(basis_values.weighted_basis.fieldTag().dataLayout().dimension(0),num_cells);
    TEST_EQUALITY(basis_values.weighted_basis.fieldTag().dataLayout().dimension(1),9);
    TEST_EQUALITY(basis_values.weighted_basis.fieldTag().dataLayout().dimension(2),num_qp);
    TEST_EQUALITY(basis_values.weighted_basis.fieldTag().name(),"prefix_weighted_basis");

    // check gradients
    TEST_EQUALITY(basis_values.grad_basis_ref.fieldTag().dataLayout().rank(),3);
    TEST_EQUALITY(basis_values.grad_basis_ref.fieldTag().dataLayout().dimension(0),9);
    TEST_EQUALITY(basis_values.grad_basis_ref.fieldTag().dataLayout().dimension(1),num_qp);
    TEST_EQUALITY(basis_values.grad_basis_ref.fieldTag().dataLayout().dimension(2),2);
    TEST_EQUALITY(basis_values.grad_basis_ref.fieldTag().name(),"prefix_grad_basis_ref");

    TEST_EQUALITY(basis_values.grad_basis.fieldTag().dataLayout().rank(),4);
    TEST_EQUALITY(basis_values.grad_basis.fieldTag().dataLayout().dimension(0),num_cells);
    TEST_EQUALITY(basis_values.grad_basis.fieldTag().dataLayout().dimension(1),9);
    TEST_EQUALITY(basis_values.grad_basis.fieldTag().dataLayout().dimension(2),num_qp);
    TEST_EQUALITY(basis_values.grad_basis.fieldTag().dataLayout().dimension(3),2);
    TEST_EQUALITY(basis_values.grad_basis.fieldTag().name(),"prefix_grad_basis");

    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().dataLayout().rank(),4);
    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().dataLayout().dimension(0),num_cells);
    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().dataLayout().dimension(1),9);
    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().dataLayout().dimension(2),num_qp);
    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().dataLayout().dimension(3),2);
    TEST_EQUALITY(basis_values.weighted_grad_basis.fieldTag().name(),"prefix_weighted_grad_basis");

    // check coordinates
    TEST_EQUALITY(basis_values.basis_coordinates_ref.fieldTag().dataLayout().rank(),2);
    TEST_EQUALITY(basis_values.basis_coordinates_ref.fieldTag().dataLayout().dimension(0),9);
    TEST_EQUALITY(basis_values.basis_coordinates_ref.fieldTag().dataLayout().dimension(1),2);
    TEST_EQUALITY(basis_values.basis_coordinates_ref.fieldTag().name(),"prefix_basis_coordinates_ref");

    TEST_EQUALITY(basis_values.basis_coordinates.fieldTag().dataLayout().rank(),3);
    TEST_EQUALITY(basis_values.basis_coordinates.fieldTag().dataLayout().dimension(0),num_cells);
    TEST_EQUALITY(basis_values.basis_coordinates.fieldTag().dataLayout().dimension(1),9);
    TEST_EQUALITY(basis_values.basis_coordinates.fieldTag().dataLayout().dimension(2),2);
    TEST_EQUALITY(basis_values.basis_coordinates.fieldTag().name(),"prefix_basis_coordinates");
  }
}
