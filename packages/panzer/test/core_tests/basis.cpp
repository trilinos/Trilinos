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
#include "Panzer_BasisIRLayout.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(basis, Q2_2D_volume)
  {
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));
    
    const int num_cells = 20;
    const int base_cell_dimension = topo->getDimension();
    const panzer::CellData cell_data(num_cells,topo);
    const int cubature_degree = 2;

    panzer::IntegrationRule int_rule(cubature_degree,cell_data);
    
    const std::string basis_type = "Q2";
    
    panzer::BasisIRLayout basis(basis_type,-1,int_rule);

    TEST_EQUALITY(basis.cardinality(), 9);
    TEST_EQUALITY(basis.numCells(), 20);
    TEST_EQUALITY(basis.numPoints(), 4);
    TEST_EQUALITY(basis.dimension(), base_cell_dimension);
    TEST_EQUALITY(basis.name(), "Q2:CubaturePoints (Degree=2,volume)");
    TEST_EQUALITY(basis.fieldName(), "Basis: Q2");
    TEST_EQUALITY(basis.fieldNameD1(), "Grad Basis: Q2");
    TEST_EQUALITY(basis.fieldNameD2(), "D2 Basis: Q2");

    Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > >
      intrepid_basis = basis.getIntrepidBasis();

    TEST_ASSERT(!Teuchos::is_null(intrepid_basis));

    const int dim = base_cell_dimension;

    TEST_EQUALITY(basis.basis_ref->size(), 9 * 4);
    TEST_EQUALITY(basis.basis->size(), num_cells * 9 * 4);
    TEST_EQUALITY(basis.basis_grad_ref->size(), 9 * 4 * dim);
    TEST_EQUALITY(basis.basis_grad->size(), num_cells * 9 * 4 * dim);
    TEST_EQUALITY(basis.basis_D2_ref->size(), 9 * 4 * dim * dim);
    TEST_EQUALITY(basis.basis_D2->size(), num_cells * 9 * 4 * dim * dim);
    TEST_EQUALITY(basis.functional->size(), num_cells * 9);
    TEST_EQUALITY(basis.functional_grad->size(), num_cells * 9 * dim);
    TEST_EQUALITY(basis.functional_D2->size(), num_cells * 9 * dim * dim);
  }

  TEUCHOS_UNIT_TEST(basis, Q2_2D_side)
  {
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));
    
    const int num_cells = 20;
    const int base_cell_dimension = topo->getDimension();
    const int cell_local_side_id = 1;
    const panzer::CellData cell_data(num_cells,cell_local_side_id,topo);
    const int cubature_degree = 2;

    panzer::IntegrationRule int_rule(cubature_degree,cell_data);
    
    const std::string basis_type = "Q2";
    
    panzer::BasisIRLayout basis(basis_type,-1,int_rule);

    TEST_EQUALITY(basis.cardinality(), 9);
    TEST_EQUALITY(basis.numCells(), 20);
    TEST_EQUALITY(basis.numPoints(), 2);
    TEST_EQUALITY(basis.dimension(), base_cell_dimension);
    TEST_EQUALITY(basis.name(), "Q2:CubaturePoints (Degree=2,side)");
    TEST_EQUALITY(basis.fieldName(), "Basis: Q2");
    TEST_EQUALITY(basis.fieldNameD1(), "Grad Basis: Q2");
    TEST_EQUALITY(basis.fieldNameD2(), "D2 Basis: Q2");

    Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > >
      intrepid_basis = basis.getIntrepidBasis();

    TEST_ASSERT(!Teuchos::is_null(intrepid_basis));

    const int dim = base_cell_dimension;

    TEST_EQUALITY(basis.basis_ref->size(), 9 * 2);
    TEST_EQUALITY(basis.basis->size(), num_cells * 9 * 2);
    TEST_EQUALITY(basis.basis_grad_ref->size(), 9 * 2 * dim);
    TEST_EQUALITY(basis.basis_grad->size(), num_cells * 9 * 2 * dim);
    TEST_EQUALITY(basis.basis_D2_ref->size(), 9 * 2 * dim * dim);
    TEST_EQUALITY(basis.basis_D2->size(), num_cells * 9 * 2 * dim * dim);
    TEST_EQUALITY(basis.functional->size(), num_cells * 9);
    TEST_EQUALITY(basis.functional_grad->size(), num_cells * 9 * dim);
    TEST_EQUALITY(basis.functional_D2->size(), num_cells * 9 * dim * dim);
  }

  TEUCHOS_UNIT_TEST(basis, TEdge1_2D_volume)
  {
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Triangle<3> >()));
    
    const int num_cells = 20;
    const int base_cell_dimension = topo->getDimension();
    const panzer::CellData cell_data(num_cells, topo);
    const std::string basis_type = "TEdge1";

    Teuchos::RCP<PureBasis> basis = Teuchos::rcp(new PureBasis(basis_type,-1,cell_data));

    TEST_EQUALITY(basis->getElementSpace(),PureBasis::HCURL);
    TEST_EQUALITY(basis->cardinality(),3);
    TEST_EQUALITY(basis->numCells(),num_cells);
    TEST_EQUALITY(basis->dimension(),base_cell_dimension);
    TEST_EQUALITY(basis->name(),basis_type);
    TEST_ASSERT(basis->getIntrepidBasis()!=Teuchos::null);
    TEST_ASSERT(basis->getCellTopology()!=Teuchos::null);
  }

  TEUCHOS_UNIT_TEST(basis, QEdge1_2D_volume)
  {
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));
    
    const int num_cells = 20;
    const int base_cell_dimension = topo->getDimension();
    const panzer::CellData cell_data(num_cells,topo);
    const std::string basis_type = "QEdge1";

    Teuchos::RCP<PureBasis> basis = Teuchos::rcp(new PureBasis(basis_type,-1,cell_data));

    TEST_EQUALITY(basis->getElementSpace(),PureBasis::HCURL);
    TEST_EQUALITY(basis->cardinality(),4);
    TEST_EQUALITY(basis->numCells(),num_cells);
    TEST_EQUALITY(basis->dimension(),base_cell_dimension);
    TEST_EQUALITY(basis->name(),basis_type);
    TEST_ASSERT(basis->getIntrepidBasis()!=Teuchos::null);
    TEST_ASSERT(basis->getCellTopology()!=Teuchos::null);
  }

  TEUCHOS_UNIT_TEST(basis, supported_bases)
  {
    const int num_cells = 20;
    Teuchos::RCP<PureBasis> basis;
      

    // Triangle
    {
      Teuchos::RCP<shards::CellTopology> topo = 
	Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Triangle<3> >()));
      const panzer::CellData cell_data(num_cells,topo);
      
      basis = Teuchos::rcp(new PureBasis("HGrad",1,cell_data));
      TEST_EQUALITY(basis->name(),"HGrad:1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("HGrad",2,cell_data));
      TEST_EQUALITY(basis->name(),"HGrad:2");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("HCurl",1,cell_data));
      TEST_EQUALITY(basis->name(),"HCurl:1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HCURL);
    }

    // Quad
    {
      Teuchos::RCP<shards::CellTopology> topo = 
	Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));
      const panzer::CellData cell_data(num_cells,topo);
      
      basis = Teuchos::rcp(new PureBasis("HGrad",1,cell_data));
      TEST_EQUALITY(basis->name(),"HGrad:1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("HGrad",2,cell_data));
      TEST_EQUALITY(basis->name(),"HGrad:2");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("HCurl",1,cell_data));
      TEST_EQUALITY(basis->name(),"HCurl:1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HCURL);
    }

    // Tet
    {
      Teuchos::RCP<shards::CellTopology> topo = 
	Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<4> >()));
      const panzer::CellData cell_data(num_cells,topo);
      
      basis = Teuchos::rcp(new PureBasis("HGrad",1,cell_data));
      TEST_EQUALITY(basis->name(),"HGrad:1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("HGrad",2,cell_data));
      TEST_EQUALITY(basis->name(),"HGrad:2");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("HCurl",1,cell_data));
      TEST_EQUALITY(basis->name(),"HCurl:1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HCURL);
    }

    // Hex
    {
      Teuchos::RCP<shards::CellTopology> topo = 
	Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<8> >()));
      const panzer::CellData cell_data(num_cells,topo);
      
      basis = Teuchos::rcp(new PureBasis("HGrad",1,cell_data));
      TEST_EQUALITY(basis->name(),"HGrad:1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("HGrad",2,cell_data));
      TEST_EQUALITY(basis->name(),"HGrad:2");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("HCurl",1,cell_data));
      TEST_EQUALITY(basis->name(),"HCurl:1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HCURL);
    }

    // Line
    {
      Teuchos::RCP<shards::CellTopology> topo = 
	Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Line<2> >()));
      const panzer::CellData cell_data(num_cells,topo);
      
      basis = Teuchos::rcp(new PureBasis("HGrad",1,cell_data));
      TEST_EQUALITY(basis->name(),"HGrad:1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
    }

  }

  TEUCHOS_UNIT_TEST(basis, deprecated_bases)
  {
    const int num_cells = 20;

    
    Teuchos::RCP<PureBasis> basis;
      

    // T1, T2
    {
      Teuchos::RCP<shards::CellTopology> topo = 
	Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Triangle<3> >()));
      const panzer::CellData cell_data(num_cells,topo);
      
      basis = Teuchos::rcp(new PureBasis("T1",0,cell_data));
      TEST_EQUALITY(basis->name(),"T1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("T2",0,cell_data));
      TEST_EQUALITY(basis->name(),"T2");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("TEdge1",0,cell_data));
      TEST_EQUALITY(basis->name(),"TEdge1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HCURL);
    }

    // Quad
    {
      Teuchos::RCP<shards::CellTopology> topo = 
	Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));
      const panzer::CellData cell_data(num_cells,topo);
      
      basis = Teuchos::rcp(new PureBasis("Q1",0,cell_data));
      TEST_EQUALITY(basis->name(),"Q1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("Q2",0,cell_data));
      TEST_EQUALITY(basis->name(),"Q2");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("QEdge1",0,cell_data));
      TEST_EQUALITY(basis->name(),"QEdge1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HCURL);
    }

    // Tet
    {
      Teuchos::RCP<shards::CellTopology> topo = 
	Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<4> >()));
      const panzer::CellData cell_data(num_cells,topo);
      
      basis = Teuchos::rcp(new PureBasis("T1",0,cell_data));
      TEST_EQUALITY(basis->name(),"T1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("T2",0,cell_data));
      TEST_EQUALITY(basis->name(),"T2");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("TEdge1",0,cell_data));
      TEST_EQUALITY(basis->name(),"TEdge1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HCURL);
    }

    // Hex
    {
      Teuchos::RCP<shards::CellTopology> topo = 
	Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<8> >()));
      const panzer::CellData cell_data(num_cells,topo);
      
      basis = Teuchos::rcp(new PureBasis("Q1",0,cell_data));
      TEST_EQUALITY(basis->name(),"Q1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("Q2",0,cell_data));
      TEST_EQUALITY(basis->name(),"Q2");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("QEdge1",0,cell_data));
      TEST_EQUALITY(basis->name(),"QEdge1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HCURL);
    }

  }

}
