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
    TEST_EQUALITY(basis.name(), "HGrad:2:CubaturePoints (Degree=2,volume)");
    TEST_EQUALITY(basis.fieldName(), "Basis: HGrad:2");
    TEST_EQUALITY(basis.fieldNameD1(), "Grad Basis: HGrad:2");
    TEST_EQUALITY(basis.fieldNameD2(), "D2 Basis: HGrad:2");

    Teuchos::RCP< Intrepid2::Basis<PHX::exec_space,double,double> >
      intrepid_basis = basis.getIntrepid2Basis();

    TEST_ASSERT(!Teuchos::is_null(intrepid_basis));

    const int dim = base_cell_dimension;

    TEST_EQUALITY(static_cast<int>(basis.basis_ref->size()), 9 * 4);
    TEST_EQUALITY(static_cast<int>(basis.basis->size()), num_cells * 9 * 4);
    TEST_EQUALITY(static_cast<int>(basis.basis_grad_ref->size()), 9 * 4 * dim);
    TEST_EQUALITY(static_cast<int>(basis.basis_grad->size()), num_cells * 9 * 4 * dim);
    TEST_EQUALITY(static_cast<int>(basis.basis_D2_ref->size()), 9 * 4 * dim * dim);
    TEST_EQUALITY(static_cast<int>(basis.basis_D2->size()), num_cells * 9 * 4 * dim * dim);
    TEST_EQUALITY(static_cast<int>(basis.functional->size()), num_cells * 9);
    TEST_EQUALITY(static_cast<int>(basis.functional_grad->size()), num_cells * 9 * dim);
    TEST_EQUALITY(static_cast<int>(basis.functional_D2->size()), num_cells * 9 * dim * dim);
    TEST_EQUALITY(static_cast<int>(basis.getBasis()->local_mat_layout->size()), num_cells * 9 * 9);

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
    TEST_EQUALITY(basis.name(), "HGrad:2:CubaturePoints (Degree=2,side)");
    TEST_EQUALITY(basis.fieldName(), "Basis: HGrad:2");
    TEST_EQUALITY(basis.fieldNameD1(), "Grad Basis: HGrad:2");
    TEST_EQUALITY(basis.fieldNameD2(), "D2 Basis: HGrad:2");

    Teuchos::RCP< Intrepid2::Basis<PHX::exec_space,double,double> >
      intrepid_basis = basis.getIntrepid2Basis();

    TEST_ASSERT(!Teuchos::is_null(intrepid_basis));

    const int dim = base_cell_dimension;

    TEST_EQUALITY(static_cast<int>(basis.basis_ref->size()), 9 * 2);
    TEST_EQUALITY(static_cast<int>(basis.basis->size()), num_cells * 9 * 2);
    TEST_EQUALITY(static_cast<int>(basis.basis_grad_ref->size()), 9 * 2 * dim);
    TEST_EQUALITY(static_cast<int>(basis.basis_grad->size()), num_cells * 9 * 2 * dim);
    TEST_EQUALITY(static_cast<int>(basis.basis_D2_ref->size()), 9 * 2 * dim * dim);
    TEST_EQUALITY(static_cast<int>(basis.basis_D2->size()), num_cells * 9 * 2 * dim * dim);
    TEST_EQUALITY(static_cast<int>(basis.functional->size()), num_cells * 9);
    TEST_EQUALITY(static_cast<int>(basis.functional_grad->size()), num_cells * 9 * dim);
    TEST_EQUALITY(static_cast<int>(basis.functional_D2->size()), num_cells * 9 * dim * dim);
    TEST_EQUALITY(static_cast<int>(basis.getBasis()->local_mat_layout->size()), num_cells * 9 * 9);

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
    TEST_EQUALITY(basis->name(),"HCurl:1");
    TEST_ASSERT(basis->getIntrepid2Basis()!=Teuchos::null);
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
    TEST_EQUALITY(basis->name(),"HCurl:1");
    TEST_ASSERT(basis->getIntrepid2Basis()!=Teuchos::null);
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
      TEST_EQUALITY(basis->name(),"HGrad:1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("T2",0,cell_data));
      TEST_EQUALITY(basis->name(),"HGrad:2");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("TEdge1",0,cell_data));
      TEST_EQUALITY(basis->name(),"HCurl:1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HCURL);
    }

    // Quad
    {
      Teuchos::RCP<shards::CellTopology> topo = 
	Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));
      const panzer::CellData cell_data(num_cells,topo);
      
      basis = Teuchos::rcp(new PureBasis("Q1",0,cell_data));
      TEST_EQUALITY(basis->name(),"HGrad:1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("Q2",0,cell_data));
      TEST_EQUALITY(basis->name(),"HGrad:2");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("QEdge1",0,cell_data));
      TEST_EQUALITY(basis->name(),"HCurl:1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HCURL);
    }

    // Tet
    {
      Teuchos::RCP<shards::CellTopology> topo = 
	Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<4> >()));
      const panzer::CellData cell_data(num_cells,topo);
      
      basis = Teuchos::rcp(new PureBasis("T1",0,cell_data));
      TEST_EQUALITY(basis->name(),"HGrad:1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("T2",0,cell_data));
      TEST_EQUALITY(basis->name(),"HGrad:2");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("TEdge1",0,cell_data));
      TEST_EQUALITY(basis->name(),"HCurl:1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HCURL);
    }

    // Hex
    {
      Teuchos::RCP<shards::CellTopology> topo = 
	Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<8> >()));
      const panzer::CellData cell_data(num_cells,topo);
      
      basis = Teuchos::rcp(new PureBasis("Q1",0,cell_data));
      TEST_EQUALITY(basis->name(),"HGrad:1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("Q2",0,cell_data));
      TEST_EQUALITY(basis->name(),"HGrad:2");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HGRAD);
      basis = Teuchos::rcp(new PureBasis("QEdge1",0,cell_data));
      TEST_EQUALITY(basis->name(),"HCurl:1");
      TEST_EQUALITY(basis->getElementSpace(),PureBasis::HCURL);
    }

  }

}
