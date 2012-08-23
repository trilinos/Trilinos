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
    const int base_cell_dimension = 2;
    const panzer::CellData cell_data(num_cells, base_cell_dimension,topo);
    const int cubature_degree = 2;

    panzer::IntegrationRule int_rule(cubature_degree, cell_data);
    
    const std::string basis_type = "Q2";
    
    panzer::BasisIRLayout basis(basis_type, int_rule);

    TEST_ASSERT(basis.getCardinality() == 9);
    TEST_ASSERT(basis.getNumCells() == 20);
    TEST_ASSERT(basis.getNumPoints() == 4);
    TEST_ASSERT(basis.getDimension() == base_cell_dimension);
    // TEST_ASSERT(basis.integrationRuleDegree() == cubature_degree);
    TEST_ASSERT(basis.name() == basis_type);
    TEST_ASSERT(basis.fieldName() == "Basis: Q2");
    TEST_ASSERT(basis.fieldNameD1() == "Grad Basis: Q2");
    TEST_ASSERT(basis.fieldNameD2() == "D2 Basis: Q2");

    Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > >
      intrepid_basis = basis.getIntrepidBasis();

    TEST_ASSERT(!Teuchos::is_null(intrepid_basis));

    const int dim = base_cell_dimension;

    TEST_ASSERT(basis.basis_ref->size() == 9 * 4);
    TEST_ASSERT(basis.basis->size() == num_cells * 9 * 4);
    TEST_ASSERT(basis.basis_grad_ref->size() == 9 * 4 * dim);
    TEST_ASSERT(basis.basis_grad->size() == num_cells * 9 * 4 * dim);
    TEST_ASSERT(basis.basis_D2_ref->size() == 9 * 4 * dim * dim);
    TEST_ASSERT(basis.basis_D2->size() == num_cells * 9 * 4 * dim * dim);
    TEST_ASSERT(basis.functional->size() == num_cells * 9);
    TEST_ASSERT(basis.functional_grad->size() == num_cells * 9 * dim);
    TEST_ASSERT(basis.functional_D2->size() == num_cells * 9 * dim * dim);
  }

  TEUCHOS_UNIT_TEST(basis, Q2_2D_side)
  {
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));
    
    const int num_cells = 20;
    const int base_cell_dimension = 2;
    const int cell_local_side_id = 1;
    const panzer::CellData cell_data(num_cells, base_cell_dimension,
				     cell_local_side_id,topo);
    const int cubature_degree = 2;

    panzer::IntegrationRule int_rule(cubature_degree, cell_data);
    
    const std::string basis_type = "Q2";
    
    panzer::BasisIRLayout basis(basis_type, int_rule);

    std::cout << "ROGER = " << basis.basis_ref->size();

    TEST_ASSERT(basis.getCardinality() == 9);
    TEST_ASSERT(basis.getNumCells() == 20);
    TEST_ASSERT(basis.getNumPoints() == 2);
    TEST_ASSERT(basis.getDimension() == base_cell_dimension);
    // TEST_ASSERT(basis.integrationRuleDegree() == cubature_degree);
    TEST_ASSERT(basis.name() == basis_type);
    TEST_ASSERT(basis.fieldName() == "Basis: Q2");
    TEST_ASSERT(basis.fieldNameD1() == "Grad Basis: Q2");
    TEST_ASSERT(basis.fieldNameD2() == "D2 Basis: Q2");

    Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > >
      intrepid_basis = basis.getIntrepidBasis();

    TEST_ASSERT(!Teuchos::is_null(intrepid_basis));

    const int dim = base_cell_dimension;

    TEST_ASSERT(basis.basis_ref->size() == 9 * 2);
    TEST_ASSERT(basis.basis->size() == num_cells * 9 * 2);
    TEST_ASSERT(basis.basis_grad_ref->size() == 9 * 2 * dim);
    TEST_ASSERT(basis.basis_grad->size() == num_cells * 9 * 2 * dim);
    TEST_ASSERT(basis.basis_D2_ref->size() == 9 * 2 * dim * dim);
    TEST_ASSERT(basis.basis_D2->size() == num_cells * 9 * 2 * dim * dim);
    TEST_ASSERT(basis.functional->size() == num_cells * 9);
    TEST_ASSERT(basis.functional_grad->size() == num_cells * 9 * dim);
    TEST_ASSERT(basis.functional_D2->size() == num_cells * 9 * dim * dim);
  }

  TEUCHOS_UNIT_TEST(basis, TEdge1_2D_volume)
  {
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Triangle<3> >()));
    
    const int num_cells = 20;
    const int base_cell_dimension = 2;
    const panzer::CellData cell_data(num_cells, base_cell_dimension,topo);
    const std::string basis_type = "TEdge1";

    Teuchos::RCP<PureBasis> basis = Teuchos::rcp(new PureBasis(basis_type,cell_data));

    TEST_EQUALITY(basis->getElementSpace(),PureBasis::HCURL);
    TEST_EQUALITY(basis->getCardinality(),3);
    TEST_EQUALITY(basis->getNumCells(),num_cells);
    TEST_EQUALITY(basis->getDimension(),base_cell_dimension);
    TEST_EQUALITY(basis->name(),basis_type);
    TEST_ASSERT(basis->getIntrepidBasis()!=Teuchos::null);
    TEST_ASSERT(basis->getCellTopology()!=Teuchos::null);
  }

  TEUCHOS_UNIT_TEST(basis, QEdge1_2D_volume)
  {
    Teuchos::RCP<shards::CellTopology> topo = 
       Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));
    
    const int num_cells = 20;
    const int base_cell_dimension = 2;
    const panzer::CellData cell_data(num_cells, base_cell_dimension,topo);
    const std::string basis_type = "QEdge1";

    Teuchos::RCP<PureBasis> basis = Teuchos::rcp(new PureBasis(basis_type,cell_data));

    TEST_EQUALITY(basis->getElementSpace(),PureBasis::HCURL);
    TEST_EQUALITY(basis->getCardinality(),4);
    TEST_EQUALITY(basis->getNumCells(),num_cells);
    TEST_EQUALITY(basis->getDimension(),base_cell_dimension);
    TEST_EQUALITY(basis->name(),basis_type);
    TEST_ASSERT(basis->getIntrepidBasis()!=Teuchos::null);
    TEST_ASSERT(basis->getCellTopology()!=Teuchos::null);
  }

}
