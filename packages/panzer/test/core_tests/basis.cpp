#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Panzer_CellData.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Basis.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(basis, Q2_2D_volume)
  {
    const int num_cells = 20;
    const int base_cell_dimension = 2;
    const panzer::CellData cell_data(num_cells, base_cell_dimension);
    const int cubature_degree = 2;
    
    panzer::IntegrationRule int_rule(cubature_degree, cell_data);
    
    const std::string basis_type = "Q2";
    
    panzer::Basis basis(basis_type, int_rule);

    TEST_ASSERT(basis.getCardinality() == 9);
    TEST_ASSERT(basis.getNumCells() == 20);
    TEST_ASSERT(basis.getNumIntPoints() == 4);
    TEST_ASSERT(basis.getDimension() == base_cell_dimension);
    TEST_ASSERT(basis.integrationRuleDegree() == cubature_degree);
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
    const int num_cells = 20;
    const int base_cell_dimension = 2;
    const int cell_local_side_id = 1;
    const panzer::CellData cell_data(num_cells, base_cell_dimension,
				     cell_local_side_id);
    const int cubature_degree = 2;
    
    panzer::IntegrationRule int_rule(cubature_degree, cell_data);
    
    const std::string basis_type = "Q2";
    
    panzer::Basis basis(basis_type, int_rule);

    std::cout << "ROGER = " << basis.basis_ref->size();

    TEST_ASSERT(basis.getCardinality() == 9);
    TEST_ASSERT(basis.getNumCells() == 20);
    TEST_ASSERT(basis.getNumIntPoints() == 2);
    TEST_ASSERT(basis.getDimension() == base_cell_dimension);
    TEST_ASSERT(basis.integrationRuleDegree() == cubature_degree);
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

}
