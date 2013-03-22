// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ArrayView.hpp"

#include "Stokhos.hpp"
#include "Stokhos_UnitTestHelpers.hpp"

namespace SmolyakBasisUtilsUnitTest {

  // Common setup for unit tests
  template <typename ordinal_type, typename value_type>
  struct UnitTestSetup {
    value_type rtol, atol, sparse_tol;
    ordinal_type p,d;
    
    UnitTestSetup() : 
      rtol(1e-12), 
      atol(1e-12), 
      sparse_tol(1e-12),
      p(5),
      d(3) {}
    
  };

  typedef int ordinal_type;
  typedef double value_type;
  typedef UnitTestSetup<ordinal_type,value_type> setup_type;
  setup_type setup;

  template <typename ordinal_type, typename value_type>
  bool testCoefficients(
    const Stokhos::ProductBasis<ordinal_type,value_type>& basis1, 
    const Stokhos::ProductBasis<ordinal_type,value_type>& basis2,
    Teuchos::FancyOStream& out) {
    bool success = true;

     // Check sizes
    TEUCHOS_TEST_EQUALITY(basis1.size(), basis2.size(), out, success);
    if (!success)
      return false;

    ordinal_type coeff_sz = basis1.size();

    // Check coefficients
    for (ordinal_type i=0; i<coeff_sz; ++i) {
      out << "coeff " << basis1.term(i) << " == " << basis2.term(i) << " : ";
      bool is_equal = basis1.term(i) == basis2.term(i);
      if (is_equal)
	out << "passed" << std::endl;
      else {
	out << "failed" << std::endl;
	success = false; 
      }
    }

    return success;
  }

  TEUCHOS_UNIT_TEST( Coefficients, IsotropicLinear ) {
    // Build isotropic Smolyak basis of dimension d and order p with 
    // linear growth -- note order of basis is determined by index set and
    // growth rule, not coordinate bases
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(setup.p, true));
    Stokhos::TotalOrderIndexSet<ordinal_type> coeff_index_set(setup.d, setup.p);
    Stokhos::SmolyakBasis<ordinal_type,value_type> smolyak_basis(
      bases, coeff_index_set);
    
    // Build isotropic total order basis of dimension d and order p
    Stokhos::TotalOrderBasis<ordinal_type,value_type> total_order_basis(bases);

    // Two basis should be identical
    success = testCoefficients(smolyak_basis, total_order_basis, out);
  }

  TEUCHOS_UNIT_TEST( Coefficients, AnisotropicLinear ) {
    // Build anisotropic Smolyak basis of dimension d with 
    // linear growth -- note order of basis is determined by index set and
    // growth rule, not coordinate bases
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    Stokhos::MultiIndex<ordinal_type> upper(setup.d);
    for (ordinal_type i=0; i<setup.d; i++) {
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(i+1, true));
      upper[i] = i+1;
    }
    Stokhos::AnisotropicTotalOrderIndexSet<ordinal_type> coeff_index_set(
      setup.d, upper);  // largest order is setup.d-1+1
    Stokhos::SmolyakBasis<ordinal_type,value_type> smolyak_basis(
      bases, coeff_index_set);
    
    // Build isotropic total order basis of dimension d and order p
    Stokhos::TotalOrderBasis<ordinal_type,value_type> total_order_basis(bases);

    // Two basis should be identical
    success = testCoefficients(smolyak_basis, total_order_basis, out);
  }

  TEUCHOS_UNIT_TEST( Sparse3Tensor, AnisotropicLinear ) {
    success = true;

    // Build anisotropic Smolyak basis of dimension d with 
    // linear growth -- note order of basis is determined by index set and
    // growth rule, not coordinate bases
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    Stokhos::MultiIndex<ordinal_type> upper(setup.d);
    for (ordinal_type i=0; i<setup.d; i++) {
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(i+1, true));
      upper[i] = i+1;
    }
    Stokhos::AnisotropicTotalOrderIndexSet<ordinal_type> coeff_index_set(
      setup.d, upper);  // largest order is setup.d-1+1
    Stokhos::SmolyakBasis<ordinal_type,value_type> smolyak_basis(
      bases, coeff_index_set);
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk1 = 
      smolyak_basis.computeTripleProductTensor();
    
    // Build isotropic total order basis of dimension d and order p
    Stokhos::TotalOrderBasis<ordinal_type,value_type> total_order_basis(bases);
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk2 = 
      total_order_basis.computeTripleProductTensor();

    // Compare Cijk tensors
    success = Stokhos::compareSparse3Tensor(*Cijk1, "Smolyak Cijk", 
					    *Cijk2, "Total Order Cijk",
					    setup.rtol, setup.atol, out);

    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk1_lin= 
      smolyak_basis.computeLinearTripleProductTensor();
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk2_lin= 
      total_order_basis.computeLinearTripleProductTensor();

    // Compare Cijk tensors
    success = success && 
      Stokhos::compareSparse3Tensor(*Cijk1_lin, "Smolyak linear Cijk", 
				    *Cijk2_lin, "Total Order linear Cijk",
				    setup.rtol, setup.atol, out);

  }

  TEUCHOS_UNIT_TEST( Sparse3Tensor, AnisotropicLinear2 ) {
    success = true;

    // Build anisotropic Smolyak basis of dimension d with 
    // linear growth -- note order of basis is determined by index set and
    // growth rule, not coordinate bases
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    Stokhos::MultiIndex<ordinal_type> upper(setup.d);
    for (ordinal_type i=0; i<setup.d; i++) {
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(i+1, true));
      upper[i] = i+1;
    }
    Stokhos::AnisotropicTotalOrderIndexSet<ordinal_type> coeff_index_set(
      setup.d, upper);  // largest order is setup.d-1+1
    Stokhos::SmolyakBasis<ordinal_type,value_type> smolyak_basis(
      bases, coeff_index_set, setup.sparse_tol);
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk = 
      smolyak_basis.computeTripleProductTensor();

    success = Stokhos::testSparse3Tensor(*Cijk, smolyak_basis,
					 setup.sparse_tol, setup.rtol, 
					 setup.atol, out);

    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk_lin = 
      smolyak_basis.computeLinearTripleProductTensor();

    success = success && 
      Stokhos::testSparse3Tensor(*Cijk_lin, smolyak_basis,
				 setup.sparse_tol, setup.rtol, 
				 setup.atol, out, true);

  }

  TEUCHOS_UNIT_TEST( Sparse3Tensor, AnisotropicModerateLinear ) {
    success = true;

    // Build anisotropic Smolyak basis of dimension d with 
    // moderate linear growth -- note order of basis is determined by index 
    // set and growth rule, not coordinate bases
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    Stokhos::MultiIndex<ordinal_type> upper(setup.d);
    for (ordinal_type i=0; i<setup.d; i++) {
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(i+1, true, Stokhos::MODERATE_GROWTH));
      upper[i] = i+1;
    }
    Stokhos::AnisotropicTotalOrderIndexSet<ordinal_type> coeff_index_set(
      setup.d, upper);  // largest order is setup.d-1+1
    Stokhos::SmolyakBasis<ordinal_type,value_type> smolyak_basis(
      bases, coeff_index_set, setup.sparse_tol);
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk = 
      smolyak_basis.computeTripleProductTensor();

    success = Stokhos::testSparse3Tensor(*Cijk, smolyak_basis,
					 setup.sparse_tol, setup.rtol, 
					 setup.atol, out);

    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk_lin = 
      smolyak_basis.computeLinearTripleProductTensor();

    success = success && 
      Stokhos::testSparse3Tensor(*Cijk_lin, smolyak_basis,
				 setup.sparse_tol, setup.rtol, 
				 setup.atol, out, true);
    
  }

#ifdef HAVE_STOKHOS_DAKOTA

  TEUCHOS_UNIT_TEST( Sparse3Tensor, AnisotropicClenshawCurtis ) {
    success = true;

    // Build anisotropic Smolyak basis of dimension d with 
    // exponential growth -- note order of basis is determined by index 
    // set and growth rule, not coordinate bases
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    Stokhos::MultiIndex<ordinal_type> upper(setup.d);
    for (ordinal_type i=0; i<setup.d; i++) {
      bases[i] = Teuchos::rcp(new Stokhos::ClenshawCurtisLegendreBasis<ordinal_type,value_type>(i+1, true));
      upper[i] = i+1;
    }
    Stokhos::AnisotropicTotalOrderIndexSet<ordinal_type> coeff_index_set(
      setup.d, upper);  // largest order is setup.d-1+1
    Stokhos::SmolyakBasis<ordinal_type,value_type> smolyak_basis(
      bases, coeff_index_set, setup.sparse_tol);
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk = 
      smolyak_basis.computeTripleProductTensor();

    success = Stokhos::testSparse3Tensor(*Cijk, smolyak_basis,
					 setup.sparse_tol, setup.rtol, 
					 setup.atol, out);

    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk_lin = 
      smolyak_basis.computeLinearTripleProductTensor();

    success = success && 
      Stokhos::testSparse3Tensor(*Cijk_lin, smolyak_basis,
				 setup.sparse_tol, setup.rtol, 
				 setup.atol, out, true);
    
  }

  TEUCHOS_UNIT_TEST( Sparse3Tensor, AnisotropicGaussPatterson ) {
    success = true;

    // Build anisotropic Smolyak basis of dimension d with 
    // exponential growth -- note order of basis is determined by index 
    // set and growth rule, not coordinate bases
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    Stokhos::MultiIndex<ordinal_type> upper(setup.d);
    for (ordinal_type i=0; i<setup.d; i++) {
      bases[i] = Teuchos::rcp(new Stokhos::GaussPattersonLegendreBasis<ordinal_type,value_type>(i+1, true));
      upper[i] = i+1;
    }
    Stokhos::AnisotropicTotalOrderIndexSet<ordinal_type> coeff_index_set(
      setup.d, upper);  // largest order is setup.d-1+1
    Stokhos::SmolyakBasis<ordinal_type,value_type> smolyak_basis(
      bases, coeff_index_set, setup.sparse_tol);
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk = 
      smolyak_basis.computeTripleProductTensor();

    success = Stokhos::testSparse3Tensor(*Cijk, smolyak_basis,
					 setup.sparse_tol, setup.rtol, 
					 setup.atol, out);

    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk_lin = 
      smolyak_basis.computeLinearTripleProductTensor();

    success = success && 
      Stokhos::testSparse3Tensor(*Cijk_lin, smolyak_basis,
				 setup.sparse_tol, setup.rtol, 
				 setup.atol, out, true);
    
  }

#endif

}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
