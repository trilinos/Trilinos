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

#include "Stokhos.hpp"
#include "Stokhos_UnitTestHelpers.hpp"

namespace TensorProductUnitTest {

  // Common setup for unit tests
  template <typename OrdinalType, typename ValueType>
  struct UnitTestSetup {
    ValueType rtol, atol, sparse_tol;
    OrdinalType p,d;
    
    UnitTestSetup() {
      rtol = 1e-12;
      atol = 1e-12;
      sparse_tol = 1e-12;
      d = 3;
      p = 4;
    }
    
  };

  typedef int ordinal_type;
  typedef double value_type;
  UnitTestSetup<ordinal_type,value_type> setup;

  TEUCHOS_UNIT_TEST( Coefficients, Isotropic ) {
    success = true;

    // Build tensor product basis of dimension d and order p
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(setup.p, true));
    Teuchos::RCP<const Stokhos::TensorProductBasis<ordinal_type,value_type> > basis = Teuchos::rcp(new Stokhos::TensorProductBasis<ordinal_type,value_type>(bases));

    // Compute expected size
    ordinal_type sz = 1;
    for (ordinal_type i=0; i<setup.d; ++i)
      sz *= setup.p+1;

    // Check sizes
    TEUCHOS_TEST_EQUALITY(sz, basis->size(), out, success);

    std::ostream_iterator<ordinal_type> out_iterator(out, " ");
    for (ordinal_type i=0; i<sz; i++) {
      const Stokhos::MultiIndex<ordinal_type>& term = basis->term(i);

      // Verify terms match
      out << "term " << term << " <= " << setup.p << " : ";
      bool is_less = true;
      for (ordinal_type j=0; j<setup.d; j++)
	is_less = is_less && term[j] <= setup.p;
      if (is_less)
	out << "passed" << std::endl;
      else {
	out << "failed" << std::endl;
	success = false; 
      }

    }

  }

  TEUCHOS_UNIT_TEST( Coeficients, Anisotropic ) {
    success = true;

    // Build anisotropic tensor product basis of dimension d
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(i+1, true));
    Teuchos::RCP<const Stokhos::TensorProductBasis<ordinal_type,value_type> > basis = Teuchos::rcp(new Stokhos::TensorProductBasis<ordinal_type,value_type>(bases));

    // Compute expected size
    ordinal_type sz = 1;
    for (ordinal_type i=0; i<setup.d; ++i)
      sz *= i+2;

    // Check sizes
    TEUCHOS_TEST_EQUALITY(sz, basis->size(), out, success);

    std::ostream_iterator<ordinal_type> out_iterator(out, " ");
    for (ordinal_type i=0; i<sz; i++) {
      const Stokhos::MultiIndex<ordinal_type>& term = basis->term(i);

      // Verify terms match
      out << "term " << term << " <= " << setup.p << " : ";
      bool is_less = true;
      for (ordinal_type j=0; j<setup.d; j++)
	is_less = is_less && term[j] <= setup.p;
      if (is_less)
	out << "passed" << std::endl;
      else {
	out << "failed" << std::endl;
	success = false; 
      }

    }

  }

  TEUCHOS_UNIT_TEST( Sparse3Tensor, Anisotropic ) {
    success = true;

    // Build anisotropic tensor product basis of dimension d
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(i+1, true));
    Teuchos::RCP<const Stokhos::ProductBasis<ordinal_type,value_type> > basis = Teuchos::rcp(new Stokhos::TensorProductBasis<ordinal_type,value_type>(bases, setup.sparse_tol));
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk = 
      basis->computeTripleProductTensor(basis->size());

    success = Stokhos::testSparse3Tensor(*Cijk, *basis, setup.sparse_tol, 
					 setup.rtol, setup.atol, out);
  }

  
}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  int res = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
  return res;
}
