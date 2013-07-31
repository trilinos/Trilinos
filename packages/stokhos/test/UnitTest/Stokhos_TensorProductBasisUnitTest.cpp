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

#include <iterator>

namespace TensorProductBasisUnitTest {

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

  TEUCHOS_UNIT_TEST( Sparse3Tensor, Anisotropic_Full ) {
    success = true;

    // Build anisotropic tensor product basis of dimension d
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(i+1, true));
    Teuchos::RCP<const Stokhos::ProductBasis<ordinal_type,value_type> > basis = Teuchos::rcp(new Stokhos::TensorProductBasis<ordinal_type,value_type>(bases, setup.sparse_tol));
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk = 
      basis->computeTripleProductTensor();

    success = Stokhos::testSparse3Tensor(*Cijk, *basis, setup.sparse_tol, 
					 setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Sparse3Tensor, Anisotropic_Linear ) {
    success = true;

    // Build anisotropic tensor product basis of dimension d
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(i+1, true));
    Teuchos::RCP<const Stokhos::ProductBasis<ordinal_type,value_type> > basis = Teuchos::rcp(new Stokhos::TensorProductBasis<ordinal_type,value_type>(bases, setup.sparse_tol));
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk = 
      basis->computeLinearTripleProductTensor();

    success = Stokhos::testSparse3Tensor(*Cijk, *basis, setup.sparse_tol, 
					 setup.rtol, setup.atol, out, true);
  }

  
}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  int res = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
  return res;
}
