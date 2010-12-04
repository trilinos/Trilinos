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

namespace Sparse3TensorUnitTest {

  typedef Stokhos::Sparse3Tensor<int,double> Cijk_type;

  // Common setup for unit tests
  template <typename OrdinalType, typename ValueType>
  struct UnitTestSetup {
    ValueType rtol, atol;
    OrdinalType sz;
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<OrdinalType,ValueType> > basis;
    Teuchos::RCP<const Stokhos::Quadrature<OrdinalType,ValueType> > quad;
    Teuchos::RCP<Cijk_type> Cijk;
    
    UnitTestSetup() {
      rtol = 1e-12;
      atol = 1e-12;
      const OrdinalType d = 3;
      const OrdinalType p = 5;
      
      // Create product basis
      Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<OrdinalType,ValueType> > > bases(d);
      for (OrdinalType i=0; i<d; i++)
	bases[i] = 
	  Teuchos::rcp(new Stokhos::LegendreBasis<OrdinalType,ValueType>(p));
      basis =
	Teuchos::rcp(new Stokhos::CompletePolynomialBasis<OrdinalType,ValueType>(bases));
      sz = basis->size();
      
      // Tensor product quadrature
      quad = 
	Teuchos::rcp(new Stokhos::TensorProductQuadrature<OrdinalType,ValueType>(basis, 3*p));
      
      
      // Triple product tensor
      Cijk = basis->computeTripleProductTensor(sz);
    }
    
  };

  UnitTestSetup<int,double> setup;

  TEUCHOS_UNIT_TEST( Stokhos_Sparse3Tensor, Values ) {
    const Teuchos::Array<double>& weights = setup.quad->getQuadWeights();
    const Teuchos::Array< Teuchos::Array<double> > & values = 
      setup.quad->getBasisAtQuadPoints();

    success = true;
    for (Cijk_type::k_iterator k_it=setup.Cijk->k_begin(); 
	 k_it!=setup.Cijk->k_end(); ++k_it) {
      int k = Stokhos::index(k_it);
      for (Cijk_type::kj_iterator j_it = setup.Cijk->j_begin(k_it); 
	   j_it != setup.Cijk->j_end(k_it); ++j_it) {
	int j = Stokhos::index(j_it);
	for (Cijk_type::kji_iterator i_it = setup.Cijk->i_begin(j_it);
	     i_it != setup.Cijk->i_end(j_it); ++i_it) {
	  int i = Stokhos::index(i_it);
	  double c = Stokhos::value(i_it);

	  double c2 = 0.0;
	  int nqp = weights.size();
	  for (int qp=0; qp<nqp; qp++)
	    c2 += weights[qp]*values[qp][i]*values[qp][j]*values[qp][k];

	  double tol = setup.atol + c*setup.rtol;
	  double err = std::abs(c-c2);
	  bool s = err < tol;
	  if (!s) {
	    out << std::endl
		<< "Check: rel_err( C(" << i << "," << j << "," << k << ") )"
		<< " = " << "rel_err( " << c << ", " << c2 << " ) = " << err 
		<< " <= " << tol << " : ";
	    if (s) out << "Passed.";
	    else 
	      out << "Failed!";
	    out << std::endl;
	  }
	  success = success && s;
	}
      }
    }
  }
  

}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
