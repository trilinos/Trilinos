// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2008) Sandia Corporation
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

  // Common setup for unit tests
  template <typename OrdinalType, typename ValueType>
  struct UnitTestSetup {
    ValueType rtol, atol;
    OrdinalType sz;
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<OrdinalType,ValueType> > basis;
    Teuchos::RCP<const Stokhos::Quadrature<OrdinalType,ValueType> > quad;
    Teuchos::RCP<const Stokhos::Sparse3Tensor<OrdinalType,ValueType> > Cijk;
    
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
      Cijk = basis->getTripleProductTensor();
    }
    
  };

  UnitTestSetup<int,double> setup;

  TEUCHOS_UNIT_TEST( Stokhos_Sparse3Tensor, Values ) {
    const Teuchos::Array<double>& weights = setup.quad->getQuadWeights();
    const Teuchos::Array< Teuchos::Array<double> > & values = 
      setup.quad->getBasisAtQuadPoints();

    success = true;
    for (int k=0; k<setup.sz; k++) {
      int nl = setup.Cijk->num_values(k);
      for (int l=0; l<nl; l++) {
	int i, j;
	double c;
	setup.Cijk->value(k,l,i,j,c);

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

  TEUCHOS_UNIT_TEST( Stokhos_Sparse3Tensor, Values2 ) {
    const Teuchos::Array<double>& weights = setup.quad->getQuadWeights();
    const Teuchos::Array< Teuchos::Array<double> > & values = 
      setup.quad->getBasisAtQuadPoints();

    success = true;
    for (int k=0; k<setup.sz; k++) {
      int nj = setup.Cijk->num_j(k);
      const Teuchos::Array<int>& j_indices = setup.Cijk->Jindices(k);
      bool ss = nj == static_cast<int>(j_indices.size());
      if (!ss) {
	out << std::endl
	    << "Check:  Cijk->num_j(" << k 
	    << ") == Cijk->Jindices(" << k << ").size() = "
	    << setup.Cijk->num_j(k) << " == " << j_indices.size() 
	    << ":  Failed!" << std::endl;
      }
      success = success && ss;
      for (int l=0; l<nj; l++) {
	int j = j_indices[l];
	const Teuchos::Array<int>& i_indices = setup.Cijk->Iindices(k,l);
	const Teuchos::Array<double>& cijk_values = setup.Cijk->values(k,l);
	int ni = i_indices.size();
	for (int m=0; m<ni; m++) {
	  int i = i_indices[m];
	  double c = cijk_values[m];

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
