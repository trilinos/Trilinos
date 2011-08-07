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

// Test the exponential of a linear Hermite expansion using an analytic
// closed form solution

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Stokhos.hpp"
#include "Stokhos_UnitTestHelpers.hpp"

TEUCHOS_UNIT_TEST( Stokhos_LogNormal, UnitTest ) {
  // Basis of dimension 3, order 5
  const int d = 3;
  const int p = 5;
  const double mean = 0.2;
  const double std_dev = 0.1;
  Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d); 
  for (int i=0; i<d; i++) {
    bases[i] = Teuchos::rcp(new Stokhos::HermiteBasis<int,double>(p));
  }
  Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
    Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
  
  // Quadrature method
  Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad = 
    Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
  
  // Triple product tensor
  Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk =
    basis->computeTripleProductTensor(basis->size());

  // Expansion method
  Stokhos::QuadOrthogPolyExpansion<int,double> expn(basis, Cijk, quad);

  // Polynomial expansions
  Stokhos::OrthogPolyApprox<int,double> u(basis), v(basis), w(basis);
  u.term(0,0) = mean;
  for (int i=0; i<d; i++)
    u.term(i,1) = std_dev / (i+1);
  
  // Compute expansion
  expn.exp(v,u);
   
  // Compute expansion analytically
  double w_mean = mean;
  for (int j=0; j<d; j++)
    w_mean += u.term(j,1)*u.term(j,1)/2.0;
  w_mean = std::exp(w_mean);
  for (int i=0; i<basis->size(); i++) {
    Teuchos::Array<int> multiIndex = basis->getTerm(i);
    double s = 1.0;
    for (int j=0; j<d; j++)
      s *= std::pow(u.term(j,1), multiIndex[j]);
    w[i] = w_mean*s/basis->norm_squared(i);
  }

  success = Stokhos::comparePCEs(v, "quad_expansion", w, "analytic", 
				 1e-12, 1e-9, out);
}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
