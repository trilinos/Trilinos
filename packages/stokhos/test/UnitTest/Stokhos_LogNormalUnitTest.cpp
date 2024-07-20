// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Test the exponential of a linear Hermite expansion using an analytic
// closed form solution

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Stokhos.hpp"
#include "Stokhos_UnitTestHelpers.hpp"

TEUCHOS_UNIT_TEST( Stokhos_LogNormal_TP, UnitTest ) {
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
    basis->computeTripleProductTensor();

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
    Stokhos::MultiIndex<int> multiIndex = basis->term(i);
    double s = 1.0;
    for (int j=0; j<d; j++)
      s *= std::pow(u.term(j,1), multiIndex[j]);
    w[i] = w_mean*s/basis->norm_squared(i);
  }

  success = Stokhos::comparePCEs(v, "quad_expansion", w, "analytic", 
				 1e-12, 1e-9, out);
}

#ifdef HAVE_STOKHOS_DAKOTA
TEUCHOS_UNIT_TEST( Stokhos_LogNormal_SG, UnitTest ) {
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
    Teuchos::rcp(new Stokhos::SparseGridQuadrature<int,double>(basis, p));
  
  // Triple product tensor
  Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk =
    basis->computeTripleProductTensor();

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
    Stokhos::MultiIndex<int> multiIndex = basis->term(i);
    double s = 1.0;
    for (int j=0; j<d; j++)
      s *= std::pow(u.term(j,1), multiIndex[j]);
    w[i] = w_mean*s/basis->norm_squared(i);
  }

  success = Stokhos::comparePCEs(v, "quad_expansion", w, "analytic", 
				 1e-12, 1e-9, out);
}
#endif

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
