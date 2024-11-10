// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// pecos_hermite_example
//
//  usage: 
//     pecos_hermite_example
//
//  output:  
//     prints the Hermite Polynomial Chaos Expansion of the simple function
//
//     v = 1/(log(u)^2+1)
//
//     where u = 1 + 0.4*H_1(x) + 0.06*H_2(x) + 0.002*H_3(x), x is a zero-mean
//     and unit-variance Gaussian random variable, and H_i(x) is the i-th
//     Hermite polynomial.
//
//     Same as hermite_example, except uses Pecos to define the Hermite basis.

#include "Stokhos.hpp"
#include "HermiteOrthogPolynomial.hpp" // from Pecos

int main(int argc, char **argv)
{
  try {

    // Basis of dimension 3, order 5
    const int d = 3;
    const int p = 5;
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d); 
    for (int i=0; i<d; i++) {
      bases[i] = Teuchos::rcp(new Stokhos::PecosOneDOrthogPolyBasis<int,double>(Teuchos::rcp(new Pecos::HermiteOrthogPolynomial), "Hermite", p));
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
    u.term(0,0) = 1.0;
    for (int i=0; i<d; i++) {
      u.term(i,1) = 0.4 / d;
      u.term(i,2) = 0.06 / d;
      u.term(i,3) = 0.002 / d;
    }

    // Compute expansion
    expn.log(v,u);
    expn.times(w,v,v);
    expn.plusEqual(w,1.0);
    expn.divide(v,1.0,w);
    //expn.times(v,u,u);

    // Print u and v
    std::cout << "v = 1.0 / (log(u)^2 + 1):" << std::endl;
    std::cout << "\tu = ";
    u.print(std::cout);
    std::cout << "\tv = ";
    v.print(std::cout);

    // Compute moments
    double mean = v.mean();
    double std_dev = v.standard_deviation();

    // Evaluate PCE and function at a point = 0.25 in each dimension
    Teuchos::Array<double> pt(d); 
    for (int i=0; i<d; i++) 
      pt[i] = 0.25;
    double up = u.evaluate(pt);
    double vp = 1.0/(std::log(up)*std::log(up)+1.0);
    double vp2 = v.evaluate(pt);
    
    // Print results
    std::cout << "\tv mean         = " << mean << std::endl;
    std::cout << "\tv std. dev.    = " << std_dev << std::endl;
    std::cout << "\tv(0.25) (true) = " << vp << std::endl;
    std::cout << "\tv(0.25) (pce)  = " << vp2 << std::endl;
    
    // Check the answer
    if (std::abs(vp - vp2) < 1e-2)
      std::cout << "\nExample Passed!" << std::endl;
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
