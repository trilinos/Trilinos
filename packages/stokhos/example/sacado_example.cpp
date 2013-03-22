// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

// pce_example
//
//  usage: 
//     pce_example
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
//     This example is the same as 
//             Trilinos/packages/stokhos/example/hermite_example.cpp
//     using the Sacado overloaded operators.

#include "Stokhos_Sacado.hpp"

// The function to compute the polynomial chaos expansion of,
// written as a template function
template <class ScalarType>
inline ScalarType simple_function(const ScalarType& u) {
  return 1.0/(std::pow(std::log(u),2.0) + 1.0);
}

// Typename of PC expansion type
typedef Sacado::ETPCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > pce_type;

int main(int argc, char **argv)
{
  try {

    // Basis of dimension 3, order 5
    const int d = 3;
    const int p = 5;
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
    Teuchos::RCP<Stokhos::QuadOrthogPolyExpansion<int,double> > expn = 
      Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<int,double>(
		     basis, Cijk, quad));

    // Polynomial expansions
    pce_type u(expn);
    u.term(0,0) = 1.0;
    for (int i=0; i<d; i++) {
      u.term(i,1) = 0.4 / d;
      u.term(i,2) = 0.06 / d;
      u.term(i,3) = 0.002 / d;
    }

    // Compute PCE expansion of function
    pce_type v = simple_function(u);

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
    double vp = simple_function(up);
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
