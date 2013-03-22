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

// hermite_example
//
//  usage: 
//     hermite_example
//
//  output:  
//     prints the Hermite Polynomial Chaos Expansion of the simple function
//
//     v = 1/(log(u)^2+1)
//
//     where u = 1 + 0.4*H_1(x) + 0.06*H_2(x) + 0.002*H_3(x), x is a zero-mean
//     and unit-variance Gaussian random variable, and H_i(x) is the i-th
//     Hermite polynomial.

#include "Stokhos.hpp"
#include "Teuchos_GlobalMPISession.hpp"

int main(int argc, char **argv)
{
  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // If applicable, set up MPI.
  Teuchos::GlobalMPISession mpiSession (&argc, &argv);

  try {

    // Basis of dimension 3, order 5
    const int d = 3;
    const int p = 5;
    Array< RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d); 
    for (int i=0; i<d; i++) {
      bases[i] = rcp(new Stokhos::HermiteBasis<int,double>(p-i));
    }
    RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
      rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));

    // Quadrature method
    RCP<const Stokhos::Quadrature<int,double> > quad = 
      rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));

    // Triple product tensor
    RCP<Stokhos::Sparse3Tensor<int,double> > Cijk =
      basis->computeTripleProductTensor();

    // Expansion method
    Stokhos::QuadOrthogPolyExpansion<int,double> expn(basis, Cijk, quad);

    // Polynomial expansions
    Stokhos::OrthogPolyApprox<int,double> u(basis), v(basis), w(basis);
    u.term(0,0) = 1.0;
    for (int i=0; i<d; i++) {
      if (bases[i]->order() >= 1)
	u.term(i,1) = 0.4 / d;
      if (bases[i]->order() >= 2)
	u.term(i,2) = 0.06 / d;
      if (bases[i]->order() >= 3)
	u.term(i,3) = 0.002 / d;
    }

    // Compute expansion
    expn.log(v,u);
    expn.times(w,v,v);
    expn.plusEqual(w,1.0);
    expn.divide(v,1.0,w);

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
    Array<double> pt(d); 
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

    Teuchos::TimeMonitor::summarize(std::cout);
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
