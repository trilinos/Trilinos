// $Id: hermite_example.cpp,v 1.13 2009/09/14 21:53:56 etphipp Exp $ 
// $Source: /space/CVS/Trilinos/packages/stokhos/example/hermite_example.cpp,v $ 
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

#include "Stokhos_KL_ExponentialRandomField.hpp"

int main(int argc, char **argv)
{
  try {

    int M = 10;
    Teuchos::ParameterList solverParams;
    solverParams.set("Number of KL Terms", M);
    solverParams.set("Mean", 1.0);
    solverParams.set("Standard Deviation", 0.1);
    int ndim = 3;
    Teuchos::Array<double> domain_upper(ndim), domain_lower(ndim), 
      correlation_length(ndim);
    for (int i=0; i<ndim; i++) {
      domain_upper[i] = 1.0;
      domain_lower[i] = 0.0;
      correlation_length[i] = 10.0;
    }
    solverParams.set("Domain Upper Bounds", domain_upper);
    solverParams.set("Domain Lower Bounds", domain_lower);
    solverParams.set("Correlation Lengths", correlation_length);
    
    Stokhos::KL::ExponentialRandomField<double> rf(solverParams);
    rf.print(std::cout);

    Teuchos::Array<double> x(ndim);
    for (int i=0; i<ndim; i++)
      x[i] = (domain_upper[i] + domain_lower[i])/2.0 + 
	0.1*(domain_upper[i] - domain_lower[i])/2.0;
    Teuchos::Array<double> rvar(M);
    for (int i=0; i<M; i++)
      rvar[i] = 1.5;
    double result = rf.evaluate(x, rvar);
    std::cout << "result = " << result << std::endl;
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
