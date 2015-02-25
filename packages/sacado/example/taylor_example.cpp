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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

// taylor_example
//
//  usage: 
//     taylor_example
//
//  output:  
//     prints the results of computing a single Taylor series expansion of
//     the solution to:
//
//           dx/dt = 1 + x^2,    x(0) = 1.0;
//
//     The exact solution is x(t) = tan(t + pi/4)

#include <iostream>

#include "Sacado_No_Kokkos.hpp"

// Function implementing RHS of ODE
template <typename ScalarT>
void func(ScalarT& f, const ScalarT& x) {
  f = 1.0 + x*x;
}

int main(int argc, char **argv)
{
  double x0 = 1.0;                      // Initial condition
  int deg = 40;                // Degree of Taylor series solution

  Sacado::Tay::Taylor<double> x = x0;   // Taylor polynomial for independent
  Sacado::Tay::Taylor<double> f;        // Taylor polynomial for dependent

  // Reserve space for degree deg coefficients
  x.reserve(deg);

  // Compute Taylor series solution to dx/dt = f(x)
  for (int k=0; k<deg; k++) {
    func(f, x);

    // Set next coefficient
    x.resize(k+1, true);

    // x_{k+1} = f_k / (k+1)
    x.fastAccessCoeff(k+1) = f.coeff(k) / (k+1);
  }

  // Print Taylor series solution
  std::cout << "Taylor series solution = " << std::endl 
	    << x << std::endl;

  // Compute Taylor series expansion of solution x(t) = tan(t+pi/4)
  double pi = std::atan(1.0)*4.0;
  Sacado::Tay::Taylor<double> t(deg);
  t.fastAccessCoeff(0) = pi/4.0;
  t.fastAccessCoeff(1) = 1.0;
  Sacado::Tay::Taylor<double> u = std::tan(t);

  // Print expansion of solution
  std::cout << "Exact expansion = " << std::endl
	    << u << std::endl;

  // Compute maximum relative error
  double max_err = 0.0;
  double err = 0.0;
  for (int k=0; k<=deg; k++) {
    err = std::fabs(x.coeff(k) - u.coeff(k)) / (1.0 + fabs(u.coeff(k)));
    if (err > max_err) max_err = err;
  }
  std::cout << "Maximum relative error = " << max_err << std::endl;

  double tol = 1.0e-12;
  if (max_err < tol){
    std::cout << "\nExample passed!" << std::endl;
    return 0;
  }
  else {
    std::cout <<"\nSomething is wrong, example failed!" << std::endl;
    return 1;
  }
}

  

