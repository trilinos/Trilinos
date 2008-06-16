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

// dfad_example
//
//  usage: 
//     dfad_example
//
//  output:  
//     prints the results of differentiating a simple function with forward
//     mode AD using the Sacado::Fad::DFad class (uses dynamic memory
//     allocation for number of derivative components).

#include <iostream>
#include <iomanip>

#include "Sacado.hpp"

// The function to differentiate
template <typename ScalarT>
ScalarT func(const ScalarT& a, const ScalarT& b, const ScalarT& c) {
  ScalarT r = std::log(b+1.)/std::sin(a);

  return r;
}

typedef Sacado::LFad::LogicalSparse<double,bool> FadType;

int main(int argc, char **argv)
{
  double pi = std::atan(1.0)*4.0;

  // Values of function arguments
  double a = pi/4;
  double b = 2.0;
  double c = 3.0;

  // Number of independent variables
  int num_deriv = 3;

  // Fad objects
  FadType afad(num_deriv, 0, a); // First (0) indep. var
  FadType bfad(num_deriv, 1, b); // Second (1) indep. var
  FadType cfad(num_deriv, 2, c); // Third (2) indep. var
  FadType rfad;                  // Result

  // Compute function
  double r = func(a, b, c);

  // Compute function and derivative with AD
  rfad = func(afad, bfad, cfad);

  std::cout << rfad << std::endl;

  // Extract value and derivatives
  double r_ad = rfad.val();     // r
  bool drda_ad = rfad.dx(0);  // dr/da
  bool drdb_ad = rfad.dx(1);  // dr/db
  bool drdc_ad = rfad.dx(2);  // dr/dc

  double tol = 1.0e-14;
  if (std::fabs(r - r_ad) < tol && drda_ad && drdb_ad && !drdc_ad) {
    std::cout << "\nExample passed!" << std::endl;
    return 0;
  }
  else {
    std::cout <<"\nSomething is wrong, example failed!" << std::endl;
    return 1;
  }
}
