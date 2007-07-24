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

// dmfad_example
//
//  usage: 
//     dmfad_example
//
//  output:  
//     prints the results of differentiating a simple function with forward
//     mode AD using the Sacado::Fad::DMFad class (uses dynamic memory
//     allocation for number of derivative components using a custom memory
//     manager).

#include <iostream>
#include <iomanip>

#include "Sacado.hpp"

template <>
Sacado::Fad::MemPool* Sacado::Fad::MemPoolStorage<double>::defaultPool_ = NULL;

// The function to differentiate
template <typename ScalarT>
ScalarT func(const ScalarT& a, const ScalarT& b, const ScalarT& c) {
  ScalarT r = c*std::log(b+1.)/std::sin(a);

  return r;
}

// The analytic derivative of func(a,b,c) with respect to a and b
void func_deriv(double a, double b, double c, double& drda, double& drdb)
{
  drda = -(c*std::log(b+1.)/std::pow(std::sin(a),2))*std::cos(a);
  drdb = c / ((b+1.)*std::sin(a));
}

int main(int argc, char **argv)
{
  double pi = std::atan(1.0)*4.0;

  // Values of function arguments
  double a = pi/4;
  double b = 2.0;
  double c = 3.0;

  // Number of independent variables
  int num_deriv = 2;

  // Memory pool & manager
  Sacado::Fad::MemPoolManager<double> poolManager(10);
  Sacado::Fad::MemPool* pool = poolManager.getMemoryPool(num_deriv);
  Sacado::Fad::DMFad<double>::setDefaultPool(pool);

  // Fad objects
  Sacado::Fad::DMFad<double> afad(num_deriv, 0, a); // First (0) indep. var
  Sacado::Fad::DMFad<double> bfad(num_deriv, 1, b); // Second (1) indep. var
  Sacado::Fad::DMFad<double> cfad(c);               // Passive variable
  Sacado::Fad::DMFad<double> rfad;                  // Result

  // Compute function
  double r = func(a, b, c);

  // Compute derivative analytically
  double drda, drdb;
  func_deriv(a, b, c, drda, drdb);

  // Compute function and derivative with AD
  rfad = func(afad, bfad, cfad);

  // Extract value and derivatives
  double r_ad = rfad.val();     // r
  double drda_ad = rfad.dx(0);  // dr/da
  double drdb_ad = rfad.dx(1);  // dr/db

  // Print the results
  int p = 4;
  int w = p+7;
  std::cout.setf(std::ios::scientific);
  std::cout.precision(p);
  std::cout << "    r =  " << r << " (original) == " << std::setw(w) << r_ad
	    << " (AD) Error = " << std::setw(w) << r - r_ad << std::endl
	    << "dr/da = " << std::setw(w) << drda << " (analytic) == " 
	    << std::setw(w) << drda_ad << " (AD) Error = " << std::setw(w) 
	    << drda - drda_ad << std::endl
	    << "dr/db = " << std::setw(w) << drdb << " (analytic) == " 
	    << std::setw(w) << drdb_ad << " (AD) Error = " << std::setw(w) 
	    << drdb - drdb_ad << std::endl;

  double tol = 1.0e-14;
  if (std::fabs(r - r_ad)       < tol &&
      std::fabs(drda - drda_ad) < tol &&
      std::fabs(drdb - drdb_ad) < tol) {
    std::cout << "\nExample passed!" << std::endl;
    return 0;
  }
  else {
    std::cout <<"\nSomething is wrong, example failed!" << std::endl;
    return 1;
  }
}
