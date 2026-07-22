// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// trad_example
//
//  usage: 
//     trad_example
//
//  output:  
//     prints the results of differentiating a simple function with reverse
//     mode AD using the Sacado::Rad::ADvar class.

#include <iostream>
#include <iomanip>

#include "Sacado_No_Kokkos.hpp"

// The function to differentiate
template <typename ScalarT>
ScalarT func(const ScalarT& a, const ScalarT& b, const ScalarT& c) {
  ScalarT r = c*std::log(b+1.)/std::sin(a);

  return r;
}

// The analytic derivative of func(a,b,c) with respect to a and b
void func_deriv(double a, double b, double c, double& drda, double& drdb)
{
  drda = -(c*std::log(b+1.)/std::pow(std::sin(a),2.))*std::cos(a);
  drdb = c / ((b+1.)*std::sin(a));
}

int main(int argc, char **argv)
{
  double pi = std::atan(1.0)*4.0;

  // Values of function arguments
  double a = pi/4;
  double b = 2.0;
  double c = 3.0;

  // Rad objects
  Sacado::Rad::ADvar<double> arad = a; 
  Sacado::Rad::ADvar<double> brad = b; 
  Sacado::Rad::ADvar<double> crad = c;              // Passive variable
  Sacado::Rad::ADvar<double> rrad;                  // Result

  // Compute function
  double r = func(a, b, c);

  // Compute derivative analytically
  double drda, drdb;
  func_deriv(a, b, c, drda, drdb);

  // Compute function and derivative with AD
  rrad = func(arad, brad, crad);

  Sacado::Rad::ADvar<double>::Gradcomp();

  // Extract value and derivatives
  double r_ad = rrad.val();     // r
  double drda_ad = arad.adj();  // dr/da
  double drdb_ad = brad.adj();  // dr/db

  // Free Rad's memory to avoid memory leaks
  Sacado::Rad::ADcontext< double >::free_all();

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
