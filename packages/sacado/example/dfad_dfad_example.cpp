// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// dfad_dfad_example
//
//  usage:
//     dfad_dfad_example
//
//  output:
//     prints the results of computing the second derivative a simple function //     with forward nested forward mode AD using the Sacado::Fad::DFad class
//     (uses dynamic memory allocation for number of derivative components).

#include <iostream>
#include <iomanip>

#include "Sacado.hpp"

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

// The analytic second derivative of func(a,b,c) with respect to a and b
void func_deriv2(double a, double b, double c, double& d2rda2, double& d2rdb2,
                 double& d2rdadb)
{
  d2rda2 = c*std::log(b+1.)/std::sin(a) + 2.*(c*std::log(b+1.)/std::pow(std::sin(a),3.))*std::pow(std::cos(a),2.);
  d2rdb2 = -c / (std::pow(b+1.,2.)*std::sin(a));
  d2rdadb = -c / ((b+1.)*std::pow(std::sin(a),2.))*std::cos(a);
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

  // Fad objects
  typedef Sacado::Fad::DFad<double> DFadType;
  Sacado::Fad::DFad<DFadType> afad(num_deriv, 0, a);
  Sacado::Fad::DFad<DFadType> bfad(num_deriv, 1, b);
  Sacado::Fad::DFad<DFadType> cfad = c;
  Sacado::Fad::DFad<DFadType> rfad;

  afad.val() = Sacado::Fad::DFad<double>(num_deriv, 0, a);
  bfad.val() = Sacado::Fad::DFad<double>(num_deriv, 1, b);

  // Compute function
  double r = func(a, b, c);

  // Compute derivative analytically
  double drda, drdb;
  func_deriv(a, b, c, drda, drdb);

  // Compute second derivative analytically
  double d2rda2, d2rdb2, d2rdadb;
  func_deriv2(a, b, c, d2rda2, d2rdb2, d2rdadb);

  // Compute function and derivative with AD
  rfad = func(afad, bfad, cfad);

  // Extract value and derivatives
  double r_ad = rfad.val().val();       // r
  double drda_ad = rfad.dx(0).val();    // dr/da
  double drdb_ad = rfad.dx(1).val();    // dr/db
  double d2rda2_ad = rfad.dx(0).dx(0);  // d^2r/da^2
  double d2rdadb_ad = rfad.dx(0).dx(1); // d^2r/dadb
  double d2rdbda_ad = rfad.dx(1).dx(0); // d^2r/dbda
  double d2rdb2_ad = rfad.dx(1).dx(1);  // d^2/db^2

  // Print the results
  int p = 4;
  int w = p+7;
  std::cout.setf(std::ios::scientific);
  std::cout.precision(p);
  std::cout << "        r = " << std::setw(w) << r << " (original) == "
            << std::setw(w) << r_ad << " (AD) Error = " << std::setw(w)
            << r - r_ad << std::endl
            << "    dr/da = " << std::setw(w) << drda << " (analytic) == "
            << std::setw(w) << drda_ad << " (AD) Error = " << std::setw(w)
            << drda - drda_ad << std::endl
            << "    dr/db = " << std::setw(w) << drdb << " (analytic) == "
            << std::setw(w) << drdb_ad << " (AD) Error = " << std::setw(w)
            << drdb - drdb_ad << std::endl
            << "d^2r/da^2 = " << std::setw(w) << d2rda2 << " (analytic) == "
            << std::setw(w) << d2rda2_ad << " (AD) Error = " << std::setw(w)
            << d2rda2 - d2rda2_ad << std::endl
            << "d^2r/db^2 = " << std::setw(w) << d2rdb2 << " (analytic) == "
            << std::setw(w) << d2rdb2_ad << " (AD) Error = " << std::setw(w)
            << d2rdb2 - d2rdb2_ad << std::endl
            << "d^2r/dadb = " << std::setw(w) << d2rdadb << " (analytic) == "
            << std::setw(w) << d2rdadb_ad << " (AD) Error = " << std::setw(w)
            << d2rdadb - d2rdadb_ad << std::endl
            << "d^2r/dbda = " << std::setw(w) << d2rdadb << " (analytic) == "
            << std::setw(w) << d2rdbda_ad << " (AD) Error = " << std::setw(w)
            << d2rdadb - d2rdbda_ad << std::endl;

  double tol = 1.0e-14;
  if (std::fabs(r - r_ad)             < tol &&
      std::fabs(drda - drda_ad)       < tol &&
      std::fabs(drdb - drdb_ad)       < tol &&
      std::fabs(d2rda2 - d2rda2_ad)   < tol &&
      std::fabs(d2rdb2 - d2rdb2_ad)   < tol &&
      std::fabs(d2rdadb - d2rdadb_ad) < tol) {
    std::cout << "\nExample passed!" << std::endl;
    return 0;
  }
  else {
    std::cout <<"\nSomething is wrong, example failed!" << std::endl;
    return 1;
  }
}
