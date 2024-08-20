// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// dfad_sfc_example
//
//  usage:
//     dfad_sfc_example
//
//  output:
//     Uses the scalar flop counter to count the flops for a derivative
//     of a simple function using DFad

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
template <typename ScalarT>
void func_deriv(const ScalarT& a, const ScalarT& b, const ScalarT& c,
                ScalarT& drda, ScalarT& drdb)
{
  drda = -(c*std::log(b+1.)/std::pow(std::sin(a),2.))*std::cos(a);
  drdb = c / ((b+1.)*std::sin(a));
}

typedef Sacado::FlopCounterPack::ScalarFlopCounter<double> SFC;
typedef Sacado::Fad::DFad<SFC> FAD_SFC;

int main(int argc, char **argv)
{
  double pi = std::atan(1.0)*4.0;

  // Values of function arguments
  double a = pi/4;
  double b = 2.0;
  double c = 3.0;

  // Number of independent variables
  int num_deriv = 2;

  // Compute function
  SFC as(a);
  SFC bs(b);
  SFC cs(c);
  SFC::resetCounters();
  SFC rs = func(as, bs, cs);
  SFC::finalizeCounters();

  std::cout << "Flop counts for function evaluation:";
  SFC::printCounters(std::cout);

  // Compute derivative analytically
  SFC drdas, drdbs;
  SFC::resetCounters();
  func_deriv(as, bs, cs, drdas, drdbs);
  SFC::finalizeCounters();

  std::cout << "\nFlop counts for analytic derivative evaluation:";
  SFC::printCounters(std::cout);

  // Compute function and derivative with AD
  FAD_SFC afad(num_deriv, 0, a);
  FAD_SFC bfad(num_deriv, 1, b);
  FAD_SFC cfad(c);
  SFC::resetCounters();
  FAD_SFC rfad = func(afad, bfad, cfad);
  SFC::finalizeCounters();

  std::cout << "\nFlop counts for AD function and derivative evaluation:";
  SFC::printCounters(std::cout);

  // Extract value and derivatives
  double r = rs.val();               // r
  double drda = drdas.val();         // dr/da
  double drdb = drdbs.val();         // dr/db

  double r_ad = rfad.val().val();     // r
  double drda_ad = rfad.dx(0).val();  // dr/da
  double drdb_ad = rfad.dx(1).val();  // dr/db

  // Print the results
  int p = 4;
  int w = p+7;
  std::cout.setf(std::ios::scientific);
  std::cout.precision(p);
  std::cout << "\nValues/derivatives of computation" << std::endl
            << "    r =  " << r << " (original) == " << std::setw(w) << r_ad
            << " (AD) Error = " << std::setw(w) << r - r_ad << std::endl
            << "dr/da = " << std::setw(w) << drda << " (analytic) == "
            << std::setw(w) << drda_ad << " (AD) Error = " << std::setw(w)
            << drda - drda_ad << std::endl
            << "dr/db = " << std::setw(w) << drdb << " (analytic) == "
            << std::setw(w) << drdb_ad << " (AD) Error = " << std::setw(w)
            << drdb - drdb_ad << std::endl;

  double tol = 1.0e-14;
  Sacado::FlopCounterPack::FlopCounts fc = SFC::getCounters();
  // The Solaris and Irix CC compilers get higher counts for operator=
  // than does g++, which avoids an extra copy when returning a function value.
  // The test on fc.totalFlopCount allows for this variation.
  if (std::fabs(r - r_ad)       < tol &&
      std::fabs(drda - drda_ad) < tol &&
      std::fabs(drdb - drdb_ad) < tol &&
      (fc.totalFlopCount == 40)) {
    std::cout << "\nExample passed!" << std::endl;
    return 0;
  }
  else {
    std::cout <<"\nSomething is wrong, example failed!" << std::endl;
    return 1;
  }
}
