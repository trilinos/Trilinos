// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// dfad_sfad_example
//
//  usage:
//     dfad_sfad_example
//
//  output:
//     prints the results of computing the second derivative times a vector
//     for a simple function with forward nested forward mode AD using the
//     Sacado::Fad::DFad and Sacado::Fad::SFad classes.

#include <iostream>
#include <iomanip>

#include "Sacado.hpp"

// The function to differentiate
template <typename ScalarT>
ScalarT func(const ScalarT& a, const ScalarT& b, const ScalarT& c) {
  ScalarT r = c*std::log(b+1.)/std::sin(a);
  return r;
}

// The analytic first and second derivative of func with respect to a and b
void analytic_deriv(double a, double b, double c,
                    double& drda, double& drdb,
                    double& d2rda2, double& d2rdb2, double& d2rdadb)
{
  drda = -(c*std::log(b+1.)/std::pow(std::sin(a),2.))*std::cos(a);
  drdb = c / ((b+1.)*std::sin(a));
  d2rda2 = c*std::log(b+1.)/std::sin(a) + 2.*(c*std::log(b+1.)/std::pow(std::sin(a),3.))*std::pow(std::cos(a),2.);
  d2rdb2 = -c / (std::pow(b+1.,2.)*std::sin(a));
  d2rdadb = -c / ((b+1.)*std::pow(std::sin(a),2.))*std::cos(a);
}

// Function that computes func and its first derivative w.r.t a & b using
// Sacado AD
template <typename ScalarT>
void func_and_deriv(const ScalarT& a, const ScalarT& b, const ScalarT& c,
                    ScalarT& r, ScalarT& drda, ScalarT& drdb) {
  typedef Sacado::Fad::DFad<ScalarT> FadType;
  FadType a_fad(2, 0, a);
  FadType b_fad(2, 1, b);
  FadType c_fad = c;

  FadType r_fad = func(a_fad, b_fad, c_fad);
  r = r_fad.val();
  drda = r_fad.dx(0);
  drdb = r_fad.dx(1);
}

// Function that computes func, its first derivative w.r.t a & b, and its
// second derivative in the direction of [v_a, v_b] with Sacado AD
//
// Define x = [a, b], v = [v_a, v_b], and y(t) = x + t*v.  Then
// df/dx*v = d/dt f(y(t)) |_{t=0}.
//
// In the code below, we differentiate with respect to t in this manner.
// Addtionally we take a short-cut and don't introduce t directly and
// compute a(t) = a + t*v_a, b(t) = b + t*v_b.  Instead we
// initialize a_fad and b_fad directly as if we had computed them in this way.
template <typename ScalarT>
void func_and_deriv2(const ScalarT& a, const ScalarT& b, const ScalarT& c,
                     const ScalarT& v_a, const ScalarT& v_b,
                     ScalarT& r, ScalarT& drda, ScalarT& drdb,
                     ScalarT& z_a, ScalarT& z_b) {
  typedef Sacado::Fad::SFad<ScalarT,1> FadType;

  // The below is equivalent to:
  // FadType t(1, 0.0); f_fad.fastAccessDx(0) = 1;
  // FadType a_fad = a + t*v_a;
  // FadType b_fad = b + t*v_b;
  FadType a_fad(1, a); a_fad.fastAccessDx(0) = v_a;
  FadType b_fad(1, b); b_fad.fastAccessDx(0) = v_b;
  FadType c_fad = c;

  FadType r_fad, drda_fad, drdb_fad;
  func_and_deriv(a_fad, b_fad, c_fad, r_fad, drda_fad, drdb_fad);
  r = r_fad.val();       // r
                         // note:  also have r_fad.dx(0) = dr/da*v_a + dr/db*v_b
  drda = drda_fad.val(); // dr/da
  drdb = drdb_fad.val(); // dr/db
  z_a = drda_fad.dx(0);  // d^2r/da^2 * v_a + d^2r/dadb * v_b
  z_b = drdb_fad.dx(0);  // d^2r/dadb * v_a + d^2r/db^2 * v_b
}

int main(int argc, char **argv)
{
  double pi = std::atan(1.0)*4.0;

  // Values of function arguments
  double a = pi/4;
  double b = 2.0;
  double c = 3.0;

  // Direction we wish to differentiate for second derivative
  double v_a = 1.5;
  double v_b = 3.6;

  // Compute derivatives via AD
  double r_ad, drda_ad, drdb_ad, z_a_ad, z_b_ad;
  func_and_deriv2(a, b, c, v_a, v_b, r_ad, drda_ad, drdb_ad, z_a_ad, z_b_ad);

  // Compute function
  double r = func(a, b, c);

  // Compute derivatives analytically
  double drda, drdb, d2rda2, d2rdb2, d2rdadb;
  analytic_deriv(a, b, c, drda, drdb, d2rda2, d2rdb2, d2rdadb);
  double z_a = d2rda2*v_a + d2rdadb*v_b;
  double z_b = d2rdadb*v_a + d2rdb2*v_b;

  // Print the results
  int p = 4;
  int w = p+7;
  std::cout.setf(std::ios::scientific);
  std::cout.precision(p);
  std::cout << "    r = " << std::setw(w) << r << " (original) == "
            << std::setw(w) << r_ad << " (AD) Error = " << std::setw(w)
            << r - r_ad << std::endl
            << "dr/da = " << std::setw(w) << drda << " (analytic) == "
            << std::setw(w) << drda_ad << " (AD) Error = " << std::setw(w)
            << drda - drda_ad << std::endl
            << "dr/db = " << std::setw(w) << drdb << " (analytic) == "
            << std::setw(w) << drdb_ad << " (AD) Error = " << std::setw(w)
            << drdb - drdb_ad << std::endl
            << "z_a   = " << std::setw(w) << z_a << " (analytic) == "
            << std::setw(w) << z_a_ad << " (AD) Error = " << std::setw(w)
            << z_a - z_a_ad << std::endl
            << "z_b   = " << std::setw(w) << z_b << " (analytic) == "
            << std::setw(w) << z_b_ad << " (AD) Error = " << std::setw(w)
            << z_b - z_b_ad << std::endl;

  double tol = 1.0e-14;
  if (std::fabs(r    - r_ad)     < tol &&
      std::fabs(drda - drda_ad)  < tol &&
      std::fabs(drdb - drdb_ad)  < tol &&
      std::fabs(z_a  - z_a_ad)   < tol &&
      std::fabs(z_b  - z_b_ad)   < tol) {
    std::cout << "\nExample passed!" << std::endl;
    return 0;
  }
  else {
    std::cout <<"\nSomething is wrong, example failed!" << std::endl;
    return 1;
  }
}
