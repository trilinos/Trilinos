// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// dfad_example
//
//  usage:
//     dfad_view_handle_example
//
//  output:
//     prints the results of differentiating a simple function with forward
//     mode AD using the Sacado::Fad::DFad class (uses dynamic memory
//     allocation for number of derivative components) and ViewFad as a
//     handle into externally stored derivative data

#include <iostream>
#include <iomanip>

#include "Sacado.hpp"

// The function to differentiate
template <typename ScalarRes, typename Scalar1, typename Scalar2>
ScalarRes func(const Scalar1& a, const Scalar1& b, const Scalar2& c) {
  ScalarRes r = c*std::log(b+1.)/std::sin(a);

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
  Kokkos::initialize();
  int ret = 0;
  {

  double pi = std::atan(1.0)*4.0;

  // Values of function arguments
  double a = pi/4;
  double b = 2.0;
  double c = 3.0;

  // View to store derivative data
  const int num_deriv = 2;
  Kokkos::View<double**,Kokkos::LayoutLeft,Kokkos::HostSpace> v( "v", 2, num_deriv );

  // Initialize derivative data
  Kokkos::deep_copy( v, 0.0 );
  v(0,0) = 1.0; // First (0) indep. var
  v(1,1) = 1.0; // Second (1) indep. var

  // The Fad type
  typedef Sacado::Fad::DFad<double> FadType;

  // View handle type -- first 0 is static length (e.g., SFad), second 0
  // is static stride, which you can make 1 if you know the View will be
  // LayoutRight (e.g., not GPU).  When values are 0, they are treated
  // dynamically
  typedef Sacado::Fad::ViewFad<double,0,0,FadType> ViewFadType;

  // Fad objects
  ViewFadType afad( &v(0,0), &a, num_deriv, v.stride_1() );
  ViewFadType bfad( &v(1,0), &b, num_deriv, v.stride_1() );
  FadType cfad(c);
  FadType rfad;

  // Compute function
  double r = func<double>(a, b, c);

  // Compute derivative analytically
  double drda, drdb;
  func_deriv(a, b, c, drda, drdb);

  // Compute function and derivative with AD
  rfad = func<FadType>(afad, bfad, cfad);

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
    ret = 0;
  }
  else {
    std::cout <<"\nSomething is wrong, example failed!" << std::endl;
    ret = 1;
  }

  }
  Kokkos::finalize();
  return ret;
}
