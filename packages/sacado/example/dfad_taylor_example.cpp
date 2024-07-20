// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
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
//           dx/dt = 1 + x^2,    x(0) = x0 = 1.0;
//
//     The exact solution is x(t) = tan(t + atan(x0)) = tan(t + pi/4)
//     Also computes the derivative of the Taylor series solution with
//     respect to x0.  The exact series derivative is 
//     dx/dx0(t) = 1/2 * sec^2(t + atan(x0))

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

  // Compute derivative w.r.t. x0
  Sacado::Fad::DFad< Sacado::Tay::Taylor<double> > x_fad(1, 0, x);
  Sacado::Fad::DFad< Sacado::Tay::Taylor<double> > f_fad;
  func(f_fad, x_fad);
  Sacado::Tay::Taylor<double> dxdx0(deg);
  dxdx0.fastAccessCoeff(0) = 1.0;
  for (int k=0; k<deg; k++) {
    dxdx0.fastAccessCoeff(k+1) = 0.0;
    for (int j=0; j<=k; j++)
      dxdx0.fastAccessCoeff(k+1) += f_fad.dx(0).coeff(k-j) * dxdx0.coeff(j);
    dxdx0.fastAccessCoeff(k+1) /= k+1;
  }

  // Print Taylor series solution
  std::cout << "Taylor series solution = " << std::endl 
	    << x << std::endl;

  // Print Taylor series solution derivative
  std::cout << "Taylor series solution derivative= " << std::endl 
	    << dxdx0 << std::endl;

  // Compute Taylor series expansion of solution x(t) = tan(t+pi/4)
  double pi = std::atan(1.0)*4.0;
  Sacado::Tay::Taylor<double> t(deg);
  t.fastAccessCoeff(0) = pi/4.0;
  t.fastAccessCoeff(1) = 1.0;
  Sacado::Tay::Taylor<double> u = std::tan(t);

  // Compute Taylor series expansion of derivative
  Sacado::Tay::Taylor<double> dudx0 = 0.5*(1.0+u*u);

  // Print expansion of solution
  std::cout << "Exact expansion = " << std::endl
	    << u << std::endl;

  // Print expansion of solution
  std::cout << "Exact expansion derivative = " << std::endl
	    << dudx0 << std::endl;

  // Compute maximum relative error
  double max_err = 0.0;
  double err = 0.0;
  for (int k=0; k<=deg; k++) {
    err = std::fabs(x.coeff(k) - u.coeff(k)) / (1.0 + fabs(u.coeff(k)));
    if (err > max_err) max_err = err;
  }
  std::cout << "Maximum relative error = " << max_err << std::endl;

  // Compute maximum derivative relative error
  double deriv_max_err = 0.0;
  double deriv_err = 0.0;
  for (int k=0; k<=deg; k++) {
    deriv_err = std::fabs(dxdx0.coeff(k) - dudx0.coeff(k)) / 
      (1.0 + fabs(dudx0.coeff(k)));
    if (deriv_err > deriv_max_err) deriv_max_err = deriv_err;
  }
  std::cout << "Maximum derivative relative error = " << deriv_max_err 
	    << std::endl;

  double tol = 1.0e-12;
  if ((max_err < tol) && (deriv_max_err < tol)) {
    std::cout << "\nExample passed!" << std::endl;
    return 0;
  }
  else {
    std::cout <<"\nSomething is wrong, example failed!" << std::endl;
    return 1;
  }

  return 0;
}

  

