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

// trad_dfad_example
//
//  usage: 
//     trad_dfad_example
//
//  output:  
//     prints the results of computing the second derivative a simple function //     with forward nested forward and reverse mode AD using the 
//     Sacado::Fad::DFad and Sacado::Rad::ADvar classes.

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
  drda = -(c*std::log(b+1.)/std::pow(std::sin(a),2))*std::cos(a);
  drdb = c / ((b+1.)*std::sin(a));
}

// The analytic second derivative of func(a,b,c) with respect to a and b
void func_deriv2(double a, double b, double c, double& d2rda2, double& d2rdb2,
		 double& d2rdadb)
{
  d2rda2 = c*std::log(b+1.)/std::sin(a) + 2.*(c*std::log(b+1.)/std::pow(std::sin(a),3))*std::pow(std::cos(a),2);
  d2rdb2 = -c / (std::pow(b+1.,2)*std::sin(a));
  d2rdadb = -c / ((b+1.)*std::pow(std::sin(a),2))*std::cos(a);
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
  Sacado::Rad::ADvar< Sacado::Fad::DFad<double> > arad = 
    Sacado::Fad::DFad<double>(num_deriv, 0, a);
  Sacado::Rad::ADvar< Sacado::Fad::DFad<double> > brad = 
    Sacado::Fad::DFad<double>(num_deriv, 1, b);
  Sacado::Rad::ADvar< Sacado::Fad::DFad<double> > crad = c;
  Sacado::Rad::ADvar< Sacado::Fad::DFad<double> > rrad;

  // Compute function
  double r = func(a, b, c);

  // Compute derivative analytically
  double drda, drdb;
  func_deriv(a, b, c, drda, drdb);

  // Compute second derivative analytically
  double d2rda2, d2rdb2, d2rdadb;
  func_deriv2(a, b, c, d2rda2, d2rdb2, d2rdadb);

  // Compute function and derivative with AD
  rrad = func(arad, brad, crad);

  Sacado::Rad::ADvar< Sacado::Fad::DFad<double> >::Gradcomp();

  // Extract value and derivatives
  double r_ad = rrad.val().val();       // r
  double drda_ad = arad.adj().val();    // dr/da
  double drdb_ad = brad.adj().val();    // dr/db
  double d2rda2_ad = arad.adj().dx(0);  // d^2r/da^2
  double d2rdadb_ad = arad.adj().dx(1); // d^2r/dadb
  double d2rdbda_ad = brad.adj().dx(0); // d^2r/dbda
  double d2rdb2_ad = brad.adj().dx(1);  // d^2/db^2

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
