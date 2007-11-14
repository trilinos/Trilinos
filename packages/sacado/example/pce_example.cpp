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

// pce_example
//
//  usage: 
//     pce_example
//
//  output:  
//     prints the Hermite Polynomial Chaos Expansion of a simple function
//     using the class Sacado::PCE::Hermite

#include <iostream>
#include <iomanip>
#include <math.h>

#include "Sacado.hpp"
#include "Sacado_PCE_Hermite.hpp"

// The function to differentiate
template <typename ScalarT>
ScalarT func(const ScalarT& a) {
  ScalarT r = std::exp(a);

  return r;
}

int main(int argc, char **argv)
{
  try {
    const unsigned int d = 30;
    
//     Sacado::PCE::HermiteEBasis<double> basis(d);
//     std::cout << basis << std::endl;
//     std::vector<double> a(d+1);
//     for (unsigned int i=0; i<=d; i++) {
//       const Sacado::PCE::StandardPoly<double>& p = basis.getBasisPoly(i);
//       basis.project(p, a);
//       std::cout << i << ":  ";
//       for (unsigned int j=0; j<=d; j++)
// 	std::cout << a[j] << " ";
//       std::cout << std::endl;
//     }

    Sacado::PCE::Hermite<double>::initWorkspace(d);
    Sacado::PCE::Hermite<double> u(d);
    Sacado::PCE::HermiteBasis<double> basis(d);
    std::vector<double> norm_squared = basis.norm_squared();

//     for (unsigned int i=0; i<=5; i++)
//       u.fastAccessCoeff(i) = 1.0*pow(10.0,-double(i+1));

    u.fastAccessCoeff(1) = 1.0;

//     u.fastAccessCoeff(0) = 1.0;
//     u.fastAccessCoeff(1) = 1.0;
//     u.fastAccessCoeff(2) = 0.1;
//     u.fastAccessCoeff(3) = 0.2;
//     u.fastAccessCoeff(4) = -0.3;

//     for (unsigned int i=0; i<=d; i++)
//       u.fastAccessCoeff(i) /= std::sqrt(norm_squared[i]);

//     Sacado::PCE::Hermite<double> v = std::exp(u);
//     Sacado::PCE::Hermite<double> w = std::log(v);

//     for (unsigned int i=0; i<=d; i++) {
//       v.fastAccessCoeff(i) *= std::sqrt(norm_squared[i]);
//       w.fastAccessCoeff(i) *= std::sqrt(norm_squared[i]);
//     }
//     std::cout << "u = " << u << std::endl;
//     std::cout << "v = " << v << std::endl;
//     std::cout << "v*v = " << v*v << std::endl;
//     std::cout << "w = " << w << std::endl;
    
 //    Sacado::PCE::Hermite<double> s = std::sin(u);
//     Sacado::PCE::Hermite<double> c = std::cos(u);
//     Sacado::PCE::Hermite<double> v = (std::exp(u) - std::exp(-u)) / 2.0;
//     Sacado::PCE::Hermite<double> w = (std::exp(u) + std::exp(-u)) / 2.0;
    
//     for (unsigned int i=0; i<=d; i++) {
//       s.fastAccessCoeff(i) *= std::sqrt(norm_squared[i]);
//       c.fastAccessCoeff(i) *= std::sqrt(norm_squared[i]);
//       v.fastAccessCoeff(i) *= std::sqrt(norm_squared[i]);
//       w.fastAccessCoeff(i) *= std::sqrt(norm_squared[i]);
//     }
//     std::cout << "u = " << u << std::endl;
//     std::cout << "s = " << s << std::endl;
//     std::cout << "v-s = " << v-s << std::endl;
//     std::cout << "c = " << c << std::endl;
//     std::cout << "w-c = " << w-c << std::endl;
//     std::cout << "c*c = " << c*c << std::endl;
//     std::cout << "s*s = " << s*s << std::endl;
//     std::cout << "c*c + s*s = " << c*c + s*s << std::endl;

    Sacado::PCE::Hermite<double> c = std::sinh(u);
    Sacado::PCE::Hermite<double> ac = std::asinh(c);
    std::cout << "c = " << c << std::endl;
    std::cout << "u = " << u << std::endl;
    std::cout << "ac = " << ac << std::endl;
    std::cout << "ac - u = " << ac - u << std::endl;

//     Sacado::PCE::Hermite<double> cc = c*c;
//     std::cout << "c*c = " << cc << std::endl;
//     Sacado::PCE::Hermite<double> ss = s*s;
//     std::cout << "s*s = " << ss << std::endl;

//     Sacado::PCE::Hermite<double> c = std::cos(u);
//     Sacado::PCE::Hermite<double> cc = std::cos(c);
//     Sacado::PCE::Hermite<double> sc = std::sin(c);
//     std::cout << c << std::endl;
//     std::cout << cc << std::endl;
//     std::cout << sc << std::endl;
//     std::cout << cc*cc + sc*sc << std::endl;
//     Sacado::PCE::Hermite<double> t = std::sin(c)/std::cos(c);
//     Sacado::PCE::Hermite<double> s = 1.0/std::cos(c);
//     std::cout << "t*t - s*s = " << t*t - s*s << std::endl;
 //    std::cout << "t*t - s*s = " << 
//       (std::sin(c)*std::sin(c)-1.0) - (std::cos(c)*std::cos(c)) << std::endl;
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
