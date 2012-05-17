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
//     prints the Hermite Polynomial Chaos Expansion of the simple function
//
//     v = 1/(log(u)^2+1)
//
//     where u = 1 + 0.4*H_1(x) + 0.06*H_2(x) + 0.002*H_3(x), x is a zero-mean
//     and unit-variance Gaussian random variable, and H_i(x) is the i-th
//     Hermite polynomial.
//
//     This example is the same as 
//             Trilinos/packages/stokhos/example/hermite_example.cpp
//     using the Sacado overloaded operators.

#include "stdio.h"
#include "Stokhos_Sacado_Kokkos.hpp"
#include "Teuchos_TimeMonitor.hpp"

// The function to compute the polynomial chaos expansion of,
// written as a template function
template <class ScalarType>
__host__ __device__
ScalarType simple_function(const ScalarType& u) {
  return 1.0/(std::pow(std::log(u),2.0) + 1.0);
}

// Typename of PC expansion type
typedef Kokkos::Cuda node_type;
typedef Sacado::MP::Vector<Stokhos::StaticFixedStorage<int,double,10,node_type>, node_type> vec_type;

__global__ void kernel() 
{
  // Polynomial expansions
  const int sz = 10;
  vec_type u(sz, 0.0), v(sz, 0.0);
  for (int i=0; i<sz; i++) {
    u.fastAccessCoeff(i) = 0.1 * (i+1);
  }

  v = simple_function(u);

  // Print u and v
  printf("u = [ ");
  for (int i=0; i<sz; i++)
    printf("%g ", u.coeff(i));
  printf("]\n");

  printf("v = [ ");
  for (int i=0; i<sz; i++)
    printf("%g ", v.coeff(i));
  printf("]\n");

  // std::cout << "v = 1.0 / (log(u)^2 + 1):" << std::endl;
  // std::cout << "\tu = " << u << std::endl;
  // std::cout << "\tv = " << v << std::endl;
  
}

int main(int argc, char **argv)
{
  try {

    {
      TEUCHOS_FUNC_TIME_MONITOR("Device PCE calculation");
      kernel<<<1,2>>>();
      cudaDeviceReset();
    }

  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
