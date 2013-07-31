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

#include <iostream>

#define KERNEL_PREFIX __device__ __host__

#include "Stokhos_Sacado.hpp"
#include "Stokhos_CUDAStorage.hpp"
#include "Stokhos_CUDAQuadOrthogPolyExpansion.hpp"
#include "Teuchos_TimeMonitor.hpp"

template <> 
Teuchos::RCP<Sacado::ETPCE::OrthogPolyImpl<float, Stokhos::CUDAStorage<int,float> >::expansion_type> Sacado::ETPCE::OrthogPolyImpl<float, Stokhos::CUDAStorage<int,float> >::const_expansion_ = Teuchos::null;

// The function to compute the polynomial chaos expansion of,
// written as a template function.
// Currently the constants don't work when run on a GPU
template <class ScalarType>
ScalarType simple_function(const ScalarType& u) {
  //return 1.0/(std::pow(std::log(u),float(2.0)) + 1.0);
  return ScalarType(1.0)/(std::pow(std::log(u),ScalarType(2.0)) + ScalarType(1.0));
}

typedef Stokhos::StandardStorage<int,float> StdStorage;
typedef Stokhos::CUDAStorage<int,float> Storage;
//typedef Stokhos::StandardStorage<int,float> Storage;
typedef Sacado::ETPCE::OrthogPoly<float,StdStorage> pce_type;
typedef Sacado::ETPCE::OrthogPoly<float,Storage> dev_pce_type;

int main(int argc, char **argv)
{
  try {

    // Basis of dimension 3, order 5
    const int d = 3;
    const int p = 5;
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,float> > > bases(d); 
    for (int i=0; i<d; i++) {
      bases[i] = Teuchos::rcp(new Stokhos::HermiteBasis<int,float>(p));
    }
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,float> > basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,float>(bases));

    // Quadrature method
    Teuchos::RCP<const Stokhos::Quadrature<int,float> > quad = 
        Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,float>(basis));

    // Triple product tensor
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,float> > Cijk =
      basis->computeTripleProductTensor();

    // Expansion method
    Teuchos::RCP<Teuchos::ParameterList> params =
      Teuchos::rcp(new Teuchos::ParameterList);
    params->set("Use Quadrature for Times", true);
    Teuchos::RCP<Stokhos::QuadOrthogPolyExpansion<int,float,StdStorage> > expn =
      Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<int,float,StdStorage>(
		     basis, Cijk, quad, params));

    // Polynomial expansions
    pce_type u(expn), v(expn);
    u.term(0,0) = 1.0;
    for (int i=0; i<d; i++) {
      u.term(i,1) = 0.4 / d;
      u.term(i,2) = 0.06 / d;
      u.term(i,3) = 0.002 / d;
    }

    // Compute expansion
    {
      Teuchos::RCP<Stokhos::QuadOrthogPolyExpansion<int,float,Storage> > 
	expn_dev =
	Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<int,float,Storage>(
		       basis, Cijk, quad, true));
      dev_pce_type u_dev(expn_dev), v_dev(expn_dev);

      {
	TEUCHOS_FUNC_TIME_MONITOR("Host->Device transfer");
	u_dev.init(u);
      }

      {
	TEUCHOS_FUNC_TIME_MONITOR("Device PCE calculation");
	v_dev = simple_function(u_dev);
      }

      {
	TEUCHOS_FUNC_TIME_MONITOR("Device->Host transfer");
	v_dev.load(v);
      }

    }

    // Print u and v
    std::cout << "v = 1.0 / (log(u)^2 + 1):" << std::endl;
    std::cout << "\tu = ";
    u.print(std::cout);
    std::cout << "\tv = ";
    v.print(std::cout);

    // Compute moments
    float mean = v.mean();
    float std_dev = v.standard_deviation();

    // Evaluate PCE and function at a point = 0.25 in each dimension
    Teuchos::Array<float> pt(d); 
    for (int i=0; i<d; i++) 
      pt[i] = 0.25;
    float up = u.evaluate(pt);
    float vp = simple_function(up);
    float vp2 = v.evaluate(pt);
    
    // Print results
    std::cout << "\tv mean         = " << mean << std::endl;
    std::cout << "\tv std. dev.    = " << std_dev << std::endl;
    std::cout << "\tv(0.25) (true) = " << vp << std::endl;
    std::cout << "\tv(0.25) (pce)  = " << vp2 << std::endl;
    
    // Check the answer
    if (std::abs(vp - vp2) < 1e-2)
      std::cout << "\nExample Passed!" << std::endl;

    Teuchos::TimeMonitor::summarize(std::cout);
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
