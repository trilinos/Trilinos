// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
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
//     where u = 1 + 0.1*x_1 + 0.05*x_2 + 0.01*x_3 x1,x2,x3 are zero-mean
//     unit-variance Gaussian random variables.

#include "Stokhos_Sacado.hpp"

// The function to compute the polynomial chaos expansion of,
// written as a template function
template <class ScalarType>
ScalarType simple_function(const ScalarType& u) {
  ScalarType z = std::log(u);
  return 1.0/(z*z + 1.0);
}

int main(int argc, char **argv)
{
  // Typename of Polynomial Chaos scalar type
  typedef Stokhos::StandardStorage<int,double> storage_type;
  typedef Sacado::ETPCE::OrthogPoly<double, storage_type> pce_type;

  // Short-hand for several classes used below
  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Stokhos::OneDOrthogPolyBasis;
  using Stokhos::HermiteBasis;
  using Stokhos::CompletePolynomialBasis;
  using Stokhos::Quadrature;
  using Stokhos::TensorProductQuadrature;
  using Stokhos::Sparse3Tensor;
  using Stokhos::QuadOrthogPolyExpansion;

  try {

    // Basis of dimension 3, order 4
    const int d = 3;
    const int p = 4;
    Array< RCP<const OneDOrthogPolyBasis<int,double> > > bases(d);
    for (int i=0; i<d; i++) {
      bases[i] = rcp(new HermiteBasis<int,double>(p));
    }
    RCP<const CompletePolynomialBasis<int,double> > basis =
      rcp(new CompletePolynomialBasis<int,double>(bases));

    // Quadrature method
    RCP<const Quadrature<int,double> > quad =
      rcp(new TensorProductQuadrature<int,double>(basis));

    // Triple product tensor
    RCP<Sparse3Tensor<int,double> > Cijk =
      basis->computeTripleProductTensor();

    // Expansion method
    RCP<QuadOrthogPolyExpansion<int,double> > expn =
      rcp(new QuadOrthogPolyExpansion<int,double>(basis, Cijk, quad));

    // Polynomial expansion of u
    pce_type u(expn);
    u.term(0,0) = 1.0;     // zeroth order term
    u.term(0,1) = 0.1;     // first order term for dimension 0
    u.term(1,1) = 0.05;    // first order term for dimension 1
    u.term(2,1) = 0.01;    // first order term for dimension 2

    // Compute PCE expansion of function
    pce_type v = simple_function(u);

    // Print u and v
    std::cout << "\tu = ";
    u.print(std::cout);
    std::cout << "\tv = ";
    v.print(std::cout);

    // Compute moments
    double mean = v.mean();
    double std_dev = v.standard_deviation();

    // Evaluate PCE and function at a point = 0.25 in each dimension
    Teuchos::Array<double> pt(d);
    for (int i=0; i<d; i++)
      pt[i] = 0.25;
    double up = u.evaluate(pt);
    double vp = simple_function(up);
    double vp2 = v.evaluate(pt);

    // Print results
    std::cout << "\tv mean         = " << mean << std::endl;
    std::cout << "\tv std. dev.    = " << std_dev << std::endl;
    std::cout << "\tv(0.25) (true) = " << vp << std::endl;
    std::cout << "\tv(0.25) (pce)  = " << vp2 << std::endl;
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
