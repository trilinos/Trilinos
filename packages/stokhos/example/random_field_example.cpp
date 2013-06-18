// $Id: hermite_example.cpp,v 1.13 2009/09/14 21:53:56 etphipp Exp $ 
// $Source: /space/CVS/Trilinos/packages/stokhos/example/hermite_example.cpp,v $ 
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

// hermite_example
//
//  usage: 
//     hermite_example
//
//  output:  
//     prints the Hermite Polynomial Chaos Expansion of the simple function
//
//     v = 1/(log(u)^2+1)
//
//     where u = 1 + 0.4*H_1(x) + 0.06*H_2(x) + 0.002*H_3(x), x is a zero-mean
//     and unit-variance Gaussian random variable, and H_i(x) is the i-th
//     Hermite polynomial.

#include "Stokhos_KL_ExponentialRandomField.hpp"

int main(int argc, char **argv)
{
  try {

    int M = 10;
    Teuchos::ParameterList solverParams;
    solverParams.set("Number of KL Terms", M);
    solverParams.set("Mean", 1.0);
    solverParams.set("Standard Deviation", 0.1);
    int ndim = 3;
    Teuchos::Array<double> domain_upper(ndim), domain_lower(ndim), 
      correlation_length(ndim);
    for (int i=0; i<ndim; i++) {
      domain_upper[i] = 1.0;
      domain_lower[i] = 0.0;
      correlation_length[i] = 10.0;
    }
    solverParams.set("Domain Upper Bounds", domain_upper);
    solverParams.set("Domain Lower Bounds", domain_lower);
    solverParams.set("Correlation Lengths", correlation_length);
    
    Stokhos::KL::ExponentialRandomField<double> rf(solverParams);
    rf.print(std::cout);

    Teuchos::Array<double> x(ndim);
    for (int i=0; i<ndim; i++)
      x[i] = (domain_upper[i] + domain_lower[i])/2.0 + 
	0.1*(domain_upper[i] - domain_lower[i])/2.0;
    Teuchos::Array<double> rvar(M);
    for (int i=0; i<M; i++)
      rvar[i] = 1.5;
    double result = rf.evaluate(x, rvar);
    std::cout << "result = " << result << std::endl;
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
