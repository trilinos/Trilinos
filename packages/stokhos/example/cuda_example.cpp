// $Id$ 
// $Source$ 
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

// cuda_example
//
//  usage: 
//     cuda_example
//
//  output:  
//     prints the Hermite Polynomial Chaos Expansion of the simple function
//
//     v = 1/(log(u)^2+1)
//
//     where u = 1 + 0.4*H_1(x) + 0.06*H_2(x) + 0.002*H_3(x), x is a zero-mean
//     and unit-variance Gaussian random variable, and H_i(x) is the i-th
//     Hermite polynomial.

#include "Stokhos.hpp"
#include "Stokhos_CUDAStorage.hpp"
#include "Stokhos_CUDAQuadOrthogPolyExpansion.hpp"
#include "Teuchos_TimeMonitor.hpp"

int main(int argc, char **argv)
{
  try {

    typedef Stokhos::StandardStorage<int,float> StdStorage;
    typedef Stokhos::CUDAStorage<int,float> Storage;
    //typedef Stokhos::StandardStorage<int,float> Storage;

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
      //basis->computeLinearTripleProductTensor();

    // Expansion method
    Stokhos::QuadOrthogPolyExpansion<int,float,Storage> expn(basis, Cijk, quad,
							     true);

    // Polynomial expansions
    Stokhos::OrthogPolyApprox<int,float,StdStorage> u(basis), v(basis);
    u.term(0,0) = 1.0;
    for (int i=0; i<d; i++) {
      u.term(i,1) = 0.4 / d;
      u.term(i,2) = 0.06 / d;
      u.term(i,3) = 0.002 / d;
    }

    // Compute expansion
    {
      Stokhos::OrthogPolyApprox<int,float,Storage> u_dev(basis), v_dev(basis), 
	w_dev(basis);
     
      {
	TEUCHOS_FUNC_TIME_MONITOR("Host->Device transfer");
	u_dev.init(u);
      }

      {
	 TEUCHOS_FUNC_TIME_MONITOR("Device PCE calculation");
	 expn.log(v_dev,u_dev);
	 expn.times(w_dev,v_dev,v_dev);
	 expn.plusEqual(w_dev,1.0);
	 expn.divide(v_dev,1.0,w_dev);
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
    float vp = 1.0/(std::log(up)*std::log(up)+1.0);
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
