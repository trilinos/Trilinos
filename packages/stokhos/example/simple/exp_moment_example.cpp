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

#include <iostream>
#include <iomanip>

#include "Stokhos.hpp"

typedef Stokhos::LegendreBasis<int,double> basis_type;

// This example uses PC expansions for computing moments of
//
// u = exp(x1 + ... + xd)
//
// where x1, ..., xd are uniform random variables on [-1,1].  The methods
// are compared to computing the "true" value via very high-order quadrature.
// Because of the structure of the exponential, the moments can easily
// be computed in one dimension.

int main(int argc, char **argv)
{
  try {

    // Compute "true" 1-D mean, std. dev using quadrature
    const unsigned int true_quad_order = 200;
    basis_type tmp_basis(1);
    Teuchos::Array<double> true_quad_points, true_quad_weights;
    Teuchos::Array< Teuchos::Array<double> > true_quad_values;
    tmp_basis.getQuadPoints(true_quad_order, true_quad_points, 
			    true_quad_weights, true_quad_values);
    double mean_1d = 0.0;
    double sd_1d = 0.0;
    for (unsigned int qp=0; qp<true_quad_points.size(); qp++) {
      double t = std::exp(true_quad_points[qp]);
      mean_1d += t*true_quad_weights[qp];
      sd_1d += t*t*true_quad_weights[qp];
    }

    const unsigned int dmin = 1;
    const unsigned int dmax = 4;
    const unsigned int pmin = 1;
    const unsigned int pmax = 5;

    // Loop over dimensions
    for (unsigned int d=dmin; d<=dmax; d++) {

      // compute "true" values
      double true_mean = std::pow(mean_1d,static_cast<double>(d));
      double true_sd = std::pow(sd_1d,static_cast<double>(d)) - 
	true_mean*true_mean;
      true_sd = std::sqrt(true_sd);
      std::cout.precision(12);
      std::cout << "true mean = " << true_mean << "\t true std. dev. = "
                << true_sd << std::endl;

      Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d); 

      // Loop over orders
      for (unsigned int p=pmin; p<=pmax; p++) {

	// Create product basis
        for (unsigned int i=0; i<d; i++)
          bases[i] = Teuchos::rcp(new basis_type(p));
	Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
	  Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));

	// Create approximation
	int sz = basis->size();
	Stokhos::OrthogPolyApprox<int,double> x(basis), u(basis);
	for (unsigned int i=0; i<d; i++) {
	  x.term(i,1) = 1.0;
	}

	// Tensor product quadrature
	Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad = 
	  Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));

	// Triple product tensor
	Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk =
	  basis->computeTripleProductTensor();

	// Quadrature expansion
	Stokhos::QuadOrthogPolyExpansion<int,double> quad_exp(basis, Cijk, 
							      quad);

	// Compute PCE via quadrature expansion
	quad_exp.exp(u,x);
	double mean = u.mean(); 
	double sd = u.standard_deviation();
        
	std::cout.precision(4);
	std::cout.setf(std::ios::scientific);
	std::cout << "d = " << d << " p = " << p
		  << " sz = " << sz
		  << "\tmean err = " 
		  << std::fabs(true_mean-mean) << "\tstd. dev. err = "
		  << std::fabs(true_sd-sd) << std::endl;
      }
      
    }
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
