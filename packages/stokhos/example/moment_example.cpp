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

// hermite_example
//
//  usage: 
//     hermite_example
//
//  output:  
//     prints the Hermite Polynomial Chaos Expansion of a simple function

#include <iostream>
#include <iomanip>

#include "Stokhos.hpp"

typedef Stokhos::LegendreBasis<int,double> basis_type;

// This example compares various PCE methods for computing moments of
//
// u = exp(x1 + ... + xd)
//
// where x1, ..., xd are uniform random variables on [-1,1].  The methods
// are compared to computing the "true" value via very high-order quadrature.
// Because of the structure of the exponential, the moments can easily
// be computed in one dimension.

// Computes first and second moments of a PCE
void pce_moments(const Stokhos::OrthogPolyApprox<int,double>& pce,
                 const Stokhos::OrthogPolyBasis<int,double>& basis,
                 double& mean, double& std_dev) {
  mean = pce[0];
  std_dev = 0.0;
  const std::vector<double> nrm2 = basis.norm_squared();
  for (int i=1; i<basis.size(); i++)
    std_dev += pce[i]*pce[i]*nrm2[i];
  std_dev = std::sqrt(std_dev);
}

int main(int argc, char **argv)
{
  try {

    // Compute "true" 1-D mean, std. dev using quadrature
    const unsigned int true_quad_order = 200;
    basis_type tmp_basis(1);
    std::vector<double> true_quad_points, true_quad_weights;
    std::vector< std::vector<double> > true_quad_values;
    tmp_basis.getQuadPoints(true_quad_order, true_quad_points, 
			    true_quad_weights, true_quad_values);

    const unsigned int nmin = 1;
    const unsigned int nmax = 10;
    const unsigned int pmin = 5;
    const unsigned int pmax = 5;
    const unsigned int d = 1;

    // Loop over moment
    for (unsigned int n=nmin; n<=nmax; n++) {

      // compute "true" values
      double moment_true = 0.0;
      for (unsigned int qp=0; qp<true_quad_points.size(); qp++)
	moment_true += std::exp(n*true_quad_points[qp])*true_quad_weights[qp];

      std::cout.precision(12);
      std::cout << "true moment = " << moment_true << std::endl;

      std::vector< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d); 

      // Loop over orders
      for (unsigned int p=pmin; p<=pmax; p++) {

	// Create product basis
        for (unsigned int i=0; i<d; i++)
          bases[i] = Teuchos::rcp(new basis_type(p));
	Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis = 
	  Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));

	// Create approximation
	int sz = basis->size();
	Stokhos::OrthogPolyApprox<int,double> x(sz), u(sz), v(sz);
	for (unsigned int i=0; i<d; i++) {
	  x.term2(*basis, i,1) = 1.0;
	}

	// Tensor product quadrature
	Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad = 
	  Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));

	// Quadrature expansion
	Stokhos::QuadOrthogPolyExpansion<int,double> quad_exp(basis, quad);

	// Compute PCE via quadrature expansion
	quad_exp.exp(u,x);

	for (unsigned int order=p; order<=3*p; order++) {
	  std::vector< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > moment_bases(d);
	  for (unsigned int i=0; i<d; i++)
	    moment_bases[i] = Teuchos::rcp(new basis_type(order));
	  Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > moment_basis = Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(moment_bases));
	  Teuchos::RCP<const Stokhos::Quadrature<int,double> > moment_quad = 
	    Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(moment_basis));
	  const std::vector<double>& moment_quad_weights = 
	    moment_quad->getQuadWeights();
	  const std::vector< std::vector<double> >& moment_quad_points = 
	    moment_quad->getQuadPoints();
	  const std::vector< std::vector<double> >& moment_quad_values =
	    moment_quad->getBasisAtQuadPoints();
	  
	  double moment_numerical = 0.0;
	  for (unsigned int i=0; i<moment_quad_weights.size(); i++) {
	    double t = u.evaluate(*moment_basis, 
				  moment_quad_points[i], 
				  moment_quad_values[i]);
	    moment_numerical += std::pow(t,static_cast<double>(n))*moment_quad_weights[i];
	  }
        
	  std::cout.precision(4);
	  std::cout.setf(std::ios::scientific);
	  std::cout << "n = " << n << " p = " << p
		    << " d = " << d << " order = " << order
		    << "\tmoment err = " 
		    << std::fabs(moment_true-moment_numerical) << std::endl;

	}
      }
      
    }
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
