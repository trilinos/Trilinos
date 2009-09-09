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

#include <iostream>
#include <iomanip>

#include "Stokhos.hpp"
#include "Stokhos_StieltjesPCEBasis.hpp"

typedef Stokhos::LegendreBasis<int,double> basis_type;

// This example compares various PCE methods for computing moments of
//
// u = exp(x1 + ... + xd)
//
// where x1, ..., xd are uniform random variables on [-1,1].  The methods
// are compared to computing the "true" value via very high-order quadrature.
// Because of the structure of the exponential, the moments can easily
// be computed in one dimension.

struct stieltjes_pce_quad_func {
  stieltjes_pce_quad_func(
   const Stokhos::OrthogPolyApprox<int,double>& pce_,
   const Stokhos::OrthogPolyBasis<int,double>& basis_) :
    pce(pce_), basis(basis_), vec(1) {}
  
  double operator() (const double& a) const {
    vec[0] = a;
    return pce.evaluate(vec);
  }
  const Stokhos::OrthogPolyApprox<int,double>& pce;
  const Stokhos::OrthogPolyBasis<int,double>& basis;
  mutable Teuchos::Array<double> vec;
};

int main(int argc, char **argv)
{
  try {

    const unsigned int dmin = 1;
    const unsigned int dmax = 1;
    const unsigned int pmin = 5;
    const unsigned int pmax = 5;

    // Loop over dimensions
    for (unsigned int d=dmin; d<=dmax; d++) {

      Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d); 

      // Loop over orders
      for (unsigned int p=pmin; p<=pmax; p++) {

	// Create product basis
        for (unsigned int i=0; i<d; i++)
          bases[i] = Teuchos::rcp(new basis_type(p));
	Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
	  Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
	std::cout << *basis << std::endl;

	// Create approximation
	Stokhos::OrthogPolyApprox<int,double> x(basis), u(basis), v(basis), 
	  v2(basis);
	for (unsigned int i=0; i<d; i++) {
	  x.term(i, 1) = 1.0;
	}

	// Tensor product quadrature
	Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad = 
	  Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis,
									4*p));
	// Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad = 
	//   Teuchos::rcp(new Stokhos::SparseGridQuadrature<int,double>(basis,
	// 							     p+1));

	// Quadrature expansion
	Stokhos::QuadOrthogPolyExpansion<int,double> quad_exp(basis, quad);

	// Compute PCE via quadrature expansion
	//quad_exp.sinh(u,x);
	quad_exp.cos(u,x);
	//quad_exp.cosh(u,x);
	//quad_exp.sin(v,u);
	quad_exp.times(v,u,u);

	// Compute Stieltjes basis
	Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > st_bases(1);
	st_bases[0] = 
	  Teuchos::rcp(new Stokhos::StieltjesPCEBasis<int,double>(p, u, *quad, 
								  true));
	Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > st_basis = 
	  Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(st_bases));
	std::cout << *st_basis << std::endl;

	Stokhos::OrthogPolyApprox<int,double>  u_st(st_basis), v_st(st_basis);
	u_st.term(0, 0) = u.mean();
	u_st.term(0, 1) = 1.0;
	
        // Tensor product quadrature
	Teuchos::RCP<const Stokhos::Quadrature<int,double> > st_quad = 
	 Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(st_basis));

	// Quadrature expansion
	Stokhos::QuadOrthogPolyExpansion<int,double> st_quad_exp(st_basis, 
								 st_quad);
	
	//st_quad_exp.sin(v_st, u_st);
	st_quad_exp.times(v_st, u_st, u_st);

	//std::cout << v_st << std::endl;

	stieltjes_pce_quad_func func(v_st, *st_basis);
	quad_exp.unary_op(func, v2, u);

	std::cout.precision(12);
	std::cout << v;
	std::cout << v2;
	std::cout << v_st;
	
	std::cout.setf(std::ios::scientific);
	std::cout << "v.mean()       = " << v.mean() << std::endl
		  << "v_st.mean()    = " << v_st.mean() << std::endl
		  << "v2.mean()      = " << v2.mean() << std::endl
		  << "v.std_dev()    = " << v.standard_deviation() << std::endl
		  << "v_st.std_dev() = " << v_st.standard_deviation() 
		  << std::endl
		  << "v2.std_dev()   = " << v2.standard_deviation()
		  << std::endl;
      }
      
    }
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
