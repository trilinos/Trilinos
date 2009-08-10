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
#include "Stokhos_StieltjesPCEBasis.hpp"
#include "Stokhos_GramSchmidtBasis.hpp"
#include "Stokhos_UserDefinedQuadrature.hpp"

typedef Stokhos::LegendreBasis<int,double> basis_type;

struct pce_quad_func {
  pce_quad_func(
   const Stokhos::OrthogPolyApprox<int,double>& pce_,
   const Stokhos::OrthogPolyBasis<int,double>& basis_) :
    pce(pce_), basis(basis_), vec(2) {}
  
  double operator() (const double& a, double& b) const {
    vec[0] = a;
    vec[1] = b;
    return pce.evaluate(basis, vec);
  }
  const Stokhos::OrthogPolyApprox<int,double>& pce;
  const Stokhos::OrthogPolyBasis<int,double>& basis;
  mutable std::vector<double> vec;
};

int main(int argc, char **argv)
{
  try {

    const unsigned int d = 1;
    const unsigned int p = 5;

    // Create product basis
    std::vector< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d);
    for (unsigned int i=0; i<d; i++)
      bases[i] = Teuchos::rcp(new basis_type(p));
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));

    // Create approximation
    int sz = basis->size();
    Stokhos::OrthogPolyApprox<int,double> x(sz), u(sz), v(sz), w(sz), w2(sz), w3(sz);
    for (unsigned int i=0; i<d; i++) {
      x.term2(*basis, i, 1) = 1.0;
    }

    // Tensor product quadrature
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad = 
      Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
    
    // Quadrature expansion
    Stokhos::QuadOrthogPolyExpansion<int,double> quad_exp(basis, quad);
    
    // Compute PCE via quadrature expansion
    quad_exp.sin(u,x);
    quad_exp.exp(v,x);
    quad_exp.times(w,v,u);
    
    // Compute tensor product Stieltjes basis for u and v
    std::vector< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > st_bases(2);
    st_bases[0] = 
      Teuchos::rcp(new Stokhos::StieltjesPCEBasis<int,double>(p, u, *basis, 
							      *quad, true));
    st_bases[1] = 
      Teuchos::rcp(new Stokhos::StieltjesPCEBasis<int,double>(p, v, *basis, 
							      *quad, true));
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > st_basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(st_bases));
    int st_sz = st_basis->size();
    Stokhos::OrthogPolyApprox<int,double>  u_st(st_sz), v_st(st_sz), w_st(st_sz);
    u_st.term2(*st_basis, 0, 0) = u.mean(*basis);
    u_st.term2(*st_basis, 0, 1) = 1.0;
    v_st.term2(*st_basis, 0, 0) = v.mean(*basis);
    v_st.term2(*st_basis, 1, 1) = 1.0;
    
    // Stieltjes quadrature expansion
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > st_quad = 
      Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(st_basis));
    Stokhos::QuadOrthogPolyExpansion<int,double> st_quad_exp(st_basis, 
							     st_quad);
    
    // Compute w_st = u_st*v_st in Stieltjes basis
    st_quad_exp.times(w_st, u_st, v_st);
    
    // Use Gram-Schmidt to orthogonalize Stieltjes basis of u and v
    std::vector<double> st_points_0;
    std::vector<double> st_weights_0;
    std::vector< std::vector<double> > st_values_0;
    st_bases[0]->getQuadPoints(p+1, st_points_0, st_weights_0, st_values_0);
    std::vector<double> st_points_1;
    std::vector<double> st_weights_1;
    std::vector< std::vector<double> > st_values_1;
    st_bases[1]->getQuadPoints(p+1, st_points_1, st_weights_1, st_values_1);
    std::vector< std::vector<double> > st_points(st_points_0.size());
    for (unsigned int i=0; i<st_points_0.size(); i++) {
      st_points[i].resize(2);
      st_points[i][0] = st_points_0[i];
      st_points[i][1] = st_points_1[i];
    }
    std::vector<double> st_weights = st_weights_0;
    
    Teuchos::RCP< Stokhos::GramSchmidtBasis<int,double> > gs_basis = 
      Teuchos::rcp(new Stokhos::GramSchmidtBasis<int,double>(st_basis,
							     st_points,
							     st_weights,
							     1e-15));
    
    // Create quadrature for Gram-Schmidt basis using quad points and 
    // and weights from original basis mapped to Stieljtes basis
    Teuchos::RCP< const std::vector< std::vector<double> > > points =
      Teuchos::rcp(&st_points,false);
    Teuchos::RCP< const std::vector<double> > weights =
      Teuchos::rcp(&st_weights,false);
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > gs_quad = 
      Teuchos::rcp(new Stokhos::UserDefinedQuadrature<int,double>(gs_basis,
								  points,
								  weights));
    
    // Gram-Schmidt quadrature expansion
    Stokhos::QuadOrthogPolyExpansion<int,double> gs_quad_exp(gs_basis, 
							     gs_quad);
    
    Stokhos::OrthogPolyApprox<int,double>  u_gs(st_sz), v_gs(st_sz), w_gs(st_sz);
    
    // Map expansion in Stieltjes basis to Gram-Schmidt basis
    gs_basis->transformCoeffs(u_st.coeff(), u_gs.coeff());
    gs_basis->transformCoeffs(v_st.coeff(), v_gs.coeff());
    
    // Compute w_gs = u_gs*v_gs in Gram-Schmidt basis
    gs_quad_exp.times(w_gs, u_gs, v_gs);
    
    // Project w_gs back to original basis
    pce_quad_func gs_func(w_gs, *gs_basis);
    quad_exp.binary_op(gs_func, w2, u, v);
    
    // Project w_st back to original basis
    pce_quad_func st_func(w_st, *st_basis);
    quad_exp.binary_op(st_func, w3, u, v);
    
    std::cout.precision(12);
    std::cout << "w = " << std::endl << w;
    std::cout << "w2 = " << std::endl << w2;
    std::cout << "w3 = " << std::endl << w3;
    std::cout << "w_gs = " << std::endl << w_gs;
    std::cout << "w_st = " << std::endl << w_st;
    
    std::cout.setf(std::ios::scientific);
    std::cout << "w.mean()       = " << w.mean(*basis) << std::endl
	      << "w2.mean()      = " << w2.mean(*basis) << std::endl
	      << "w3.mean()      = " << w3.mean(*basis) << std::endl
	      << "w_gs.mean()    = " << w_gs.mean(*gs_basis) << std::endl
	      << "w_st.mean()    = " << w_st.mean(*st_basis) << std::endl
	      << "w.std_dev()    = " << w.standard_deviation(*basis) 
	      << std::endl
	      << "w2.std_dev()   = " << w2.standard_deviation(*basis)
	      << std::endl
	      << "w3.std_dev()   = " << w3.standard_deviation(*basis)
	      << std::endl
	      << "w_gs.std_dev() = " << w_gs.standard_deviation(*gs_basis)
	      << std::endl
	      << "w_st.std_dev() = " << w_st.standard_deviation(*st_basis)
	      << std::endl;
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
