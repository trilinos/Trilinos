// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
  
  double operator() (const double& a, const double& b) const {
    vec[0] = a;
    vec[1] = b;
    return pce.evaluate(vec);
  }
  const Stokhos::OrthogPolyApprox<int,double>& pce;
  const Stokhos::OrthogPolyBasis<int,double>& basis;
  mutable Teuchos::Array<double> vec;
};

int main(int argc, char **argv)
{
  try {

    const unsigned int d = 2;
    const unsigned int p = 5;

    // Create product basis
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d);
    for (unsigned int i=0; i<d; i++)
      bases[i] = Teuchos::rcp(new basis_type(p));
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));

    // Create approximation
    Stokhos::OrthogPolyApprox<int,double> x(basis), u(basis), v(basis), 
      w(basis), w2(basis), w3(basis);
    for (unsigned int i=0; i<d; i++) {
      x.term(i, 1) = 1.0;
    }

    // Tensor product quadrature
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad = 
      Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));

    // Triple product tensor
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk =
      basis->computeTripleProductTensor();
    
    // Quadrature expansion
    Stokhos::QuadOrthogPolyExpansion<int,double> quad_exp(basis, Cijk, quad);
    
    // Compute PCE via quadrature expansion
    quad_exp.sin(u,x);
    quad_exp.exp(v,x);
    quad_exp.times(w,v,u);
    
    // Compute tensor product Stieltjes basis for u and v
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > st_bases(2);
    st_bases[0] = 
      Teuchos::rcp(new Stokhos::StieltjesPCEBasis<int,double>(
		     p, Teuchos::rcp(&u,false), quad, true));
    st_bases[1] = 
      Teuchos::rcp(new Stokhos::StieltjesPCEBasis<int,double>(
		     p, Teuchos::rcp(&v,false), quad, true));
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > st_basis =
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(st_bases));
    Stokhos::OrthogPolyApprox<int,double>  u_st(st_basis), v_st(st_basis), 
      w_st(st_basis);
    u_st.term(0, 0) = u.mean();
    u_st.term(0, 1) = 1.0;
    v_st.term(0, 0) = v.mean();
    v_st.term(1, 1) = 1.0;
    
    // Use Gram-Schmidt to orthogonalize Stieltjes basis of u and v
    Teuchos::Array<double> st_points_0;
    Teuchos::Array<double> st_weights_0;
    Teuchos::Array< Teuchos::Array<double> > st_values_0;
    st_bases[0]->getQuadPoints(p+1, st_points_0, st_weights_0, st_values_0);
    Teuchos::Array<double> st_points_1;
    Teuchos::Array<double> st_weights_1;
    Teuchos::Array< Teuchos::Array<double> > st_values_1;
    st_bases[1]->getQuadPoints(p+1, st_points_1, st_weights_1, st_values_1);
    Teuchos::Array< Teuchos::Array<double> > st_points(st_points_0.size());
    for (int i=0; i<st_points_0.size(); i++) {
      st_points[i].resize(2);
      st_points[i][0] = st_points_0[i];
      st_points[i][1] = st_points_1[i];
    }
    Teuchos::Array<double> st_weights = st_weights_0;
    
    Teuchos::RCP< Stokhos::GramSchmidtBasis<int,double> > gs_basis = 
      Teuchos::rcp(new Stokhos::GramSchmidtBasis<int,double>(st_basis,
							     st_points,
							     st_weights,
							     1e-15));
    
    // Create quadrature for Gram-Schmidt basis using quad points and 
    // and weights from original basis mapped to Stieljtes basis
    Teuchos::RCP< const Teuchos::Array< Teuchos::Array<double> > > points =
      Teuchos::rcp(&st_points,false);
    Teuchos::RCP< const Teuchos::Array<double> > weights =
      Teuchos::rcp(&st_weights,false);
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > gs_quad = 
      Teuchos::rcp(new Stokhos::UserDefinedQuadrature<int,double>(gs_basis,
								  points,
								  weights));

    // Triple product tensor
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > gs_Cijk =
      gs_basis->computeTripleProductTensor();
    
    // Gram-Schmidt quadrature expansion
    Stokhos::QuadOrthogPolyExpansion<int,double> gs_quad_exp(gs_basis, 
							     gs_Cijk,
							     gs_quad);
    
    Stokhos::OrthogPolyApprox<int,double>  u_gs(gs_basis), v_gs(gs_basis), 
      w_gs(gs_basis);
    
    // Map expansion in Stieltjes basis to Gram-Schmidt basis
    gs_basis->transformCoeffs(u_st.coeff(), u_gs.coeff());
    gs_basis->transformCoeffs(v_st.coeff(), v_gs.coeff());
    
    // Compute w_gs = u_gs*v_gs in Gram-Schmidt basis
    gs_quad_exp.times(w_gs, u_gs, v_gs);
    
    // Project w_gs back to original basis
    pce_quad_func gs_func(w_gs, *gs_basis);
    quad_exp.binary_op(gs_func, w2, u, v);

    // Triple product tensor
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > st_Cijk =
      st_basis->computeTripleProductTensor();

    // Stieltjes quadrature expansion
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > st_quad = 
      Teuchos::rcp(new Stokhos::UserDefinedQuadrature<int,double>(st_basis,
								  points,
								  weights));
    Stokhos::QuadOrthogPolyExpansion<int,double> st_quad_exp(st_basis,
							     st_Cijk,
							     st_quad);
    
    // Compute w_st = u_st*v_st in Stieltjes basis
    st_quad_exp.times(w_st, u_st, v_st);
    
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
    std::cout << "w.mean()       = " << w.mean() << std::endl
	      << "w2.mean()      = " << w2.mean() << std::endl
	      << "w3.mean()      = " << w3.mean() << std::endl
	      << "w_gs.mean()    = " << w_gs.mean() << std::endl
	      << "w_st.mean()    = " << w_st.mean() << std::endl
	      << "w.std_dev()    = " << w.standard_deviation() << std::endl
	      << "w2.std_dev()   = " << w2.standard_deviation() << std::endl
	      << "w3.std_dev()   = " << w3.standard_deviation() << std::endl
	      << "w_gs.std_dev() = " << w_gs.standard_deviation() << std::endl
	      << "w_st.std_dev() = " << w_st.standard_deviation() << std::endl;
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
