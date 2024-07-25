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
#include "Stokhos_ReducedBasisFactory.hpp"

#include "Teuchos_CommandLineProcessor.hpp"

typedef Stokhos::LegendreBasis<int,double> basis_type;

int main(int argc, char **argv)
{
  try {

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This example runs a Gram-Schmidt-based dimension reduction example.\n");
 
    int d = 2;
    CLP.setOption("d", &d, "Stochastic dimension");

    int p = 5;
    CLP.setOption("p", &p, "Polynomial order");

    CLP.parse( argc, argv );

    // Create product basis
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d);
    for (int i=0; i<d; i++)
      bases[i] = Teuchos::rcp(new basis_type(p, true));
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));

    std::cout << "original basis size = " << basis->size() << std::endl;

    // Create approximation
    Stokhos::OrthogPolyApprox<int,double> x(basis), u(basis), v(basis), 
      w(basis), w2(basis);
    for (int i=0; i<d; i++) {
      x.term(i, 1) = 1.0;
    }

    // Tensor product quadrature
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad = 
      Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
    // Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad = 
    //   Teuchos::rcp(new Stokhos::SparseGridQuadrature<int,double>(basis,
    // 								 p+1));

    std::cout << "original quadrature size = " << quad->size() << std::endl;

    // Triple product tensor
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk =
      basis->computeTripleProductTensor();
    
    // Quadrature expansion
    Stokhos::QuadOrthogPolyExpansion<int,double> quad_exp(basis, Cijk, quad);
    
    // Compute PCE via quadrature expansion
    quad_exp.sin(u,x);
    quad_exp.exp(v,x);
    quad_exp.times(w,v,u);
    
    // Create new basis from u and v
    Teuchos::Array< Stokhos::OrthogPolyApprox<int,double> > pces(2);
    pces[0] = u;
    pces[1] = v;
    Teuchos::ParameterList params;
    params.set("Reduced Basis Method", "Monomial Proj Gram-Schmidt");
    Stokhos::ReducedBasisFactory<int,double> factory(params);
    Teuchos::RCP< Stokhos::ReducedPCEBasis<int,double> > gs_basis = 
      factory.createReducedBasis(p, pces, quad, Cijk);
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > gs_quad =
      gs_basis->getReducedQuadrature();
    Stokhos::OrthogPolyApprox<int,double>  u_gs(gs_basis), v_gs(gs_basis), 
      w_gs(gs_basis);
    gs_basis->transformFromOriginalBasis(u.coeff(), u_gs.coeff());
    gs_basis->transformFromOriginalBasis(v.coeff(), v_gs.coeff());

    std::cout << "reduced basis size = " << gs_basis->size() << std::endl;
    std::cout << "reduced quadrature size = " << gs_quad->size() << std::endl;

    Stokhos::OrthogPolyApprox<int,double> u2(basis), v2(basis); 
    gs_basis->transformToOriginalBasis(u_gs.coeff(), u2.coeff());
    gs_basis->transformToOriginalBasis(v_gs.coeff(), v2.coeff());
    double err_u = 0.0;
    double err_v = 0.0;
    for (int i=0; i<basis->size(); i++) {
      double eu = std::abs(u[i]-u2[i]);
      double ev = std::abs(v[i]-v2[i]);
      if (eu > err_u) err_u = eu;
      if (ev > err_v) err_v = ev;
    }
    std::cout << "error in u transformation = " << err_u << std::endl;
    std::cout << "error in v transformation = " << err_v << std::endl;
    
    // Triple product tensor
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > gs_Cijk =
      Teuchos::null;
    
    // Gram-Schmidt quadrature expansion
    Teuchos::RCP< Teuchos::ParameterList > gs_exp_params = 
      Teuchos::rcp(new Teuchos::ParameterList);
    gs_exp_params->set("Use Quadrature for Times", true);
    Stokhos::QuadOrthogPolyExpansion<int,double> gs_quad_exp(gs_basis, 
							     gs_Cijk,
							     gs_quad,
							     gs_exp_params);
    
    // Compute w_gs = u_gs*v_gs in Gram-Schmidt basis
    gs_quad_exp.times(w_gs, u_gs, v_gs);
    
    // Project w_gs back to original basis
    gs_basis->transformToOriginalBasis(w_gs.coeff(), w2.coeff());

    std::cout.precision(12);
    std::cout << "w = " << std::endl << w;
    std::cout << "w2 = " << std::endl << w2;
    std::cout << "w_gs = " << std::endl << w_gs;

    double err_w = 0.0;
    for (int i=0; i<basis->size(); i++) {
      double ew = std::abs(w[i]-w2[i]);
      if (ew > err_w) err_w = ew;
    }
    
    std::cout.setf(std::ios::scientific);
    std::cout << "w.mean()       = " << w.mean() << std::endl
	      << "w2.mean()      = " << w2.mean() << std::endl
	      << "mean error     = " 
	      << std::abs(w.mean()-w2.mean()) << std::endl
	      << "w.std_dev()    = " << w.standard_deviation() << std::endl
	      << "w2.std_dev()   = " << w2.standard_deviation() << std::endl
	      << "std_dev error  = " 
	      << std::abs(w.standard_deviation()-w2.standard_deviation()) 
	      << std::endl
	      << "w coeff error  = " << err_w << std::endl;
    
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
