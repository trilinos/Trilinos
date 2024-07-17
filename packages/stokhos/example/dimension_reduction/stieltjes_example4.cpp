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
#include "Stokhos_LanczosProjPCEBasis.hpp"
#include "Stokhos_LanczosPCEBasis.hpp"
#include "Stokhos_UserDefinedQuadrature.hpp"

typedef Stokhos::LegendreBasis<int,double> basis_type;

template <int Num>
struct s_quad_func {
  static const int N = Num;
  double a;
  
  s_quad_func(double a_) : a(a_) {}
  
  double operator() (double x[Num]) const {
    double prod = 1;
    for (int i=0; i<Num; i++)
      prod *= x[i];
    return 1.0 / (prod - a);
  }
};

struct pce_quad_func {
  pce_quad_func(
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

double rel_err(double a, double b) {
  return std::abs(a-b)/std::abs(b);
}

int main(int argc, char **argv)
{
  try {

    const int d = 5;
    const int pmin = 1;
    const int pmax = 10;
    const int np = pmax-pmin+1;
    const double a = 1.5;
    bool use_pce_quad_points = false;
    bool normalize = false;
    bool project_integrals = false;
    bool lanczos = true;
    bool sparse_grid = true;
#ifndef HAVE_STOKHOS_DAKOTA
    sparse_grid = false;
#endif
    Teuchos::Array<double> mean(np), mean_st(np), std_dev(np), std_dev_st(np);
    Teuchos::Array<double> pt(np), pt_st(np);

    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d); 

    Teuchos::Array<double> eval_pt(d, 0.56789);
    double pt_true;
    
    // Loop over orders
    int n = 0;
    for (int p=pmin; p<=pmax; p++) {

      std::cout << "p = " << p << std::endl;
      
      // Create product basis
      for (int i=0; i<d; i++)
	bases[i] = Teuchos::rcp(new basis_type(p));
      Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
	Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
      
      // Create approximation
      Stokhos::OrthogPolyApprox<int,double> s(basis), t1(basis), t2(basis);
      Teuchos::Array< Stokhos::OrthogPolyApprox<int,double>* > xi(d);
      for (int i=0; i<d; i++) {
	xi[i] = new Stokhos::OrthogPolyApprox<int,double>(basis);
	xi[i]->term(i, 1) = 1.0;
      }
      
      // Quadrature
      Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad;
#ifdef HAVE_STOKHOS_DAKOTA
      if (sparse_grid)
      	quad = 
	  Teuchos::rcp(new Stokhos::SparseGridQuadrature<int,double>(basis, p));
#endif
      if (!sparse_grid)
	quad = 
	  Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));

      // Triple product tensor
      Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk =
	basis->computeTripleProductTensor();
      
      // Quadrature expansion
      Stokhos::QuadOrthogPolyExpansion<int,double> quad_exp(basis, Cijk, quad);
      
      // Compute PCE via quadrature expansion
      s_quad_func<d> s_func(a);
      const Stokhos::OrthogPolyApprox<int,double> **xip = 
	new const Stokhos::OrthogPolyApprox<int,double>*[d];
      for (int i=0; i<d; i++)
	xip[i] = xi[i];
      quad_exp.nary_op(s_func,s,xip);
      quad_exp.divide(t1,1.0,s);
      delete [] xip;

      // compute true point
      Teuchos::Array<double> xx(d);
      for (int i=0; i<d; i++)
	xx[i] = xi[i]->evaluate(eval_pt);
      pt_true = s_func(&(xx[0]));
      pt_true = 1.0/pt_true;
	
      // Compute Stieltjes basis
      Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > st_bases(1);
      Teuchos::RCP< Stokhos::LanczosProjPCEBasis<int,double> > stp_basis_s;
      Teuchos::RCP< Stokhos::LanczosPCEBasis<int,double> > st_basis_s;
      if (lanczos) {
        if (project_integrals) {
	  stp_basis_s = 
	    Teuchos::rcp(new Stokhos::LanczosProjPCEBasis<int,double>(
	  		 p, Teuchos::rcp(&s,false), Cijk, normalize, true));
	  st_bases[0] = stp_basis_s;
        }
        else {
  	  st_basis_s = 
	    Teuchos::rcp(new Stokhos::LanczosPCEBasis<int,double>(
			 p, Teuchos::rcp(&s,false), quad, normalize, true));
	  st_bases[0] = st_basis_s;
        }
      }
      else {
        st_bases[0] =
          Teuchos::rcp(new Stokhos::StieltjesPCEBasis<int,double>(
                       p, Teuchos::rcp(&s,false), quad, use_pce_quad_points,
                       normalize, project_integrals, Cijk));
      }
      
      Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > 
	st_basis = 
	Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(st_bases));
      //std::cout << *st_basis << std::endl;

      Stokhos::OrthogPolyApprox<int,double>  s_st(st_basis), t_st(st_basis);
      if (lanczos) {
        if (project_integrals) {
	  s_st.term(0, 0) = stp_basis_s->getNewCoeffs(0);
	  s_st.term(0, 1) = stp_basis_s->getNewCoeffs(1);
        }
        else {
	  s_st.term(0, 0) = st_basis_s->getNewCoeffs(0);
	  s_st.term(0, 1) = st_basis_s->getNewCoeffs(1);
        }
      }
      else {
        s_st.term(0, 0) = s.mean();
        s_st.term(0, 1) = 1.0;
      }

      // Triple product tensor
      Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > st_Cijk =
	st_basis->computeTripleProductTensor();
	
      // Tensor product quadrature
      Teuchos::RCP<const Stokhos::Quadrature<int,double> > st_quad;
#ifdef HAVE_STOKHOS_DAKOTA
      if (sparse_grid)
	st_quad = Teuchos::rcp(new Stokhos::SparseGridQuadrature<int,double>(st_basis, p));
#endif
      if (!sparse_grid)
	st_quad = Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(st_basis));
      
      // Quadrature expansion
      Stokhos::QuadOrthogPolyExpansion<int,double> st_quad_exp(st_basis, 
							       st_Cijk,
							       st_quad);
      
      // Compute t_st = 1/s_st in Stieltjes basis
      st_quad_exp.divide(t_st, 1.0, s_st);
      
      // Project t_st back to original basis
      pce_quad_func st_func(t_st, *st_basis);
      quad_exp.unary_op(st_func, t2, s);

      // std::cout.precision(12);
      // std::cout << w;
      // std::cout << w2;
      // std::cout << w_st;
      mean[n] = t1.mean();
      mean_st[n] = t2.mean();
      std_dev[n] = t1.standard_deviation();
      std_dev_st[n] = t2.standard_deviation();
      pt[n] = t1.evaluate(eval_pt);
      pt_st[n] = t2.evaluate(eval_pt);
      n++;

      for (int i=0; i<d; i++)
	delete xi[i];
    }

    n = 0;
    int wi=10;
    std::cout << "Statistical error:" << std::endl;
    std::cout << "p  " 
	      << std::setw(wi) << "mean" << "  " 
	      << std::setw(wi) << "mean_st" << "  "
	      << std::setw(wi) << "std_dev" << "  "
	      << std::setw(wi) << "std_dev_st" << "  "
	      << std::setw(wi) << "point" << "  "
	      << std::setw(wi) << "point_st" << std::endl;
    for (int p=pmin; p<pmax; p++) {
      std::cout.precision(3);
      std::cout.setf(std::ios::scientific);
      std::cout << p << "  " 
		<< std::setw(wi) << rel_err(mean[n], mean[np-1]) << "  "
		<< std::setw(wi) << rel_err(mean_st[n], mean[np-1]) << "  "
		<< std::setw(wi) << rel_err(std_dev[n], std_dev[np-1]) << "  "
		<< std::setw(wi) << rel_err(std_dev_st[n], std_dev[np-1]) 
		<< "  "
		<< std::setw(wi) << rel_err(pt[n], pt_true) << "  "
		<< std::setw(wi) << rel_err(pt_st[n], pt_true) 
		<< std::endl;
      n++;
    }
      
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
