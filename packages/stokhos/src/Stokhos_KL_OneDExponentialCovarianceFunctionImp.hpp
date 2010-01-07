// $Id: Stokhos_Quadrature.hpp,v 1.4 2009/09/14 18:35:48 etphipp Exp $ 
// $Source: /space/CVS/Trilinos/packages/stokhos/src/Stokhos_Quadrature.hpp,v $ 
// @HEADER
// ***********************************************************************
// 
//                     Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_TestForException.hpp"

template <typename value_type>
Stokhos::KL::OneDExponentialCovarianceFunction<value_type>::
OneDExponentialCovarianceFunction(int M, 
				  const value_type& a, 
				  const value_type& b, 
				  const value_type& L_,
				  const std::string& dim_name,
				  Teuchos::ParameterList& solverParams) :
  L(L_),
  eig_pair(M)
{
  // Get parameters with default values
  magnitude_type eps = solverParams.get("Bound Perturbation Size", 1e-6);
  magnitude_type tol = solverParams.get("Nonlinear Solver Tolerance", 1e-10);
  int max_it = solverParams.get("Maximum Nonlinear Solver Iterations", 100);

  value_type aa, alpha, omega, lambda;
  int i=0;
  double pi = 4.0*std::atan(1.0);
  int idx = 0;
  
  aa = (b-a)/2.0;
  while (i < M-1) {
    alpha = aa/L;
    omega = bisection(EigFuncCos(alpha), idx*pi, idx*pi+pi/2.0-eps, 
		      tol, max_it) / aa;
    lambda = 2.0*L/(L*L*omega*omega + 1.0);
    eig_pair[i].eig_val = lambda;
    eig_pair[i].eig_func = Teuchos::rcp(new 
      ExponentialOneDEigenFunction<value_type>(
	ExponentialOneDEigenFunction<value_type>::COS, a, b, omega, dim_name)
      );
    i++;

    omega = bisection(EigFuncSin(alpha), idx*pi+pi/2.0+eps, (idx+1)*pi,
		      tol, max_it) / aa;
    lambda = 2.0*L/(L*L*omega*omega + 1.0);
    eig_pair[i].eig_val = lambda;
    eig_pair[i].eig_func = Teuchos::rcp(new
      ExponentialOneDEigenFunction<value_type>(
	ExponentialOneDEigenFunction<value_type>::SIN, a, b, omega, dim_name)
      );
    i++;

    idx++;
  }
  if (i < M) {
    omega = bisection(EigFuncCos(alpha), idx*pi, idx*pi+pi/2.0-eps, 
		      tol, max_it) / aa;
    lambda = 2.0*L/(L*L*omega*omega + 1.0);
    eig_pair[i].eig_val = lambda;
    eig_pair[i].eig_func = Teuchos::rcp(new
      ExponentialOneDEigenFunction<value_type>(
	ExponentialOneDEigenFunction<value_type>::COS, a, b, omega, dim_name)
      );
  }
}

template <typename value_type>
value_type
Stokhos::KL::OneDExponentialCovarianceFunction<value_type>::
evaluateCovariance(const value_type& x, const value_type& xp) const
{
  return std::exp(-std::abs(x-xp)/L);
}

template <typename value_type>
const Teuchos::Array< Stokhos::KL::OneDEigenPair<value_type> >&
Stokhos::KL::OneDExponentialCovarianceFunction<value_type>::
getEigenPairs() const
{
  return eig_pair;
}

template <typename value_type>
template <class Func>
value_type
Stokhos::KL::OneDExponentialCovarianceFunction<value_type>::
newton(const Func& func, const value_type& a, const value_type& b, 
       magnitude_type tol, int max_num_its)
{
  value_type u = (a+b)/2.0;
  value_type f = func.eval(u);
  int nit = 0;
  while (Teuchos::ScalarTraits<value_type>::magnitude(f) > tol && 
	 nit < max_num_its) {
    std::cout << "u = " << u << " f = " << f << std::endl;
    value_type dfdu = func.deriv(u);
    u -= f / dfdu;
    f = func.eval(u);
    ++nit;
  }
  TEST_FOR_EXCEPTION(nit >= max_num_its, std::logic_error,
		     "Nonlinear solver did not converge!" << std::endl);

  return u;
}

template <typename value_type>
template <class Func>
value_type
Stokhos::KL::OneDExponentialCovarianceFunction<value_type>::
bisection(const Func& func, const value_type& a, const value_type& b, 
	  magnitude_type tol, int max_num_its)
{
  value_type low, hi;
  value_type fa = func.eval(a);
  value_type fb = func.eval(b);
  TEST_FOR_EXCEPTION(fa*fb > 0, std::logic_error,
		     "Bounds must bracket the root!" << std::endl)

  if (fa <= 0.0) {
    low = a;
    hi = b;
  }
  else {
    low = b;
    hi = a;
  }

  int nit = 0;
  value_type u = low + (hi - low)/2.0;
  value_type f = func.eval(u);
  while ((Teuchos::ScalarTraits<value_type>::magnitude(hi - low) > 2.0*tol || 
	  Teuchos::ScalarTraits<value_type>::magnitude(f) > tol) && 
	  nit < max_num_its) {
    //std::cout << "u = " << u << " f = " << f << std::endl;
    if (f <= 0.0)
      low = u;
    else
      hi = u;
    u = low + (hi - low)/2.0;
    f = func.eval(u);
    ++nit;
  }
  TEST_FOR_EXCEPTION(nit >= max_num_its, std::logic_error,
		     "Nonlinear solver did not converge!" << std::endl);

  return u;
}
