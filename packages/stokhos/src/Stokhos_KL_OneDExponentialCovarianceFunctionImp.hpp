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

#include "Teuchos_Assert.hpp"

template <typename value_type>
Stokhos::KL::OneDExponentialCovarianceFunction<value_type>::
OneDExponentialCovarianceFunction(int M,
                                  const value_type& a,
                                  const value_type& b,
                                  const value_type& L_,
                                  const int dim_name,
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
    eig_pair[i].eig_func = ExponentialOneDEigenFunction<value_type>(
      ExponentialOneDEigenFunction<value_type>::COS, a, b, omega, dim_name);
    i++;

    omega = bisection(EigFuncSin(alpha), idx*pi+pi/2.0+eps, (idx+1)*pi,
                      tol, max_it) / aa;
    lambda = 2.0*L/(L*L*omega*omega + 1.0);
    eig_pair[i].eig_val = lambda;
    eig_pair[i].eig_func = ExponentialOneDEigenFunction<value_type>(
      ExponentialOneDEigenFunction<value_type>::SIN, a, b, omega, dim_name);
    i++;

    idx++;
  }
  if (i < M) {
    omega = bisection(EigFuncCos(alpha), idx*pi, idx*pi+pi/2.0-eps,
                      tol, max_it) / aa;
    lambda = 2.0*L/(L*L*omega*omega + 1.0);
    eig_pair[i].eig_val = lambda;
    eig_pair[i].eig_func = ExponentialOneDEigenFunction<value_type>(
      ExponentialOneDEigenFunction<value_type>::COS, a, b, omega, dim_name);
  }
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
  TEUCHOS_TEST_FOR_EXCEPTION(nit >= max_num_its, std::logic_error,
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
  TEUCHOS_TEST_FOR_EXCEPTION(fa*fb > value_type(0.0), std::logic_error,
    "Bounds [" << a << "," << b << "] must bracket the root!" << std::endl <<
    "f(a) = " << fa << ", f(b) = " << fb << std::endl)

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
  TEUCHOS_TEST_FOR_EXCEPTION(nit >= max_num_its, std::logic_error,
                     "Nonlinear solver did not converge!" << std::endl);

  return u;
}
