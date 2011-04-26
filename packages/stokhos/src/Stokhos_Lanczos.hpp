// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_LANCZOS_HPP
#define STOKHOS_LANCZOS_HPP

#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_Array.hpp"

namespace Stokhos {

  /*! 
   * \brief Applies Lanczos procedure to a given matrix
   */
  /*!
   * Given a matrix \f$A\f$, a starting vector \f$u_0\f$, and integer \f$k>0\f$,
   * applies the Lanczos procedure to compute an orthogonal basis for the
   * Krylov subspace 
   * \f[
   *   \mathcal{K}(A, u_0, k) = \{ u_0, A_u_0, A^2 u_0, \dots, A^k u_0 \}.
   * \f]
   * The basis vectors are given by 
   * \f[
   *   u_{i+1} = Au_i - \alpha_i u_i - \beta_i u_{i-1}, \;\; i=0,\dots,k
   * \f]
   * where $u_{-1} = 0$ and
   * \f[
   *   alpha_i = \frac{u_i^T W (Au_i)}{u_i^T W u_i}, \;\; 
   *   beta_i = \frac{u_i^T W u_i}{u_{i-1}^T W u_{i-1}}.
   * \f]
   * Here \f$W\f$ is a diagonal weighting matrix.
   */
  template <typename ordinal_type, typename value_type>
  class Lanczos {
  public:

    typedef Teuchos::SerialDenseMatrix<ordinal_type, value_type> matrix_type;
    typedef Teuchos::SerialDenseVector<ordinal_type, value_type> vector_type;

    //! Compute Lanczos basis
    static void compute(ordinal_type n,
			ordinal_type k, 
			const matrix_type& A,
			const Teuchos::Array<value_type>& w,
			const vector_type& u0,
			Teuchos::Array<vector_type>& u,
			Teuchos::Array<value_type>& alpha,
			Teuchos::Array<value_type>& beta,
			Teuchos::Array<value_type>& nrm_sqrd) {
      u.resize(k+1);
      alpha.resize(k);
      beta.resize(k);
      nrm_sqrd.resize(k);

      beta[0] = 1.0;
      u[0] = u0;

      value_type nrm;
      vector_type v(n);
      for (ordinal_type i=0; i<k; i++) {

	// Compute u_i^T*W*u_i
	nrm = value_type(0);
	for (ordinal_type j=0; j<n; j++)
	  nrm += w[j]*u[i][j]*u[i][j];
	nrm_sqrd[i] = nrm;

	// Compute v = A*u_i
	v.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 
		   value_type(1), A, u[i], value_type(0));

	// Compute v^T*W*u_i
	nrm = value_type(0);
	for (ordinal_type j=0; j<n; j++)
	  nrm += w[j]*v[j]*u[i][j];

	// Compute alpha = v^T*W*u_i / u_i^T*W*u_i
	alpha[i] = nrm / nrm_sqrd[i];

	// Compute beta = u_i^T*W*u_i / u_{i-1}^T*W*u_{i-1}
	if (i > 0)
	  beta[i] = nrm_sqrd[i] / nrm_sqrd[i-1];

	// Compute u_{i+1} = v - alpha_i*u_i - beta_i*u_{i-1}
	u[i+1].resize(n);
	if (i == 0) 
	  for (ordinal_type j=0; j<n; j++)
	    u[i+1][j] = v[j] - alpha[i]*u[i][j];
	else
	  for (ordinal_type j=0; j<n; j++)
	    u[i+1][j] = v[j] - alpha[i]*u[i][j] - beta[i]*u[i-1][j];

	std::cout << "i = " << i 
		  << " alpha = " << alpha[i] << " beta = " << beta[i]
		  << " nrm = " << nrm_sqrd[i] << std::endl;
    TEST_FOR_EXCEPTION(nrm_sqrd[i] < 0, std::logic_error,
		       "Stokhos::LanczosProjPCEBasis::lanczos():  "
		       << " Polynomial " << i << " out of " << k
		       << " has norm " << nrm_sqrd[i]
		       << "!  Try increasing number of quadrature points");

      }
    }

    //! Compute Lanczos basis
    static void computeDiag(ordinal_type n,
			    ordinal_type k, 
			    const Teuchos::Array<value_type>& A,
			    const Teuchos::Array<value_type>& w,
			    const Teuchos::Array<value_type>& u0,
			    Teuchos::Array< Teuchos::Array<value_type> >& u,
			    Teuchos::Array<value_type>& alpha,
			    Teuchos::Array<value_type>& beta,
			    Teuchos::Array<value_type>& nrm_sqrd) {
      u.resize(k+1);
      // alpha.resize(k);
      // beta.resize(k);
      // nrm_sqrd.resize(k);

      beta[0] = 1.0;
      u[0] = u0;

      value_type nrm;
      Teuchos::Array<value_type> v(n);
      for (ordinal_type i=0; i<k; i++) {

	// Compute u_i^T*W*u_i
	nrm = value_type(0);
	for (ordinal_type j=0; j<n; j++)
	  nrm += w[j]*u[i][j]*u[i][j];
	nrm_sqrd[i] = nrm;

	// Compute v = A*u_i
	for (ordinal_type j=0; j<n; j++)
	  v[j] = A[j]*u[i][j];

	// Compute v^T*W*u_i
	nrm = value_type(0);
	for (ordinal_type j=0; j<n; j++)
	  nrm += w[j]*v[j]*u[i][j];

	// Compute alpha = v^T*W*u_i / u_i^T*W*u_i
	alpha[i] = nrm / nrm_sqrd[i];

	// Compute beta = u_i^T*W*u_i / u_{i-1}^T*W*u_{i-1}
	if (i > 0)
	  beta[i] = nrm_sqrd[i] / nrm_sqrd[i-1];

	// Compute u_{i+1} = v - alpha_i*u_i - beta_i*u_{i-1}
	u[i+1].resize(n);
	if (i == 0) 
	  for (ordinal_type j=0; j<n; j++)
	    u[i+1][j] = v[j] - alpha[i]*u[i][j];
	else
	  for (ordinal_type j=0; j<n; j++)
	    u[i+1][j] = v[j] - alpha[i]*u[i][j] - beta[i]*u[i-1][j];

	std::cout << "i = " << i 
		  << " alpha = " << alpha[i] << " beta = " << beta[i]
		  << " nrm = " << nrm_sqrd[i] << std::endl;
    TEST_FOR_EXCEPTION(nrm_sqrd[i] < 0, std::logic_error,
		       "Stokhos::LanczosProjPCEBasis::lanczos():  "
		       << " Polynomial " << i << " out of " << k
		       << " has norm " << nrm_sqrd[i]
		       << "!  Try increasing number of quadrature points");

      }
    }
  };
}

#endif // STOKHOS_LANCZOS_HPP

