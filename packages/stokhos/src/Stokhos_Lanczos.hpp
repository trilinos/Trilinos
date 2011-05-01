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

#include "Teuchos_Array.hpp"

namespace Stokhos {

  template <typename ord_type, typename val_type>
  class WeightedVectorSpace {
  public:
    typedef ord_type ordinal_type;
    typedef val_type value_type;
    typedef Teuchos::Array<value_type> vector_type;

    WeightedVectorSpace(const vector_type& weights) :
      w(weights),
      n(weights.size())
    {
    }

    ~WeightedVectorSpace() {}

    vector_type create_vector() const { return vector_type(n); }
    
    value_type 
    inner_product(const vector_type& u, const vector_type& v) const {
                  value_type val = value_type(0);
      for (ordinal_type j=0; j<n; j++)
	val += w[j]*u[j]*v[j];
      return val;
    }

    void
    add2(const value_type& a, const vector_type& u1,
         const value_type& b, const vector_type& u2, vector_type& v) const {
      for (ordinal_type j=0; j<n; j++)
        v[j] = a*u1[j] + b*u2[j];
    }

    void
    add3(const value_type& a, const vector_type& u1,
         const value_type& b, const vector_type& u2, 
         const value_type& c, const vector_type& u3, vector_type& v) const {
      for (ordinal_type j=0; j<n; j++)
        v[j] = a*u1[j] + b*u2[j] +c*u3[j];
    }

  protected:

    const vector_type& w;
    ordinal_type n;

  };

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
  template <typename vectorspace_type, typename matrix_type> 
  class Lanczos {
  public:

    typedef typename matrix_type::ordinal_type ordinal_type;
    typedef typename matrix_type::value_type value_type;
    typedef typename vectorspace_type::vector_type vector_type;

    //! Compute Lanczos basis
    static void compute(ordinal_type k, 
                        const vectorspace_type& vs,
			const matrix_type& A,
			const vector_type& u0,
			Teuchos::Array<vector_type>& u,
			Teuchos::Array<value_type>& alpha,
			Teuchos::Array<value_type>& beta,
			Teuchos::Array<value_type>& nrm_sqrd) {
      beta[0] = 1.0;
      u[0] = u0;

      value_type nrm;
      vector_type v = vs.create_vector();
      for (ordinal_type i=0; i<k; i++) {

	// Compute (u_i,u_i)
	nrm_sqrd[i] = vs.inner_product(u[i], u[i]);

	// Compute v = A*u_i
        A.apply(u[i], v);

	// Compute (v,u_i)
	nrm = vs.inner_product(u[i], v);

	// Compute alpha = (v,u_i) / (u_i,u_i)
	alpha[i] = nrm / nrm_sqrd[i];

	// Compute beta = (u_i,u_i) / (u_{i-1}.u_{i-1})
	if (i > 0)
	  beta[i] = nrm_sqrd[i] / nrm_sqrd[i-1];

	// Compute u_{i+1} = v - alpha_i*u_i - beta_i*u_{i-1}
	if (i == 0) 
          vs.add2(value_type(1), v, -alpha[i], u[i], u[i+1]);
	else
          vs.add3(value_type(1), v, -alpha[i], u[i], -beta[i], u[i-1], u[i+1]);

	std::cout << "i = " << i 
		  << " alpha = " << alpha[i] << " beta = " << beta[i]
		  << " nrm = " << nrm_sqrd[i] << std::endl;
        TEST_FOR_EXCEPTION(nrm_sqrd[i] < 0, std::logic_error,
	       	           "Stokhos::LanczosProjPCEBasis::lanczos():  "
		           << " Polynomial " << i << " out of " << k
		           << " has norm " << nrm_sqrd[i] << "!");

      }
    }

  };
}

#endif // STOKHOS_LANCZOS_HPP

