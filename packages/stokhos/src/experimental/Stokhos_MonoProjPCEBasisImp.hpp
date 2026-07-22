// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Assert.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::MonoProjPCEBasis<ordinal_type, value_type>::
MonoProjPCEBasis(
   ordinal_type p,
   const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& pce,
   const Stokhos::Quadrature<ordinal_type, value_type>& quad,
   const Stokhos::Sparse3Tensor<ordinal_type, value_type>& Cijk,
   bool limit_integration_order_) :
  RecurrenceBasis<ordinal_type, value_type>("Monomial Projection", p, true),
  limit_integration_order(limit_integration_order_),
  pce_sz(pce.basis()->size()),
  pce_norms(pce.basis()->norm_squared()),
  a(pce_sz), 
  b(pce_sz),
  basis_vecs(pce_sz, p+1),
  new_pce(p+1)
{
  // If the original basis is normalized, we can use the standard QR
  // factorization.  For simplicity, we renormalize the PCE coefficients
  // for a normalized basis
  Stokhos::OrthogPolyApprox<ordinal_type, value_type> normalized_pce(pce);
  for (ordinal_type i=0; i<pce_sz; i++) {
    pce_norms[i] = std::sqrt(pce_norms[i]);
    normalized_pce[i] *= pce_norms[i];
  }

  // Evaluate PCE at quad points
  ordinal_type nqp = quad.size();
  Teuchos::Array<value_type> pce_vals(nqp);
  const Teuchos::Array<value_type>& weights = quad.getQuadWeights();
  const Teuchos::Array< Teuchos::Array<value_type> >& quad_points =
    quad.getQuadPoints();
  const Teuchos::Array< Teuchos::Array<value_type> >& basis_values =
    quad.getBasisAtQuadPoints();
  for (ordinal_type i=0; i<nqp; i++) {
    pce_vals[i] = normalized_pce.evaluate(quad_points[i], basis_values[i]);
  }

  // Form Kylov matrix up to order pce_sz
  matrix_type K(pce_sz, pce_sz);

  // Compute matrix
  matrix_type A(pce_sz, pce_sz);
  typedef Stokhos::Sparse3Tensor<ordinal_type, value_type> Cijk_type;
  for (typename Cijk_type::k_iterator k_it = Cijk.k_begin();
       k_it != Cijk.k_end(); ++k_it) {
    ordinal_type k = index(k_it);
    for (typename Cijk_type::kj_iterator j_it = Cijk.j_begin(k_it); 
	 j_it != Cijk.j_end(k_it); ++j_it) {
      ordinal_type j = index(j_it);
      value_type val = 0;
      for (typename Cijk_type::kji_iterator i_it = Cijk.i_begin(j_it);
	   i_it != Cijk.i_end(j_it); ++i_it) {
	ordinal_type i = index(i_it);
	value_type c = value(i_it) / (pce_norms[j]*pce_norms[k]);
	val += pce[i]*c;
      }
      A(k,j) = val;
    }
  }

  // Each column i is given by projection of the i-th order monomial 
  // onto original basis
  vector_type u0 = Teuchos::getCol(Teuchos::View, K, 0);
  u0(0) = 1.0;
  for (ordinal_type i=1; i<pce_sz; i++)
    u0(i) = 0.0;
  for (ordinal_type k=1; k<pce_sz; k++) {
    vector_type u = Teuchos::getCol(Teuchos::View, K, k);
    vector_type up = Teuchos::getCol(Teuchos::View, K, k-1);
    u.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, up, 0.0);
  }
  /*
  for (ordinal_type j=0; j<pce_sz; j++) {
    for (ordinal_type i=0; i<pce_sz; i++) {
      value_type val = 0.0;
      for (ordinal_type k=0; k<nqp; k++)
	val += weights[k]*std::pow(pce_vals[k],j)*basis_values[k][i];
      K(i,j) = val;
    }
  }
  */

  std::cout << K << std::endl << std::endl;

  // Compute QR factorization of K
  ordinal_type ws_size, info;
  value_type ws_size_query;
  Teuchos::Array<value_type> tau(pce_sz);
  Teuchos::LAPACK<ordinal_type,value_type> lapack;
  lapack.GEQRF(pce_sz, pce_sz, K.values(), K.stride(), &tau[0], 
	       &ws_size_query, -1, &info);
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error, 
		     "GEQRF returned value " << info);
  ws_size = static_cast<ordinal_type>(ws_size_query);
  Teuchos::Array<value_type> work(ws_size);
  lapack.GEQRF(pce_sz, pce_sz, K.values(), K.stride(), &tau[0], 
	       &work[0], ws_size, &info);
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error, 
		     "GEQRF returned value " << info);
  
  // Get Q
  lapack.ORGQR(pce_sz, pce_sz, pce_sz, K.values(), K.stride(), &tau[0], 
	       &ws_size_query, -1, &info);
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error, 
		     "ORGQR returned value " << info);
  ws_size = static_cast<ordinal_type>(ws_size_query);
  work.resize(ws_size);
  lapack.ORGQR(pce_sz, pce_sz, pce_sz, K.values(), K.stride(), &tau[0], 
	       &work[0], ws_size, &info);
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error, 
		     "ORGQR returned value " << info);

  // Get basis vectors
  for (ordinal_type j=0; j<p+1; j++)
    for (ordinal_type i=0; i<pce_sz; i++)
      basis_vecs(i,j) = K(i,j);

  
  // Compute T = Q'*A*Q
  matrix_type prod(pce_sz,pce_sz);
  prod.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, K, A, 0.0);
  matrix_type T(pce_sz,pce_sz);
  T.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, prod, K, 0.0);

  //std::cout << T << std::endl;

  // Recursion coefficients are diagonal and super diagonal
  b[0] = 1.0;
  for (ordinal_type i=0; i<pce_sz-1; i++) {
    a[i] = T(i,i);
    b[i+1] = T(i,i+1);
  }
  a[pce_sz-1] = T(pce_sz-1,pce_sz-1);

  // Setup rest of basis
  this->setup();

  // Project original PCE into the new basis
  vector_type u(pce_sz);
  for (ordinal_type i=0; i<pce_sz; i++)
    u[i] = normalized_pce[i];
  new_pce.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, basis_vecs, u, 
		   0.0);
  for (ordinal_type i=0; i<=p; i++)
    new_pce[i] /= this->norms[i];
}

template <typename ordinal_type, typename value_type>
Stokhos::MonoProjPCEBasis<ordinal_type, value_type>::
~MonoProjPCEBasis()
{
}

template <typename ordinal_type, typename value_type>
void
Stokhos::MonoProjPCEBasis<ordinal_type, value_type>::
getQuadPoints(ordinal_type quad_order,
	      Teuchos::Array<value_type>& quad_points,
	      Teuchos::Array<value_type>& quad_weights,
	      Teuchos::Array< Teuchos::Array<value_type> >& quad_values) const
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::MonoProjPCEBasis -- compute Gauss points");
#endif

  // Call base class
  ordinal_type num_points = 
    static_cast<ordinal_type>(std::ceil((quad_order+1)/2.0));

  // We can't always reliably generate quadrature points of order > 2*p
  // when using sparse grids for the underlying quadrature
  if (limit_integration_order && quad_order > 2*this->p)
    quad_order = 2*this->p;
  Stokhos::RecurrenceBasis<ordinal_type,value_type>::getQuadPoints(quad_order, 
								   quad_points, 
								   quad_weights,
								   quad_values);

  // Fill in the rest of the points with zero weight
  if (quad_weights.size() < num_points) {
    ordinal_type old_size = quad_weights.size();
    quad_weights.resize(num_points);
    quad_points.resize(num_points);
    quad_values.resize(num_points);
    for (ordinal_type i=old_size; i<num_points; i++) {
      quad_weights[i] = value_type(0);
      quad_points[i] = quad_points[0];
      quad_values[i].resize(this->p+1);
      evaluateBases(quad_points[i], quad_values[i]);
    }
  }
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > 
Stokhos::MonoProjPCEBasis<ordinal_type, value_type>::
cloneWithOrder(ordinal_type p) const
{
   return Teuchos::rcp(new Stokhos::MonoProjPCEBasis<ordinal_type,value_type>(
			 p,*this));
}

template <typename ordinal_type, typename value_type>
value_type
Stokhos::MonoProjPCEBasis<ordinal_type, value_type>::
getNewCoeffs(ordinal_type i) const
{
  return new_pce[i];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::MonoProjPCEBasis<ordinal_type, value_type>::
transformCoeffs(const value_type *in, value_type *out) const
{
  // Transform coefficients to normalized basis
  Teuchos::BLAS<ordinal_type, value_type> blas;
  blas.GEMV(Teuchos::NO_TRANS, pce_sz, this->p+1, 
	    value_type(1.0), basis_vecs.values(), pce_sz, 
	    in, ordinal_type(1), value_type(0.0), out, ordinal_type(1));

  // Transform from normalized to original
  for (ordinal_type i=0; i<pce_sz; i++)
    out[i] /= pce_norms[i];
}

template <typename ordinal_type, typename value_type>
bool
Stokhos::MonoProjPCEBasis<ordinal_type, value_type>::
computeRecurrenceCoefficients(ordinal_type n,
			      Teuchos::Array<value_type>& alpha,
			      Teuchos::Array<value_type>& beta,
			      Teuchos::Array<value_type>& delta,
			      Teuchos::Array<value_type>& gamma) const
{
  // Get recurrence coefficients from the full set we already computed
  for (ordinal_type i=0; i<n; i++) {
    alpha[i] = a[i];
    beta[i] = b[i];
    delta[i] = value_type(1.0);
    gamma[i] = b[i];

    std::cout << "i = " << i << " alpha = " << alpha[i] << " beta = " << beta[i]
	      << " gamma = " << gamma[i] << std::endl;
  }

  return true;
}

template <typename ordinal_type, typename value_type> 
Stokhos::MonoProjPCEBasis<ordinal_type, value_type>::
MonoProjPCEBasis(ordinal_type p, const MonoProjPCEBasis& basis) :
  RecurrenceBasis<ordinal_type, value_type>("Lanczos PCE", p, false),
  limit_integration_order(basis.limit_integration_order),
  pce_sz(basis.pce_sz),
  pce_norms(basis.pce_norms),
  a(basis.a),
  b(basis.b),
  basis_vecs(basis.basis_vecs),
  new_pce(basis.new_pce)
{
  this->setup();
}
