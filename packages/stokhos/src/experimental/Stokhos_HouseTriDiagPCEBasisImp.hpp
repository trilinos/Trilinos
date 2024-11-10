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
#include "Teuchos_TimeMonitor.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::HouseTriDiagPCEBasis<ordinal_type, value_type>::
HouseTriDiagPCEBasis(
  ordinal_type p,
  const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& pce,
  const Stokhos::Sparse3Tensor<ordinal_type, value_type>& Cijk,
  bool limit_integration_order_) :
  RecurrenceBasis<ordinal_type, value_type>("Householder Tridiagonalization PCE", p, true),
  limit_integration_order(limit_integration_order_),
  pce_sz(pce.basis()->size()),
  a(pce_sz+1), 
  b(pce_sz),
  basis_vecs(pce_sz, p+1),
  new_pce(p+1)
{
  pce_norms = pce.basis()->norm_squared();
  
  // Compute matrix -- rescale basis to unit norm
  matrix_type A(pce_sz+1, pce_sz+1);
  A(0,0) = 1.0;
  for (ordinal_type i=0; i<pce_sz; i++) {
    A(0,i+1) = 1.0/pce_sz;
    A(i+1,0) = 1.0/pce_sz;
  }
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
	value_type c = value(i_it) / std::sqrt(pce_norms[i]*pce_norms[j]*pce_norms[k]);
	val += std::sqrt(pce_norms[i])*pce[i]*c;
      }
      A(k+1,j+1) = val;
    }
  }

  // Call LAPACK Householder tridiagonalization algorithm
  // Householder vectors are stored in lower-triangular part
  ordinal_type ws_size, info;
  value_type ws_size_query;
  Teuchos::Array<value_type> tau(pce_sz-1);
  lapack.SYTRD('L', pce_sz+1, A.values(), A.stride(), &a[0], &b[0], &tau[0], 
	       &ws_size_query, -1, &info);
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error, 
		     "SYTRD returned value " << info);
  ws_size = static_cast<ordinal_type>(ws_size_query);
  Teuchos::Array<value_type> work(ws_size);
  lapack.SYTRD('L', pce_sz+1, A.values(), A.stride(), &a[0], &b[0], &tau[0], 
	       &work[0], ws_size, &info);
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error, 
		     "SYTRD returned value " << info);

  // Set sub-diagonal terms to zero
  for (ordinal_type j=0; j<pce_sz; j++)
    A(j+1,j) = 0.0;
  
  // Now compute orthogonal matrix
  lapack.ORGQR(pce_sz+1, pce_sz+1, pce_sz-1, A.values(), A.stride(), &tau[0], 
	       &ws_size_query, -1, &info);
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error, 
		     "ORGQR returned value " << info);
  ws_size = static_cast<ordinal_type>(ws_size_query);
  work.resize(ws_size);
  lapack.ORGQR(pce_sz+1, pce_sz+1, pce_sz-1, A.values(), A.stride(), &tau[0], 
	       &work[0], ws_size, &info);
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error, 
		     "ORGQR returned value " << info);

  // Get basis vectors
  for (ordinal_type j=0; j<p+1; j++)
    for (ordinal_type i=0; i<pce_sz; i++)
      basis_vecs(i,j) = A(i+1,j+1);

  // Setup of rest of recurrence basis
  this->setup();

  // Project original PCE into the new basis
  vector_type u(pce_sz);
  for (ordinal_type i=0; i<pce_sz; i++)
    u[i] = pce[i]*std::sqrt(pce_norms[i]);
  new_pce.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, basis_vecs, u, 
		   0.0);
  for (ordinal_type i=0; i<=p; i++)
    new_pce[i] /= this->norms[i];

  std::cout << new_pce << std::endl;

  // Test orthogonality of basis vectors
  matrix_type prod(p+1,p+1);
  prod.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, basis_vecs, basis_vecs, 
		0.0);
  for (ordinal_type j=0; j<p+1; j++)
    prod(j,j) -= 1.0;
  std::cout << "orthogonalization error = " << prod.normInf() << std::endl;
  //std::cout << prod << std::endl;
}

template <typename ordinal_type, typename value_type>
Stokhos::HouseTriDiagPCEBasis<ordinal_type, value_type>::
~HouseTriDiagPCEBasis()
{
}

template <typename ordinal_type, typename value_type>
void
Stokhos::HouseTriDiagPCEBasis<ordinal_type, value_type>::
getQuadPoints(ordinal_type quad_order,
	      Teuchos::Array<value_type>& quad_points,
	      Teuchos::Array<value_type>& quad_weights,
	      Teuchos::Array< Teuchos::Array<value_type> >& quad_values) const
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::LanczosPCEBasis -- compute Gauss points");
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
Stokhos::HouseTriDiagPCEBasis<ordinal_type, value_type>::
cloneWithOrder(ordinal_type p) const
{
   return 
     Teuchos::rcp(new Stokhos::HouseTriDiagPCEBasis<ordinal_type,value_type>(
		    p,*this));
}

template <typename ordinal_type, typename value_type>
value_type
Stokhos::HouseTriDiagPCEBasis<ordinal_type, value_type>::
getNewCoeffs(ordinal_type i) const
{
  return new_pce[i];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::HouseTriDiagPCEBasis<ordinal_type, value_type>::
transformCoeffsFromHouse(const value_type *in, value_type *out) const
{
  blas.GEMV(Teuchos::NO_TRANS, pce_sz, this->p+1, 
	    value_type(1.0), basis_vecs.values(), pce_sz, 
	    in, ordinal_type(1), value_type(0.0), out, ordinal_type(1));
  for (ordinal_type i=0; i<pce_sz; i++)
    out[i] /= pce_norms[i];
}

template <typename ordinal_type, typename value_type>
bool
Stokhos::HouseTriDiagPCEBasis<ordinal_type, value_type>::
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
Stokhos::HouseTriDiagPCEBasis<ordinal_type, value_type>::
HouseTriDiagPCEBasis(ordinal_type p, const HouseTriDiagPCEBasis& basis) :
  RecurrenceBasis<ordinal_type, value_type>("Householder Tridiagonalization PCE", p, false),
  limit_integration_order(basis.limit_integration_order),
  pce_sz(basis.pce_sz),
  a(basis.a),
  b(basis.b),
  basis_vecs(basis.basis_vecs),
  new_pce(basis.new_pce)
{
  this->setup();
}
