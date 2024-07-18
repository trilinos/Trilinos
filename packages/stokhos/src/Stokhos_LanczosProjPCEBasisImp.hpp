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
Stokhos::LanczosProjPCEBasis<ordinal_type, value_type>::
LanczosProjPCEBasis(
  ordinal_type p,
  const Teuchos::RCP< const Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pce_,
  const Teuchos::RCP< const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk,
  bool normalize,
  bool limit_integration_order_) :
  RecurrenceBasis<ordinal_type, value_type>("Lanczos-proj PCE", p, normalize),
  pce(pce_),
  limit_integration_order(limit_integration_order_),
  pce_sz(pce->basis()->size()),
  Cijk_matrix(pce_sz,pce_sz),
  weights(Teuchos::Copy, 
	  const_cast<value_type*>(pce->basis()->norm_squared().getRawPtr()), 
	  pce_sz),
  u0(pce_sz),
  lanczos_vecs(pce_sz, p+1),
  new_pce(p+1)
{
  u0[0] = value_type(1);

  pce_norms = pce->basis()->norm_squared();
  for (ordinal_type i=0; i<pce_sz; i++) {
    pce_norms[i] = std::sqrt(pce_norms[i]);
    weights[i] = value_type(1);
  }
  
  // Compute matrix -- For the matrix to be symmetric, the original basis
  // must be normalized.  However we don't want to require this, so we
  // rescale the pce coefficients for a normalized basis
  typedef Stokhos::Sparse3Tensor<ordinal_type, value_type> Cijk_type;
  for (typename Cijk_type::k_iterator k_it = Cijk->k_begin();
       k_it != Cijk->k_end(); ++k_it) {
    ordinal_type k = index(k_it);
    for (typename Cijk_type::kj_iterator j_it = Cijk->j_begin(k_it); 
	 j_it != Cijk->j_end(k_it); ++j_it) {
      ordinal_type j = index(j_it);
      value_type val = 0;
      for (typename Cijk_type::kji_iterator i_it = Cijk->i_begin(j_it);
	   i_it != Cijk->i_end(j_it); ++i_it) {
	ordinal_type i = index(i_it);
	value_type c = value(i_it);
	val += (*pce)[i]*c / (pce_norms[j]*pce_norms[k]);
      }
      Cijk_matrix(k,j) = val;
    }
  }
  
  // Setup of rest of recurrence basis
  this->setup();

  
}

template <typename ordinal_type, typename value_type>
Stokhos::LanczosProjPCEBasis<ordinal_type, value_type>::
~LanczosProjPCEBasis()
{
}

template <typename ordinal_type, typename value_type>
void
Stokhos::LanczosProjPCEBasis<ordinal_type, value_type>::
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
      this->evaluateBases(quad_points[i], quad_values[i]);
    }
  }
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > 
Stokhos::LanczosProjPCEBasis<ordinal_type, value_type>::
cloneWithOrder(ordinal_type p) const
{
   return 
     Teuchos::rcp(new Stokhos::LanczosProjPCEBasis<ordinal_type,value_type>(
		    p,*this));
}

template <typename ordinal_type, typename value_type>
value_type
Stokhos::LanczosProjPCEBasis<ordinal_type, value_type>::
getNewCoeffs(ordinal_type i) const
{
  return new_pce[i];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::LanczosProjPCEBasis<ordinal_type, value_type>::
transformCoeffsFromLanczos(const value_type *in, value_type *out) const
{
  // Transform coefficients to normalized basis
  Teuchos::BLAS<ordinal_type, value_type> blas;
  blas.GEMV(Teuchos::NO_TRANS, pce_sz, this->p+1, 
	    value_type(1.0), lanczos_vecs.values(), pce_sz, 
	    in, ordinal_type(1), value_type(0.0), out, ordinal_type(1));

  // Transform from normalized to original
  for (ordinal_type i=0; i<pce_sz; i++)
    out[i] /= pce_norms[i];
}

template <typename ordinal_type, typename value_type>
bool
Stokhos::LanczosProjPCEBasis<ordinal_type, value_type>::
computeRecurrenceCoefficients(ordinal_type n,
			      Teuchos::Array<value_type>& alpha,
			      Teuchos::Array<value_type>& beta,
			      Teuchos::Array<value_type>& delta,
			      Teuchos::Array<value_type>& gamma) const
{
  Teuchos::Array<value_type> nrm(n);
  vectorspace_type vs(weights);
  operator_type A(Cijk_matrix);

  // Create space to store lanczos vectors -- use lanczos_vecs if 
  // we are requesting p+1 vectors
  Teuchos::RCP<matrix_type> lv;
  if (n == this->p+1)
    lv = Teuchos::rcp(&lanczos_vecs, false);
  else
    lv = Teuchos::rcp(new matrix_type(pce_sz,n));

  if (this->normalize)
    lanczos_type::computeNormalized(n, vs, A, u0, *lv, alpha, beta, nrm);
  else
    lanczos_type::compute(n, vs, A, u0, *lv, alpha, beta, nrm);

  for (ordinal_type i=0; i<n; i++) {
    delta[i] = value_type(1.0);
  }
  if (this->normalize)
    gamma = beta;
  else
    for (ordinal_type i=0; i<n; i++)
      gamma[i] = value_type(1.0);

  /*
  matrix_type slv(pce_sz, n);
  for (ordinal_type j=0; j<n; j++)
    for (ordinal_type i=0; i<pce_sz; i++)
      slv(i,j) = (*lv)(i,j) * weights[i];
  matrix_type prod(n,n);
  prod.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, *lv, slv, 0.0);
  for (ordinal_type j=0; j<n; j++) {
    for (ordinal_type i=0; i<n; i++)
      prod(i,j) /= std::sqrt(nrm[i]*nrm[j]);
    prod(j,j) -= 1.0;
  }
  std::cout << "orthogonalization error = " << prod.normInf() << std::endl;
  */

  return this->normalize;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::LanczosProjPCEBasis<ordinal_type, value_type>::
setup() 
{
  RecurrenceBasis<ordinal_type,value_type>::setup();

  // Project original PCE into the new basis
  vector_type u(pce_sz);
  for (ordinal_type i=0; i<pce_sz; i++)
    u[i] = (*pce)[i]*pce_norms[i];
  new_pce.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, lanczos_vecs, u, 
		   0.0);
  for (ordinal_type i=0; i<=this->p; i++)
    new_pce[i] /= this->norms[i];
}

template <typename ordinal_type, typename value_type> 
Stokhos::LanczosProjPCEBasis<ordinal_type, value_type>::
LanczosProjPCEBasis(ordinal_type p, const LanczosProjPCEBasis& basis) :
  RecurrenceBasis<ordinal_type, value_type>("Lanczos-proj PCE", p, false),
  pce(basis.pce),
  limit_integration_order(basis.limit_integration_order),
  pce_sz(basis.pce_sz),
  pce_norms(basis.pce_norms),
  Cijk_matrix(basis.Cijk_matrix),
  weights(basis.weights),
  u0(basis.u0),
  lanczos_vecs(pce_sz, p+1),
  new_pce()
{
  this->setup();
}
