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
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
LanczosPCEBasis(
   ordinal_type p,
   const Teuchos::RCP< const Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pce_,
   const Teuchos::RCP< const Stokhos::Quadrature<ordinal_type, value_type> >& quad_,
   bool normalize,
   bool limit_integration_order_) :
  RecurrenceBasis<ordinal_type, value_type>("Lanczos PCE", p, normalize),
  pce(pce_),
  quad(quad_),
  limit_integration_order(limit_integration_order_),
  nqp(quad->size()),
  pce_weights(Teuchos::Copy, 
	      const_cast<value_type*>(quad->getQuadWeights().getRawPtr()), 
	      nqp),
  pce_vals(nqp),
  u0(nqp),
  lanczos_vecs(nqp, p+1),
  fromStieltjesMat(),
  new_pce()
{
  // Evaluate PCE at quad points
  const Teuchos::Array< Teuchos::Array<value_type> >& quad_points =
    quad->getQuadPoints();
  const Teuchos::Array< Teuchos::Array<value_type> >& basis_values =
    quad->getBasisAtQuadPoints();
  for (ordinal_type i=0; i<nqp; i++) {
    pce_vals[i] = pce->evaluate(quad_points[i], basis_values[i]);
    u0[i] = value_type(1);
  }

  // Setup rest of basis
  this->setup();
}

template <typename ordinal_type, typename value_type>
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
~LanczosPCEBasis()
{
}

template <typename ordinal_type, typename value_type>
void
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
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
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
cloneWithOrder(ordinal_type p) const
{
   return Teuchos::rcp(new Stokhos::LanczosPCEBasis<ordinal_type,value_type>(
			 p,*this));
}

template <typename ordinal_type, typename value_type>
value_type
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
getNewCoeffs(ordinal_type i) const
{
  return new_pce[i];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
transformCoeffsFromLanczos(const value_type *in, value_type *out) const
{
  Teuchos::BLAS<ordinal_type, value_type> blas;
  ordinal_type sz = fromStieltjesMat.numRows();
  blas.GEMV(Teuchos::NO_TRANS, sz, this->p+1, 
	    value_type(1.0), fromStieltjesMat.values(), sz, 
	    in, ordinal_type(1), value_type(0.0), out, ordinal_type(1));
}

template <typename ordinal_type, typename value_type>
bool
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
computeRecurrenceCoefficients(ordinal_type n,
			      Teuchos::Array<value_type>& alpha,
			      Teuchos::Array<value_type>& beta,
			      Teuchos::Array<value_type>& delta,
			      Teuchos::Array<value_type>& gamma) const
{
  Teuchos::Array<value_type> nrm(n);
  vectorspace_type vs(pce_weights);
  operator_type A(pce_vals);
  
  // Create space to store lanczos vectors -- use lanczos_vecs if 
  // we are requesting p+1 vectors
  Teuchos::RCP<matrix_type> lv;
  if (n == this->p+1)
    lv = Teuchos::rcp(&lanczos_vecs, false);
  else
    lv = Teuchos::rcp(new matrix_type(nqp,n));

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

  return this->normalize;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
setup() 
{
  RecurrenceBasis<ordinal_type,value_type>::setup();

  // Compute transformation matrix back to original basis
  ordinal_type sz = pce->size();
  fromStieltjesMat.shape(sz, this->p+1);
  fromStieltjesMat.putScalar(0.0);
  const Teuchos::Array< Teuchos::Array<value_type> >& basis_values =
    quad->getBasisAtQuadPoints();
  for (ordinal_type i=0; i<sz; i++) {
    for (ordinal_type j=0; j<=this->p; j++) {
      for (ordinal_type k=0; k<nqp; k++)
	fromStieltjesMat(i,j) += 
	  pce_weights[k]*lanczos_vecs(k,j)*basis_values[k][i];
      fromStieltjesMat(i,j) /= pce->basis()->norm_squared(i);
    }
  }

  // Project original PCE into the new basis
  new_pce.resize(this->p+1);
  vector_type u(sz);
  for (ordinal_type i=0; i<sz; i++)
    u[i] = (*pce)[i]*pce->basis()->norm_squared(i);
  new_pce.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, fromStieltjesMat, u, 
		   0.0);
  for (ordinal_type i=0; i<=this->p; i++)
    new_pce[i] /= this->norms[i];
}

template <typename ordinal_type, typename value_type> 
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
LanczosPCEBasis(ordinal_type p, const LanczosPCEBasis& basis) :
  RecurrenceBasis<ordinal_type, value_type>(p, basis),
  pce(basis.pce),
  quad(basis.quad),
  limit_integration_order(basis.limit_integration_order),
  nqp(basis.nqp),
  pce_weights(basis.pce_weights),
  pce_vals(basis.pce_vals),
  u0(basis.u0),
  lanczos_vecs(nqp, p+1),
  fromStieltjesMat(),
  new_pce()
{
  this->setup();
}
