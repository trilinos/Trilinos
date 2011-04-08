// $Id$ 
// $Source$ 
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

#include "Teuchos_TestForException.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_TimeMonitor.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
StieltjesPCEBasis(
   ordinal_type p,
   const Teuchos::RCP<const Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pce_,
   const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad_,
   bool use_pce_quad_points_,
   bool normalize,
   bool project_integrals_,
   const Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk_) :
  RecurrenceBasis<ordinal_type, value_type>("Stieltjes PCE", p, normalize),
  pce(pce_),
  quad(quad_),
  pce_weights(quad->getQuadWeights()),
  basis_values(quad->getBasisAtQuadPoints()),
  pce_vals(),
  phi_vals(),
  use_pce_quad_points(use_pce_quad_points_),
  fromStieltjesMat(p+1,pce->size()),
  project_integrals(project_integrals_),
  basis(pce->basis()),
  Cijk(Cijk_),
  phi_pce_coeffs()
{
  // Evaluate PCE at quad points
  const Teuchos::Array< Teuchos::Array<value_type> >& quad_points =
    quad->getQuadPoints();
  ordinal_type nqp = pce_weights.size();
  pce_vals.resize(nqp);
  phi_vals.resize(nqp);
  double mean = pce->mean();
  double std_dev = pce->standard_deviation();
  std::cout << "mean = " << mean << " std_dev = " << std_dev << std::endl;
  for (ordinal_type i=0; i<nqp; i++) {
    // pce_vals[i] = 
    //   (pce->evaluate(quad_points[i], basis_values[i])-mean)/std_dev;
    pce_vals[i] = 
      (pce->evaluate(quad_points[i], basis_values[i]));
    phi_vals[i].resize(p+1);
  }

  if (project_integrals)
    phi_pce_coeffs.resize(basis->size());
  
  // Compute coefficients via Stieltjes
  stieltjes(0, p+1, pce_weights, pce_vals, this->alpha, this->beta, this->norms,
	    phi_vals);
  for (ordinal_type i=0; i<=p; i++)
    this->delta[i] = value_type(1.0);

  ordinal_type sz = pce->size();
  fromStieltjesMat.putScalar(0.0);
  for (ordinal_type i=0; i<=p; i++) {
    for (ordinal_type j=0; j<sz; j++) {
      for (ordinal_type k=0; k<nqp; k++)
	fromStieltjesMat(i,j) += 
	  pce_weights[k]*phi_vals[k][i]*basis_values[k][j];
      fromStieltjesMat(i,j) /= basis->norm_squared(j);
    }
  }

  // Setup rest of recurrence basis
  //this->setup();
  this->gamma[0] = value_type(1);
  for (ordinal_type k=1; k<=p; k++) {
    this->gamma[k] = value_type(1);
  } 

  //If you want normalized polynomials, set gamma and reset the norms to 1.
  if( normalize ) {
    for (ordinal_type k=0; k<=p; k++) {
      this->gamma[k] = value_type(1)/std::sqrt(this->norms[k]);
      this->norms[k] = value_type(1);
    }
  }
}

template <typename ordinal_type, typename value_type>
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
~StieltjesPCEBasis()
{
}

template <typename ordinal_type, typename value_type>
void
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
getQuadPoints(ordinal_type quad_order,
	      Teuchos::Array<value_type>& quad_points,
	      Teuchos::Array<value_type>& quad_weights,
	      Teuchos::Array< Teuchos::Array<value_type> >& quad_values) const
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::StieltjesPCEBasis -- compute Gauss points");
#endif

  // Use underlying pce's quad points, weights, values
  if (use_pce_quad_points) {
    quad_points = pce_vals;
    quad_weights = pce_weights;
    quad_values = phi_vals;
    return;
  }

  // Call base class
  ordinal_type num_points = 
    static_cast<ordinal_type>(std::ceil((quad_order+1)/2.0));

  // We can't reliably generate quadrature points of order > 2*p
  //std::cout << "quad_order = " << quad_order << ", 2*p = " << 2*this->p << std::endl;
  // if (quad_order > 2*this->p)
  //   quad_order = 2*this->p;
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
void
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
computeRecurrenceCoefficients(ordinal_type n,
			      Teuchos::Array<value_type>& alpha,
			      Teuchos::Array<value_type>& beta,
			      Teuchos::Array<value_type>& delta) const
{
  ordinal_type nqp = phi_vals.size();
  Teuchos::Array<value_type> nrm(n);
  Teuchos::Array< Teuchos::Array<value_type> > vals(nqp);
  for (ordinal_type i=0; i<nqp; i++)
    vals[i].resize(n);
  stieltjes(0, n, pce_weights, pce_vals, alpha, beta, nrm, vals);
  for (ordinal_type i=0; i<n; i++)
    delta[i] = value_type(1.0);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
stieltjes(ordinal_type nstart,
	  ordinal_type nfinish,
	  const Teuchos::Array<value_type>& weights,
	  const Teuchos::Array<value_type>& points,
	  Teuchos::Array<value_type>& a,
	  Teuchos::Array<value_type>& b,
	  Teuchos::Array<value_type>& nrm,
	  Teuchos::Array< Teuchos::Array<value_type> >& phi_vals) const
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::StieltjesPCEBasis -- Discretized Stieltjes Procedure");
#endif

  value_type val1, val2;   
  ordinal_type start = nstart;
  if (nstart == 0) {
    if (project_integrals)
      integrateBasisSquaredProj(0, a, b, weights, points, phi_vals, val1, val2);
    else
      integrateBasisSquared(0, a, b, weights, points, phi_vals, val1, val2);
    nrm[0] = val1;
    a[0] = val2/val1;
    b[0] = value_type(1);
    start = 1;
  }
  for (ordinal_type i=start; i<nfinish; i++) {
    if (project_integrals)
      integrateBasisSquaredProj(i, a, b, weights, points, phi_vals, val1, val2);
    else
      integrateBasisSquared(i, a, b, weights, points, phi_vals, val1, val2);
    // std::cout << "i = " << i << " val1 = " << val1 << " val2 = " << val2
    // 	      << std::endl;
    TEST_FOR_EXCEPTION(val1 < 0, std::logic_error,
		     "Stokhos::StieltjesPCEBasis::stieltjes():  "
		       << " Polynomial " << i << " out of " << nfinish 
		       << " has norm " << val1 
		       << "!  Try increasing number of quadrature points");
    nrm[i] = val1;
    a[i] = val2/val1;
    b[i] = nrm[i]/nrm[i-1];
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
integrateBasisSquared(ordinal_type k, 
		      const Teuchos::Array<value_type>& a,
		      const Teuchos::Array<value_type>& b,
		      const Teuchos::Array<value_type>& weights,
		      const Teuchos::Array<value_type>& points,
		      Teuchos::Array< Teuchos::Array<value_type> >& phi_vals,
		      value_type& val1, value_type& val2) const
{
  evaluateRecurrence(k, a, b, points, phi_vals);
  ordinal_type nqp = weights.size();
  val1 = value_type(0);
  val2 = value_type(0);
  for (ordinal_type i=0; i<nqp; i++) {
    val1 += weights[i]*phi_vals[i][k]*phi_vals[i][k];
    val2 += weights[i]*phi_vals[i][k]*phi_vals[i][k]*points[i];
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
evaluateRecurrence(ordinal_type k,
		   const Teuchos::Array<value_type>& a,
		   const Teuchos::Array<value_type>& b,
		   const Teuchos::Array<value_type>& points,
		   Teuchos::Array< Teuchos::Array<value_type> >& values) const
{
  ordinal_type np = points.size();
  if (k == 0)
    for (ordinal_type i=0; i<np; i++)
      values[i][k] = 1.0;
  else if (k == 1)
    for (ordinal_type i=0; i<np; i++)
      values[i][k] = points[i] - a[k-1];
  else
    for (ordinal_type i=0; i<np; i++)
      values[i][k] = 
	(points[i] - a[k-1])*values[i][k-1] - b[k-1]*values[i][k-2];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
integrateBasisSquaredProj(
  ordinal_type k, 
  const Teuchos::Array<value_type>& a,
  const Teuchos::Array<value_type>& b,
  const Teuchos::Array<value_type>& weights,
  const Teuchos::Array<value_type>& points,
  Teuchos::Array< Teuchos::Array<value_type> >& phi_vals,
  value_type& val1, value_type& val2) const
{
  ordinal_type nqp = weights.size();
  ordinal_type npc = basis->size();
  const Teuchos::Array<value_type>& norms = basis->norm_squared();

  // Compute PC expansion of phi_k in original basis
  evaluateRecurrence(k, a, b, points, phi_vals);
  for (ordinal_type j=0; j<npc; j++) {
    value_type c = value_type(0);
    for (ordinal_type i=0; i<nqp; i++)
      c += weights[i]*phi_vals[i][k]*basis_values[i][j];
    c /= norms[j];
    phi_pce_coeffs[j] = c;
  }

  // Compute \int phi_k^2(\eta) d\eta
  val1 = value_type(0);
  for (ordinal_type j=0; j<npc; j++)
    val1 += phi_pce_coeffs[j]*phi_pce_coeffs[j]*norms[j];

  // Compute \int \eta phi_k^2(\eta) d\eta
  val2 = value_type(0);
  for (typename Cijk_type::k_iterator k_it = Cijk->k_begin();
       k_it != Cijk->k_end(); ++k_it) {
    ordinal_type l = index(k_it);
    for (typename Cijk_type::kj_iterator j_it = Cijk->j_begin(k_it); 
	 j_it != Cijk->j_end(k_it); ++j_it) {
      ordinal_type j = index(j_it);
      for (typename Cijk_type::kji_iterator i_it = Cijk->i_begin(j_it);
	   i_it != Cijk->i_end(j_it); ++i_it) {
	ordinal_type i = index(i_it);
	value_type c = value(i_it);
	val2 += phi_pce_coeffs[i]*phi_pce_coeffs[j]*(*pce)[l]*c;
      }
    }
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
transformCoeffsFromStieltjes(const value_type *in, value_type *out) const
{
  Teuchos::BLAS<ordinal_type, value_type> blas;
  blas.GEMV(Teuchos::TRANS, fromStieltjesMat.numRows(), 
	    fromStieltjesMat.numCols(), 1.0, fromStieltjesMat.values(), 
	    fromStieltjesMat.numRows(), in, 1, 0.0, out, 1);
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > 
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
cloneWithOrder(ordinal_type p) const
{
   return Teuchos::rcp(new Stokhos::StieltjesPCEBasis<ordinal_type,value_type>(p,pce,quad,use_pce_quad_points,this->normalize,project_integrals,Cijk));
}
