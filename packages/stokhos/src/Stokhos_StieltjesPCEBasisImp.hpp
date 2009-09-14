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

template <typename ordinal_type, typename value_type>
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
StieltjesPCEBasis(
   ordinal_type p,
   const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& pce,
   const Stokhos::Quadrature<ordinal_type, value_type>& quad,
   bool use_pce_quad_points_,
   bool normalize) :
  RecurrenceBasis<ordinal_type, value_type>("Stieltjes PCE", p, normalize),
  pce_vals(),
  phi_vals(),
  use_pce_quad_points(use_pce_quad_points_),
  fromStieltjesMat(p+1,pce.size())
{
  // Evaluate PCE at quad points
  const Teuchos::Array< Teuchos::Array<value_type> >& quad_points =
    quad.getQuadPoints();
  pce_weights = quad.getQuadWeights();
  const Teuchos::Array< Teuchos::Array<value_type> >& basis_values =
    quad.getBasisAtQuadPoints();
  ordinal_type nqp = pce_weights.size();
  pce_vals.resize(nqp);
  phi_vals.resize(nqp);
  for (ordinal_type i=0; i<nqp; i++) {
    pce_vals[i] = pce.evaluate(quad_points[i], basis_values[i]);
    phi_vals[i].resize(p+1);
  }
  
  // Compute coefficients via Stieltjes
  stieltjes(0, p+1, pce_weights, pce_vals, this->alpha, this->beta, this->norms,
	    phi_vals);
  for (ordinal_type i=0; i<=p; i++)
    this->delta[i] = value_type(1.0);

  ordinal_type sz = pce.size();
  fromStieltjesMat.putScalar(0.0);
  for (ordinal_type i=0; i<=p; i++) {
    for (ordinal_type j=0; j<sz; j++) {
      for (ordinal_type k=0; k<nqp; k++)
	fromStieltjesMat(i,j) += 
	  pce_weights[k]*phi_vals[k][i]*basis_values[k][j];
      fromStieltjesMat(i,j) /= pce.basis()->norm_squared(j);
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
  // Use underlying pce's quad points, weights, values
  if (use_pce_quad_points) {
    quad_points = pce_vals;
    quad_weights = pce_weights;
    quad_values = phi_vals;
    return;
  }

  // Call base class
  Stokhos::RecurrenceBasis<ordinal_type,value_type>::getQuadPoints(quad_order, 
								   quad_points, 
								   quad_weights,
								   quad_values);
}

template <typename ordinal_type, typename value_type>
ordinal_type 
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
getRule() const
{
  return 10;
}

template <typename ordinal_type, typename value_type>
value_type 
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
getQuadWeightFactor() const
{
  return 1.0;
}

template <typename ordinal_type, typename value_type>
value_type 
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
getQuadPointFactor() const
{
  return 1.0;
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
  value_type val1, val2;   
  ordinal_type start = nstart;
  if (nstart == 0) {
    integrateBasisSquared(0, a, b, weights, points, phi_vals, val1, val2);
    nrm[0] = val1;
    a[0] = val2/val1;
    b[0] = value_type(1);
    start = 1;
  }
  for (ordinal_type i=start; i<nfinish; i++) {
    integrateBasisSquared(i, a, b, weights, points, phi_vals, val1, val2);
    //std::cout << "i = " << i << " val1 = " << val1 << " val2 = " << val2
    //	      << std::endl;
    TEST_FOR_EXCEPTION(val1 < 1.0e-14, std::logic_error,
		     "Stokhos::StieltjesPCEBasis::stieltjes():  "
		     << " Polynomial is zero!  Try increasing number of quadrature points");
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
transformCoeffsFromStieltjes(const value_type *in, value_type *out) const
{
  Teuchos::BLAS<ordinal_type, value_type> blas;
  blas.GEMV(Teuchos::TRANS, fromStieltjesMat.numRows(), 
	    fromStieltjesMat.numCols(), 1.0, fromStieltjesMat.values(), 
	    fromStieltjesMat.numRows(), in, 1, 0.0, out, 1);
}
