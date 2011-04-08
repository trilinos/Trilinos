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
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
LanczosPCEBasis(
   ordinal_type p,
   const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& pce,
   const Stokhos::Quadrature<ordinal_type, value_type>& quad) :
  RecurrenceBasis<ordinal_type, value_type>("Lanczos PCE", p, false)
{
  // Evaluate PCE at quad points
  pce_weights = quad.getQuadWeights();
  ordinal_type nqp = pce_weights.size();
  pce_vals.resize(nqp);
  const Teuchos::Array< Teuchos::Array<value_type> >& quad_points =
    quad.getQuadPoints();
  const Teuchos::Array< Teuchos::Array<value_type> >& basis_values =
    quad.getBasisAtQuadPoints();
  for (ordinal_type i=0; i<nqp; i++)
    pce_vals[i] = pce.evaluate(quad_points[i], basis_values[i]);
  
  // Compute coefficients via Lanczos
  for (ordinal_type i=0; i<nqp; i++)
    pce_weights[i] = std::sqrt(pce_weights[i]);
  Teuchos::Array<value_type> g(p+1);
  lanczos(nqp, p, pce_vals, pce_weights, this->alpha, g);
  for (ordinal_type i=0; i<p; i++) {
    this->alpha[i] = this->alpha[i]/std::sqrt(g[i+1]);
    this->beta[i] = std::sqrt(g[i])/std::sqrt(g[i+1]);
    this->delta[i] = value_type(1.0)/std::sqrt(g[i+1]);
    this->gamma[i] = value_type(1);
  }
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
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
computeRecurrenceCoefficients(ordinal_type n,
			      Teuchos::Array<value_type>& alpha,
			      Teuchos::Array<value_type>& beta,
			      Teuchos::Array<value_type>& delta) const
{
  Teuchos::Array<value_type> g(n+1);
  lanczos(nqp, n, pce_vals, pce_weights, alpha, g);
  for (ordinal_type i=0; i<n; i++) {
    alpha[i] = alpha[i]/std::sqrt(g[i+1]);
    beta[i] = std::sqrt(g[i])/std::sqrt(g[i+1]);
    delta[i] = value_type(1.0)/std::sqrt(g[i+1]);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::LanczosPCEBasis<ordinal_type, value_type>::
lanczos(ordinal_type n,
	ordinal_type nsteps,
	const Teuchos::Array<value_type>& A_diag,
	const Teuchos::Array<value_type>& h0,
	Teuchos::Array<value_type>& a,
	Teuchos::Array<value_type>& g) const
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::LanczosPCEBasis -- Lanczos Procedure");
#endif

  Teuchos::Array< Teuchos::Array<value_type> > h(nsteps+1);
  for (ordinal_type i=0; i<nsteps+1; i++)
    h[i].resize(n);
  h[0] = h0;
  g[0] = 1.0;

  for (ordinal_type i=0; i<nsteps; i++) {
    // compute alpha = h_i^T*A*h_i
    value_type alpha = 0.0;
    for (ordinal_type j=0; j<n; j++)
      alpha += h[i][j]*A_diag[j]*h[i][j];
    a[i] = alpha;

    // compute
    for (ordinal_type j=0; j<n; j++)
      h[i+1][j] = (A_diag[j] - alpha)*h[i][j] - gamma[i]*h[i-1][j];

    // compute gamma = ||h_{i+1}||
    value_type gamma = 0.0;
    for (ordinal_type j=0; j<n; j++)
      gamma += h[i+1][j]*h[i+1][j];
    g[i+1] = gamma;

    // compute h_{i+1}/gamma
    for (ordinal_type j=0; j<n; j++)
      h[i+1][j] /= gamma;

    std::cout << "i = " << i << " alpha = " << alpha << " gamma = " << gamma
    	      << std::endl;
    TEST_FOR_EXCEPTION(gamma < 0, std::logic_error,
		     "Stokhos::LanczosPCEBasis::lanczos():  "
		       << " Polynomial " << i << " out of " << nsteps
		       << " has norm " << gamma
		       << "!  Try increasing number of quadrature points");
  }
}
