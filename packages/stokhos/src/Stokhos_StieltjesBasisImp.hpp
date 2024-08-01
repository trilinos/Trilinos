// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Assert.hpp"
#include "Teuchos_TimeMonitor.hpp"

template <typename ordinal_type, typename value_type, typename func_type>
Stokhos::StieltjesBasis<ordinal_type, value_type, func_type>::
StieltjesBasis(
   ordinal_type p,
   const Teuchos::RCP<const func_type>& func_,
   const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad_,
   bool use_pce_quad_points_,
   bool normalize) :
  RecurrenceBasis<ordinal_type, value_type>("Stieltjes PCE", p, normalize),
  func(func_),
  quad(quad_),
  pce_weights(quad->getQuadWeights()),
  basis_values(quad->getBasisAtQuadPoints()),
  func_vals(),
  phi_vals(),
  use_pce_quad_points(use_pce_quad_points_)
{
  // Setup rest of recurrence basis
  this->setup();
}

template <typename ordinal_type, typename value_type, typename func_type>
Stokhos::StieltjesBasis<ordinal_type, value_type, func_type>::
~StieltjesBasis()
{
}

template <typename ordinal_type, typename value_type, typename func_type>
void
Stokhos::StieltjesBasis<ordinal_type, value_type, func_type>::
getQuadPoints(ordinal_type quad_order,
	      Teuchos::Array<value_type>& quad_points,
	      Teuchos::Array<value_type>& quad_weights,
	      Teuchos::Array< Teuchos::Array<value_type> >& quad_values) const
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::StieltjesBasis -- compute Gauss points");
#endif

  // Use underlying pce's quad points, weights, values
  if (use_pce_quad_points) {
    quad_points = func_vals;
    quad_weights = pce_weights;
    quad_values = phi_vals;
    return;
  }

  // Call base class
  ordinal_type num_points = 
    static_cast<ordinal_type>(std::ceil((quad_order+1)/2.0));

  // We can't reliably generate quadrature points of order > 2*p
  if (quad_order > 2*this->p)
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

template <typename ordinal_type, typename value_type, typename func_type>
bool
Stokhos::StieltjesBasis<ordinal_type, value_type, func_type>::
computeRecurrenceCoefficients(ordinal_type n,
			      Teuchos::Array<value_type>& alpha,
			      Teuchos::Array<value_type>& beta,
			      Teuchos::Array<value_type>& delta,
			      Teuchos::Array<value_type>& gamma) const
{
  ordinal_type nqp = phi_vals.size();
  Teuchos::Array<value_type> nrm(n);
  Teuchos::Array< Teuchos::Array<value_type> > vals(nqp);
  for (ordinal_type i=0; i<nqp; i++)
    vals[i].resize(n);
  stieltjes(0, n, pce_weights, func_vals, alpha, beta, nrm, vals);
  for (ordinal_type i=0; i<n; i++) {
    delta[i] = value_type(1.0);
    gamma[i] = value_type(1.0);
  }

  // Save basis functions at quad point values
  if (n == this->p+1)
    phi_vals = vals;

  return false;
}

template <typename ordinal_type, typename value_type, typename func_type>
void
Stokhos::StieltjesBasis<ordinal_type, value_type, func_type>::
setup() 
{
  // Evaluate func at quad points
  const Teuchos::Array< Teuchos::Array<value_type> >& quad_points =
    quad->getQuadPoints();
  ordinal_type nqp = pce_weights.size();
  func_vals.resize(nqp);
  phi_vals.resize(nqp);
  for (ordinal_type i=0; i<nqp; i++) {
    func_vals[i] = (*func)(quad_points[i]);
    phi_vals[i].resize(this->p+1);
  }

  RecurrenceBasis<ordinal_type,value_type>::setup();
}

template <typename ordinal_type, typename value_type, typename func_type>
void
Stokhos::StieltjesBasis<ordinal_type, value_type, func_type>::
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
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos::StieltjesBasis -- Discretized Stieltjes Procedure");
#endif

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
    // std::cout << "i = " << i << " val1 = " << val1 << " val2 = " << val2
    // 	      << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(val1 < 0.0, std::logic_error,
		     "Stokhos::StieltjesBasis::stieltjes():  "
		       << " Polynomial " << i << " out of " << nfinish 
		       << " has norm " << val1 
		       << "!  Try increasing number of quadrature points");
    nrm[i] = val1;
    a[i] = val2/val1;
    b[i] = nrm[i]/nrm[i-1];
  }
}

template <typename ordinal_type, typename value_type, typename func_type>
void
Stokhos::StieltjesBasis<ordinal_type, value_type, func_type>::
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

template <typename ordinal_type, typename value_type, typename func_type>
void
Stokhos::StieltjesBasis<ordinal_type, value_type, func_type>::
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

template <typename ordinal_type, typename value_type, typename func_type>
Teuchos::RCP<Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > 
Stokhos::StieltjesBasis<ordinal_type, value_type, func_type>::
cloneWithOrder(ordinal_type p) const
{
   return Teuchos::rcp(new Stokhos::StieltjesBasis<ordinal_type,value_type,func_type>(p,*this));
}

template <typename ordinal_type, typename value_type, typename func_type>
Stokhos::StieltjesBasis<ordinal_type, value_type, func_type>::
StieltjesBasis(ordinal_type p, const StieltjesBasis& basis) :
  RecurrenceBasis<ordinal_type, value_type>(p, basis),
  func(basis.func),
  quad(basis.quad),
  pce_weights(quad->getQuadWeights()),
  basis_values(quad->getBasisAtQuadPoints()),
  func_vals(),
  phi_vals(),
  use_pce_quad_points(basis.use_pce_quad_points)
{
  this->setup();
}
