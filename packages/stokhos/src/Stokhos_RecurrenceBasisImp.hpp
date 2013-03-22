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

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::RecurrenceBasis<ordinal_type, value_type>::
RecurrenceBasis(const std::string& name_, ordinal_type p_, bool normalize_,
		Stokhos::GrowthPolicy growth_) :
  name(name_),
  p(p_),
  normalize(normalize_),
  growth(growth_),
  quad_zero_tol(1.0e-14),
#ifdef HAVE_STOKHOS_DAKOTA
  sparse_grid_growth_rule(webbur::level_to_order_linear_nn),
#else
  sparse_grid_growth_rule(NULL),
#endif
  alpha(p+1, value_type(0.0)),
  beta(p+1, value_type(0.0)),
  delta(p+1, value_type(0.0)),
  gamma(p+1, value_type(0.0)),
  norms(p+1, value_type(0.0))
{
}

template <typename ordinal_type, typename value_type>
Stokhos::RecurrenceBasis<ordinal_type, value_type>::
RecurrenceBasis(ordinal_type p_, const RecurrenceBasis& basis) :
  name(basis.name),
  p(p_),
  normalize(basis.normalize),
  growth(basis.growth),
  quad_zero_tol(basis.quad_zero_tol),
  sparse_grid_growth_rule(basis.sparse_grid_growth_rule),
  alpha(p+1, value_type(0.0)),
  beta(p+1, value_type(0.0)),
  delta(p+1, value_type(0.0)),
  gamma(p+1, value_type(0.0)),
  norms(p+1, value_type(0.0))
{
}

template <typename ordinal_type, typename value_type>
void
Stokhos::RecurrenceBasis<ordinal_type, value_type>::
setup()
{
  bool is_normalized = 
    computeRecurrenceCoefficients(p+1, alpha, beta, delta, gamma);
  if (normalize && !is_normalized) {
    normalizeRecurrenceCoefficients(alpha, beta, delta, gamma);
  }
  

  // Compute norms -- when polynomials are normalized, this should work
  // out to one (norms[0] == 1, delta[k] == 1, beta[k] == gamma[k]
  norms[0] = beta[0]/(gamma[0]*gamma[0]);
  for (ordinal_type k=1; k<=p; k++) {
    norms[k] = (beta[k]/gamma[k])*(delta[k-1]/delta[k])*norms[k-1];
  } 
}

template <typename ordinal_type, typename value_type>
Stokhos::RecurrenceBasis<ordinal_type, value_type>::
~RecurrenceBasis()
{
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::RecurrenceBasis<ordinal_type, value_type>::
order() const
{
  return p;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::RecurrenceBasis<ordinal_type, value_type>::
size() const
{
  return p+1;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array<value_type>&
Stokhos::RecurrenceBasis<ordinal_type, value_type>::
norm_squared() const
{
  return norms;
}

template <typename ordinal_type, typename value_type>
const value_type&
Stokhos::RecurrenceBasis<ordinal_type, value_type>::
norm_squared(ordinal_type i) const
{
  return norms[i];
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Stokhos::Dense3Tensor<ordinal_type, value_type> >
Stokhos::RecurrenceBasis<ordinal_type, value_type>::
computeTripleProductTensor() const
{
  // Compute Cijk = < \Psi_i \Psi_j \Psi_k >
  ordinal_type sz = size();
  Teuchos::RCP< Stokhos::Dense3Tensor<ordinal_type, value_type> > Cijk = 
    Teuchos::rcp(new Dense3Tensor<ordinal_type, value_type>(sz));
  Teuchos::Array<value_type> points, weights;
  Teuchos::Array< Teuchos::Array<value_type> > values;
  getQuadPoints(3*p, points, weights, values);

  for (ordinal_type i=0; i<sz; i++) {
    for (ordinal_type j=0; j<sz; j++) {
      for (ordinal_type k=0; k<sz; k++) {
	value_type triple_product = 0;
	for (ordinal_type l=0; l<static_cast<ordinal_type>(points.size());
	     l++){
	  triple_product += 
	    weights[l]*(values[l][i])*(values[l][j])*(values[l][k]);
	}
	(*Cijk)(i,j,k) = triple_product;
      }
    }
  }

  return Cijk;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::RecurrenceBasis<ordinal_type, value_type>::
computeSparseTripleProductTensor(ordinal_type order) const
{
  // Compute Cijk = < \Psi_i \Psi_j \Psi_k >
  value_type sparse_tol = 1.0e-15;
  ordinal_type sz = size();
  Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk = 
    Teuchos::rcp(new Sparse3Tensor<ordinal_type, value_type>());
  Teuchos::Array<value_type> points, weights;
  Teuchos::Array< Teuchos::Array<value_type> > values;
  getQuadPoints(3*p, points, weights, values);

  for (ordinal_type i=0; i<sz; i++) {
    for (ordinal_type j=0; j<sz; j++) {
      for (ordinal_type k=0; k<order; k++) {
	value_type triple_product = 0;
	for (ordinal_type l=0; l<static_cast<ordinal_type>(points.size());
	     l++){
	  triple_product += 
	    weights[l]*(values[l][i])*(values[l][j])*(values[l][k]);
	}
	if (std::abs(triple_product/norms[i]) > sparse_tol)
	  Cijk->add_term(i,j,k,triple_product);
      }
    }
  }
  Cijk->fillComplete();

  return Cijk;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type, value_type> >
Stokhos::RecurrenceBasis<ordinal_type, value_type>::
computeDerivDoubleProductTensor() const
{
  // Compute Bij = < \Psi_i' \Psi_j >
  Teuchos::Array<value_type> points, weights;
  Teuchos::Array< Teuchos::Array<value_type> > values, derivs;
  getQuadPoints(2*p, points, weights, values);
  ordinal_type nqp = weights.size();
  derivs.resize(nqp);
  ordinal_type sz = size();
  for (ordinal_type i=0; i<nqp; i++) {
    derivs[i].resize(sz);
    evaluateBasesAndDerivatives(points[i], values[i], derivs[i]);
  }
  Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type, value_type> > Bij = 
    Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type, value_type>(sz,sz));
  for (ordinal_type i=0; i<sz; i++) {
    for (ordinal_type j=0; j<sz; j++) {
      value_type b = value_type(0.0);
      for (int qp=0; qp<nqp; qp++)
	b += weights[qp]*derivs[qp][i]*values[qp][j];
      (*Bij)(i,j) = b;
    }
  }

  return Bij;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::RecurrenceBasis<ordinal_type,value_type>::
evaluateBases(const value_type& x, Teuchos::Array<value_type>& basis_pts) const
{
  // Evaluate basis polynomials P(x) using 3 term recurrence
  // P_0(x) = 1/gamma[0], P_-1 = 0
  // P_i(x) = [(delta[i-1]*x-alpha[i-1])*P_{i-1}(x) - 
  //           beta[i-1]*P_{i-2}(x)]/gamma[i], 
  // i=2,3,...
  basis_pts[0] = value_type(1)/gamma[0];
  if (p >= 1)
    basis_pts[1] = (delta[0]*x-alpha[0])*basis_pts[0]/gamma[1];
  for (ordinal_type i=2; i<=p; i++)
    basis_pts[i] = ((delta[i-1]*x-alpha[i-1])*basis_pts[i-1] - 
		    beta[i-1]*basis_pts[i-2])/gamma[i];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::RecurrenceBasis<ordinal_type,value_type>::
evaluateBasesAndDerivatives(const value_type& x, 
			    Teuchos::Array<value_type>& vals,
			    Teuchos::Array<value_type>& derivs) const
{
  evaluateBases(x, vals);
  derivs[0] = 0.0;
  if (p >= 1)
    derivs[1] = delta[0]/(gamma[0]*gamma[1]);
  for (ordinal_type i=2; i<=p; i++)
    derivs[i] = (delta[i-1]*vals[i-1] + (delta[i-1]*x-alpha[i-1])*derivs[i-1] -
		 beta[i-1]*derivs[i-2])/gamma[i];
}

template <typename ordinal_type, typename value_type>
value_type
Stokhos::RecurrenceBasis<ordinal_type,value_type>::
evaluate(const value_type& x, ordinal_type k) const
{
  // Evaluate basis polynomials P(x) using 3 term recurrence
  // P_0(x) = 1/gamma[0], P_-1 = 0
  // P_i(x) = [(delta[i-1]*x-alpha[i-1])*P_{i-1}(x) - 
  //           beta[i-1]*P_{i-2}(x)]/gamma[i], 
  // i=2,3,...

  value_type v0 = value_type(1.0)/gamma[0];
  if (k == 0)
    return v0;

  value_type v1 = (x*delta[0]-alpha[0])*v0/gamma[1];
  if (k == 1)
    return v1;

  value_type v2 = value_type(0.0);
  for (ordinal_type i=2; i<=k; i++) {
    v2 = ((delta[i-1]*x-alpha[i-1])*v1 - beta[i-1]*v0)/gamma[i];
    v0 = v1;
    v1 = v2;
  }

  return v2;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::RecurrenceBasis<ordinal_type, value_type>::
print(std::ostream& os) const
{
  os << name << " basis of order " << p << "." << std::endl;

  os << "Alpha recurrence coefficients:\n\t";
  for (ordinal_type i=0; i<=p; i++)
    os << alpha[i] << " ";
  os << std::endl;

  os << "Beta recurrence coefficients:\n\t";
  for (ordinal_type i=0; i<=p; i++)
    os << beta[i] << " ";
  os << std::endl;

  os << "Delta recurrence coefficients:\n\t";
  for (ordinal_type i=0; i<=p; i++)
    os << delta[i] << " ";
  os << std::endl;

  os << "Gamma recurrence coefficients:\n\t";
  for (ordinal_type i=0; i<=p; i++)
    os << gamma[i] << " ";
  os << std::endl;
       
  os << "Basis polynomial norms (squared):\n\t";
  for (ordinal_type i=0; i<=p; i++)
    os << norms[i] << " ";
  os << std::endl;
}

template <typename ordinal_type, typename value_type>
const std::string&
Stokhos::RecurrenceBasis<ordinal_type, value_type>::
getName() const
{
  return name;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::RecurrenceBasis<ordinal_type,value_type>::
getQuadPoints(ordinal_type quad_order,
	      Teuchos::Array<value_type>& quad_points,
	      Teuchos::Array<value_type>& quad_weights,
	      Teuchos::Array< Teuchos::Array<value_type> >& quad_values) const
{
  
  //This is a transposition into C++ of Gautschi's code for taking the first 
  // N recurrance coefficient and generating a N point quadrature rule.  
  // The MATLAB version is available at
  // http://www.cs.purdue.edu/archives/2002/wxg/codes/gauss.m
  ordinal_type num_points = 
    static_cast<ordinal_type>(std::ceil((quad_order+1)/2.0));
  Teuchos::Array<value_type> a(num_points,0);
  Teuchos::Array<value_type> b(num_points,0);
  Teuchos::Array<value_type> c(num_points,0);
  Teuchos::Array<value_type> d(num_points,0);
  
  // If we don't have enough recurrance coefficients, get some more.
  if(num_points > p+1){
    bool is_normalized = computeRecurrenceCoefficients(num_points, a, b, c, d);
    if (!is_normalized)
      normalizeRecurrenceCoefficients(a, b, c, d);
  }
  else {  //else just take the ones we already have.
    for(ordinal_type n = 0; n<num_points; n++){
      a[n] = alpha[n];
      b[n] = beta[n];
      c[n] = delta[n];
      d[n] = gamma[n];
    }
    if (!normalize)
      normalizeRecurrenceCoefficients(a, b, c, d);
  }
  
  // With normalized coefficients, A is symmetric and tri-diagonal, with
  // diagonal = a, and off-diagonal = b, starting at b[1]
  Teuchos::SerialDenseMatrix<ordinal_type,value_type> eig_vectors(num_points,
								  num_points);
  Teuchos::Array<value_type> workspace(2*num_points);
  ordinal_type info_flag;
  Teuchos::LAPACK<ordinal_type,value_type> my_lapack;
  
  // compute the eigenvalues (stored in a) and right eigenvectors.
  if (num_points == 1)
    my_lapack.STEQR('I', num_points, &a[0], &b[0], eig_vectors.values(), 
		    num_points, &workspace[0], &info_flag);
  else
    my_lapack.STEQR('I', num_points, &a[0], &b[1], eig_vectors.values(), 
		    num_points, &workspace[0], &info_flag);
  
  // eigenvalues are sorted by STEQR
  quad_points.resize(num_points);
  quad_weights.resize(num_points);
  for (ordinal_type i=0; i<num_points; i++) {
    quad_points[i] = a[i];
    if (std::abs(quad_points[i]) < quad_zero_tol)
      quad_points[i] = 0.0;
    quad_weights[i] = beta[0]*eig_vectors[i][0]*eig_vectors[i][0];
  }
  
  // Evalute basis at gauss points
  quad_values.resize(num_points);
  for (ordinal_type i=0; i<num_points; i++) {
    quad_values[i].resize(p+1);
    evaluateBases(quad_points[i], quad_values[i]);
  }
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::RecurrenceBasis<ordinal_type,value_type>::
quadDegreeOfExactness(ordinal_type n) const
{
  return ordinal_type(2)*n-ordinal_type(1);
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::RecurrenceBasis<ordinal_type,value_type>::
coefficientGrowth(ordinal_type n) const
{
  if (growth == SLOW_GROWTH)
    return n;

  // else moderate
  return ordinal_type(2)*n;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::RecurrenceBasis<ordinal_type,value_type>::
pointGrowth(ordinal_type n) const
{
  if (growth == SLOW_GROWTH) {
    if (n % ordinal_type(2) == ordinal_type(1))
      return n + ordinal_type(1);
    else
      return n;
  }

  // else moderate
  return n;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::RecurrenceBasis<ordinal_type,value_type>::
getRecurrenceCoefficients(Teuchos::Array<value_type>& a,
			  Teuchos::Array<value_type>& b,
			  Teuchos::Array<value_type>& c,
			  Teuchos::Array<value_type>& g) const
{
  a = alpha;
  b = beta;
  c = delta;
  g = gamma;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::RecurrenceBasis<ordinal_type,value_type>::
normalizeRecurrenceCoefficients(
  Teuchos::Array<value_type>& a,
  Teuchos::Array<value_type>& b,
  Teuchos::Array<value_type>& c,
  Teuchos::Array<value_type>& g) const
{
  ordinal_type n = a.size();
  a[0] = a[0] / c[0];
  b[0] = std::sqrt(b[0]);
  g[0] = b[0];
  for (ordinal_type k=1; k<n; k++) {
    a[k] = a[k] / c[k];
    b[k] = std::sqrt((b[k]*g[k])/(c[k]*c[k-1]));
    g[k] = b[k];
  }
  for (ordinal_type k=0; k<n; k++)
    c[k] = value_type(1);
}
