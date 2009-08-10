// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2008) Sandia Corporation
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

template <typename ordinal_type, typename value_type>
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
StieltjesPCEBasis(
   ordinal_type p_,
   const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& pce,
   const Stokhos::OrthogPolyBasis<ordinal_type, value_type>& pce_basis,
   const Stokhos::Quadrature<ordinal_type, value_type>& quad,
   bool use_pce_quad_points_) :
  name("Stieltjes PCE"),
  p(p_),
  norms(p+1),
  alpha(p+1),
  beta(p+1),
  pce_vals(),
  phi_vals(),
  Cijk(),
  use_pce_quad_points(use_pce_quad_points_)
{
  // Evaluate PCE at quad points
  const std::vector< std::vector<value_type> >& quad_points =
    quad.getQuadPoints();
  pce_weights = quad.getQuadWeights();
  const std::vector< std::vector<value_type> >& basis_values =
    quad.getBasisAtQuadPoints();
  ordinal_type nqp = pce_weights.size();
  pce_vals.resize(nqp);
  phi_vals.resize(nqp);
  for (ordinal_type i=0; i<nqp; i++) {
    pce_vals[i] = pce.evaluate(pce_basis, quad_points[i], basis_values[i]);
    phi_vals[i].resize(p+1);
  }
  
  // Compute coefficients via Stieltjes
  stieltjes(0, p+1, pce_weights, pce_vals, alpha, beta, norms, phi_vals);

  /*
  std::cout << "pce_weights = ";
  for (ordinal_type i=0; i<nqp; i++)
    std::cout << pce_weights[i] << " ";
  std::cout << std::endl;
 
  std::cout << "pce_vals = ";
  for (ordinal_type i=0; i<nqp; i++)
    std::cout << pce_vals[i] << " ";
  std::cout << std::endl;

  for (ordinal_type j=0; j<nqp; j++) {
    std::cout << "phi_vals[" << j << "] = ";
    for (ordinal_type i=0; i<p+1; i++)
      std::cout << phi_vals[j][i] << " ";
    std::cout << std::endl;
  }

  std::cout << "alpha = ";
  for (ordinal_type i=0; i<=p; i++)
    std::cout << alpha[i] << " ";
  std::cout << std::endl;

  std::cout << "beta = ";
  for (ordinal_type i=0; i<=p; i++)
    std::cout << beta[i] << " ";
  std::cout << std::endl;

  std::cout << "norms = ";
  for (ordinal_type i=0; i<=p; i++)
    std::cout << norms[i] << " ";
  std::cout << std::endl;
  */
}

template <typename ordinal_type, typename value_type>
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
~StieltjesPCEBasis()
{
}

template <typename ordinal_type, typename value_type>
ordinal_type 
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
order() const
{
  return p;
}

template <typename ordinal_type, typename value_type>
ordinal_type 
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
dimension() const
{
  return 1;
}

template <typename ordinal_type, typename value_type>
ordinal_type 
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
size() const
{
  return p+1;
}

template <typename ordinal_type, typename value_type>
const std::vector<value_type>&
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
norm_squared() const
{
  return norms;
}

template <typename ordinal_type, typename value_type>
const value_type&
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
norm_squared(ordinal_type i) const
{
  return norms[i];
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< const Stokhos::Dense3Tensor<ordinal_type, value_type> >
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
getTripleProductTensor() const
{
  ordinal_type sz = size();
  
  // Compute Cijk = < \Psi_i \Psi_j \Psi_k >
  if (Cijk == Teuchos::null) {
    std::vector<value_type> points, weights;
    std::vector< std::vector<value_type> > values;
    getQuadPoints(3*p, points, weights, values);
    Cijk = Teuchos::rcp(new Dense3Tensor<ordinal_type, value_type>(sz));
    
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
  }

  return Cijk;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< const Teuchos::SerialDenseMatrix<ordinal_type, value_type> >
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
getDerivDoubleProductTensor() const
{
  // TEST_FOR_EXCEPTION(true, std::logic_error,
  // 		     "Stokhos::StieltjesPCEBasis::getDerivDoubleProductTensor():  "
  // 		     << " Method not implemented!");
  return Teuchos::null;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
projectPoly(const Stokhos::Polynomial<value_type>& x, 
	    std::vector<value_type>& coeffs) const
{
  TEST_FOR_EXCEPTION(true, std::logic_error,
		     "Stokhos::StieltjesPCEBasis::projectPoly():  "
		     << " Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
projectProduct(ordinal_type i, ordinal_type j,
	       std::vector<value_type>& coeffs) const
{
  TEST_FOR_EXCEPTION(true, std::logic_error,
		     "Stokhos::StieltjesPCEBasis::projectProduct():  "
		     << " Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
projectDerivative(ordinal_type i, std::vector<value_type>& coeffs) const
{
  TEST_FOR_EXCEPTION(true, std::logic_error,
		     "Stokhos::StieltjesPCEBasis::projectDerivative():  "
		     << " Method not implemented!");
}

template <typename ordinal_type, typename value_type>
Stokhos::Polynomial<value_type>
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
toStandardBasis(const value_type coeffs[], ordinal_type n) const
{
  TEST_FOR_EXCEPTION(true, std::logic_error,
		     "Stokhos::StieltjesPCEBasis::toStandardBasis():  "
		     << " Method not implemented!");
}

template <typename ordinal_type, typename value_type>
value_type
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
evaluateZero(ordinal_type i) const
{
  TEST_FOR_EXCEPTION(true, std::logic_error,
		     "Stokhos::StieltjesPCEBasis::evaluateZero():  "
		     << " Method not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
evaluateBases(const value_type& x, std::vector<value_type>& basis_pts) const
{
  // Evaluate basis polynomials P(x) using 3 term recurrence
  basis_pts[0] = 1.0;
  if (p >= 1)
    basis_pts[1] = x - alpha[0];
  for (ordinal_type k=2; k<=p; k++)
    basis_pts[k] = (x - alpha[k-1])*basis_pts[k-1] - beta[k-1]*basis_pts[k-2];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
print(std::ostream& os) const
{
  os << name << " basis of order " << p << ".\n";
  os << "Basis vector norms (squared):\n\t";
  for (ordinal_type i=0; i<static_cast<ordinal_type>(norms.size()); i++)
    os << norms[i] << " ";
  os << "\n";
}

template <typename ordinal_type, typename value_type>
const std::string&
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
getName() const
{
  return name;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
getQuadPoints(ordinal_type quad_order,
	      std::vector<value_type>& quad_points,
	      std::vector<value_type>& quad_weights,
	      std::vector< std::vector<value_type> >& quad_values) const
{
  // Use underlying pce's quad points, weights, values
  if (use_pce_quad_points) {
    quad_points = pce_vals;
    quad_weights = pce_weights;
    quad_values = phi_vals;
    return;
  }

  // Compute gauss points, weights
  ordinal_type num_points = 
    static_cast<ordinal_type>(std::ceil((quad_order+1)/2.0));
  std::vector<double> x(num_points), w(num_points);
  
  //This is a transposition into C++ of Gautschi's code for taking the first N recurrance coefficients
  //and generating a N point quadrature rule.  The MATLAB version is available at
  // http://www.cs.purdue.edu/archives/2002/wxg/codes/gauss.m
  
  //If we don't have enough recurrance coefficients, get some more.
  std::vector<value_type> a(num_points), b(num_points);
  if (num_points > p+1) {
    ordinal_type nqp = phi_vals.size();
    std::vector<value_type> nrm(num_points);
    std::vector< std::vector<value_type> > vals(nqp);
    for (ordinal_type i=0; i<nqp; i++)
      vals[i].resize(num_points);
    stieltjes(0, num_points, pce_weights, pce_vals, a, b, nrm, vals);
  }
  else 
    for (ordinal_type i=0; i<num_points; i++) {
      a[i] = alpha[i];
      b[i] = beta[i];
    }
    
  
  // A is symmetric and tridiagonal with A[i][i] = alpha[i] 
  // and A[i-1][i] = A[i][i-1] = sqrt(beta[i]) 
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> A(num_points,num_points);
  for (ordinal_type i=0; i<num_points; i++) 
    A(i,i) = a[i];
  for (ordinal_type i=1; i<num_points; i++) { 
    value_type t = std::sqrt(b[i]);
    A(i,i-1) = t;
    A(i-1,i) = t;
  }

  Teuchos::SerialDenseMatrix<ordinal_type,value_type> eig_vectors(num_points,
								  num_points);
  
  Teuchos::SerialDenseVector<ordinal_type,value_type> workspace(8*num_points);
  ordinal_type info_flag;
  Teuchos::LAPACK<ordinal_type,value_type> my_lapack;
  Teuchos::SerialDenseVector<ordinal_type,value_type> eig_real(num_points);
  Teuchos::SerialDenseVector<ordinal_type,value_type> eig_imag(num_points);
  
  //compute the eigenvalues and right eigenvectors.
  my_lapack.GEEV('N', 'V', num_points, A.values(), num_points,
		 eig_real.values(), eig_imag.values(), A.values(), num_points,
		 eig_vectors.values(), num_points, workspace.values(), 
		 8*num_points, &info_flag);
  
  quad_points.resize(num_points);
  quad_weights.resize(num_points);
  
  //the eigenvalues are the quadrature points, sort these and keep track of the indices.
  std::vector< ordinal_type > idx(num_points,0);
  value_type temp1,temp2;
  for(ordinal_type i = 0; i< num_points; i++) 
    idx[i] = i;
  for (ordinal_type i = 0; i< num_points; i++) {
    for (ordinal_type j = 0; j< num_points; j++) {
      if (eig_real[i]< eig_real[j]) {
        temp1 = eig_real[j];
        temp2 = idx[j];
        eig_real[j] = eig_real[i];
        idx[j] = idx[i];
        eig_real[i] = temp1;
        idx[i] = temp2;
      }
    }
  }
  
  
  for (ordinal_type i = 0; i< num_points; i++) {
    quad_points[i] = eig_real[i];
    quad_weights[i] = b[0]*pow(eig_vectors[idx[i]][0],2);
  }

  // Evalute basis at gauss points
  quad_values.resize(num_points);
  for (ordinal_type i=0; i<num_points; i++) {
    quad_values[i].resize(num_points);
    evaluateBases(quad_points[i], quad_values[i]);
  }
  
//   std::cout << "StieltjesPCE quadrature points, weights, values = " << std::endl;
//   for (int i=0; i<n; i++) {
//     std::cout << "\t" << quad_points[i] 
//               << "\t" << quad_weights[i];
//     for (ordinal_type j=0; j<p+1; j++)
//       std::cout << "\t" << quad_values[i][j];
//     cout << std::endl;
//   }
}

template <typename ordinal_type, typename value_type>
ordinal_type 
Stokhos::StieltjesPCEBasis<ordinal_type, value_type>::
getRule() const
{
  return 0;
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
stieltjes(ordinal_type nstart,
	  ordinal_type nfinish,
	  const std::vector<value_type>& weights,
	  const std::vector<value_type>& points,
	  std::vector<value_type>& a,
	  std::vector<value_type>& b,
	  std::vector<value_type>& nrm,
	  std::vector< std::vector<value_type> >& phi_vals) const
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
		      const std::vector<value_type>& a,
		      const std::vector<value_type>& b,
		      const std::vector<value_type>& weights,
		      const std::vector<value_type>& points,
		      std::vector< std::vector<value_type> >& phi_vals,
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
		   const std::vector<value_type>& a,
		   const std::vector<value_type>& b,
		   const std::vector<value_type>& points,
		   std::vector< std::vector<value_type> >& values) const
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
