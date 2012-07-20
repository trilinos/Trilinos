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

#include "Teuchos_SerialDenseHelpers.hpp"

#define DGEQP3_F77  F77_BLAS_MANGLE(dgeqp3,DGEQP3)
extern "C" {
void DGEQP3_F77(int*, int*, double*, int*, int*, double*, double*, int*, int*);
}

#ifdef HAVE_STOKHOS_CLP
#include "coin/ClpSimplex.hpp"
#endif

#ifdef HAVE_STOKHOS_GLPK
extern "C" {
#include "glpk.h"
}
#endif

template <typename ordinal_type, typename value_type>
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
MonomialGramSchmidtSimplexPCEBasis(
  ordinal_type p_,
  const Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pce,
  const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad_,
  const Teuchos::ParameterList& params_) :
  name("Monomial Gram Schmidt Simplex PCE Basis"),
  quad(quad_),
  params(params_),
  pce_sz(pce[0].size()),
  p(p_),
  d(pce.size()),
  verbose(params.get("Verbose", false)),
  reduction_tol(params.get("Reduction Tolerance", 1.0e-12)),
  orthogonalization_method(params.get("Orthogonalization Method", 
				      "Classical Gram-Schmidt"))
{
  
  
  // Compute basis terms -- 2-D array giving powers for each linear index
  compute_terms(p, d, sz, terms, num_terms);

  // Basis is orthonormal by construction
  norms.resize(sz, 1.0);

  // Get quadrature data
  const Teuchos::Array<value_type>& weights = quad->getQuadWeights();
  const Teuchos::Array< Teuchos::Array<value_type> >& points = 
    quad->getQuadPoints(); 
  const Teuchos::Array< Teuchos::Array<value_type> >& basis_values = 
    quad->getBasisAtQuadPoints();
  ordinal_type nqp = weights.size();

  // Original basis at quadrature points -- needed to transform expansions
  // in this basis back to original
  A.reshape(nqp, pce_sz);
  for (ordinal_type i=0; i<nqp; i++)
    for (ordinal_type j=0; j<pce_sz; j++)
      A(i,j) = basis_values[i][j];

  // Compute F matrix -- PCEs evaluated at all quadrature points
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> F(nqp, d);
  Teuchos::Array< Teuchos::Array<value_type> > values(nqp);
  for (ordinal_type i=0; i<nqp; i++)
    for (ordinal_type j=0; j<d; j++)
      F(i,j) = pce[j].evaluate(points[i], basis_values[i]);

  // Compute B matrix -- monomials in F
  // for i=0,...,nqp-1
  //   for j=0,...,sz2-1
  //      B(i,j) = F(i,1)^terms[j][1] * ... * F(i,d)^terms[j][d]
  // where sz2 is the total size of a basis up to order 2*p and terms[j] 
  // is an array of powers for each term in the total-order basis
  Teuchos::Array< Teuchos::Array<ordinal_type> > terms2;
  Teuchos::Array<ordinal_type> num_terms2;
  ordinal_type sz2;
  compute_terms(2*p, d, sz2, terms2, num_terms2);
  if (verbose)
    std::cout << "sz2 = " << sz2 << std::endl;
  //TEUCHOS_ASSERT(sz2 <= nqp);
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> B2(nqp, sz2);
  for (ordinal_type i=0; i<nqp; i++) {
    for (ordinal_type j=0; j<sz2; j++) {
      B2(i,j) = 1.0;
      for (ordinal_type k=0; k<d; k++)
	B2(i,j) *= std::pow(F(i,k), terms2[j][k]);
    }
  }

  // Get the first sz columns of B to define our basis
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> B(
    Teuchos::View, B2, nqp, sz, 0, 0);

  // Compute QR factorization of B using Gram-Schmidt
  // Q defines our new basis
  Q.shape(nqp, sz);
  R.shape(sz, sz);
  if (orthogonalization_method == "Classical Gram-Schmidt")
    computeQR_CGS(sz, B, weights, Q, R);
  else if (orthogonalization_method == "Modified Gram-Schmidt")
    computeQR_MGS(sz, B, weights, Q, R);
  else
    TEUCHOS_TEST_FOR_EXCEPTION(
	true, std::logic_error, 
	"Invalid orthogonalization method " << orthogonalization_method);

  Teuchos::RCP< Teuchos::Array<value_type> > reduced_weights;
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > > reduced_points;
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > > reduced_values;

  // Compute reduced quadrature rule
  std::string reduction_method = 
    params.get("Reduced Quadrature Method", "L1 Minimization");
  if (reduction_method == "Column-Pivoted QR")
    reducedQuadrature_QRCP(B2, Q, F, weights, 
			   reduced_weights, reduced_points, reduced_values);
  else if (reduction_method == "L1 Minimization") {
    std::string solver = params.get("LP Solver", "GLPK");
    if (solver == "GLPK")
       reducedQuadrature_GLPK(B2, Q, F, weights, 
			      reduced_weights, reduced_points, reduced_values);
    else
      TEUCHOS_TEST_FOR_EXCEPTION(
	true, std::logic_error, "Invalid LP solver method " << solver);
  }
  else if (reduction_method == "None") {
    reduced_weights =
      Teuchos::rcp(new Teuchos::Array<value_type>(weights));
    reduced_points =
      Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(points));
    reduced_values =
      Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(nqp));
    for (ordinal_type i=0; i<nqp; i++) {
      (*reduced_values)[i].resize(sz);
      for (ordinal_type j=0; j<sz; j++)
	(*reduced_values)[i][j] = Q(i,j);
    }
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(
	true, std::logic_error, 
	"Invalid dimension reduction method " << reduction_method);
 
  
  // Build reduced quadrature object
  Teuchos::RCP< const Teuchos::Array<value_type> > creduced_weights =
    reduced_weights;
  Teuchos::RCP< const Teuchos::Array< Teuchos::Array<value_type> > > creduced_points = reduced_points;
  Teuchos::RCP< const Teuchos::Array< Teuchos::Array<value_type> > > creduced_values = reduced_values;
  reduced_quad =
    Teuchos::rcp(new UserDefinedQuadrature<ordinal_type,value_type>(
		   creduced_points,
		   creduced_weights,
		   creduced_values));
  
}

template <typename ordinal_type, typename value_type>
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
~MonomialGramSchmidtSimplexPCEBasis()
{
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
order() const
{
  return p;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
dimension() const
{
  return d;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
size() const
{
  return sz;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array<value_type>&
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
norm_squared() const
{
  return norms;
}

template <typename ordinal_type, typename value_type>
const value_type&
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
norm_squared(ordinal_type i) const
{
  return norms[i];
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
computeTripleProductTensor(ordinal_type order) const

{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented!");
}

template <typename ordinal_type, typename value_type>
value_type
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
evaluateZero(ordinal_type i) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
evaluateBases(const Teuchos::Array<value_type>& point,
	      Teuchos::Array<value_type>& basis_vals) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
print(std::ostream& os) const
{
  os << "Gram-Schmidt basis of order " << p << ", dimension " << d 
     << ", and size " << sz << ".  Matrix coefficients:\n";
  os << Q << std::endl;
  os << "Basis vector norms (squared):\n\t";
  for (ordinal_type i=0; i<sz; i++)
    os << norms[i] << " ";
  os << "\n";
}

template <typename ordinal_type, typename value_type>
const std::string&
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
getName() const
{
  return name;
}

template <typename ordinal_type, typename value_type>
Teuchos::Array<ordinal_type>
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
getTerm(ordinal_type i) const
{
  return terms[i];
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
getIndex(const Teuchos::Array<ordinal_type>& term) const
{
  return compute_index(term);
}

template <typename ordinal_type, typename value_type>
Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type, value_type> > >
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
getCoordinateBases() const
{
  return Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type, value_type> > >();
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
getReducedQuadrature() const
{
  return reduced_quad;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
transformToOriginalBasis(const value_type *in, value_type *out) const
{
  Teuchos::SerialDenseVector<ordinal_type, value_type> zbar(
    Teuchos::View, const_cast<value_type*>(in), sz);
  Teuchos::SerialDenseVector<ordinal_type, value_type> z(
    Teuchos::View, out, pce_sz);
  const Teuchos::Array<value_type>& weights = quad->getQuadWeights();
  ordinal_type nqp = weights.size();
  Teuchos::SerialDenseVector<ordinal_type, value_type> tmp(nqp);

  // Compute z = A^T*W*Q*zbar
  tmp.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Q, zbar, 0.0);
  for (ordinal_type i=0; i<nqp; i++)
    tmp[i] *= weights[i];

  z.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, A, tmp, 0.0);
  
}

template <typename ordinal_type, typename value_type>
void
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
computeTransformedPCE(
  ordinal_type i,
  Stokhos::OrthogPolyApprox<ordinal_type,value_type>& pce) const
{
  // a <- R*a
  Teuchos::SerialDenseVector<ordinal_type, value_type> a(
    Teuchos::View, pce.coeff(), sz);
  a[i+1] = 1.0;
  blas.TRMM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
	    Teuchos::NON_UNIT_DIAG, sz, 1, 1.0, R.values(), sz, a.values(),
	    sz);
}

template <typename ordinal_type, typename value_type>
void 
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
getBasisAtOriginalQuadraturePoints(
  Teuchos::Array< Teuchos::Array<double> >& red_basis_vals) const
{
  ordinal_type nqp = Q.numRows();
  red_basis_vals.resize(nqp);
  for (ordinal_type i=0; i<nqp; i++) {
    red_basis_vals[i].resize(sz);
    for (ordinal_type j=0; j<sz; j++)
      red_basis_vals[i][j] = Q(i,j);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
compute_terms(ordinal_type p, ordinal_type d, ordinal_type& sz,
	      Teuchos::Array< Teuchos::Array<ordinal_type> >& terms,
	      Teuchos::Array<ordinal_type>& num_terms)
{
  // The approach here for ordering the terms is inductive on the total
  // order p.  We get the terms of total order p from the terms of total
  // order p-1 by incrementing the orders of the first dimension by 1.
  // We then increment the orders of the second dimension by 1 for all of the
  // terms whose first dimension order is 0.  We then repeat for the third
  // dimension whose first and second dimension orders are 0, and so on.
  // How this is done is most easily illustrated by an example of dimension 3:
  //
  // Order  terms   cnt  Order  terms   cnt
  //   0    0 0 0          4    4 0 0  15 5 1
  //                            3 1 0
  //   1    1 0 0  3 2 1        3 0 1
  //        0 1 0               2 2 0
  //        0 0 1               2 1 1
  //                            2 0 2
  //   2    2 0 0  6 3 1        1 3 0
  //        1 1 0               1 2 1
  //        1 0 1               1 1 2
  //        0 2 0               1 0 3
  //        0 1 1               0 4 0
  //        0 0 2               0 3 1
  //                            0 2 2
  //   3    3 0 0  10 4 1       0 1 3
  //        2 1 0               0 0 4
  //        2 0 1
  //        1 2 0
  //        1 1 1
  //        1 0 2
  //        0 3 0
  //        0 2 1
  //        0 1 2
  //        0 0 3

  // Temporary array of terms grouped in terms of same order
  Teuchos::Array< Teuchos::Array< Teuchos::Array<ordinal_type> > > terms_order(p+1);

  // Store number of terms up to each order
  num_terms.resize(p+2, ordinal_type(0));

  // Set order 0
  terms_order[0].resize(1);
  terms_order[0][0].resize(d, ordinal_type(0));
  num_terms[0] = 1;

  Teuchos::Array<ordinal_type> basis_orders(d, p);

  // The array "cnt" stores the number of terms we need to increment for each
  // dimension.  
  Teuchos::Array<ordinal_type> cnt(d), cnt_next(d), term(d);
  for (ordinal_type j=0; j<d; j++) {
    if (basis_orders[j] >= 1)
      cnt[j] = 1;
    else
      cnt[j] = 0;
    cnt_next[j] = 0;
  }

  sz = 1;
  // Loop over orders
  for (ordinal_type k=1; k<=p; k++) {

    num_terms[k] = num_terms[k-1];

    // Stores the index of the term we copying
    ordinal_type prev = 0;

    // Loop over dimensions
    for (ordinal_type j=0; j<d; j++) {

      // Increment orders of cnt[j] terms for dimension j
      for (ordinal_type i=0; i<cnt[j]; i++) {
	if (terms_order[k-1][prev+i][j] < basis_orders[j]) {
	  term = terms_order[k-1][prev+i];
	  ++term[j];
	  terms_order[k].push_back(term);
	  ++sz;
	  num_terms[k]++;
	  for (ordinal_type l=0; l<=j; l++)
	    ++cnt_next[l];
	}
      }

      // Move forward to where all orders for dimension j are 0
      if (j < d-1)
	prev += cnt[j] - cnt[j+1];

    }

    // Update the number of terms we must increment for the new order
    for (ordinal_type j=0; j<d; j++) {
      cnt[j] = cnt_next[j];
      cnt_next[j] = 0;
    }

  }

  num_terms[p+1] = sz;

  // Copy into final terms array
  terms.resize(sz);
  ordinal_type i = 0;
  for (ordinal_type k=0; k<=p; k++) {
    ordinal_type num_k = terms_order[k].size();
    for (ordinal_type j=0; j<num_k; j++)
      terms[i++] = terms_order[k][j];
  }

  /*
  std::cout << "sz = " << sz << std::endl;
  for (ordinal_type i=0; i<sz; i++) {
    std::cout << i << ":  ";
    for (ordinal_type j=0; j<d; j++)
      std::cout << terms[i][j] << " ";
    std::cout << std::endl;
  }
  std::cout << "num_terms = ";
  for (ordinal_type i=0; i<=p; i++)
    std::cout << num_terms[i] << " ";
  std::cout << std::endl;
  */
}

template <typename ordinal_type, typename value_type>
ordinal_type 
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
compute_index(const Teuchos::Array<ordinal_type>& term) const
{
  // The approach here for computing the index is to find the order block
  // corresponding to this term by adding up the component orders.  We then
  // do a linear search through the terms_order array for this order

  // First compute order of term
  ordinal_type ord = 0;
  for (ordinal_type i=0; i<d; i++)
    ord += term[i];
  TEUCHOS_TEST_FOR_EXCEPTION(ord < 0 || ord > p, std::logic_error,
		     "Stokhos::CompletePolynomialBasis::compute_index(): " <<
		     "Term has invalid order " << ord);

  // Now search through terms of that order to find a match
  ordinal_type k;
  if (ord == 0)
    k = 0;
  else
    k = num_terms[ord-1];
  ordinal_type k_max=num_terms[ord];
  bool found = false;
  while (k < k_max && !found) {
    bool found_term = true;
    for (ordinal_type j=0; j<d; j++) {
      found_term = found_term && (term[j] == terms[k][j]);
      if (!found_term)
	break;
    }
    found = found_term;
    ++k;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(k >= k_max && !found, std::logic_error,
		     "Stokhos::CompletePolynomialBasis::compute_index(): " <<
		     "Could not find specified term.");

  return k-1;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
computeQR_CGS(
  ordinal_type k,
  Teuchos::SerialDenseMatrix<ordinal_type, value_type>& A,
  const Teuchos::Array<value_type>& w,
  Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
  Teuchos::SerialDenseMatrix<ordinal_type, value_type>& R)
{
  // Compute thin QR factorization using classical Gram-Schmidt
  using Teuchos::getCol;
  typedef Teuchos::SerialDenseVector<ordinal_type,value_type> SDV;
  
  for (ordinal_type j=0; j<k; j++) {
    SDV Aj = getCol(Teuchos::View, A, j);
    SDV Qj = getCol(Teuchos::View, Q, j);
    Qj.assign(Aj);
    for (ordinal_type i=0; i<j; i++) {
      SDV Qi = getCol(Teuchos::View, Q, i);
      R(i,j) = weighted_inner_product(Qi, Aj, w);
      saxpy(1.0, Qj, -R(i,j), Qi);  // Q(:,j) = 1.0*Q(:,j) - R(i,j)*Q(:,i)
    }
    R(j,j) = std::sqrt(weighted_inner_product(Qj, Qj, w));
    Qj.scale(1.0/R(j,j));
  }

  value_type r_max = std::abs(R(0,0));
  value_type r_min = std::abs(R(0,0));
  for (ordinal_type i=1; i<k; i++) {
    if (std::abs(R(i,i)) > r_max)
      r_max = R(i,i);
    if (std::abs(R(i,i)) < r_min)
      r_min = R(i,i);
  }
  value_type cond_r = r_max / r_min;
  std::cout << "Condition number of R = " << cond_r << std::endl;

  // Check Q^T*W*Q = I
  if (verbose) {
    ordinal_type m = Q.numRows();
    ordinal_type n = Q.numCols();
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> Qt(m,n);
    for (ordinal_type i=0; i<m; i++)
      for (ordinal_type j=0; j<n; j++)
	Qt(i,j) = w[i]*Q(i,j);
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> err1(n,n);
    err1.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Q, Qt, 0.0);
    for (ordinal_type i=0; i<n; i++)
      err1(i,i) -= 1.0;
    std::cout << "||Q^T*W*Q - I||_infty = " << err1.normInf() << std::endl;
  }
  
  // Check A = QR
  if (verbose) {
    ordinal_type m = Q.numRows();
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> AA(
      Teuchos::View, A, m, k, 0, 0);
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> err2(m,k);
    err2.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Q, R, 0.0);
    err2 -= AA;
    std::cout << "||QR-A||_infty = " << err2.normInf() << std::endl;
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
computeQR_MGS(
  ordinal_type k,
  Teuchos::SerialDenseMatrix<ordinal_type, value_type>& A,
  const Teuchos::Array<value_type>& w,
  Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
  Teuchos::SerialDenseMatrix<ordinal_type, value_type>& R)
{
  // Compute thin QR factorization using classical Gram-Schmidt
  using Teuchos::getCol;
  typedef Teuchos::SerialDenseVector<ordinal_type,value_type> SDV;

  /*
  for (ordinal_type i=0; i<k; i++) {
    SDV Ai = getCol(Teuchos::View, A, i);
    SDV Qi = getCol(Teuchos::View, Q, i);
    Qi.assign(Ai);
  }
  for (ordinal_type i=0; i<k; i++) {
    SDV Qi = getCol(Teuchos::View, Q, i);
    R(i,i) = std::sqrt(weighted_inner_product(Qi, Qi, w));
    Qi.scale(1.0/R(i,i));
    for (ordinal_type j=i+1; j<k; j++) {
      SDV Aj = getCol(Teuchos::View, A, j);
      SDV Qj = getCol(Teuchos::View, Q, j);
      R(i,j) = weighted_inner_product(Qi, Aj, w);
      saxpy(1.0, Qj, -R(i,j), Qi);  // Q(:,j) = 1.0*Q(:,j) - R(i,j)*Q(:,i)
    }
  }
  */
  
  for (ordinal_type j=0; j<k; j++) {
    SDV Aj = getCol(Teuchos::View, A, j);
    SDV Qj = getCol(Teuchos::View, Q, j);
    Qj.assign(Aj);
    for (ordinal_type i=0; i<j; i++) {
      SDV Qi = getCol(Teuchos::View, Q, i);
      R(i,j) = weighted_inner_product(Qi, Qj, w);
      saxpy(1.0, Qj, -R(i,j), Qi);  // Q(:,j) = 1.0*Q(:,j) - R(i,j)*Q(:,i)
    }
    R(j,j) = std::sqrt(weighted_inner_product(Qj, Qj, w));
    Qj.scale(1.0/R(j,j));
  }

  value_type r_max = std::abs(R(0,0));
  value_type r_min = std::abs(R(0,0));
  for (ordinal_type i=1; i<k; i++) {
    if (std::abs(R(i,i)) > r_max)
      r_max = R(i,i);
    if (std::abs(R(i,i)) < r_min)
      r_min = R(i,i);
  }
  value_type cond_r = r_max / r_min;
  std::cout << "Condition number of R = " << cond_r << std::endl;

  // Check Q^T*W*Q = I
  if (verbose) {
    ordinal_type m = Q.numRows();
    ordinal_type n = Q.numCols();
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> Qt(m,n);
    for (ordinal_type i=0; i<m; i++)
      for (ordinal_type j=0; j<n; j++)
	Qt(i,j) = w[i]*Q(i,j);
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> err1(n,n);
    err1.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Q, Qt, 0.0);
    for (ordinal_type i=0; i<n; i++)
      err1(i,i) -= 1.0;
    std::cout << "||Q^T*W*Q - I||_infty = " << err1.normInf() << std::endl;
  }
  
  // Check A = QR
  if (verbose) {
    ordinal_type m = Q.numRows();
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> AA(
      Teuchos::View, A, m, k, 0, 0);
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> err2(m,k);
    err2.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Q, R, 0.0);
    err2 -= AA;
    std::cout << "||QR-A||_infty = " << err2.normInf() << std::endl;
  }
}
      
template <typename ordinal_type, typename value_type>
value_type
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
weighted_inner_product(
  const Teuchos::SerialDenseVector<ordinal_type,value_type>& x,
  const Teuchos::SerialDenseVector<ordinal_type,value_type>& y,
  const Teuchos::Array<value_type>& w)
{
  ordinal_type n = x.length();
  value_type t = 0;
  for (ordinal_type i=0; i<n; i++)
    t += x[i]*w[i]*y[i];

  return t;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
saxpy(
  const value_type& alpha,
  Teuchos::SerialDenseVector<ordinal_type, value_type>& x,
  const value_type& beta,
  Teuchos::SerialDenseVector<ordinal_type, value_type>& y)
{
  ordinal_type n = x.length();
  for (ordinal_type i=0; i<n; i++)
    x[i] = alpha*x[i] + beta*y[i];
}

template <typename ordinal_type, typename value_type>
void
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
reducedQuadrature_QRCP(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& B2,
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
  const Teuchos::Array<value_type>& weights,
  Teuchos::RCP< Teuchos::Array<value_type> >& reduced_weights,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_points,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_values)
{
  //
  // Find reduced quadrature weights by applying column-pivoted QR to
  // problem B2^T*w = e_1
  //
  ordinal_type nqp = B2.numRows();
  ordinal_type sz2 = B2.numCols();
  ordinal_type m = nqp;
  if (sz2 < nqp)
    m = sz2;
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> Bt(sz2, nqp);
  Teuchos::Array<ordinal_type> piv(nqp);
  Teuchos::Array<value_type> tau(m);
  ordinal_type info;
  for (ordinal_type i=0; i<nqp; i++)
    for (ordinal_type j=0; j<sz2; j++)
      Bt(j,i) = B2(i,j);

  // Workspace query
  Teuchos::Array<value_type> work(1);
  ordinal_type lwork = -1;
  DGEQP3_F77(&sz2, &nqp, Bt.values(), &sz2, &piv[0], &tau[0], &work[0], &lwork, 
	     &info);
  TEUCHOS_TEST_FOR_EXCEPTION(
    info < 0, std::logic_error, "dgeqp3 returned info = " << info);

  // Column pivoted QR
  lwork = work[0];
  work.resize(lwork);
  DGEQP3_F77(&sz2, &nqp, Bt.values(), &sz2, &piv[0], &tau[0], &work[0], &lwork, 
	     &info);
  TEUCHOS_TEST_FOR_EXCEPTION(
    info < 0, std::logic_error, "dgeqp3 returned info = " << info);

  // Determine rank
  ordinal_type rank;
  for (rank=0; rank<m; rank++)
    if (std::abs(Bt(rank,rank)) < reduction_tol)
      break;

  if (verbose)
    std::cout << "rank = " << rank << std::endl;

  // Apply b = Q^T*B2^T*w
  Teuchos::SerialDenseVector<ordinal_type,value_type> b(sz2);
  // Teuchos::SerialDenseVector<ordinal_type,value_type> w(
  //   Teuchos::View, const_cast<value_type*>(&weights[0]), nqp);
  // b.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, B2, w, 1.0);
  b.putScalar(0.0);
  b[0] = 1.0;

  //std::cout << "B2^T*w = " << b << std::endl;

  lapack.ORMQR('L', 'T', sz2, 1, sz2, Bt.values(), sz2, &tau[0], 
	       b.values(), sz2, &work[0], lwork, &info);
  TEUCHOS_TEST_FOR_EXCEPTION(
    info < 0, std::logic_error, "dormqr returned info = " << info);

  // Compute [R11^{-1}*b1; 0] where b1 is the first r rows of b
  blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, 
	    Teuchos::NON_UNIT_DIAG, rank, 1, 1.0, Bt.values(), sz2, b.values(),
	    sz2);
  for (ordinal_type i=rank; i<sz2; i++)
    b[i] = 0.0;

  // Get reduced weights, points and values
  Teuchos::SerialDenseVector<ordinal_type,value_type> wt(nqp);
  wt.putScalar(0.0);
  for (ordinal_type i=0; i<rank; i++)
    wt[piv[i]-1] = b[i];

  if (verbose)
    std::cout << "reduced weights = " << wt << std::endl;

  reduced_weights =
    Teuchos::rcp(new Teuchos::Array<value_type>(rank));
  reduced_points =
    Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(rank));
  reduced_values =
    Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(rank));
  ordinal_type idx = 0;
  for (ordinal_type i=0; i<nqp; i++) {
    if (wt[i] != 0.0) {
      (*reduced_weights)[idx] = wt[i];
      (*reduced_points)[idx].resize(d);
      for (ordinal_type j=0; j<d; j++)
	(*reduced_points)[idx][j] = F(i,j);
      (*reduced_values)[idx].resize(sz);
      for (ordinal_type j=0; j<sz; j++)
	(*reduced_values)[idx][j] = Q(i,j);
      idx++;
    }
  }
  TEUCHOS_ASSERT(idx == rank);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::MonomialGramSchmidtSimplexPCEBasis<ordinal_type, value_type>::
reducedQuadrature_GLPK(
  Teuchos::SerialDenseMatrix<ordinal_type, value_type>& B2,
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
  const Teuchos::Array<value_type>& weights,
  Teuchos::RCP< Teuchos::Array<value_type> >& reduced_weights,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_points,
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_values)
{
#ifdef HAVE_STOKHOS_GLPK
  //
  // Find reduced quadrature weights by solving linear program
  // min b^T*u s.t. Q2^T*u = e_1, u >= 0 where B2 = Q2*R2
  //
  ordinal_type nqp = B2.numRows();
  ordinal_type sz2 = B2.numCols();

  // Compute QR factorization of B2 using Gram-Schmidt
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> Q2(nqp, sz2);
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> R2(sz2, sz2);
  if (orthogonalization_method == "Classical Gram-Schmidt")
    computeQR_CGS(sz2, B2, weights, Q2, R2);
  else if (orthogonalization_method == "Modified Gram-Schmidt")
    computeQR_MGS(sz2, B2, weights, Q2, R2);
  else
    TEUCHOS_TEST_FOR_EXCEPTION(
	true, std::logic_error, 
	"Invalid orthogonalization method " << orthogonalization_method);

  // Setup linear program
  LPX *lp = lpx_create_prob();
  lpx_set_prob_name(lp, "Monomial PCE Reduction");
  lpx_set_obj_dir(lp, LPX_MIN);
  lpx_add_rows(lp, sz2);
  lpx_add_cols(lp, nqp);

  // Compute e_1 = Q2^T*w which is our constraint
  Teuchos::SerialDenseVector<ordinal_type,value_type> e1(sz2);
  // e1.putScalar(0.0);
  // e1[0] = 1.0;
  Teuchos::SerialDenseVector<ordinal_type,value_type> w(
    Teuchos::View, const_cast<value_type*>(&weights[0]), nqp);
  e1.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Q2, w, 0.0);
  for (ordinal_type i=0; i<sz2; i++)
    lpx_set_row_bnds(lp, i+1, LPX_FX, e1[i], e1[i]);

  // Set columns bounds and object coefficients
  for (ordinal_type j=0; j<nqp; j++) {
    lpx_set_col_bnds(lp, j+1, LPX_LO, 0.0, 0.0);
    lpx_set_obj_coef(lp, j+1, 1.0);
  }

  // Set constraint matrix = Q2^T
  int **cols = new int*[sz2];
  double **vals = new double*[sz2];
  for (ordinal_type i=0; i<sz2; i++) {
    cols[i] = new int[nqp+1];
    vals[i] = new double[nqp+1];
    for (ordinal_type j=0; j<nqp; j++) {
      cols[i][j+1] = j+1;
      vals[i][j+1] = Q2(j,i);
    }
    lpx_set_mat_row(lp, i+1, nqp, cols[i], vals[i]);
  }
  
  // Solve linear program
  lpx_simplex(lp);
  int status = lpx_get_status(lp);
  if (verbose) {
    std::cout << "glpk status = " << status << std::endl;
    double Z = lpx_get_obj_val(lp);
    std::cout << "glpk objective = " << Z << std::endl;
  }
  Teuchos::SerialDenseVector<ordinal_type,value_type> u(nqp);
  for (ordinal_type i=0; i<nqp; i++)
    u[i] = lpx_get_col_prim(lp, i+1);
  if (verbose)
    std::cout << "reduced weights = " << u << std::endl;

  // Clean up linear program
  lpx_delete_prob(lp);
  for (ordinal_type i=0; i<sz2; i++) {
    delete [] cols[i];
    delete [] vals[i];
  }
  delete [] cols;
  delete [] vals;
  
  ordinal_type rank = 0;
  for (ordinal_type i=0; i<nqp; i++)
    if (std::abs(u[i]) > reduction_tol) ++rank;

  // Get reduced weights, points and values
  reduced_weights =
    Teuchos::rcp(new Teuchos::Array<value_type>(rank));
  reduced_points =
    Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(rank));
  reduced_values =
    Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(rank));
  ordinal_type idx = 0;
  for (ordinal_type i=0; i<nqp; i++) {
    if (std::abs(u[i]) > reduction_tol) {
      (*reduced_weights)[idx] = u[i];
      (*reduced_points)[idx].resize(d);
      for (ordinal_type j=0; j<d; j++)
	(*reduced_points)[idx][j] = F(i,j);
      (*reduced_values)[idx].resize(sz);
      for (ordinal_type j=0; j<sz; j++)
	(*reduced_values)[idx][j] = Q(i,j);
      idx++;
    }
  }
  TEUCHOS_ASSERT(idx == rank);
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
			     "GLPK solver called but not enabled!");
#endif
}

