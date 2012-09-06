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

#include "Stokhos_SDMUtils.hpp"
#include "Stokhos_ReducedQuadratureFactory.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::MonomialGramSchmidtPCEBasis<ordinal_type, value_type>::
MonomialGramSchmidtPCEBasis(
  ordinal_type max_p,
  const Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pce,
  const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad,
  const Teuchos::ParameterList& params_) :
  name("Monomial Gram Schmidt  PCE Basis"),
  params(params_),
  pce_basis(pce[0].basis()),
  pce_sz(pce_basis->size()),
  p(max_p),
  d(pce.size()),
  verbose(params.get("Verbose", false)),
  rank_threshold(params.get("Rank Threshold", 1.0e-12)),
  orthogonalization_method(params.get("Orthogonalization Method", 
				      "Modified Gram-Schmidt"))
{
  // Check for pce's that are constant and don't represent true random
  // dimensions
  Teuchos::Array< const Stokhos::OrthogPolyApprox<ordinal_type, value_type>* > pce2;
  for (ordinal_type i=0; i<pce.size(); i++) {
    if (pce[i].standard_deviation() > 1.0e-15)
      pce2.push_back(&pce[i]);
  }
  d = pce2.size();

  // Get quadrature data
  const Teuchos::Array<value_type>& weights = quad->getQuadWeights();
  const Teuchos::Array< Teuchos::Array<value_type> >& points = 
    quad->getQuadPoints(); 
  const Teuchos::Array< Teuchos::Array<value_type> >& basis_values = 
    quad->getBasisAtQuadPoints();
  ordinal_type nqp = weights.size();

  // Original basis at quadrature points -- needed to transform expansions
  // in this basis back to original
  SDM A(nqp, pce_sz);
  for (ordinal_type i=0; i<nqp; i++)
    for (ordinal_type j=0; j<pce_sz; j++)
      A(i,j) = basis_values[i][j];

  // Compute norms of each pce for rescaling
  Teuchos::Array<value_type> pce_norms(d, 0.0);
  for (ordinal_type j=0; j<d; j++) {
    for (ordinal_type i=0; i<d; i++)
      pce_norms[j] += (*pce2[j])[i]*(*pce2[j])[i]*pce_basis->norm_squared(i);
    pce_norms[j] = std::sqrt(pce_norms[j]);
  }

  // Compute F matrix -- PCEs evaluated at all quadrature points
  // Since F is used in the reduced quadrature below as the quadrature points
  // for this reduced basis, does scaling by the pce_norms mess up the points?
  // No -- F essentially defines the random variables this basis is a function
  // of, and thus they can be scaled in any way we want.  Because we don't 
  // explicitly write the basis in terms of F, the scaling is implicit.
  SDM F(nqp, d);
  Teuchos::Array< Teuchos::Array<value_type> > values(nqp);
  for (ordinal_type i=0; i<nqp; i++) 
    for (ordinal_type j=0; j<d; j++)
      F(i,j) = pce2[j]->evaluate(points[i], basis_values[i]);

  // Build the reduced basis
  sz = buildReducedBasis(max_p, A, F, weights, terms, num_terms, Qp, Q);

  // Compute reduced quadrature rule
  Teuchos::ParameterList quad_params = params.sublist("Reduced Quadrature");
  Stokhos::ReducedQuadratureFactory<ordinal_type,value_type> quad_factory(
    quad_params);
  SDM Q2;
  if (quad_params.isParameter("Reduced Quadrature Method") &&
      quad_params.get<std::string>("Reduced Quadrature Method") == "Q2") {
    Teuchos::Array< Teuchos::Array<ordinal_type> > terms2;
    Teuchos::Array<ordinal_type> num_terms2;
    SDM Qp2;
    //ordinal_type sz2 = 
    buildReducedBasis(2*max_p, A, F, weights, terms2, num_terms2, Qp2, Q2);
  }
  reduced_quad = quad_factory.createReducedQuadrature(Q, Q2, F, weights);

  // Basis is orthonormal by construction
  norms.resize(sz, 1.0);
}

template <typename ordinal_type, typename value_type>
Stokhos::MonomialGramSchmidtPCEBasis<ordinal_type, value_type>::
~MonomialGramSchmidtPCEBasis()
{
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::MonomialGramSchmidtPCEBasis<ordinal_type, value_type>::
order() const
{
  return p;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::MonomialGramSchmidtPCEBasis<ordinal_type, value_type>::
dimension() const
{
  return d;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::MonomialGramSchmidtPCEBasis<ordinal_type, value_type>::
size() const
{
  return sz;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array<value_type>&
Stokhos::MonomialGramSchmidtPCEBasis<ordinal_type, value_type>::
norm_squared() const
{
  return norms;
}

template <typename ordinal_type, typename value_type>
const value_type&
Stokhos::MonomialGramSchmidtPCEBasis<ordinal_type, value_type>::
norm_squared(ordinal_type i) const
{
  return norms[i];
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::MonomialGramSchmidtPCEBasis<ordinal_type, value_type>::
computeTripleProductTensor(ordinal_type order) const

{
  return Teuchos::null;
}

template <typename ordinal_type, typename value_type>
value_type
Stokhos::MonomialGramSchmidtPCEBasis<ordinal_type, value_type>::
evaluateZero(ordinal_type i) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::MonomialGramSchmidtPCEBasis<ordinal_type, value_type>::
evaluateBases(const Teuchos::Array<value_type>& point,
	      Teuchos::Array<value_type>& basis_vals) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::MonomialGramSchmidtPCEBasis<ordinal_type, value_type>::
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
Stokhos::MonomialGramSchmidtPCEBasis<ordinal_type, value_type>::
getName() const
{
  return name;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::MonomialGramSchmidtPCEBasis<ordinal_type, value_type>::
transformToOriginalBasis(const value_type *in, value_type *out,
			 ordinal_type ncol, bool transpose) const
{
  if (transpose) {
    SDM zbar(Teuchos::View, const_cast<value_type*>(in), ncol, ncol, sz);
    SDM z(Teuchos::View, out, ncol, ncol, pce_sz);
    ordinal_type ret = 
      z.multiply(Teuchos::NO_TRANS, Teuchos::TRANS, 1.0, zbar, Qp, 0.0);
    TEUCHOS_ASSERT(ret == 0);
  }
  else {
    SDM zbar(Teuchos::View, const_cast<value_type*>(in), sz, sz, ncol);
    SDM z(Teuchos::View, out, pce_sz, pce_sz, ncol);
    ordinal_type ret = 
      z.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Qp, zbar, 0.0);
    TEUCHOS_ASSERT(ret == 0);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::MonomialGramSchmidtPCEBasis<ordinal_type, value_type>::
transformFromOriginalBasis(const value_type *in, value_type *out,
			 ordinal_type ncol, bool transpose) const
{
  if (transpose) {
    SDM z(Teuchos::View, const_cast<value_type*>(in), ncol, ncol, pce_sz);
    SDM zbar(Teuchos::View, out, ncol, ncol, sz);
    ordinal_type ret = 
      zbar.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, z, Qp, 0.0);
    TEUCHOS_ASSERT(ret == 0);
  }
  else {
    SDM z(Teuchos::View, const_cast<value_type*>(in), pce_sz, pce_sz, ncol);
    SDM zbar(Teuchos::View, out, sz, sz, ncol);
    ordinal_type ret = 
      zbar.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Qp, z, 0.0);
    TEUCHOS_ASSERT(ret == 0);
  }
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >
Stokhos::MonomialGramSchmidtPCEBasis<ordinal_type, value_type>::
getReducedQuadrature() const
{
  return reduced_quad;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::MonomialGramSchmidtPCEBasis<ordinal_type, value_type>::
buildReducedBasis(ordinal_type max_p, const SDM& A, const SDM& F,
		  const Teuchos::Array<value_type>& weights, 
		  Teuchos::Array< Teuchos::Array<ordinal_type> >& terms_,
		  Teuchos::Array<ordinal_type>& num_terms_,
		  SDM& Qp_, SDM& Q_)
{
  // Compute basis terms -- 2-D array giving powers for each linear index
  ordinal_type max_sz;
  CPBUtils::compute_terms(max_p, d, max_sz, terms_, num_terms_);

  // Compute B matrix -- monomials in F
  // for i=0,...,nqp-1
  //   for j=0,...,sz-1
  //      B(i,j) = F(i,1)^terms_[j][1] * ... * F(i,d)^terms_[j][d]
  // where sz is the total size of a basis up to order p and terms_[j] 
  // is an array of powers for each term in the total-order basis
  ordinal_type nqp = weights.size();
  SDM B(nqp, max_sz);
  for (ordinal_type i=0; i<nqp; i++) {
    for (ordinal_type j=0; j<max_sz; j++) {
      B(i,j) = 1.0;
      for (ordinal_type k=0; k<d; k++)
	B(i,j) *= std::pow(F(i,k), terms_[j][k]);
    }
  }

  // Rescale columns of B to have unit norm
  for (ordinal_type j=0; j<max_sz; j++) {
    value_type nrm = 0.0;
    for (ordinal_type i=0; i<nqp; i++)
      nrm += B(i,j)*B(i,j)*weights[i];
    nrm = std::sqrt(nrm);
    for (ordinal_type i=0; i<nqp; i++)
      B(i,j) /= nrm;
  }

  // Compute our new basis -- each column of Q is the new basis evaluated
  // at the original quadrature points
  ordinal_type sz_;
  
  // Compute QR factorization of B using column-pivoted QR
  // By setting the first d+1 entries of piv, we enforce that they are
  // permuted to the front of B*P
  // "Q" in the QR factorization defines the new basis
  SDM R;
  Teuchos::Array<ordinal_type> piv(max_sz);
  for (int i=0; i<d+1; i++)
    piv[i] = 1;
  if (orthogonalization_method == "Modified Gram-Schmidt")
    sz_ = CPQR_MGS_threshold(rank_threshold, B, weights, Q_, R, piv);
  else if (orthogonalization_method == "Classical Gram-Schmidt")
    sz_ = CPQR_CGS_threshold(rank_threshold, B, weights, Q_, R, piv);
  else
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, 
      "Invalid orthogonalization method " << orthogonalization_method);
  
  if (verbose) {
    std::cout << "piv = [";
    for (ordinal_type i=0; i<sz_; i++)
      std::cout << piv[i] << " ";
    std::cout << "]" << std::endl;
    
    std::cout << "diag(R) = [ ";
    for (ordinal_type i=0; i<sz_; i++)
      std::cout << R(i,i) << " ";
    std::cout << "]" << std::endl;
    
    std::cout << "rank = " << sz_ << std::endl;
    
    // Check B = Q*R
    std::cout << "||A*P-Q*R||_infty = " 
	      << Stokhos::residualCPQRError(B,Q_,R,piv) << std::endl;
    
    // Check Q_^T*W*Q_ = I
    std::cout << "||I - Q^T*W*Q||_infty = " 
	      << weightedQROrthogonalizationError(Q_, weights) << std::endl;
  }

  // Compute Qp = A^T*W*Q
  SDM tmp(nqp, sz_);
  Qp_.reshape(pce_sz, sz_);
  for (ordinal_type i=0; i<nqp; i++)
    for (ordinal_type j=0; j<sz_; j++)
      tmp(i,j) = Q_(i,j)*weights[i];
  ordinal_type ret = 
    Qp_.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, A, tmp, 0.0);
  TEUCHOS_ASSERT(ret == 0);

  // It isn't clear that Qp is orthogonal, but if you derive the projection
  // matrix from the original space to the reduced, you end up with 
  // Q^T*W*A = Qp^T

  return sz_;
}
