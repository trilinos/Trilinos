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

template <typename ordinal_type, typename value_type>
Stokhos::MonomialGramSchmidtPCEBasis<ordinal_type, value_type>::
MonomialGramSchmidtPCEBasis(
  ordinal_type max_p,
  const Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pce,
  const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad,
  const Teuchos::ParameterList& params) :
  GSReducedPCEBasisBase<ordinal_type,value_type>(max_p, pce, quad, params),
  name("Monomial Gram Schmidt PCE Basis")
{
  this->setup(max_p, pce, quad);
}

template <typename ordinal_type, typename value_type>
Stokhos::MonomialGramSchmidtPCEBasis<ordinal_type, value_type>::
~MonomialGramSchmidtPCEBasis()
{
}

template <typename ordinal_type, typename value_type>
const std::string&
Stokhos::MonomialGramSchmidtPCEBasis<ordinal_type, value_type>::
getName() const
{
  return name;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::MonomialGramSchmidtPCEBasis<ordinal_type, value_type>::
buildReducedBasis(
  ordinal_type max_p, 
  const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& A, 
  const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& F,
  const Teuchos::Array<value_type>& weights, 
  Teuchos::Array< Teuchos::Array<ordinal_type> >& terms_,
  Teuchos::Array<ordinal_type>& num_terms_,
  Teuchos::SerialDenseMatrix<ordinal_type,value_type>& Qp_, 
  Teuchos::SerialDenseMatrix<ordinal_type,value_type>& Q_)
{
  // Compute basis terms -- 2-D array giving powers for each linear index
  ordinal_type max_sz;
  CPBUtils::compute_terms(max_p, this->d, max_sz, terms_, num_terms_);

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
      for (ordinal_type k=0; k<this->d; k++)
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
  for (int i=0; i<this->d+1; i++)
    piv[i] = 1;
  if (this->orthogonalization_method == "Modified Gram-Schmidt")
    sz_ = CPQR_MGS_threshold(this->rank_threshold, B, weights, Q_, R, piv);
  else if (this->orthogonalization_method == "Classical Gram-Schmidt")
    sz_ = CPQR_CGS_threshold(this->rank_threshold, B, weights, Q_, R, piv);
  else
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, 
      "Invalid orthogonalization method " << this->orthogonalization_method);
  
  if (this->verbose) {
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
  Qp_.reshape(this->pce_sz, sz_);
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
