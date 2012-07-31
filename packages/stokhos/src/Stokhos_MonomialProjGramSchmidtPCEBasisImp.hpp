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
Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type, value_type>::
MonomialProjGramSchmidtPCEBasis(
  ordinal_type max_p,
  const Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pce,
  const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad,
  const Teuchos::ParameterList& params_) :
  name("MonomialProj Gram Schmidt  PCE Basis"),
  params(params_),
  pce_basis(pce[0].basis()),
  pce_sz(pce_basis->size()),
  p(max_p),
  d(pce.size()),
  verbose(params.get("Verbose", false)),
  rank_threshold(params.get("Rank Threshold", 1.0e-12)),
  basis_reduction_method(params.get("Basis Reduction Method", 
				    "Column-pivoted QR")),
  orthogonalization_method(params.get("Orthogonalization Method", 
				      "Householder"))
{
  // Compute basis terms -- 2-D array giving powers for each linear index
  ordinal_type max_sz;
  CPBUtils::compute_terms(max_p, d, max_sz, terms, num_terms);

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
    for (ordinal_type i=0; i<pce_sz; i++)
      pce_norms[j] += pce[j][i]*pce[j][i]*pce_basis->norm_squared(i);
    pce_norms[j] = std::sqrt(pce_norms[j]);
  }

  // Compute F matrix -- PCEs evaluated at all quadrature points
  // Since F is used in the reduced quadrature below as the quadrature points
  // for this reduced basis, does scaling by the pce_norms mess up the points?
  SDM F(nqp, d);
  Teuchos::Array< Teuchos::Array<value_type> > values(nqp);
  for (ordinal_type i=0; i<nqp; i++)
    for (ordinal_type j=0; j<d; j++)
      //F(i,j) = pce[j].evaluate(points[i], basis_values[i]) / pce_norms[j];
      F(i,j) = pce[j].evaluate(points[i], basis_values[i]);

  // Compute B matrix -- monomials in F
  // for i=0,...,nqp-1
  //   for j=0,...,sz-1
  //      B(i,j) = F(i,1)^terms[j][1] * ... * F(i,d)^terms[j][d]
  // where sz is the total size of a basis up to order p and terms[j] 
  // is an array of powers for each term in the total-order basis
  SDM B(nqp, max_sz);
  for (ordinal_type i=0; i<nqp; i++) {
    for (ordinal_type j=0; j<max_sz; j++) {
      B(i,j) = 1.0;
      for (ordinal_type k=0; k<d; k++)
	B(i,j) *= std::pow(F(i,k) / pce_norms[k], terms[j][k]);
    }
  }

  // Project B into original basis -- should use SPAM for this
  SDM Bp(pce_sz, max_sz);
  const Teuchos::Array<value_type>& basis_norms = pce_basis->norm_squared();
  for (ordinal_type i=0; i<pce_sz; i++) {
    for (ordinal_type j=0; j<max_sz; j++) {
      Bp(i,j) = 0.0;
      for (ordinal_type k=0; k<nqp; k++)
	Bp(i,j) += weights[k]*B(k,j)*A(k,i);
      Bp(i,j) /= basis_norms[i];
    }
  }

  // Rescale columns of Bp to have unit norm
  for (ordinal_type j=0; j<max_sz; j++) {
    value_type nrm = 0.0;
    for (ordinal_type i=0; i<pce_sz; i++)
      nrm += Bp(i,j)*Bp(i,j)*basis_norms[i];
    nrm = std::sqrt(nrm);
    for (ordinal_type i=0; i<pce_sz; i++)
      Bp(i,j) /= nrm;
  }

  // Compute our new basis -- each column of Qp is the coefficients of the
  // new basis in the original basis
  if (basis_reduction_method == "Column-pivoted QR") {
    // Compute QR factorization of Bp using column-pivoted QR
    // By setting the first d+1 entries of piv, we enforce that they are
    // permuted to the front of Bp*P
    // "Q" in the QR factorization defines the new basis
    Teuchos::Array<value_type> w(pce_sz, 1.0);
    SDM R;
    Teuchos::Array<ordinal_type> piv(max_sz);
    for (int i=0; i<d+1; i++)
    //for (int i=0; i<max_sz; i++)
      piv[i] = 1;
    if (orthogonalization_method == "Householder")
      sz = CPQR_Householder_threshold(rank_threshold, Bp, w, Qp, R, piv);
    else if (orthogonalization_method == "Modified Gram-Schmidt")
      sz = CPQR_MGS_threshold(rank_threshold, Bp, w, Qp, R, piv);
    else if (orthogonalization_method == "Classical Gram-Schmidt")
      sz = CPQR_CGS_threshold(rank_threshold, Bp, w, Qp, R, piv);
    else
      TEUCHOS_TEST_FOR_EXCEPTION(
	true, std::logic_error, 
	"Invalid orthogonalization method " << orthogonalization_method);

    if (verbose) {
      std::cout << "piv = [";
      for (ordinal_type i=0; i<sz; i++)
	std::cout << piv[i] << " ";
      std::cout << "]" << std::endl;
    
      std::cout << "diag(R) = [ ";
      for (ordinal_type i=0; i<sz; i++)
	std::cout << R(i,i) << " ";
      std::cout << "]" << std::endl;
      
      std::cout << "rank = " << sz << std::endl;

      // Check Bpp = Qp*R
      std::cout << "||A*P-Q*R||_infty = " 
		<< Stokhos::residualCPQRError(Bp,Qp,R,piv) << std::endl;
      
      // Check Qp^T*Qp = I
      std::cout << "||I - Q^T*Q||_infty = " 
		<< QROrthogonalizationError(Qp) << std::endl;
    }
  }
  else if (basis_reduction_method == "SVD") {
    // Compute SVD of Bp using standard SVD algorithm
    // "U" in the SVD defines the new basis
    Teuchos::Array<value_type> sigma;
    SDM Vt;
    sz = svd_threshold(rank_threshold, Bp, sigma, Qp, Vt);

    if (verbose) {
      std::cout << "diag(sigma) = [ ";
      for (ordinal_type i=0; i<sz; i++)
	std::cout << sigma[i] << " ";
      std::cout << "]" << std::endl;
      
      std::cout << "rank = " << sz << std::endl;
    }
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(
	true, std::logic_error, 
	"Invalid basis reduction method " << basis_reduction_method);

  // Evaluate new basis at original quadrature points
  Q.reshape(nqp, sz);
  ordinal_type ret = 
    Q.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, Qp, 0.0);
  TEUCHOS_ASSERT(ret == 0);

  // Compute reduced quadrature rule
  Stokhos::ReducedQuadratureFactory<ordinal_type,value_type> quad_factory(
    params.sublist("Reduced Quadrature"));
  reduced_quad = quad_factory.createReducedQuadrature(Q, F, weights);

  // Basis is orthonormal by construction
  norms.resize(sz, 1.0);
}

template <typename ordinal_type, typename value_type>
Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type, value_type>::
~MonomialProjGramSchmidtPCEBasis()
{
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type, value_type>::
order() const
{
  return p;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type, value_type>::
dimension() const
{
  return d;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type, value_type>::
size() const
{
  return sz;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array<value_type>&
Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type, value_type>::
norm_squared() const
{
  return norms;
}

template <typename ordinal_type, typename value_type>
const value_type&
Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type, value_type>::
norm_squared(ordinal_type i) const
{
  return norms[i];
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type, value_type>::
computeTripleProductTensor(ordinal_type order) const

{
  return Teuchos::null;
}

template <typename ordinal_type, typename value_type>
value_type
Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type, value_type>::
evaluateZero(ordinal_type i) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type, value_type>::
evaluateBases(const Teuchos::Array<value_type>& point,
	      Teuchos::Array<value_type>& basis_vals) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type, value_type>::
print(std::ostream& os) const
{
  os << "Gram-Schmidt basis of order " << p << ", dimension " << d 
     << ", and size " << sz << ".  Matrix coefficients:\n";
  os << Qp << std::endl;
  os << "Basis vector norms (squared):\n\t";
  for (ordinal_type i=0; i<sz; i++)
    os << norms[i] << " ";
  os << "\n";
}

template <typename ordinal_type, typename value_type>
const std::string&
Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type, value_type>::
getName() const
{
  return name;
}

template <typename ordinal_type, typename value_type>
Teuchos::Array<ordinal_type>
Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type, value_type>::
getTerm(ordinal_type i) const
{
  return terms[i];
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type, value_type>::
getIndex(const Teuchos::Array<ordinal_type>& term) const
{
  return CPBUtils::compute_index(term, terms, num_terms, p);
}

template <typename ordinal_type, typename value_type>
Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type, value_type> > >
Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type, value_type>::
getCoordinateBases() const
{
  return Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type, value_type> > >();
}

template <typename ordinal_type, typename value_type>
void
Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type, value_type>::
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
Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type, value_type>::
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
Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type, value_type>::
getReducedQuadrature() const
{
  return reduced_quad;
}

template <typename ordinal_type, typename value_type>
void 
Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type, value_type>::
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
