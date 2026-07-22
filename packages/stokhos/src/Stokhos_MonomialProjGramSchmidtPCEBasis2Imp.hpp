// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_ReducedQuadratureFactory.hpp"
#include "Stokhos_BasisFactory.hpp"
#include "Stokhos_QuadratureFactory.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"
#include "Stokhos_OrthogonalizationFactory.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::MonomialProjGramSchmidtPCEBasis2<ordinal_type, value_type>::
MonomialProjGramSchmidtPCEBasis2(
  ordinal_type max_p,
  const Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pce,
  const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad,
  const Teuchos::ParameterList& params_) :
  name("Monomial Proj Gram Schmidt PCE Basis"),
  params(params_),
  pce_basis(pce[0].basis()),
  pce_sz(pce_basis->size()),
  p(max_p),
  d(pce.size()),
  verbose(params.get("Verbose", false)),
  rank_threshold(params.get("Rank Threshold", 1.0e-12)),
  orthogonalization_method(params.get("Orthogonalization Method", 
				      "Householder"))
{
  // Check for pce's that are constant and don't represent true random
  // dimensions
  Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> > non_const_pce;
  for (ordinal_type i=0; i<pce.size(); i++) {
    if (pce[i].standard_deviation() > 1.0e-15)
      non_const_pce.push_back(pce[i]);
  }
  d = non_const_pce.size();

  // Build Q, Qp matrices
  SDM A, F;
  sz = buildQ(max_p, rank_threshold, non_const_pce, quad, terms, num_terms, 
	      Qp, A, F);
  Q.reshape(A.numRows(), sz);
  ordinal_type ret = 
    Q.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, Qp, 0.0);
  TEUCHOS_ASSERT(ret == 0);

//print_matlab(std::cout << "Qp = ", Qp);

  // Compute reduced quadrature rule
  Teuchos::ParameterList quad_params = params.sublist("Reduced Quadrature");
  Stokhos::ReducedQuadratureFactory<ordinal_type,value_type> quad_factory(
    quad_params);
  SDM Q2;
  if (quad_params.isParameter("Reduced Quadrature Method") &&
      quad_params.get<std::string>("Reduced Quadrature Method") == "Q2") {
    Teuchos::Array< Stokhos::MultiIndex<ordinal_type> > terms2;
    Teuchos::Array<ordinal_type> num_terms2;
    value_type rank_threshold2 = quad_params.get("Q2 Rank Threshold", 
						 rank_threshold);
    SDM Qp2, A2, F2;

    // Build basis, quadrature of order 2*max_p
    Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> > pce2(non_const_pce);
    Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> > quad2 = quad;
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type,value_type> > basis2 = pce_basis;
    if (2*max_p > pce_basis->order()) {
      
      // Basis
      Teuchos::RCP<const Stokhos::ProductBasis<ordinal_type,value_type> > prod_basis = Teuchos::rcp_dynamic_cast<const Stokhos::ProductBasis<ordinal_type,value_type> >(pce_basis);
      ordinal_type dim = prod_basis->dimension();
      Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(dim);
      for (ordinal_type i=0; i<dim; i++)
	bases[i] = prod_basis->getCoordinateBases()[i]->cloneWithOrder(2*max_p);
      Teuchos::RCP< const Stokhos::CompletePolynomialBasis<ordinal_type,value_type> > cp_basis2 = Teuchos::rcp(new Stokhos::CompletePolynomialBasis<ordinal_type,value_type>(bases));
      basis2 = cp_basis2;
      //quad_params.sublist("Basis").set("Stochastic Galerkin Basis", basis2);
      std::cout << "built new basis of dimension " << basis2->dimension() 
		<< " and order " << basis2->order() 
		<< " with total size " << basis2->size() << std::endl;

      // Quadrature
      // quad_params.sublist("Quadrature").set("Quadrature Order", 2*max_p);
      // quad2 = 
      // 	Stokhos::QuadratureFactory<ordinal_type,value_type>::create(quad_params);

      quad2 = Teuchos::rcp(new Stokhos::TensorProductQuadrature<ordinal_type,value_type>(cp_basis2));
      std::cout << "built new quadrature with total size " << quad2->size()
      		<< std::endl;

      // Project pce to new basis
      for (ordinal_type i=0; i<d; i++) {
      	pce2[i].reset(basis2); // this keeps lower order coeffs and sets 
      	                       // higher order ones to 0
      }
    }

    // Build Q matrix of order 2*max_p
    ordinal_type sz2 = 
      buildQ(2*max_p, rank_threshold2, pce2, quad2, terms2, num_terms2, 
	     Qp2, A2, F2);
    //print_matlab(std::cout << "Qp2 = ", Qp2);

    // Get quadrature data
    const Teuchos::Array<value_type>& weights = quad->getQuadWeights();
    const Teuchos::Array< Teuchos::Array<value_type> >& points = 
      quad->getQuadPoints(); 
    ordinal_type nqp = weights.size();

    // Original basis at quadrature points -- needed to transform expansions
    // in this basis back to original
    ordinal_type pce_sz2 = basis2->size();
    SDM AA(nqp, pce_sz2);
    Teuchos::Array<value_type> basis_vals(pce_sz2);
    for (ordinal_type i=0; i<nqp; i++) {
      basis2->evaluateBases(points[i], basis_vals);
      for (ordinal_type j=0; j<pce_sz2; j++)
    	AA(i,j) = basis_vals[j];
    }
    Q2.reshape(nqp, sz2);
    ret = Q2.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, AA, Qp2, 0.0);
    TEUCHOS_ASSERT(ret == 0);
    reduced_quad = quad_factory.createReducedQuadrature(Q, Q2, F, weights);

    // // Get quadrature data
    // const Teuchos::Array<value_type>& weights2 = quad2->getQuadWeights();
    // ordinal_type nqp2 = weights2.size();
    // Q2.reshape(nqp2, sz2);
    // ret = Q2.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A2, Qp2, 0.0);
    // TEUCHOS_ASSERT(ret == 0);
    // reduced_quad = quad_factory.createReducedQuadrature(Q, Q2, F2, weights2);
  } 
  else {
    const Teuchos::Array<value_type>& weights = quad->getQuadWeights();
    reduced_quad = quad_factory.createReducedQuadrature(Q, Q2, F, weights);
  }

  // Basis is orthonormal by construction
  norms.resize(sz, 1.0);
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::MonomialProjGramSchmidtPCEBasis2<ordinal_type, value_type>::
buildQ(
  ordinal_type max_p,
  value_type threshold, 
  const Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pce,
  const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad,
  Teuchos::Array< Stokhos::MultiIndex<ordinal_type> >& terms_,
  Teuchos::Array<ordinal_type>& num_terms_,
  SDM& Qp_, SDM& A_, SDM& F_)
{
  Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type,value_type> > pce_basis_ = pce[0].basis();
  ordinal_type pce_sz_ = pce_basis_->size();

  // Get quadrature data
  const Teuchos::Array<value_type>& weights = quad->getQuadWeights();
  const Teuchos::Array< Teuchos::Array<value_type> >& points = 
    quad->getQuadPoints(); 
  const Teuchos::Array< Teuchos::Array<value_type> >& basis_values = 
    quad->getBasisAtQuadPoints();
  ordinal_type nqp = weights.size();

  // Original basis at quadrature points -- needed to transform expansions
  // in this basis back to original
  A_.reshape(nqp, pce_sz_);
  for (ordinal_type i=0; i<nqp; i++)
    for (ordinal_type j=0; j<pce_sz_; j++)
      A_(i,j) = basis_values[i][j];

  // Compute norms of each pce for rescaling
  Teuchos::Array<value_type> pce_norms(d, 0.0);
  for (ordinal_type j=0; j<d; j++) {
    for (ordinal_type i=0; i<pce_sz_; i++)
      pce_norms[j] += (pce[j])[i]*(pce[j])[i]*pce_basis_->norm_squared(i);
    pce_norms[j] = std::sqrt(pce_norms[j]);
  }

  // Compute F matrix -- PCEs evaluated at all quadrature points
  // Since F is used in the reduced quadrature below as the quadrature points
  // for this reduced basis, does scaling by the pce_norms mess up the points?
  // No -- F essentially defines the random variables this basis is a function
  // of, and thus they can be scaled in any way we want.  Because we don't 
  // explicitly write the basis in terms of F, the scaling is implicit.
  F_.reshape(nqp, d);
  Teuchos::Array< Teuchos::Array<value_type> > values(nqp);
  for (ordinal_type i=0; i<nqp; i++) 
    for (ordinal_type j=0; j<d; j++)
      F_(i,j) = pce[j].evaluate(points[i], basis_values[i]);

  // Build the reduced basis
  // Compute basis terms -- 2-D array giving powers for each linear index
  ordinal_type max_sz;
  CPBUtils::compute_terms(max_p, d, max_sz, terms_, num_terms_);

  // Compute B matrix -- monomials in F
  // for i=0,...,nqp-1
  //   for j=0,...,sz-1
  //      B(i,j) = F(i,1)^terms_[j][1] * ... * F(i,d)^terms_[j][d]
  // where sz is the total size of a basis up to order p and terms_[j] 
  // is an array of powers for each term in the total-order basis
  SDM B(nqp, max_sz);
  for (ordinal_type i=0; i<nqp; i++) {
    for (ordinal_type j=0; j<max_sz; j++) {
      B(i,j) = 1.0;
      for (ordinal_type k=0; k<d; k++)
	B(i,j) *= std::pow(F_(i,k), terms_[j][k]);
    }
  }

  // Project B into original basis -- should use SPAM for this
  SDM Bp(pce_sz_, max_sz);
  const Teuchos::Array<value_type>& basis_norms = 
    pce_basis_->norm_squared();
  for (ordinal_type i=0; i<pce_sz_; i++) {
    for (ordinal_type j=0; j<max_sz; j++) {
      Bp(i,j) = 0.0;
      for (ordinal_type k=0; k<nqp; k++)
	Bp(i,j) += weights[k]*B(k,j)*A_(k,i);
      Bp(i,j) /= basis_norms[i];
    }
  }

  // Rescale columns of Bp to have unit norm
  for (ordinal_type j=0; j<max_sz; j++) {
    value_type nrm = 0.0;
    for (ordinal_type i=0; i<pce_sz_; i++)
      nrm += Bp(i,j)*Bp(i,j)*basis_norms[i];
    nrm = std::sqrt(nrm);
    for (ordinal_type i=0; i<pce_sz_; i++)
      Bp(i,j) /= nrm;
  }

  // Compute our new basis -- each column of Qp is the coefficients of the
  // new basis in the original basis.  Constraint pivoting so first d+1
  // columns and included in Qp.
  Teuchos::Array<value_type> w(pce_sz_, 1.0);
  SDM R;
  Teuchos::Array<ordinal_type> piv(max_sz);
  for (int i=0; i<d+1; i++)
    piv[i] = 1;
  typedef Stokhos::OrthogonalizationFactory<ordinal_type,value_type> SOF;
  ordinal_type sz_ = SOF::createOrthogonalBasis(
    orthogonalization_method, threshold, verbose, Bp, w, Qp_, R, piv);
  
  return sz_;
}

template <typename ordinal_type, typename value_type>
Stokhos::MonomialProjGramSchmidtPCEBasis2<ordinal_type, value_type>::
~MonomialProjGramSchmidtPCEBasis2()
{
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::MonomialProjGramSchmidtPCEBasis2<ordinal_type, value_type>::
order() const
{
  return p;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::MonomialProjGramSchmidtPCEBasis2<ordinal_type, value_type>::
dimension() const
{
  return d;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::MonomialProjGramSchmidtPCEBasis2<ordinal_type, value_type>::
size() const
{
  return sz;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array<value_type>&
Stokhos::MonomialProjGramSchmidtPCEBasis2<ordinal_type, value_type>::
norm_squared() const
{
  return norms;
}

template <typename ordinal_type, typename value_type>
const value_type&
Stokhos::MonomialProjGramSchmidtPCEBasis2<ordinal_type, value_type>::
norm_squared(ordinal_type i) const
{
  return norms[i];
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::MonomialProjGramSchmidtPCEBasis2<ordinal_type, value_type>::
computeTripleProductTensor() const

{
  return Teuchos::null;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::MonomialProjGramSchmidtPCEBasis2<ordinal_type, value_type>::
computeLinearTripleProductTensor() const

{
  return Teuchos::null;
}

template <typename ordinal_type, typename value_type>
value_type
Stokhos::MonomialProjGramSchmidtPCEBasis2<ordinal_type, value_type>::
evaluateZero(ordinal_type i) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::MonomialProjGramSchmidtPCEBasis2<ordinal_type, value_type>::
evaluateBases(const Teuchos::ArrayView<const value_type>& point,
	      Teuchos::Array<value_type>& basis_vals) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented!");
}

template <typename ordinal_type, typename value_type>
const std::string&
Stokhos::MonomialProjGramSchmidtPCEBasis2<ordinal_type, value_type>::
getName() const
{
  return name;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::MonomialProjGramSchmidtPCEBasis2<ordinal_type, value_type>::
print(std::ostream& os) const
{
  os << "Gram-Schmidt basis of order " << p << ", dimension " << d 
     << ", and size " << sz << ".  Matrix coefficients:\n";
  os << printMat(Qp) << std::endl;
  os << "Basis vector norms (squared):\n\t";
  for (ordinal_type i=0; i<sz; i++)
    os << norms[i] << " ";
  os << "\n";
}

template <typename ordinal_type, typename value_type>
void
Stokhos::MonomialProjGramSchmidtPCEBasis2<ordinal_type, value_type>::
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
Stokhos::MonomialProjGramSchmidtPCEBasis2<ordinal_type, value_type>::
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
Stokhos::MonomialProjGramSchmidtPCEBasis2<ordinal_type, value_type>::
getReducedQuadrature() const
{
  return reduced_quad;
}
