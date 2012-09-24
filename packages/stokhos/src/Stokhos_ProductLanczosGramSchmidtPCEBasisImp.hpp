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
#include "Stokhos_StieltjesPCEBasis.hpp"
#include "Stokhos_LanczosPCEBasis.hpp"
#include "Stokhos_LanczosProjPCEBasis.hpp"
#include "Stokhos_QuadratureFactory.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_OrthogonalizationFactory.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::ProductLanczosGramSchmidtPCEBasis<ordinal_type, value_type>::
ProductLanczosGramSchmidtPCEBasis(
  ordinal_type max_p,
  const Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pce,
  const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad,
  const Teuchos::RCP< const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk,
  const Teuchos::ParameterList& params_) :
  name("Product Lanczos Gram-Schmidt PCE Basis"),
  params(params_),
  p(max_p)
{
  Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type,value_type> > pce_basis = pce[0].basis();
  pce_sz = pce_basis->size();

  // Check if basis is a product basis
  Teuchos::RCP<const Stokhos::ProductBasis<ordinal_type,value_type> > prod_basis = Teuchos::rcp_dynamic_cast<const Stokhos::ProductBasis<ordinal_type,value_type> >(pce_basis);
  Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type,value_type> > > coord_bases;
  if (prod_basis != Teuchos::null)
    coord_bases = prod_basis->getCoordinateBases();

  // Build Lanczos basis for each pce
  bool project = params.get("Project", true);
  bool normalize = params.get("Normalize", true);
  bool limit_integration_order = params.get("Limit Integration Order", false);
  bool use_stieltjes = params.get("Use Old Stieltjes Method", false);
  Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double > > > coordinate_bases;
  Teuchos::Array<int> is_invariant(pce.size(),-2);
  for (ordinal_type i=0; i<pce.size(); i++) {

    // Check for pce's lying in invariant subspaces, which are pce's that
    // depend on only a single dimension.  In this case use the corresponding
    // original coordinate basis.  Convention is:  -2 -- not invariant, -1 --
    // constant, i >= 0 pce depends only on dimension i
    if (prod_basis != Teuchos::null)
      is_invariant[i] = isInvariant(pce[i]);
    if (is_invariant[i] >= 0) {
      coordinate_bases.push_back(coord_bases[is_invariant[i]]);
    }

    // Exclude constant pce's from the basis since they don't represent
    // stochastic dimensions
    else if (is_invariant[i] != -1) {
      if (use_stieltjes) {
	coordinate_bases.push_back(
	  Teuchos::rcp(
	    new Stokhos::StieltjesPCEBasis<ordinal_type,value_type>(
	      p, Teuchos::rcp(&(pce[i]),false), quad, false,
	      normalize, project, Cijk)));
      }
      else {
	if (project) 
	  coordinate_bases.push_back(
	    Teuchos::rcp(
	      new Stokhos::LanczosProjPCEBasis<ordinal_type,value_type>(
		p, Teuchos::rcp(&(pce[i]),false), Cijk,
		normalize, limit_integration_order)));
	else
	  coordinate_bases.push_back(
	    Teuchos::rcp(
	      new Stokhos::LanczosPCEBasis<ordinal_type,value_type>(
		p, Teuchos::rcp(&(pce[i]),false), quad,
		normalize, limit_integration_order)));
      }
    }
  }
  d = coordinate_bases.size();

  // Build tensor product basis
  tensor_lanczos_basis = 
    Teuchos::rcp(
      new Stokhos::CompletePolynomialBasis<ordinal_type,value_type>(
	coordinate_bases,
	params.get("Cijk Drop Tolerance", 1.0e-15),
	params.get("Use Old Cijk Algorithm", false)));

  // Build Psi matrix -- Psi_ij = Psi_i(x^j)*w_j/<Psi_i^2>
  const Teuchos::Array<value_type>& weights = quad->getQuadWeights();
  const Teuchos::Array< Teuchos::Array<value_type> >& points = 
    quad->getQuadPoints(); 
  const Teuchos::Array< Teuchos::Array<value_type> >& basis_vals = 
    quad->getBasisAtQuadPoints();
  ordinal_type nqp = weights.size();
  SDM Psi(pce_sz, nqp);
  for (ordinal_type i=0; i<pce_sz; i++)
    for (ordinal_type k=0; k<nqp; k++)
      Psi(i,k) = basis_vals[k][i]*weights[k]/pce_basis->norm_squared(i);

  // Build Phi matrix -- Phi_ij = Phi_i(y(x^j))
  sz = tensor_lanczos_basis->size();
  Teuchos::Array<value_type> red_basis_vals(sz);
  Teuchos::Array<value_type> pce_vals(d);
  SDM Phi(sz, nqp);
  SDM F(nqp, d);
  for (int k=0; k<nqp; k++) {
    ordinal_type jdx = 0;
    for (int j=0; j<pce.size(); j++) {

      // Exclude constant pce's
      if (is_invariant[j] != -1) {

	// Use the identity mapping for invariant subspaces
	if (is_invariant[j] >= 0)
	  pce_vals[jdx] = points[k][is_invariant[j]];
	else
	  pce_vals[jdx] = pce[j].evaluate(points[k], basis_vals[k]);
	F(k,jdx) = pce_vals[jdx];
	jdx++;

      }

    }
    tensor_lanczos_basis->evaluateBases(pce_vals, red_basis_vals);
    for (int i=0; i<sz; i++)
      Phi(i,k) = red_basis_vals[i];
  }

  bool verbose = params.get("Verbose", false);
 
  // Compute matrix A mapping reduced space to original
  SDM A(pce_sz, sz);
  ordinal_type ret = 
    A.multiply(Teuchos::NO_TRANS, Teuchos::TRANS, 1.0, Psi, Phi, 0.0);
  TEUCHOS_ASSERT(ret == 0);

  // Rescale columns of A to have unit norm
  const Teuchos::Array<value_type>& basis_norms = pce_basis->norm_squared();
  for (ordinal_type j=0; j<sz; j++) {
    value_type nrm = 0.0;
    for (ordinal_type i=0; i<pce_sz; i++)
      nrm += A(i,j)*A(i,j)*basis_norms[i];
    nrm = std::sqrt(nrm);
    for (ordinal_type i=0; i<pce_sz; i++)
      A(i,j) /= nrm;
  }

  // Compute our new basis -- each column of Qp is the coefficients of the
  // new basis in the original basis. Constraint pivoting so first d+1
  // columns and included in Qp.
  value_type rank_threshold = params.get("Rank Threshold", 1.0e-12);
  std::string  orthogonalization_method = 
    params.get("Orthogonalization Method", "Householder");
  Teuchos::Array<value_type> w(pce_sz, 1.0);
  SDM R;
  Teuchos::Array<ordinal_type> piv(sz);
  for (int i=0; i<d+1; i++)
    piv[i] = 1;
  typedef Stokhos::OrthogonalizationFactory<ordinal_type,value_type> SOF;
  sz = SOF::createOrthogonalBasis(
    orthogonalization_method, rank_threshold, verbose, A, w, Qp, R, piv);

  // Original basis at quadrature points -- needed to transform expansions
  // in this basis back to original
  SDM B(nqp, pce_sz);
  for (ordinal_type i=0; i<nqp; i++)
    for (ordinal_type j=0; j<pce_sz; j++)
      B(i,j) = basis_vals[i][j];

  // Evaluate new basis at original quadrature points
  Q.reshape(nqp, sz);
  ret = Q.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, B, Qp, 0.0);
  TEUCHOS_ASSERT(ret == 0);

  // Compute reduced quadrature rule
  Stokhos::ReducedQuadratureFactory<ordinal_type,value_type> quad_factory(
    params.sublist("Reduced Quadrature"));
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> Q2;
  reduced_quad = quad_factory.createReducedQuadrature(Q, Q2, F, weights);

  // Basis is orthonormal by construction
  norms.resize(sz, 1.0);
}

template <typename ordinal_type, typename value_type>
Stokhos::ProductLanczosGramSchmidtPCEBasis<ordinal_type, value_type>::
~ProductLanczosGramSchmidtPCEBasis()
{
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::ProductLanczosGramSchmidtPCEBasis<ordinal_type, value_type>::
order() const
{
  return p;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::ProductLanczosGramSchmidtPCEBasis<ordinal_type, value_type>::
dimension() const
{
  return d;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::ProductLanczosGramSchmidtPCEBasis<ordinal_type, value_type>::
size() const
{
  return sz;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array<value_type>&
Stokhos::ProductLanczosGramSchmidtPCEBasis<ordinal_type, value_type>::
norm_squared() const
{
  return norms;
}

template <typename ordinal_type, typename value_type>
const value_type&
Stokhos::ProductLanczosGramSchmidtPCEBasis<ordinal_type, value_type>::
norm_squared(ordinal_type i) const
{
  return norms[i];
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::ProductLanczosGramSchmidtPCEBasis<ordinal_type, value_type>::
computeTripleProductTensor(ordinal_type order) const

{
  return Teuchos::null;
}

template <typename ordinal_type, typename value_type>
value_type
Stokhos::ProductLanczosGramSchmidtPCEBasis<ordinal_type, value_type>::
evaluateZero(ordinal_type i) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ProductLanczosGramSchmidtPCEBasis<ordinal_type, value_type>::
evaluateBases(const Teuchos::Array<value_type>& point,
	      Teuchos::Array<value_type>& basis_vals) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented!");
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ProductLanczosGramSchmidtPCEBasis<ordinal_type, value_type>::
print(std::ostream& os) const
{
  os << "Product Lanczos Gram-Schmidt basis of order " << p << ", dimension " 
     << d 
     << ", and size " << sz << ".  Matrix coefficients:\n";
  os << Qp << std::endl;
  os << "Basis vector norms (squared):\n\t";
  for (ordinal_type i=0; i<sz; i++)
    os << norms[i] << " ";
  os << "\n";
}

template <typename ordinal_type, typename value_type>
const std::string&
Stokhos::ProductLanczosGramSchmidtPCEBasis<ordinal_type, value_type>::
getName() const
{
  return name;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ProductLanczosGramSchmidtPCEBasis<ordinal_type, value_type>::
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
Stokhos::ProductLanczosGramSchmidtPCEBasis<ordinal_type, value_type>::
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
Stokhos::ProductLanczosGramSchmidtPCEBasis<ordinal_type, value_type>::
getReducedQuadrature() const
{
  return reduced_quad;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::ProductLanczosGramSchmidtPCEBasis<ordinal_type, value_type>::
isInvariant(const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& pce) const
{
  Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type,value_type> > basis =
    pce.basis();
  ordinal_type dim = basis->dimension();
  value_type tol = 1.0e-15;

  // Check if basis is a product basis
  Teuchos::RCP<const Stokhos::ProductBasis<ordinal_type,value_type> > prod_basis = Teuchos::rcp_dynamic_cast<const Stokhos::ProductBasis<ordinal_type,value_type> >(basis);
  if (prod_basis == Teuchos::null)
    return -2;

  // Build list of dimensions pce depends on by looping over each dimension, 
  // computing norm of pce with just that dimension -- note we don't include
  // the constant term
  Teuchos::Array<ordinal_type> dependent_dims;
  tmp_pce.reset(basis);
  for (ordinal_type i=0; i<dim; i++) {
    ordinal_type p = prod_basis->getCoordinateBases()[i]->order();
    tmp_pce.init(0.0);
    for (ordinal_type j=1; j<=p; j++)
      tmp_pce.term(i,j) = pce.term(i,j);
    value_type nrm = tmp_pce.two_norm();
    if (nrm > tol) dependent_dims.push_back(i);
  }

  // If dependent_dims has length 1, pce a function of a single variable,
  // which is an invariant subspace
  if (dependent_dims.size() == 1)
    return dependent_dims[0];

  // If dependent_dims has length 0, pce is constant
  else if (dependent_dims.size() == 0)
    return -1;

  // Otherwise pce depends on more than one variable
  return -2;
}
