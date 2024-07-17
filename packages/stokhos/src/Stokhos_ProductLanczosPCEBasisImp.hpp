// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_SDMUtils.hpp"
#include "Stokhos_StieltjesPCEBasis.hpp"
#include "Stokhos_LanczosPCEBasis.hpp"
#include "Stokhos_LanczosProjPCEBasis.hpp"
#include "Stokhos_QuadratureFactory.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::ProductLanczosPCEBasis<ordinal_type, value_type>::
ProductLanczosPCEBasis(
  ordinal_type p,
  const Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pce,
  const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad,
  const Teuchos::RCP< const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk,
  const Teuchos::ParameterList& params_) :
  name("Product Lanczos PCE Basis"),
  params(params_)
{
  Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type,value_type> > pce_basis = pce[0].basis();
  ordinal_type pce_sz = pce_basis->size();

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
  ordinal_type d = coordinate_bases.size();

  // Build tensor product basis
  tensor_lanczos_basis = 
    Teuchos::rcp(
      new Stokhos::CompletePolynomialBasis<ordinal_type,value_type>(
	coordinate_bases,
	params.get("Cijk Drop Tolerance", 1.0e-15),
	params.get("Use Old Cijk Algorithm", false)));

  // Build reduced quadrature
  Teuchos::ParameterList sg_params;
  sg_params.sublist("Basis").set< Teuchos::RCP< const Stokhos::OrthogPolyBasis<ordinal_type,value_type> > >("Stochastic Galerkin Basis", tensor_lanczos_basis);
  sg_params.sublist("Quadrature") = params.sublist("Reduced Quadrature");
  reduced_quad = 
    Stokhos::QuadratureFactory<ordinal_type,value_type>::create(sg_params);

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
  ordinal_type sz = tensor_lanczos_basis->size();
  Teuchos::Array<value_type> red_basis_vals(sz);
  Teuchos::Array<value_type> pce_vals(d);
  Phi.shape(sz, nqp);
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
	jdx++;

      }

    }
    tensor_lanczos_basis->evaluateBases(pce_vals, red_basis_vals);
    for (int i=0; i<sz; i++)
      Phi(i,k) = red_basis_vals[i];
  }

  bool verbose = params.get("Verbose", false);
 
  // Compute matrix A mapping reduced space to original
  A.shape(pce_sz, sz);
  ordinal_type ret = 
    A.multiply(Teuchos::NO_TRANS, Teuchos::TRANS, 1.0, Psi, Phi, 0.0);
  TEUCHOS_ASSERT(ret == 0);
  //print_matlab(std::cout << "A = " << std::endl, A);

  // Compute pseudo-inverse of A mapping original space to reduced
  // A = U*S*V^T -> A^+ = V*S^+*U^T = (S^+*V^T)^T*U^T where 
  // S^+ is a diagonal matrix comprised of the inverse of the diagonal of S
  // for each nonzero, and zero otherwise
  Teuchos::Array<value_type> sigma;
  SDM U, Vt;
  value_type rank_threshold = params.get("Rank Threshold", 1.0e-12);
  ordinal_type rank = svd_threshold(rank_threshold, A, sigma, U, Vt);
  Ainv.shape(sz, pce_sz);
  TEUCHOS_ASSERT(rank == Vt.numRows());
  for (ordinal_type i=0; i<Vt.numRows(); i++)
    for (ordinal_type j=0; j<Vt.numCols(); j++)
      Vt(i,j) = Vt(i,j) / sigma[i];
  ret = Ainv.multiply(Teuchos::TRANS, Teuchos::TRANS, 1.0, Vt, U, 0.0);
  TEUCHOS_ASSERT(ret == 0);
  //print_matlab(std::cout << "Ainv = " << std::endl, Ainv);

  if (verbose) {
    std::cout << "rank = " << rank << std::endl;
    
    std::cout << "diag(S) = [";
    for (ordinal_type i=0; i<rank; i++)
      std::cout << sigma[i] << " ";
    std::cout << "]" << std::endl;

    // Check A = U*S*V^T
    SDM SVt(rank, Vt.numCols());
    for (ordinal_type i=0; i<Vt.numRows(); i++)
      for (ordinal_type j=0; j<Vt.numCols(); j++)
	SVt(i,j) = Vt(i,j) * sigma[i] * sigma[i];  // since we divide by sigma 
                                                   // above
    SDM err_A(pce_sz,sz);
    err_A.assign(A);
    ret = err_A.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, -1.0, U, SVt, 
			 1.0);
    TEUCHOS_ASSERT(ret == 0);
    std::cout << "||A - U*S*V^T||_infty = " << err_A.normInf() << std::endl;
    //print_matlab(std::cout << "A - U*S*V^T = " << std::endl, err_A);
 
    // Check Ainv*A == I
    SDM err(sz,sz);
    err.putScalar(0.0);
    for (ordinal_type i=0; i<sz; i++)
      err(i,i) = 1.0;
    ret = err.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Ainv, A, 
		       -1.0);
    TEUCHOS_ASSERT(ret == 0);
    std::cout << "||Ainv*A - I||_infty = " << err.normInf() << std::endl;
    //print_matlab(std::cout << "Ainv*A-I = " << std::endl, err);

    // Check A*Ainv == I
    SDM err2(pce_sz,pce_sz);
    err2.putScalar(0.0);
    for (ordinal_type i=0; i<pce_sz; i++)
      err2(i,i) = 1.0;
    ret = err2.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, Ainv, -1.0);
    TEUCHOS_ASSERT(ret == 0);
    std::cout << "||A*Ainv - I||_infty = " << err2.normInf() << std::endl;
    //print_matlab(std::cout << "A*Ainv-I = " << std::endl, err2);
  }
}

template <typename ordinal_type, typename value_type>
Stokhos::ProductLanczosPCEBasis<ordinal_type, value_type>::
~ProductLanczosPCEBasis()
{
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::ProductLanczosPCEBasis<ordinal_type, value_type>::
order() const
{
  return tensor_lanczos_basis->order();
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::ProductLanczosPCEBasis<ordinal_type, value_type>::
dimension() const
{
  return tensor_lanczos_basis->dimension();
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::ProductLanczosPCEBasis<ordinal_type, value_type>::
size() const
{
  return tensor_lanczos_basis->size();
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array<value_type>&
Stokhos::ProductLanczosPCEBasis<ordinal_type, value_type>::
norm_squared() const
{
  return tensor_lanczos_basis->norm_squared();
}

template <typename ordinal_type, typename value_type>
const value_type&
Stokhos::ProductLanczosPCEBasis<ordinal_type, value_type>::
norm_squared(ordinal_type i) const
{
  return tensor_lanczos_basis->norm_squared(i);
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::ProductLanczosPCEBasis<ordinal_type, value_type>::
computeTripleProductTensor() const

{
  Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk =
    tensor_lanczos_basis->computeTripleProductTensor();
  //std::cout << *Cijk << std::endl;
  return Cijk;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::ProductLanczosPCEBasis<ordinal_type, value_type>::
computeLinearTripleProductTensor() const

{
  Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk =
    tensor_lanczos_basis->computeLinearTripleProductTensor();
  //std::cout << *Cijk << std::endl;
  return Cijk;
}

template <typename ordinal_type, typename value_type>
value_type
Stokhos::ProductLanczosPCEBasis<ordinal_type, value_type>::
evaluateZero(ordinal_type i) const
{
  return tensor_lanczos_basis->evaluateZero(i);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ProductLanczosPCEBasis<ordinal_type, value_type>::
evaluateBases(const Teuchos::ArrayView<const value_type>& point,
	      Teuchos::Array<value_type>& basis_vals) const
{
  return tensor_lanczos_basis->evaluateBases(point, basis_vals);
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ProductLanczosPCEBasis<ordinal_type, value_type>::
print(std::ostream& os) const
{
  tensor_lanczos_basis->print(os);
}

template <typename ordinal_type, typename value_type>
const std::string&
Stokhos::ProductLanczosPCEBasis<ordinal_type, value_type>::
getName() const
{
  return name;
}

template <typename ordinal_type, typename value_type>
const Stokhos::MultiIndex<ordinal_type>&
Stokhos::ProductLanczosPCEBasis<ordinal_type, value_type>::
term(ordinal_type i) const
{
  return tensor_lanczos_basis->term(i);
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::ProductLanczosPCEBasis<ordinal_type, value_type>::
index(const Stokhos::MultiIndex<ordinal_type>& term) const
{
  return tensor_lanczos_basis->index(term);
}

template <typename ordinal_type, typename value_type>
Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type, value_type> > >
Stokhos::ProductLanczosPCEBasis<ordinal_type, value_type>::
getCoordinateBases() const
{
  return tensor_lanczos_basis->getCoordinateBases();
}

template <typename ordinal_type, typename value_type>
Stokhos::MultiIndex<ordinal_type>
Stokhos::ProductLanczosPCEBasis<ordinal_type, value_type>::
getMaxOrders() const
{
  return tensor_lanczos_basis->getMaxOrders();
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ProductLanczosPCEBasis<ordinal_type, value_type>::
transformToOriginalBasis(const value_type *in, value_type *out,
			 ordinal_type ncol, bool transpose) const
{
  ordinal_type pce_sz = A.numRows();
  ordinal_type sz = A.numCols();
  if (transpose) {
    SDM zbar(Teuchos::View, const_cast<value_type*>(in), ncol, ncol, sz);
    SDM z(Teuchos::View, out, ncol, ncol, pce_sz);
    ordinal_type ret = 
      z.multiply(Teuchos::NO_TRANS, Teuchos::TRANS, 1.0, zbar, A, 0.0);
    TEUCHOS_ASSERT(ret == 0);
  }
  else {
    SDM zbar(Teuchos::View, const_cast<value_type*>(in), sz, sz, ncol);
    SDM z(Teuchos::View, out, pce_sz, pce_sz, ncol);
    ordinal_type ret = 
      z.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, zbar, 0.0);
    TEUCHOS_ASSERT(ret == 0);
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::ProductLanczosPCEBasis<ordinal_type, value_type>::
transformFromOriginalBasis(const value_type *in, value_type *out,
			 ordinal_type ncol, bool transpose) const
{
  ordinal_type pce_sz = A.numRows();
  ordinal_type sz = A.numCols();
  if (transpose) {
    SDM z(Teuchos::View, const_cast<value_type*>(in), ncol, ncol, pce_sz);
    SDM zbar(Teuchos::View, out, ncol, ncol, sz);
    ordinal_type ret = 
      zbar.multiply(Teuchos::NO_TRANS, Teuchos::TRANS, 1.0, z, Ainv, 0.0);
    TEUCHOS_ASSERT(ret == 0);
  }
  else {
    SDM z(Teuchos::View, const_cast<value_type*>(in), pce_sz, pce_sz, ncol);
    SDM zbar(Teuchos::View, out, sz, sz, ncol);
    ordinal_type ret = 
      zbar.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Ainv, z, 0.0);
    TEUCHOS_ASSERT(ret == 0);
  }
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >
Stokhos::ProductLanczosPCEBasis<ordinal_type, value_type>::
getReducedQuadrature() const
{
  return reduced_quad;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::ProductLanczosPCEBasis<ordinal_type, value_type>::
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
