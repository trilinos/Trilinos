// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_OneDOrthogPolyBasis.hpp"
#include "Stokhos_StieltjesPCEBasis.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::StieltjesGramSchmidtBuilder<ordinal_type,value_type>::
StieltjesGramSchmidtBuilder(
  const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad_,
  const Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pces,
  ordinal_type new_order,  bool use_pce_qp, bool normalize) :
  quad(quad_)
{
  // Create array to store new coordinate bases
  ordinal_type new_dim = pces.size();
  Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type, value_type> > > new_coordinate_bases(new_dim);

  // Create Stieltjes basis for each pce
  for (ordinal_type k=0; k<new_dim; k++) {
    new_coordinate_bases[k] = Teuchos::rcp(
      new StieltjesPCEBasis<ordinal_type,value_type>(
	new_order, Teuchos::rcp(&(pces[k]),false), quad, use_pce_qp, 
	normalize)
      );
  }
  
  // Create tensor product basis from coordinate bases
  tensor_basis = Teuchos::rcp(
    new CompletePolynomialBasis<ordinal_type,value_type>(new_coordinate_bases)
    );
  
  // Use Gram-Schmidt to orthogonalize tensor product bases
  const Teuchos::Array<value_type>& weights = quad->getQuadWeights();
  const Teuchos::Array< Teuchos::Array<value_type> >& points = 
    quad->getQuadPoints();
  ordinal_type nqp = points.size();
  Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > > new_points =
    Teuchos::rcp(new Teuchos::Array< Teuchos::Array<value_type> >(nqp));
  Teuchos::RCP< Teuchos::Array<value_type> > new_weights =
    Teuchos::rcp(new Teuchos::Array<value_type>(weights));
  for (ordinal_type i=0; i<nqp; i++)
    (*new_points)[i].resize(new_dim);
  for (ordinal_type k=0; k<new_dim; k++) {
    Teuchos::Array<value_type> st_points;
    Teuchos::Array<value_type> st_weights;
    Teuchos::Array< Teuchos::Array<value_type> > st_values;
    new_coordinate_bases[k]->getQuadPoints(new_order+1, st_points, st_weights,
					   st_values);
    
    for (ordinal_type i=0; i<nqp; i++)
      (*new_points)[i][k] = st_points[i];
  }
  gs_basis = Teuchos::rcp(
    new GramSchmidtBasis<ordinal_type,value_type>(tensor_basis, 
						  *new_points, 
						  *new_weights, 
						  1e-15)
    );
  
  // Create new quadrature object
  Teuchos::RCP<const OrthogPolyBasis<ordinal_type,value_type> > new_basis =
    gs_basis;
  gs_quad = Teuchos::rcp(
    new UserDefinedQuadrature<ordinal_type,value_type>(new_basis, 
						       new_points, 
						       new_weights)
    );

}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >
Stokhos::StieltjesGramSchmidtBuilder<ordinal_type,value_type>::
getReducedBasis() const
{
  return gs_basis;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<Stokhos::Quadrature<ordinal_type, value_type> >
Stokhos::StieltjesGramSchmidtBuilder<ordinal_type,value_type>::
getReducedQuadrature() const
{
  return gs_quad;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::StieltjesGramSchmidtBuilder<ordinal_type,value_type>::
computeReducedPCEs(
  const Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pces,
  Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& new_pces)
{
  // Map pce coefficients to tensor basis to Gram-Schmidt basis
  ordinal_type dim = pces.size();
  if (new_pces.size() != pces.size())
    new_pces.resize(dim);
  for (ordinal_type k=0; k<dim; k++) {
    OrthogPolyApprox<ordinal_type,value_type> p_tensor(tensor_basis);
    p_tensor.term(k, 0) = pces[k].mean();
    p_tensor.term(k, 1) = 1.0; 
    new_pces[k].reset(gs_basis);
    gs_basis->transformCoeffs(p_tensor.coeff(), new_pces[k].coeff());
  }
}
