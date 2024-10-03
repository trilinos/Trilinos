// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_TestForException.hpp"
#include "Stokhos_MonomialProjGramSchmidtPCEBasis.hpp"
#include "Stokhos_MonomialProjGramSchmidtPCEBasis2.hpp"
#include "Stokhos_MonomialGramSchmidtPCEBasis.hpp"
#include "Stokhos_ProductLanczosPCEBasis.hpp"
#include "Stokhos_ProductLanczosGramSchmidtPCEBasis.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::ReducedBasisFactory<ordinal_type, value_type>::
ReducedBasisFactory(
  const Teuchos::ParameterList& params_) :
  params(params_),
  reduction_method(params.get("Reduced Basis Method", 
			      "Monomial Proj Gram-Schmidt"))
{
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP<Stokhos::ReducedPCEBasis<ordinal_type, value_type> >
Stokhos::ReducedBasisFactory<ordinal_type, value_type>::
createReducedBasis(
  ordinal_type p,
  const Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pce,
  const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad,
  const Teuchos::RCP< const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk) const
{
  // Compute reduced basis
  Teuchos::RCP<Stokhos::ReducedPCEBasis<ordinal_type, value_type> > red_basis;

  if (reduction_method == "Monomial Proj Gram-Schmidt")
    red_basis = Teuchos::rcp(new Stokhos::MonomialProjGramSchmidtPCEBasis<ordinal_type,value_type>(p, pce, quad, params));

  else if (reduction_method == "Monomial Proj Gram-Schmidt2")
    red_basis = Teuchos::rcp(new Stokhos::MonomialProjGramSchmidtPCEBasis2<ordinal_type,value_type>(p, pce, quad, params));

  else if (reduction_method == "Monomial Gram-Schmidt")
    red_basis = Teuchos::rcp(new Stokhos::MonomialGramSchmidtPCEBasis<ordinal_type,value_type>(p, pce, quad, params));

  else if (reduction_method == "Product Lanczos")
    red_basis = Teuchos::rcp(new Stokhos::ProductLanczosPCEBasis<ordinal_type,value_type>(p, pce, quad, Cijk, params));

  else if (reduction_method == "Product Lanczos Gram-Schmidt")
    red_basis = Teuchos::rcp(new Stokhos::ProductLanczosGramSchmidtPCEBasis<ordinal_type,value_type>(p, pce, quad, Cijk, params));

  else
    TEUCHOS_TEST_FOR_EXCEPTION(
	true, std::logic_error, 
	"Invalid reduced basis method " << reduction_method);

  return red_basis;
}
