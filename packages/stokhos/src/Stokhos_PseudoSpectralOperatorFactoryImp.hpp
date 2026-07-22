// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_BasisFactory.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_TensorProductPseudoSpectralOperator.hpp"
#include "Stokhos_SmolyakPseudoSpectralOperator.hpp"
#include "Stokhos_QuadraturePseudoSpectralOperator.hpp"
#include "Stokhos_QuadratureFactory.hpp"

template <typename ordinal_type, typename value_type>
Teuchos::RCP<const typename Stokhos::PseudoSpectralOperatorFactory<ordinal_type, value_type>::psop_type>
Stokhos::PseudoSpectralOperatorFactory<ordinal_type, value_type>::
create(Teuchos::ParameterList& sgParams)
{
  // Check if operator is already there
  Teuchos::ParameterList& psopParams = 
    sgParams.sublist("Pseudospectral Operator");
  Teuchos::RCP<const psop_type> psop = 
    psopParams.template get< Teuchos::RCP<const psop_type> >("Stochastic Galerkin Pseudospectral Operator", Teuchos::null);
  if (psop != Teuchos::null)
    return psop;

  // Get basis
  Teuchos::ParameterList& basisParams = sgParams.sublist("Basis");
  Teuchos::RCP< const OrthogPolyBasis<ordinal_type,value_type> > basis;
  if (basisParams.template isType< Teuchos::RCP< const OrthogPolyBasis<ordinal_type,value_type> > >("Stochastic Galerkin Basis"))
    basis = basisParams.template get< Teuchos::RCP<const OrthogPolyBasis<ordinal_type,value_type> > >("Stochastic Galerkin Basis");
  else
    basis = BasisFactory<ordinal_type,value_type>::create(sgParams);

  // Create operator
  std::string type = psopParams.get("Type", "Tensor Product");

  if (type == "Tensor Product") {
    bool use_pst = psopParams.get("Use PST", false);
    Teuchos::RCP<const ProductBasis<ordinal_type,value_type> > product_basis = Teuchos::rcp_dynamic_cast<const ProductBasis<ordinal_type,value_type> >(basis, true);
    psop = 
      Teuchos::rcp(new TensorProductPseudoSpectralOperator<ordinal_type,value_type>(*product_basis, use_pst));
  }

  else if (type == "Smolyak") {
    bool use_pst = psopParams.get("Use PST", true);
    bool use_smolyak = psopParams.get("Use Smolyak Apply", true);
    Teuchos::RCP<const SmolyakBasis<ordinal_type,value_type> > smolyak_basis = Teuchos::rcp_dynamic_cast<const SmolyakBasis<ordinal_type,value_type> >(basis, true);
    psop = 
      Teuchos::rcp(new SmolyakPseudoSpectralOperator<ordinal_type,value_type>(
		     *smolyak_basis, use_smolyak, use_pst));
  }

  else if (type == "Quadrature") {
    Teuchos::ParameterList& quadParams = sgParams.sublist("Quadrature");
    Teuchos::RCP<const Stokhos::Quadrature<ordinal_type,value_type> > quad;
    if (quadParams.template isType<Teuchos::RCP<const Stokhos::Quadrature<ordinal_type,value_type> > >("Stochastic Galerkin Quadrature"))
      quad = quadParams.template get<Teuchos::RCP<const Stokhos::Quadrature<ordinal_type,value_type> > >("Stochastic Galerkin Quadrature");
    else {
      quad = 
	Stokhos::QuadratureFactory<ordinal_type,value_type>::create(sgParams);
      quadParams.set("Stochastic Galerkin Quadrature", quad);
    }
    psop = 
      Teuchos::rcp(new QuadraturePseudoSpectralOperator<ordinal_type,value_type>(*basis, *quad));
  }

  else
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, Teuchos::Exceptions::InvalidParameter,
      std::endl << 
      "Invalid pseudospectral operator type  " << type << std::endl);

  psopParams.set("Stochastic Galerkin Pseudospectral Operator", psop);
  return psop;
}
