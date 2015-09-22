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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
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
