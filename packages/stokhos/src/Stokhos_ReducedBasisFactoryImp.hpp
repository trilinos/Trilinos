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

#include "Teuchos_TestForException.hpp"
#include "Stokhos_MonomialProjGramSchmidtPCEBasis.hpp"
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
